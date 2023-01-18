

# Function for calculating Hardy Weinberg P values



# Create sample data for troubleshooting
# snp_mat<-matrix( sample(c(0,1,2), replace = T, size = 5000*50), nrow=5000, ncol=50)
# colnames(snp_mat)<-paste("id_", 1:50)
# rownames(snp_mat)<-paste("snp_", 1:5000)
# snp_mat<-as.data.frame(snp_mat)%>%rownames_to_column("snp_id")

require(dplyr)
require(pbapply)
require(HardyWeinberg)


calcHWp<-function(snp_mat, snp_col="snpID", method="chisq"){


  
# Convert snp data.frame to matrix of genotype counts 
# Note: this could/should be sped up with data.table notation
print("Converting SNP matrix to genotype counts.")
  
geno_count_matrix<-
  snp_mat %>%
  dplyr::rename("snpID" = snp_col)%>%
  tidyr::pivot_longer(-c(snpID)) %>% # pivot to long format
  dplyr::mutate(value=factor(value, levels=c("0", "1", "2")))%>% # code genotype as factor for counting
  dplyr::group_by(snpID)%>% # perform tallies by SNP
  dplyr::count(value, .drop = FALSE) %>% # count genotypes, including empties
  dplyr::filter(!is.na(value))%>% # remove missing genotypes
  dplyr::mutate(value = dplyr::recode(value, `0` ="AA", `1`="AB", `2`="BB"))%>% # recode for better column names
  tidyr::pivot_wider(names_from = value, values_from = n)%>% # pivot wider
  dplyr::mutate_at(c('AA', 'AB', 'BB'), as.numeric)%>%
  tibble::column_to_rownames("snpID")%>%
  as.matrix()%>%
  suppressWarnings()

# Calculate Hardy-Weinberg
print("Calculating Hardy-Weinberg Chisq P Vals.")
if(method=="chisq"){
snp_HWp<-
  as.data.frame(geno_count_matrix)%>%
  tibble::rownames_to_column(snp_col)%>%
  dplyr::mutate(HW_p = pbapply(geno_count_matrix, 1, function(x) paste(HWChisq(x, verbose = F)[c("chisq", "D", "pval")], collapse="v")))%>% # Calculate HW metrics with progress bar
  tidyr::separate(HW_p, into = c("HW_chisq", "HW_D", "HW_p"), sep="v") # split metrics
}else if(method == "exact"){
  
snp_HWp<-
  as.data.frame(geno_count_matrix)%>%
  tibble::rownames_to_column(snp_col)%>%
  dplyr::mutate(HW_p = pbapply(geno_count_matrix, 1, function(x) HWExact(x, verbose = F)$pval))%>% # Calculate HW metrics with progress bar
  tidyr::separate(HW_p, into = c("HW_stat",  "HW_p"), sep="v") # split metrics
}

return(snp_HWp)
}


# Try to rewrite code using data.table

# geno_count_matrix_DT<-
#   snp_mat %>%
#   dplyr::rename("snpID" = snp_col)%>% # rename snpID column for convenience
#   setDT()%>% # convert to data.table for faster operations
#   melt(id.vars = "snpID", # "pivot" table to long format
#        measure.vars = colnames(snp_mat)[!colnames(snp_mat) == snp_col],
#        variable.name = "name",
#        value.name    = "value")
# 
# geno_count_matrix_DT[ , value := factor(value, levels=c("0", "1", "2"))] # convert genotypes to factor for counting
# 
# geno_count_matrix_DT[, .N, by=.(snpID, value)]  
#   
#   # pivot_longer(-c(snpID)) %>% # pivot to long format
#   # mutate(value=factor(value, levels=c("0", "1", "2")))%>% # code genotype as factor for counting
#   group_by(snpID)%>% # perform tallies by SNP
#   count(value, .drop = FALSE) %>% # count genotypes, including empties
#   filter(!is.na(value))%>% # remove missing genotypes
#   mutate(value = dplyr::recode(value, `0` ="AA", `1`="AB", `2`="BB"))%>% # recode for better column names
#   pivot_wider(names_from = value, values_from = n)%>% # pivot wider
#   mutate_at(c('AA', 'AB', 'BB'), as.numeric)%>%
#   column_to_rownames("snpID")%>%
#   as.matrix()
