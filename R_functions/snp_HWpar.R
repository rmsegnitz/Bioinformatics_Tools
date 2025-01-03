######################################
# Hardy-Weinberg Equilibrium Testing # 
######################################

# The assumption of Hardy-Weinberg equilibrium (HWE) was tested for each SNP in the dataset using allele frequencies to calculate the
# chi-square test statistic with continuity correction. The p values for the chi-square statistic were computed as the probability of 
# observing a chi-square statistic equal to or greater than the calculated value, under the null hypothesis that the population is in 
# Hardy-Weinberg equilibrium.  Parallel processing was implemented using the parallel package in R to improve computational efficiency.

# Inputs: The expected in put is a matrix with snps/features in rows and individuals in columns. snp data should be encoded as minor 
# allele counts (0, 1, 2). 
# Returns: Returns a vector of HWE p values of length equal to input rows. 


HWEmat <- function(genotype_matrix, parallel = TRUE, cores = 2, chunk_size = 1e5) {
  # Ensure the input matrix is integer for memory efficiency
  if (!is.integer(genotype_matrix)) {
    genotype_matrix <- as.integer(genotype_matrix)
  }
  
  # Internal function to calculate HWE p-values for a chunk of SNPs
  calc_hwe_chunk <- function(chunk) {
    # Observed counts for all SNPs
    obs_counts <- sapply(0:2, function(x) rowSums(chunk == x))
    
    # Total individuals and allele frequencies
    total_individuals <- rowSums(obs_counts)
    total_alleles <- 2 * total_individuals
    p <- (2 * obs_counts[, 1] + obs_counts[, 2]) / total_alleles
    q <- 1 - p
    
    # Expected counts
    exp_counts <- cbind(
      total_individuals * p^2,
      total_individuals * 2 * p * q,
      total_individuals * q^2
    )
    
    # Chi-square test statistics with continuity correction
    chi_sq <- rowSums(((abs(obs_counts - exp_counts) - 0.5)^2) / exp_counts, na.rm = TRUE)
    
    # P-values
    pchisq(chi_sq, df = 1, lower.tail = FALSE)
  }
  
  # Split SNPs into chunks
  num_snps <- nrow(genotype_matrix)
  snp_indices <- split(seq_len(num_snps), ceiling(seq_len(num_snps) / chunk_size))
  
  # Parallel or sequential processing
  if (parallel) {
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = c("calc_hwe_chunk"), envir = environment())
    pvals <- unlist(parLapply(cl, snp_indices, function(indices) {
      calc_hwe_chunk(genotype_matrix[indices, , drop = FALSE])
    }))
    stopCluster(cl)
  } else {
    pvals <- unlist(lapply(snp_indices, function(indices) {
      calc_hwe_chunk(genotype_matrix[indices, , drop = FALSE])
    }))
  }
  
  return(pvals)
}
