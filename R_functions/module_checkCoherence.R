
#---------------------------
# CHECK MODULE COHERENCE
#--------------------------

# Adapted from Matt Altman's code

module_checkCoherence<-function(module_gene_sets, voom_obj, geneSet = "geneSet",module_set="STUDY", sample_set="STUDY", remove_sets = NULL){

# Load libraries
require(gtools)
require(tidyverse) 
`%notin%`<-Negate(`%in%`)
  
moduleSets<- module_gene_sets
voomData<-voom_obj

#-----------------------------  
# Calculate the gene set means (module expression)
#-----------------------------
uniGS<-moduleSets$geneSet%>%unique()%>%as.character()%>%mixedsort() # Pull unique modules & order correctly

GSsub<-c()
GSabsent<-list()

for (i in uniGS){
  curGS<-i
  curIDs<-as.character(moduleSets$ensemblID[which(moduleSets$geneSet==curGS)])
  matchIndex<-match(curIDs, rownames(voomData))
  if (any(is.na(matchIndex)))
    matchIndex<-matchIndex[-which(is.na(matchIndex))]
    GSabsent[[i]]<-curIDs[which(curIDs %notin% rownames(voomData))]
  if (length(curIDs) > 1)
    GSsub<-rbind(GSsub, apply(voomData$E[matchIndex,], 2, mean))
  else
    GSsub<-rbind(GSsub, voomData$E[matchIndex,])
}
rownames(GSsub)<-uniGS

if(length(unlist(GSabsent))>0){
print(paste0(length(unique(unlist(GSabsent))), " Module Genes absent in voom object. Missing genes accesible in .$GSabsent"))}

#----------------------------
#  Set function for correlation p values
#----------------------------

# Correlation p value function 
cor2pvalue = function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}


# remove gene sets that are zero
if(!is.null(remove_sets)){uniGS<-uniGS[-which(uniGS %in% remove_sets)]}

# Setup storage
SubGeneCorDF<-c()
SubGeneCorMedian<-c()
SubGenePMedian<-c()

for (i in uniGS){
  curSet<-i
  curIDs<-as.character(moduleSets$ensemblID[which(moduleSets$geneSet==curSet)])
  matchIndex<-match(curIDs, rownames(voomData))
  if (any(is.na(matchIndex))){
    matchIndex<-matchIndex[-which(is.na(matchIndex))]}
  
  setCor<-cor(t(voomData$E[matchIndex,])) # Calculate pairwise gene correlations from expression data
  n<-t(!is.na(t(voomData$E[matchIndex,]))) %*% (!is.na(t(voomData$E[matchIndex,]))) # create matrix with sample size to calculate correlation p (dim=curIDs*curIDs, value = sample size)
  
  allCorInfo<-cor2pvalue(setCor, n)
  setP<-allCorInfo$p
  allCorVals<-setCor[upper.tri(setCor)] # Extract upper triangle correlation values
  allPVals<-setP[upper.tri(setP)] # Extract upper triangle correlation p values
  
  SubGeneCorDF<-rbind(SubGeneCorDF, cbind(allCorVals, allPVals, rep(curSet, length(allCorVals))))
}



# Assemble and format results dataframe
SubGeneCorDF<-as.data.frame(SubGeneCorDF)%>%
  rename(Cor = allCorVals, P = allPVals, Set = V3)%>%
  mutate(Cor = as.numeric(as.character(Cor)))%>%
  mutate(P = as.numeric(as.character(P)))%>%
  mutate(Set = factor(Set, levels =uniGS))%>%
  dplyr::select(Set, Cor, P)

SubGeneCorDF_summary<-
  SubGeneCorDF%>%
  group_by(Set)%>%
  summarise(median_mod_correlation = median(Cor),
            median_mod_p = median(P))


#Boxplot of all pairwise correlations per module
coherence_boxplot_cor<-
SubGeneCorDF%>%
  ggplot(aes(y=Cor, x=Set))+
  geom_hline(yintercept =0, color="black")+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  ylab("Correlation Strength")+
  xlab(paste0(module_set, " Modules"))+
  labs(title = paste0(module_set, " Module Coherence in ", sample_set, " Samples"),
       subtitle = "Inter-Gene PearsonCorrelation")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5))

coherence_boxplot_p<-
  SubGeneCorDF%>%
  ggplot(aes(y=-log10(P), x=Set))+
  geom_hline(yintercept =0, color="black")+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = -log10(0.05), color="red", linetype="dashed")+
  ylab("-log10(Correlation P Value)")+
  xlab(paste0(module_set, " Modules"))+
  labs(title = paste0(module_set, " Module Coherence in ", sample_set, " Samples"),
       subtitle = "Inter-Gene PearsonCorrelation")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5))

# Or with facets
coherence_boxplot_faceted<-
  SubGeneCorDF%>%
  mutate(negLogP = -log10(P))%>%
  pivot_longer(cols = c(Cor, negLogP))%>%
  mutate(h_line = ifelse(name == "negLogP", -log10(0.05), NA))%>%
  ggplot(aes(y=value, x=Set))+
  geom_hline(yintercept =0, color="black")+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(aes(yintercept = h_line), color = "red", linetype="dashed")+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  #ylab("Inter-Gene PearsonCorrelation")+
  xlab(paste0(module_set, " Modules"))+
  #ggtitle(paste0(module_set, " Module Coherence in ", sample_set, " Samples"))+
  labs(title = paste0(module_set, " Module Coherence in ", sample_set, " Samples"),
       subtitle = "Inter-Gene PearsonCorrelation")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5))+
  facet_wrap(~name, scales = "free_y", strip.position = "left", labeller = labeller(name = c("Cor" = "Correlation Strength", "negLogP" = "-log10(Correlation P Value)")))


print(coherence_boxplot_faceted)

return(list(coherence_boxplot_combined = coherence_boxplot_faceted, 
            coherence_boxplot_cor = coherence_boxplot_cor, 
            coherence_boxplot_p = coherence_boxplot_p,  
            subgene_correlation_df = SubGeneCorDF))

}
