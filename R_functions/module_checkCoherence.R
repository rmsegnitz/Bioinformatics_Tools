
#---------------------------
# CHECK MODULE COHERENCE
#--------------------------

# Adapted from Matt Altman's code

module_checkCoherence<-function(module_gene_sets, voom_obj, geneSet = "geneSet",module_set="STUDY", sample_set="STUDY", remove_sets = NULL){

moduleSets<- module_gene_sets
voomData<-voom_obj
  
# calculate the gene set means
uniGS<-as.character(unique(moduleSets$geneSet))
uniGS<-sort(uniGS)
nasalGSsub<-c()
for (i in 1:length(uniGS)){
  curGS<-uniGS[i]
  curIDs<-as.character(moduleSets$ensemblID[which(moduleSets$geneSet==curGS)])
  matchIndex<-match(curIDs, rownames(voomData))
  if (any(is.na(matchIndex)))
    matchIndex<-matchIndex[-which(is.na(matchIndex))]
  if (length(curIDs) > 1)
    nasalGSsub<-rbind(nasalGSsub, apply(voomData$E[matchIndex,], 2, mean))
  else
    nasalGSsub<-rbind(nasalGSsub, voomData$E[matchIndex,])
}
rownames(nasalGSsub)<-uniGS





# Correlation p value function 
cor2pvalue = function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}

uniNasal<-sort(unique(as.character(moduleSets$geneSet)))
# remove gene sets that are zero
if(!is.null(remove_sets)){remSets<-remove_sets
uniNasal<-uniNasal[-which(uniNasal %in% remSets)]}

# Setup storage
nasalSubGeneCorDF<-c()
nasalSubGeneCorMedian<-c()
nasalSubGenePMedian<-c()

for (i in 1:length(uniNasal)){
  curSet<-uniNasal[i]
  curIDs<-as.character(moduleSets$ensemblID[which(moduleSets$geneSet==curSet)])
  matchIndex<-match(curIDs, rownames(voomData))
  if (any(is.na(matchIndex)))
    matchIndex<-matchIndex[-which(is.na(matchIndex))]
  
  setCor<-cor(t(voomData$E[matchIndex,]))
  n<-t(!is.na(t(voomData$E[matchIndex,]))) %*% (!is.na(t(voomData$E[matchIndex,])))
  
  allCorInfo<-cor2pvalue(setCor, n)
  setP<-allCorInfo$p
  allCorVals<-setCor[upper.tri(setCor)]
  allPVals<-setP[upper.tri(setP)]
  
  nasalSubGeneCorDF<-rbind(nasalSubGeneCorDF, cbind(allCorVals, rep(curSet, length(allCorVals))))
  nasalSubGeneCorMedian<-c(nasalSubGeneCorMedian, median(allCorVals))
  nasalSubGenePMedian<-c(nasalSubGenePMedian, median(allPVals))
  #  print(i) #52 modules
}


nasalSubGeneCorDF<-as.data.frame(nasalSubGeneCorDF)
colnames(nasalSubGeneCorDF)<-c("Cor","Set")
nasalSubGeneCorDF$Cor<-as.numeric(as.character(nasalSubGeneCorDF$Cor))
#Reorder sets
nasalSubGeneCorDF$Set<-factor(as.character(nasalSubGeneCorDF$Set), levels=c(uniNasal[c(1:22,33,44,47:52,23:32,34:43,45:46)]))

#Boxplot of all pairwise correlations per module

coherence_boxplot<-
nasalSubGeneCorDF%>%
  mutate(Set = factor(Set, levels=unique(nasalSubGeneCorDF$Set)[order(as.numeric(as.character(unique(nasalSubGeneCorDF$Set))))]))%>%
  ggplot(aes(y=Cor, x=Set))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  ylab("Correlation")+
  xlab(paste0(module_set, " Modules"))+
  ggtitle(paste0(module_set, " Module Coherence in ", sample_set, " Samples"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), plot.title = element_text(hjust=0.5))

print(coherence_boxplot)

return(list(coherence_boxplot = coherence_boxplot, subgene_correlation_df = nasalSubGeneCorDF))

}
