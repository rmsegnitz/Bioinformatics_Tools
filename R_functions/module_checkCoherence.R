
# Testing of Module Coherence Across Datasets 
# Author: Max Segnitz, msegnitz@uw.edu
# Started January 2021
#

# © Richard M Segnitz 2020 (Adapted and expanded from code provided by Matt Altman)
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#


# DESCRIPTION:
# Contains function to test the coherence of modules built from Dataset A in data from Dataset B.
#
# module_checkCoherence() calculates module expression values for the new dataset, generates within-module pairwise gene correlations,
# and p values.
# 
# Outputs:  (subgene_correlation_df) a df of all inter-gene correlation values and p values, 
#           (subgene_correlation_summary) a summary with median correlation values and p values for each module.
#           (coherence_boxplot_cor) a ggplot of correlation values by module.
#           (coherence_boxplot_p) a ggplot of correlation -log10(p) values by module
#           (coherence_boxplot_combined) a combined plot with both correlation and p values. 

#---------------------------
# 
#--------------------------

###############################
#   module_checkCoherence()   #
###############################

# REQUIRED

# module_gene_sets  = (dataframe) dataframe giving module membership. Should contain column "ensemblID" with ensembl IDs.
# voom_obj = (voom object) voom object containign RNASeq dataset in which module coherence is to be tested.

# OPTIONAL
# geneSet (character string) = name of the column in module_gene_sets which defines modules. Defaults to "geneSet"
# module_set (character string) = name the study from which modules were built. This is used simply for labeling of outputs. Default = "STUDY"
# sample_set (character string) = name of the study from which the data come. This is used simply for labeling of outputs. Default = "STUDY"
# remove_sets (vector) = a vector of character strings naming modules which you want removed from the analysis (eg. "0").
# R_cutoff (numeric) = a vector or single numeric value at which to draw red cutoff lines in correlation (R) plot. Default = 0.3
# P_cutoff (numeric) = a vector or single numeric value at which to draw red cutoff lines in significance (P) plot. Default = 0.01

# EXAMPLE USAGE

# StudyA_in_StudyB_coherence<-
#   module_checkCoherence(module_gene_sets = StudyA_mods, 
#                         voom_obj = StudyB_voom, 
#                         module_set = "StudyA", 
#                         sample_set = "StudyB",
#                         R_cutoff = 0.3, P_cutoff = 0.01)

########### DEFINE INPUTS ###############

module_checkCoherence<-function(module_gene_sets, voom_obj, geneSet = "geneSet",module_set="STUDY", sample_set="STUDY", remove_sets = NULL, R_cutoff = 0.3, P_cutoff = 0.01){

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
  geom_hline(yintercept = R_cutoff, linetype = "dashed", color = "red")+
  scale_y_continuous(breaks=c(seq(-1,1,0.2), R_cutoff))+
  ylab("Correlation Strength")+
  xlab(paste0(module_set, " Modules"))+
  labs(title = paste0(module_set, " Module Coherence in ", sample_set, " Samples"),
       subtitle = "Inter-Gene PearsonCorrelation")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5))
  
labels_df_P<-data.frame(name=rep("negLogP", length(P_cutoff)),
                        h_pos = -log10(P_cutoff)-0.2,
                        h_label = paste("p=",P_cutoff,sep=""))
  
coherence_boxplot_p<-
  SubGeneCorDF%>%
  ggplot(aes(y=-log10(P), x=Set))+
  geom_hline(yintercept =0, color="black")+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = -log10(P_cutoff), color="red", linetype="dashed")+
  geom_text(data=labels_df_P, aes(label = h_label, x=3, y=h_pos), color="red")+
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
labels_df<-data.frame(name=c(rep("Cor", length(R_cutoff)), rep("negLogP", length(P_cutoff))),
                        h_line = c(R_cutoff, -log10(P_cutoff)),
                        h_pos = c(R_cutoff-0.02, -log10(P_cutoff)-0.2),
                        h_label = c(paste("r=",R_cutoff,sep=""), paste("p=",P_cutoff,sep="")))

coherence_boxplot_faceted<-
  SubGeneCorDF%>%
  mutate(negLogP = -log10(P))%>%
  pivot_longer(cols = c(Cor, negLogP))%>%
  ggplot(aes(y=value, x=Set))+
  geom_hline(yintercept =0, color="black")+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(data=labels_df, aes(yintercept = h_line), color = "red", linetype="dashed")+
  geom_text(data=labels_df, aes(label = h_label, x=3, y=h_pos), color="red")+
  xlab(paste0(module_set, " Modules"))+
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
            subgene_correlation_df = SubGeneCorDF,
            subgene_correlation_summary = SubGeneCorDF_summary))

}
