 # FRAMEWORK FOR WITHIN SAMPLE GENE CORRELATION #


# 1. FOR EACH SAMPLE, FOR EACH MODULE:
#    a) define all pairwise gene combinations
#    b) plot all gene1 x gene2
#    c) calculate correlation

# IS THIS MEANINGFUL?

sample<-"lib28555"
module<-"eos5"
mod.genes<-filter(m1.mod, geneSet==module)%>%dplyr::select(ensemblID)%>%unlist()

  sample_expression<-dat.voom$E[mod.genes ,sample]
  
  genes_pairwise<-
    gtools::combinations(n=length(names(sample_expression)),
                         v = names(sample_expression), r = 2,
                         repeats.allowed = F)%>%
    as.data.frame()%>%
    dplyr::rename(gene1 = V1, gene2 = V2)%>%
    mutate(gene1_ex = sample_expression[gene1], 
           gene2_ex =  sample_expression[gene2])
  
  genes_pairwise_plot<-
    genes_pairwise%>%
    ggplot(aes(x=gene1_ex, y=gene2_ex))+
    geom_point()+
    geom_smooth(method="lm", se=F)
  