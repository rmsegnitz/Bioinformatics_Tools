

module_calcModExpression<-function(module_gene_sets, voom_obj, gene_set = "geneSet", names = NULL){
  
  # Load libraries
  require(gtools)
  require(tidyverse) 
  `%notin%`<-Negate(`%in%`)
  
  moduleSets<- module_gene_sets
  voomData<-voom_obj
  
  #------------
  # Formatting
  #------------
  moduleSets<-
    moduleSets%>%
    mutate(geneSet = get(gene_set))%>%
    mutate(geneSet = as.character(geneSet))
  
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
  
  if(!is.null(names)){
    GSsub<-GSsub%>%
      as.data.frame()%>%
      rownames_to_column("geneSet")%>%
      mutate(geneSet = as.character(geneSet))%>%
      left_join(dplyr::select(moduleSets, geneSet, names), by = geneSet)%>%
      distinct()%>%
      dplyr::select(-geneSet)%>%
      column_to_rownames(names)%>%
      as.matrix()}
  
  #-----------------
  # Return output
  #-----------------
  return(GSsub)
  
}
  