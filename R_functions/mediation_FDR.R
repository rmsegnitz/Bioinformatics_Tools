


mediation_FDR_correct = function(mediation_output, method = "BH"){
  
  mediation.contrasts <-mediation_output$contrasts
  
  # get contrast names in df named format
  contrast_names<-list()
  
  for(i in 1:length(mediation.contrasts)){
    contrast_names[i]<-paste("_T-", mediation.contrasts[[i]][1], "_C-", mediation.contrasts[[i]][2], sep = "")
  }
  
  
  # list dfs for each contrast
  contrst1_df_list<-
    names(mediation.analysis.mod2N$summary)[which(grepl(contrast1, names(mediation.analysis.mod2N$summary)))]
  
  contrst2_df_list<-
    names(mediation.analysis.mod2N$summary)[which(grepl(contrast2, names(mediation.analysis.mod2N$summary)))]
  
  # List genes for each df list (should be indentical)
  contrst1_df_list_genes<-
    unlist(lapply(contrst1_df_list, function(x){str_split(x, "_")[[1]][1]}))
  contrst2_df_list_genes<-
    unlist(lapply(contrst2_df_list, function(x){str_split(x, "_")[[1]][1]}))
  
  # Assign names to genes
  names(contrst1_df_list_genes)<-contrst1_df_list
  names(contrst2_df_list_genes)<-contrst2_df_list
  
  # concatenate lists
  contrst_df_list_genes<-c(contrst1_df_list_genes, contrst2_df_list_genes)
  
  # List contrast for each df list (should be indentical)
  contrst1_df_list_contrast<-rep(contrast1, length(contrst1_df_list))
  contrst2_df_list_contrast<-rep(contrast2, length(contrst2_df_list))
  
  # Assign names to genes
  names(contrst1_df_list_contrast)<-contrst1_df_list
  names(contrst2_df_list_contrast)<-contrst2_df_list
  
  # concatenate lists
  contrast_df_list_contrast<-c(contrst1_df_list_contrast, contrst2_df_list_contrast)
  
  ## convert summaries to wide format
  summaries_transformed.mod2N<-
    lapply(mediation.analysis.mod2N$summary, summary_to_wide)
  
  # Assign gene names into each wide summary
  for(i in names(summaries_transformed.mod2N)){
    summaries_transformed.mod2N[[i]]$gene<-contrst_df_list_genes[[i]]
  }
  
  # Assign contrast names into each wide summary
  for(i in names(summaries_transformed.mod2N)){
    summaries_transformed.mod2N[[i]]$contrast <- contrast_df_list_contrast[[i]]
  }
  
  # Bind summaries
  summaries_transformed.mod2N_all<-
    summaries_transformed.mod2N%>%
    bind_rows()%>%
    dplyr::select(gene, contrast, everything())
  
  
  ## Run FDR correction over p-values
  colnames(summaries_transformed.mod2N_all)[grepl(x = colnames(summaries_transformed.mod2N_all), pattern = "p-value")]
  
  
  summaries_transformed.mod2N_all_FDR<-
    summaries_transformed.mod2N_all%>%
    group_by(contrast)%>%
    mutate(`ACME (average)_p-value_adjusted` = p.adjust(`ACME (average)_p-value`, method = "BH"))%>%
    mutate(`ACME (control)_p-value_adjusted` = p.adjust(`ACME (control)_p-value`, method = "BH"))%>%
    mutate(`ACME (treated)_p-value_adjusted` = p.adjust(`ACME (treated)_p-value`, method = "BH"))%>%
    mutate(`ADE (average)_p-value_adjusted` = p.adjust(`ADE (average)_p-value`, method = "BH"))%>%
    mutate(`ADE (control)_p-value_adjusted` = p.adjust(`ADE (control)_p-value`, method = "BH"))%>%
    mutate(`ADE (treated)_p-value_adjusted` = p.adjust(`ADE (treated)_p-value`, method = "BH"))%>%
    mutate(`Prop. Mediated (average)_p-value_adjusted` = p.adjust(`Prop. Mediated (control)_p-value`,
                                                                  method = "BH"))%>%
    mutate(`Prop. Mediated (control)_p-value_adjusted` = p.adjust(`Prop. Mediated (control)_p-value`, 
                                                                  method = "BH"))%>%
    mutate(`Prop. Mediated (treated)_p-value_adjusted` = p.adjust(`Prop. Mediated (treated)_p-value`, 
                                                                  method = "BH"))%>%
    mutate(`Total Effect_p-value_adjusted` = p.adjust(`Total Effect_p-value`, method = "BH"))
  
  
  
}