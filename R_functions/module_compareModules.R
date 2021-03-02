
#---------------------------
# COMPARE MODULES
#--------------------------

# This script is designed to compare module composition across two sets of gene modules derived via WGCNA.
# The script compares module composition and calculates similarity indexes for the purposes of identifying 
# modules that may have similar composition across studies or across iterations of WGCNA within study.

# Input: A dataframe of modules from multiple studies/runs. 
# These dataframes must contain columns for "ensemblID", "module", "study".
# 
# module_gene_sets = (data.frame) A data frame of module composition including the following:
#                     module_gene_sets$module = (character) Names of modules.
#                     module_gene_sets$ensemblID = (character) Genes defining each module. This could be based on other identifiers as well (eg. HGNC symbol, but the function looks for "ensemblID")
#                     module_gene_sets$study = (character) Name of study or iteration that distinguishes the module sets. Must distinguish module sets in order to generate relevent comparisons.
#
# p.adjustment = (character) Method by which to fdr adjust hypergeometric p values. Defaults to more aggressive "bonferonni" correction.

module_compareModules<-function(module_gene_sets, p.adjustment="bonferroni"){
  
  # Load libraries & tools
  require(gtools)
  require(patchwork)
  require(tidyverse) 
  require(fuzzySim)
  library(spaa)
  `%notin%`<-Negate(`%in%`)
  
  moduleSets<- module_gene_sets%>% # assign modules dataframe 
    unite(module_unique, c(module, study), sep="_", remove = F) # Give each module a unique identifier
  
  # Check that the dataframe is sensible
  if(length(unique(moduleSets$study))<2){stop("Multiple module sets not present. Did you distinguish among studies/runs?")}
  
  # Convert module gene lists to presence absense matrix.
  all_mods_pres_abs<-
    moduleSets%>%
    dplyr::select(module_unique, ensemblID)%>%
    fuzzySim::splist2presabs(sites.col = "module_unique", sp.col = "ensemblID", data = .)%>%
    column_to_rownames("module_unique")%>%
    as.matrix()
  
  # Calculate Sorensen's Similarity
  
  # Create Dataframe to hold Module Comparisons
  sorensen_shared_genes<-matrix(nrow = length(unique(moduleSets$module_unique)), ncol = length(unique(moduleSets$module_unique)), dimnames = list(unique(moduleSets$module_unique), unique(moduleSets$module_unique)))
  
  #Find intersection of module pairs & calculate Sorensen similarity
  for(r in rownames(sorensen_shared_genes)){
    for(c in colnames(sorensen_shared_genes)){
      sorensen_shared_genes[r,c]<-sorensen(filter(moduleSets, module_unique==r)$ensemblID%>%unlist()%>%as.character(), 
                                           filter(moduleSets, module_unique==c)$ensemblID%>%unlist()%>%as.character())
    }
  }
  
  # Convert matrix to pairwise
  sorensen_sim_pairwise<-
    data.frame(expand.grid(rownames(sorensen_shared_genes), rownames(sorensen_shared_genes)), 
               value = stack(as.data.frame(sorensen_shared_genes))$values)%>%
    dplyr::rename(module_unique_1 = Var1, module_unique_2 = Var2)%>%
    left_join(dplyr::select(moduleSets, module_unique, module, study), by=c("module_unique_1"="module_unique"))%>%
    distinct()%>%
    dplyr::rename(module_1 = module, study_1 = study)%>%
    left_join(dplyr::select(moduleSets, module_unique, module, study), by=c("module_unique_2"="module_unique"))%>%
    distinct()%>%
    dplyr::rename(module_2 = module, study_2 = study)
  
  # Reduce to meaningul comparisons (no self comparisons, no comparisons within study)
  sorensen_sim_pairwise_meaningful<-
    sorensen_sim_pairwise%>%
    filter(module_unique_1 != module_unique_2)%>% # Remove self comparisons
    filter(study_1 != study_2) # remove comparisons made within study
  
  
  # remove redundant comparisons (pairwise redundancy)
  delRows = NULL # the rows to be removed
  for(i in 1:nrow(sorensen_sim_pairwise_meaningful)){
    j = which(sorensen_sim_pairwise_meaningful$module_unique_1 == sorensen_sim_pairwise_meaningful$module_unique_2[i] & sorensen_sim_pairwise_meaningful$module_unique_2 == sorensen_sim_pairwise_meaningful$module_unique_1[i])
    j = j [j > i]
    if (length(j) > 0){
      delRows = c(delRows, j)
    }
  }
  
  # remove redundant rows
  sorensen_sim_pairwise_meaningful = sorensen_sim_pairwise_meaningful[-delRows,]%>%
    rename(sorensen_index = value)
  
  # Show distribution of similarities
  sorenson_boxplot<-
  sorensen_sim_pairwise_meaningful%>%
    ggplot(aes(y=sorensen_index))+
    geom_boxplot(outlier.shape = 1, outlier.size = 3)+
    geom_hline(yintercept = 0.05, color = "red", linetype ="dashed")+
    ylab("Sorensen - Dice Similarity")+
    xlab("")+
    ggtitle("Pairwise Module \nSimilarity")+
    labs(subtitle = paste(unique(moduleSets$study), collapse = " : "))+
    theme_bw()+
    theme(aspect.ratio = 4, axis.text.x = element_blank(),
          panel.grid = element_blank(), 
          plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))
  
  
  
  #---------------------------------------
  #   CALCULATE HYPERGEOMETRIC P VALUES
  #---------------------------------------
  
  # Calculate Hypergeometric P Values
  
  total_genes<-length(unique(moduleSets$ensemblID))
  
  
  # Calculate Hypergeometric p value
  
  # All Module Comparisons
  phyper_shared_genes<-matrix(nrow = length(unique(moduleSets$module_unique)), ncol = length(unique(moduleSets$module_unique)), dimnames = list(unique(moduleSets$module_unique), unique(moduleSets$module_unique)))
  
  #Find intersection of module pairs & calculate hypergeometric p
  for(r in rownames(phyper_shared_genes)){
    for(c in colnames(phyper_shared_genes)){
      
      geneSet_r<-filter(moduleSets, module_unique==r)$ensemblID%>%unlist()%>%as.character()
      geneSet_c<-filter(moduleSets, module_unique==c)$ensemblID%>%unlist()%>%as.character()
      geneSet_overlap<-length(intersect(geneSet_r, geneSet_c))
      
      phyper_shared_genes[r,c]<-phyper(geneSet_overlap-1, length(geneSet_c), total_genes-length(geneSet_c), length(geneSet_r), lower.tail= FALSE)
    }
  }
  
  # Convert matrix to pairwise
  phyper_pairwise<-
    data.frame(expand.grid(rownames(phyper_shared_genes), rownames(phyper_shared_genes)), 
               value = stack(as.data.frame(phyper_shared_genes))$values)%>%
    dplyr::rename(module_unique_1 = Var1, module_unique_2 = Var2)%>%
    left_join(dplyr::select(moduleSets, module_unique, module, study), by=c("module_unique_1"="module_unique"))%>%
    distinct()%>%
    dplyr::rename(module_1 = module, study_1 = study)%>%
    left_join(dplyr::select(moduleSets, module_unique, module, study), by=c("module_unique_2"="module_unique"))%>%
    distinct()%>%
    dplyr::rename(module_2 = module, study_2 = study)
  
  
  phyper_pairwise_meaningful<-
    phyper_pairwise%>%
    filter(module_unique_1 != module_unique_2)%>% # Remove self comparisons
    filter(study_1 != study_2) # remove comparisons made within study
  
  
  # remove redundant comparisons (pairwise redundancy)
  delRows = NULL # the rows to be removed
  for(i in 1:nrow(phyper_pairwise_meaningful)){
    j = which(phyper_pairwise_meaningful$module_1 == phyper_pairwise_meaningful$module_2[i] & phyper_pairwise_meaningful$module_2 == phyper_pairwise_meaningful$module_1[i])
    j = j [j > i]
    if (length(j) > 0){
      delRows = c(delRows, j)
    }
  }
  # drop redundant rows & fdr correct
  phyper_pairwise_meaningful = phyper_pairwise_meaningful[-delRows,]%>%
    rename(p_hyper = value)%>%
    mutate(p_hyper_adj=p.adjust(p_hyper, method = p.adjustment))
  
  # Show distribution of similarities
  p_val_boxplot<-
  phyper_pairwise_meaningful%>%
    ggplot(aes(y=p_hyper_adj))+
    geom_boxplot(outlier.shape = 1, outlier.size = 3)+
    geom_hline(yintercept = 0.05, color = "red", linetype ="dashed")+
    ylab("Hypergeometric P Value (FDR)")+
    xlab("")+
    ggtitle("Pairwise Module \nSimilarity P")+
    labs(subtitle = paste(unique(moduleSets$study), collapse = " : "))+
    theme_bw()+
    theme(aspect.ratio = 4, axis.text.x = element_blank(),
          panel.grid = element_blank(), 
          plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))
  
  # Combine with Sorensen similarity
  combined_pairwise_meaningful<-
    sorensen_sim_pairwise_meaningful%>%
    full_join(phyper_pairwise_meaningful)%>%
    dplyr::select(study_1, module_1, module_unique_1, study_2, module_2, module_unique_2, sorensen_index, p_hyper, p_hyper_adj)%>%
    mutate(p_adjust_method = p.adjustment)%>%
    mutate(module_1 = factor(module_1, levels = mixedsort(unique(module_1))),
           module_2 = factor(module_2, levels = mixedsort(unique(module_2))))
  
  
  # Plot values together
  plot_labels<-
    combined_pairwise_meaningful%>%
    mutate(label = ifelse(p_hyper_adj<0.05 & sorensen_index>0.1, as.character(module_1), ""),
           study = study_1)%>%
    bind_rows(mutate(combined_pairwise_meaningful,
                     label = ifelse(p_hyper_adj<0.05 & sorensen_index>0.1, as.character(module_2), ""),
                     study = study_2))%>%
    dplyr::select(module_1, study_1, module_2, study_2, sorensen_index, p_hyper, p_hyper_adj, label, study)%>%
    mutate(module_1_label = ifelse(p_hyper_adj<0.05 & sorensen_index>0.1, as.character(module_1), ""),
           module_2_label =  ifelse(p_hyper_adj<0.05 & sorensen_index>0.1, as.character(module_2), ""))
  
  module_comparison_plot<-
  combined_pairwise_meaningful%>%
    ggplot(aes(y=-log10(p_hyper_adj), x=sorensen_index))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
    geom_vline(xintercept = 0.1, color = "red", linetype = "dashed")+
    ggrepel::geom_label_repel(data = plot_labels,
                             aes(label=label, x =sorensen_index , y = -log10(p_hyper_adj),  fill = study),
                             max.overlaps = Inf,
                             min.segment.length = 0, color = "white", segment.color = "black",
                             fontface = 'bold',
                             box.padding = unit(0.35, "lines"),
                             point.padding = unit(0.5, "lines"))+
    

  
    ylab(" -log10(Adjusted Hypergeometric P Value)")+
    xlab("Sorensen-Dice Similarity Index")+
    labs(fill = "Study")+
    theme_bw()
  
  
  #-------------------
  # CREATE MULTI PLOTS
  #-------------------

# Setup storage for plots 
multi_plot_list<-list()
    
    for(i in unique(moduleSets$study)){ # make plot centric to each study included in the geneSet list 
      
      if (i %in% unique(plot_labels$study_1)){ # Determine whether to draw comparisons from $module_1 or $module_2
        plot_labels_sub<-filter(plot_labels, study_1==i)
        multi_sub<-
          combined_pairwise_meaningful%>%
          filter(study_1==i)%>%
          ggplot(aes(y=-log10(p_hyper_adj), x=sorensen_index))+
          geom_point()+
          geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
          geom_vline(xintercept = 0.1, color = "red", linetype = "dashed")+
          ggrepel::geom_text_repel(data = dplyr::select(plot_labels, -c(study, label))%>%distinct(),
                                   aes(label=module_2_label, x =sorensen_index , y = -log10(p_hyper_adj)),
                                   max.overlaps = Inf,
                                   min.segment.length = 0, segment.color = "black",
                                   fontface = 'bold',
                                   box.padding = unit(0.35, "lines"),
                                   point.padding = unit(0.5, "lines"))+
          
          
          
          ylab(" -log10(Adjusted Hypergeometric P Value)")+
          xlab("Sorensen-Dice Similarity Index")+
          ggtitle(paste("Focus on ", i))+
          labs(fill = "Study")+
          theme_bw()+
          theme(plot.title = element_text(hjust=0.5))+
          facet_wrap(~module_1)
          
        multi_plot_list[[i]]<-multi_sub
        
        } else if (i %notin% plot_labels$study_1 & i %in% plot_labels$study_2){
        plot_labels_sub<-filter(plot_labels, study_2==i)
        multi_sub<-
          combined_pairwise_meaningful%>%
          filter(study_2==i)%>%
          ggplot(aes(y=-log10(p_hyper_adj), x=sorensen_index))+
          geom_point()+
          geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
          geom_vline(xintercept = 0.1, color = "red", linetype = "dashed")+
          ggrepel::geom_text_repel(data = dplyr::select(plot_labels, -c(study, label))%>%distinct(),
                                   aes(label=module_1_label, x =sorensen_index , y = -log10(p_hyper_adj)),
                                   max.overlaps = Inf,
                                   min.segment.length = 0, segment.color = "black",
                                   fontface = 'bold',
                                   box.padding = unit(0.35, "lines"),
                                   point.padding = unit(0.5, "lines"))+
          
          
          
          ylab(" -log10(Adjusted Hypergeometric P Value)")+
          xlab("Sorensen-Dice Similarity Index")+
          ggtitle(paste("Focus on ", i))+
          labs(fill = "Study")+
          theme_bw()+
          theme(plot.title = element_text(hjust=0.5))+
          facet_wrap(~module_2)
        
        multi_plot_list[[i]]<-multi_sub
      }}
  
# Assemble plots into multipanel plots.
module_comparison_plot_multi<-eval(parse(text= paste(paste0("multi_plot_list$", names(multi_plot_list)), collapse = " + ")))
 
  

  return(list(module_comp_df =combined_pairwise_meaningful , module_comp_plot = module_comparison_plot, module_comp_plot_multi=module_comparison_plot_multi))
}
  
  #-------------------------------------------------
  #
  #  CALCULATE SORENSEN-DICE SIMILARITY
  #
  #-------------------------------------------------
  
  # Note: This function, while simple, was lifted from the {OmicsMarkeR} package by Charles E. Determan Jr. 
  
  # Calculates Dice-Sorensen's index between two vectors of 
  # features.  In brief, the closer to 1 the more similar the vectors.  
  # The two vectors may have an arbitrary cardinality (i.e. don't need 
  # same length).  Very similar to the Jaccard Index \code{\link{jaccard}} 
  # but Dice-Sorensen is the harmonic mean of the ratio.

  sorensen <- function(x,y){
    if(!is.character(x)){stop("Both input vectors must be class (character).")}
    if(!is.character(y)){stop("Both input vectors must be class (character).")}
    
    index <- 
      2*(length(intersect(x,y)))/(2*(length(intersect(x,y)))+
                                    length(setdiff(x,y))+
                                    length(setdiff(y,x)))
    return(index)
  }
  