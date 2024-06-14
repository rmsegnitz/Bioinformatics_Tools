# Partitioned Enrichment via the ptEnrich() Function

# Author: Max Segnitz, msegnitz@uw.edu

# December October 2023

# © R MaxSegnitz 2023
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).

# Description:
# The purpose of this function is to run enrichment analysis on subdivisions of a pre-defined gene module or geneset.
# Possible use cases of this analysis include understanding how subsets of a module may be behaving with respect to
# changes in the internal correlation structure. 

# The function takes as input a define parent module/geneset and 2 subset gene lists. These gene lists can be generated
# in different ways, but was developed using quartiles defined by changes in within-module correlation. Ideally these
# gene lists are of different sizes, though this is not necessary. 

# First hypergeometric enrichment of each subset is run, and subsequently random permutations are run in which randomly
# partitioned subsets are analyzed. Permutation-derived significance is determined, reflecting how frequently we find 
# differences in enrichment as large or large as what is observed. Optionally, a baseline enrichment of the full module 
# is also run and returned for comparison.

# Dependencies include {fenr} for preparation of the annotations as well as running fast implimentation of the hypergeometric
# enrichment calculations. A slight modification of {fenr}'s core enrichment function is defined below and used in ptEnrich()


###########################################################
#                                                         #
#     1) Prepare reference annotations for analysis.      #
#                                                         #
###########################################################

# This runs some preparatory formatting of reference genesets. 
# Gere sets can be supplied either as:
# 1) a list of named vectors containing hugo gene symbols.
# 2) a dataframe of geneset assignments containing at least 2 columns to identify the pathways names and gene IDs.

# Note that a list of 2 dataframes is returned. After this step, one can modify the term names to be more usable if desired. 


ptPrep<-function(pathway_list=NULL, pathway_df=NULL, 
                 term_id="term_name", gene_id="gene_symbol"){
  require(dplyr)
  if(!is.null(pathway_list)){
    
    suppressWarnings(
      pt_prepped<-list(
      terms=data.frame(term_id=names(pathway_list))%>%
        dplyr::mutate(term_name=term_id),
      
      mapping=(as.data.frame(do.call(rbind, pathway_list))%>%
                 tibble::rownames_to_column("term_id")%>%
                 tidyr::pivot_longer(contains("V"), names_to = NULL, values_to = "gene_symbol")%>%
                 dplyr::distinct())))
    
    return(pt_prepped)
    
  } else if(!is.null(pathway_df)){
    
    suppressWarnings(
      pt_prepped<-list(
      terms=data.frame(term_id=unique(pathway_df%>%dplyr::pull(term_id)))%>%
        dplyr::mutate(term_name=term_id),
      
      mapping=pathway_df%>%dplyr::rename("gene_symbol"=gene_id, "term_id"=term_id)
    ))
    return(pt_prepped)
    
  } #else {print("Error; reference pathways must be provided as either a list or dataframe; see documentation.")}
  
}
  
  
###########################################################
#                                                         #
#     2) Run Partitioned Enrichment                       #
#                                                         #
###########################################################
  

# Running the main enrichment comparison involves several steps. 

#  1) Prepare the reference genesets for use by fenr.
#  2) Run enrichment on the full parent module if desired. 
#  3) Run enrichment on the provided partitions
#  4) Run random permutations to test significance of any differences.
 
ptEnrich<-function(parent_geneset, geneset_sub1, geneset_sub2, background, pt_prepped,
                   enrich_parent=F, permutations=NULL, ncores=NULL, random_seed=32){
  require(dplyr)
  require(fenr)
  require(doParallel)
  require(foreach)
  require(doRNG)
  #require(doSNOW)
  require(assertthat)
  
  # Define partitions to include in results
  partitions<-list(partition1=geneset_sub1, partition2=geneset_sub2)
  
  # Define slight modification of the core hypergeomtric enrichment function from "fenr" 
  # such that pathways with no representation in the queried gene list are still included in 
  # results output & FDR calculation. 
  
  functional_enrichment_mod<-
    function (feat_all, feat_sel, term_data, feat2name = NULL) 
    {
      require(assertthat)
      
      N_with <- n_with_sel <- n_expect <- enrichment <- odds_ratio <- NULL
      desc <- p_value <- p_adjust <- NULL
      assertthat::assert_that(is(term_data, "fenr_terms"))
      if (!any(feat_sel %in% feat_all)) 
        return(NULL)
      
      #our_terms <- unique(unlist(purrr::map(feat_sel, ~term_data$feature2term[[.x]])))
      our_terms <- unique(names(term_data$term2name))
      
      N_sel <- length(feat_sel)
      N_tot <- length(feat_all)
      res <- purrr::map_dfr(our_terms, function(term_id) {
        tfeats <- term_data$term2feature[[term_id]]
        tfeats <- tfeats[tfeats %in% feat_all]
        tfeats_sel <- tfeats[tfeats %in% feat_sel]
        N_with <- length(tfeats)
        N_without <- N_tot - N_with
        n_with_sel <- length(tfeats_sel)
        n_without_sel <- N_sel - n_with_sel
        n_with_nsel <- N_with - n_with_sel
        n_without_nsel <- N_tot - (n_with_sel + n_without_sel + 
                                     n_with_nsel)
        # if (n_with_sel < 2) 
        #   return(NULL)
        n_expect <- N_with * N_sel/N_tot
        odds_ratio <- (n_with_sel/n_without_sel)/(n_with_nsel/n_without_nsel)
        p <- 1 - stats::phyper(n_with_sel - 1, N_with, N_without, 
                               N_sel)
        if (!is.null(feat2name)) 
          tfeats_sel <- unname(feat2name[as.character(tfeats_sel)])
        term_name <- term_data$term2name[[term_id]]
        if (is.null(term_name)) 
          term_name <- NA_character_
        c(term_id = term_id, term_name = term_name, 
          term_in_background = N_with, 
          query_in_background=N_sel,
          query_in_term = n_with_sel,
          n_expect = n_expect, 
          enrichment = n_with_sel/n_expect, 
          odds_ratio = odds_ratio, ids = paste(tfeats_sel, 
                                               collapse = ", "), p_value = p)
      })
      if (nrow(res) == 0) {
        res <- NULL
      }
      else {
        res <- dplyr::arrange(
          dplyr::mutate(res, 
                        dplyr::across(c(term_in_background, query_in_term), as.integer), 
                        dplyr::across(c(n_expect, enrichment, odds_ratio, p_value), as.numeric), 
                        p_adjust = stats::p.adjust(p_value,method = "BH"), 
                        dplyr::across(c(enrichment, odds_ratio, p_value, p_adjust), ~signif(.x, 3)), 
                        n_expect = round(n_expect,2)), desc(odds_ratio))
      }
      return(res)
    }
  
  # Prep reference pathways for fast enrichment
  pt_prepped_f <-
    fenr::prepare_for_enrichment(terms = pt_prepped$terms, mapping = pt_prepped$mapping, all_features = background,
                                                                       feature_name = "gene_symbol")
  # Run enrichment of parent module if desired
  if(enrich_parent){
    parent_enrichment=
      functional_enrichment_mod(feat_all = background, 
                                feat_sel = parent_geneset,
                                term_data =  pt_prepped_f)%>%
      arrange(p_adjust)
  }
  
  # Run enrichment on gene subsets
  subset_sizes<-c(length(geneset_sub1), length(geneset_sub2))
  
  sub1_enrichment<-functional_enrichment_mod(feat_all = background, 
                                             feat_sel = geneset_sub1,
                                             term_data =  pt_prepped_f)%>%
    dplyr::mutate(partition="pt1", 
                  partition_size=subset_sizes[1],
                  pct_pres_sel=query_in_term/partition_size)
  
  sub2_enrichment<-functional_enrichment_mod(feat_all = background, 
                                             feat_sel = geneset_sub2,
                                             term_data =  pt_prepped_f)%>%
    dplyr::mutate(partition="pt2", 
                  partition_size=subset_sizes[2], 
                  pct_pres_sel=query_in_term/partition_size)
  
  # Compile Comparison of enrichments
  pt1_pt2_compdf<-
    dplyr::bind_rows(sub1_enrichment, sub2_enrichment)%>%
    tidyr::pivot_wider(id_cols = c(term_id, term_name, term_in_background),
                names_from=partition,
                values_from = c(partition_size, query_in_term, n_expect, pct_pres_sel, 
                                enrichment, odds_ratio, ids, p_value, p_adjust))%>%
    dplyr::mutate(partition_dPct=abs(pct_pres_sel_pt1-pct_pres_sel_pt2),
                  partition_dOdds=abs(odds_ratio_pt1-odds_ratio_pt2),
                  partition_dEnrichment=abs(enrichment_pt1-enrichment_pt2),
                  partition_dpval=abs(-log10(p_value_pt1)+log10(p_value_pt2)))
  
  # pt_comp_deltaPlot1<-
  #   pt1_pt2_compdf%>%
  #   arrange(term_in_background)%>%
  #   ggplot(aes(y=partition_dpval, x=partition_dEnrichment))+
  #   geom_jitter(width = 0.002, aes(color=term_in_background), alpha=0.5)+
  #   scale_color_viridis_c()+
  #   labs(x="Diff. in Enrichment", 
  #        y="Diff in -log10(enrichment pval)", 
  #        color="Pathway Genes Present\nIn Background")+
  #   theme_bw()
  
  
  # Return if not running permutations
  terms_not_in_background<- filter(pt1_pt2_compdf, term_in_background==0)$term_id
  
  
  
  if(is.null(permutations)){
    if(enrich_parent){
      parent_enrichment<-filter(parent_enrichment, term_in_background>0)
      return(list(parent_enrichment=arrange(parent_enrichment, p_value), 
                  partition_enrichment=arrange(pt1_pt2_compdf, delta_enrichment_empirical_pval), 
                  partitions=partitions,
                  terms_not_in_background=terms_not_in_background))
    } else {
      return(list(partition_enrichment=arrange(pt1_pt2_compdf, delta_enrichment_empirical_pval), 
                  partitions=partitions,
                  terms_not_in_background=terms_not_in_background))
    }
  } else {
    
    ########################
    #   RUN PERMUTATIONS
    ######################## 
    
    # Create storage matrices
    perm_mat_rownames<-unique(pt_prepped$terms$term_id)
    # permutation_dP<-matrix(nrow = length(unique(pt_prepped$terms$term_id)))
    # rownames(permutation_dP)<-perm_mat_rownames
    # permutation_dPct<- permutation_dOdds<-permutation_dP

    # Establish parallel computing specifications
    if(is.null(ncores)){ncores=detectCores()-4}
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    # Set up progress bar
    pb <- txtProgressBar(max = permutations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Run calculations
    permutation_matrix<-
      foreach(i = 1:permutations, .options.snow = opts, .options.RNG = random_seed)%dorng%{
        require(dplyr)
        `%notin%`<-Negate(`%in%`)
        print(paste0("Running ", i, " of ", permutations, "on ", ncores,  " cores."))
        # randomly sample gene subsets
        genes_perm_pt1<-sample(parent_geneset, 
                               subset_sizes[1], replace = F)
        
        genes_perm_pt2<-sample(parent_geneset[which(parent_geneset %notin% c(genes_perm_pt1))], 
                               subset_sizes[2], replace = F)
        
        # Run enrichment on permuted partitions
        
        enrich_pt1_perm_i<-functional_enrichment_mod(feat_all = background, 
                                  feat_sel = genes_perm_pt1,
                                  term_data =  pt_prepped_f)%>%
          dplyr::mutate(partition="pt1", 
                        partition_size=subset_sizes[1], 
                        pct_pres_sel=query_in_term/partition_size)
        
        
        enrich_pt2_perm_i<-functional_enrichment_mod(feat_all = background, 
                                                     feat_sel = genes_perm_pt2,
                                                     term_data =  pt_prepped_f)%>%
          dplyr::mutate(partition="pt2", 
                        partition_size=subset_sizes[2], 
                        pct_pres_sel=query_in_term/partition_size)
        
        
        pt1_pt2_compdf_i<-
          dplyr::bind_rows(enrich_pt1_perm_i, enrich_pt2_perm_i)%>%
          tidyr::pivot_wider(id_cols = c(term_id, term_name, term_in_background),
                      names_from=partition,
                      values_from = c(partition_size, query_in_term, n_expect, pct_pres_sel, 
                                      enrichment, odds_ratio, ids, p_value, p_adjust))%>%
          dplyr::mutate(partition_dPct=abs(pct_pres_sel_pt1-pct_pres_sel_pt2),
                        partition_dOdds=abs(odds_ratio_pt1-odds_ratio_pt2),
                        partition_dEnrichment=abs(enrichment_pt1-enrichment_pt2))%>%
          tibble::column_to_rownames("term_id")
        
        
        # Save permutations
        
        permutation_output<-data.frame(perm_dPct=pt1_pt2_compdf_i[perm_mat_rownames, "partition_dPct"], 
                           perm_dOdds=pt1_pt2_compdf_i[perm_mat_rownames, "partition_dOdds"], 
                           perm_dEnrichment=pt1_pt2_compdf_i[perm_mat_rownames, "partition_dEnrichment"], 
                           row.names = perm_mat_rownames)
        
        # return permuted calculations
        return(permutation_output)
        
      }
    
    
    names(permutation_matrix)<-paste("perm", 1:permutations)# Assign names to permutation list
    parallel::stopCluster(cl)
    
    
    # Assemble permutation matrices
    
    dEnrichment_mat<-
      dplyr::bind_cols(lapply(permutation_matrix, FUN = function(x){x[,"perm_dEnrichment"]}))%>%
      as.matrix()
    rownames(dEnrichment_mat)<-rownames(permutation_matrix[[1]])
    
    dPct_mat<-
      dplyr::bind_cols(lapply(permutation_matrix, FUN = function(x){x[,"perm_dPct"]}))%>%
      as.matrix()
    rownames(dPct_mat)<-rownames(permutation_matrix[[1]])
    
    dOdds_mat<-
      dplyr::bind_cols(lapply(permutation_matrix, FUN = function(x){x[,"perm_dOdds"]}))%>%
      as.matrix()
    rownames(dOdds_mat)<-rownames(permutation_matrix[[1]])
    
    
    
    # Calculate empirical pvals from permutations

    dEnrichment_ref<-
      tibble::column_to_rownames(pt1_pt2_compdf, "term_id")[perm_mat_rownames,"partition_dEnrichment"]
    
    dEnrichment_Empirical<-
      (dEnrichment_mat>=dEnrichment_ref)%>% # determine if permutation is >= observation
      apply(., MARG=1, FUN=function(x){sum(x)/length(x)})
    
    dPct_ref<-
      tibble::column_to_rownames(pt1_pt2_compdf, "term_id")[perm_mat_rownames,"partition_dPct"]
    
    dPct_Empirical<-
      (dPct_mat>=dPct_ref)%>% # determine if permutation is >= observation
      apply(., MARG=1, FUN=function(x){sum(x)/length(x)})
    
    dOdds_ref<-
      tibble::column_to_rownames(pt1_pt2_compdf, "term_id")[perm_mat_rownames,"partition_dOdds"]
    
    dOdds_Empirical<-
      (dOdds_mat>=dOdds_ref)%>% # determine if permutation is >= observation
      apply(., MARG=1, FUN=function(x){sum(x)/length(x)})
    
    
    # combine results
    if(all(names(dEnrichment_Empirical)==names(dPct_Empirical))){
      permutation_tests<-
        data.frame(term_id=names(dEnrichment_Empirical), 
                   delta_enrichment_observed=dEnrichment_ref,
                   delta_enrichment_empirical_pval=dEnrichment_Empirical, 
                   delta_pct_observed=dPct_ref,
                   delta_pct_empirical_pval=dPct_Empirical,
                   delta_oddsRatio_observed=dOdds_ref,
                   delta_oddsRatio_empirical_pval=dOdds_Empirical,
                   nperm=permutations)}
    
    
  }
  
  pt1_pt2_compdf_wPerm<-
    dplyr::left_join(pt1_pt2_compdf, permutation_tests)
 
  #Final filtering to remove unrepresented pathways
  pt1_pt2_compdf_wPerm<-filter(pt1_pt2_compdf_wPerm, term_in_background>0)
  
  return(list(partition_enrichment=arrange(pt1_pt2_compdf_wPerm, delta_enrichment_empirical_pval), 
              parent_enrichment=arrange(parent_enrichment, p_value), 
              partitions=partitions,
              terms_not_in_background=terms_not_in_background))
  
}



###########################################################
#                                                         #
#     2) Plot Partitioned Enrichment                      #
#                                                         #
###########################################################

# Combine partitioned enrichment with parent enrichment results and plot comparisons. 
# Most likely needs manual curation but serves to provide overview and returns dataframe for plotting further.

# ptEntich_obj = object returned by ptEnrich(); sig_cutoff = significance label p value cutoff; term_cleanup = string to remove from term_id
plot_ptEnrich=function(ptEntich_obj, sig_cutoff_p=0.01, term_cleanup=NULL, plot_label=NULL, limit_term_labs=NULL, numbered.key=F){
  
  require(ggtext)
  
  if(!all(c("partition_enrichment", "parent_enrichment") %in% names(ptEntich_obj))){
    stop("Error: Must contain both partitioned enrichment & parent enrichment. Re-run ptEnrich with enrich_parent=T.")
  }
  
  # Format dataframe for comparison
  enrichment_comparison<-
    ptEntich_obj$partition_enrichment%>%
    
    left_join(ptEntich_obj$parent_enrichment%>%dplyr::select(term_name, enrichment, p_value, p_adjust))%>%
    mutate(delta_enrichment_empirical_fdr=p.adjust(delta_enrichment_empirical_pval, method="BH"))%>%
    mutate(partition_result_category=ifelse(delta_enrichment_empirical_pval<sig_cutoff_p, 
                                            case_when(p_value_pt1<sig_cutoff_p & p_value_pt2>=sig_cutoff_p ~ "partition 1 only",
                                                      p_value_pt2<sig_cutoff_p & p_value_pt1>=sig_cutoff_p ~ "partition 2 only",
                                                      p_value_pt2<sig_cutoff_p & p_value_pt1<sig_cutoff_p ~ "partition 1 & 2", 
                                                      .default = "ns"), "ns"))%>%
    mutate(partition_trend=ifelse(delta_enrichment_empirical_pval<sig_cutoff_p, 
                                  case_when(enrichment_pt1>enrichment_pt2 ~"partition1",
                                            enrichment_pt1<enrichment_pt2 ~"partition2"), "ns"))%>%
    rowwise()%>%
    mutate(label_all=ifelse(!is.null(!!term_cleanup),
             str_to_sentence(str_replace_all(str_remove_all(term_name, !!term_cleanup), "_", " ")),
             str_to_sentence(str_replace_all(term_name, "_", " "))))%>%
    mutate(label_sig=ifelse(delta_enrichment_empirical_pval<sig_cutoff_p, label_all, NA))%>%
    ungroup()%>%
    arrange(delta_enrichment_empirical_pval, -delta_enrichment_observed)%>%
    mutate(label_sig_numeric=row_number())%>%
    mutate(label_sig_numeric=ifelse(is.na(label_sig), NA, label_sig_numeric))
  
  # Assemble plots to highlight significantly different enrichments
  nperm=unique(enrichment_comparison$nperm)
  
  plot1_labs<-filter(enrichment_comparison, !is.na(label_sig))
  
  if(!is.null(limit_term_labs)){plot1_labs<-arrange(plot1_labs, delta_enrichment_empirical_pval)%>%head(limit_term_labs)}
  if(is.null(limit_term_labs) & dim(plot1_labs)[1]>20){
    plot1_labs<-plot1_labs%>%mutate(label_sig_og=label_sig, label_sig=label_sig_numeric)} 
  
  plot1<-
    enrichment_comparison%>%
    ggplot(aes(y=-log10(delta_enrichment_empirical_pval+1/(nperm+100)), x=-log10(p_value+1/(nperm+100))))+
    geom_point(aes(color=partition_result_category, size=delta_enrichment_observed, 
                   shape=delta_enrichment_empirical_pval<sig_cutoff_p))+
    geom_vline(xintercept = -log10(0.05),  linetype="dashed")+
    geom_vline(xintercept = -log10(0.01),  linetype="dotted")+
    geom_hline(yintercept = -log10(0.05), linetype="dashed")+
    geom_hline(yintercept = -log10(0.01), linetype="dotted")+
    geom_label(x=0.1, y=-log10(0.05), label="p=0.05")+
    geom_label(x=0.1, y=-log10(0.01), label="p=0.01")+
    ggrepel::geom_text_repel(data=plot1_labs, 
                              aes(label=label_sig, color=partition_result_category), max.overlaps = Inf)+
    scale_color_manual(values=c("partition 1 only"="dodgerblue",  
                                "partition 2 only"="red3", 
                                "partition 1 & 2"="purple3"))+
    scale_shape_manual(values=c("TRUE"=16, "FALSE"=1))+
    xlim(c(0,NA))+
    ylim(c(0,NA))+
    theme_bw()+theme(panel.grid = element_blank())+
    labs(y="Partitioned Enrichment -log10(pval)", 
         x="Full Module Enrichment \n-log10(pval)", 
         shape=paste0("Partitioned Enrichment \nPermuatation p<", sig_cutoff_p),
         color=paste0("Individual Partition \nEnrichment p<", sig_cutoff_p),
         size="Partitions Difference \nin Enrichment Score")
  
  
  plot2_df<-
    enrichment_comparison%>%
    filter((p_value)<sig_cutoff_p | p_value_pt1<sig_cutoff_p | p_value_pt2<sig_cutoff_p | delta_enrichment_empirical_pval<sig_cutoff_p)%>%
    dplyr::select(term_name, label_all, label_sig, label_sig_numeric, term_in_background,partition_result_category, p_value, p_value_pt1, p_value_pt2, 
                  enrichment_pt1, enrichment_pt2, delta_enrichment_empirical_pval)%>%
    pivot_longer(c(p_value_pt1, p_value_pt2, enrichment_pt1, enrichment_pt2), 
                   names_to = "metric", values_to = "value")%>%
    separate(metric, into = c("metric", "partition"), sep = "_pt")%>%
    mutate(metric=case_when(metric== "p_value"~"partition_p", metric=="enrichment" ~ "partition_enrichment"))%>%
    pivot_wider(names_from = "metric", values_from = "value")
 
  plot2_df_labs<-
    filter(plot2_df, !is.na(label_sig))%>%
    group_by(label_sig,label_sig_numeric, partition_result_category, delta_enrichment_empirical_pval)%>%summarize(p_value=unique(p_value), min_pt_p=min(partition_p))
  if(is.null(limit_term_labs) & dim(plot2_df_labs)[1]>20){
    plot2_df_labs<-plot2_df_labs%>%mutate(label_sig_og=label_sig, label_sig=label_sig_numeric)} 
  if(!is.null(limit_term_labs)){plot2_df_labs<-arrange(plot2_df_labs, delta_enrichment_empirical_pval)%>%head(limit_term_labs)}
   
  
  plot2_df_segs<-
    plot2_df%>%
    dplyr::select(term_name, p_value, partition, partition_p, delta_enrichment_empirical_pval)%>%
    pivot_wider(values_from = partition_p, names_from = partition, names_prefix = "p_pt")
    
  plot2<-
    plot2_df%>%
    ggplot(aes(x=-log10(p_value+(1/(nperm+100)))))+
    geom_abline(intercept = 0, slope=1, linetype="dashed")+
    
    geom_segment(data=plot2_df_segs,
                aes(x=-log10(p_value+(1/(nperm+100))), xend=-log10(p_value+(1/(nperm+100))), 
                     y=-log10(p_pt1+(1/(nperm+100))), yend=-log10(p_pt2+(1/(nperm+100)))), color="lightgrey")+
    
    geom_segment(data=filter(plot2_df_segs, delta_enrichment_empirical_pval<sig_cutoff_p) ,
                 aes(x=-log10(p_value+(1/(nperm+100))), xend=-log10(p_value+(1/(nperm+100))), 
                     y=-log10(p_pt1+(1/(nperm+100))), yend=-log10(p_pt2+(1/(nperm+100)))), color="black")+
    
    geom_point(aes(y=-log10(partition_p+1/(nperm+100)), size=partition_enrichment), color="white")+

    geom_point(aes(y=-log10(partition_p+1/(nperm+100)), size=partition_enrichment, 
                   shape=delta_enrichment_empirical_pval<sig_cutoff_p, color=partition))+
    
    ggrepel::geom_text_repel(data= plot2_df_labs,
                               aes(label=label_sig, x=-log10(p_value+(1/(nperm+100))), y=-log10(min_pt_p+(1/(nperm+100))),
                                   color=partition_result_category),
                               max.overlaps = Inf)+
    ggtext::geom_richtext(label="x=y", x=0.25, y=0.25, angle=45)+
    scale_color_manual(values=c("1"="dodgerblue", "2"="red3", "ns"="grey",
                                "partition 1 only"="dodgerblue", "partition 2 only"="red3", "partition 1 & 2"="purple3"),
                       breaks = c("1", "2"))+
    scale_shape_manual(values=c("TRUE"=16, "FALSE"=1))+
    xlim(c(-0.1,NA))+
    ylim(c(-0.1,NA))+
    labs(y="Individual Partition \nEnrichment -log10(pval)", 
         x="Full Module Enrichment \n-log10(pval)", 
         shape=paste0("Partitioned Enrichment \nPermuatation p<", sig_cutoff_p),
         color="Partition",
         size="IndividualPartition Enrichment Score")+
    theme_bw()+
    theme(panel.grid = element_blank())+
    coord_fixed()
  
  if((is.null(limit_term_labs) & dim(plot2_df_labs)[2]>20) | numbered.key){
    plot2_df_labs2<-
      plot2_df_labs%>%
      mutate(facet_group=ceiling(label_sig/50))%>%
      group_by(facet_group)%>%
      mutate(group_size=length(facet_group))%>%
      mutate(facet_group=ifelse(group_size<15, facet_group-1, facet_group))%>%
      ungroup()%>%
      arrange(-label_sig)%>%
      mutate(label_sig=factor(label_sig, levels=.$label_sig))
    
    plot_key<-
      plot2_df_labs2%>%
      ggplot(aes(y=label_sig, x=0))+
      geom_text(aes(label=label_sig_og, color=partition_result_category), hjust=0)+
      scale_color_manual(values=c("1"="dodgerblue", "2"="red3", "ns"="grey40",
                                  "partition 1 only"="dodgerblue", "partition 2 only"="red3", "partition 1 & 2"="purple3"),
                         breaks = c("1", "2"), guide="none")+
      scale_x_continuous(limits=c(0,20))+
      labs(y="", x="")+
      facet_wrap(~facet_group, nrow=1, scales = "free")+
      theme_bw()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x =  element_blank(),
            strip.text = element_blank(), 
            strip.background = element_rect(fill="white"))
    
  }
  
  
  if(!is.null(plot_label)){
    plot1<-plot1+ggtitle(plot_label)
    plot2<-plot2+ggtitle(plot_label)
  }
  
  if(exists("plot_key", inherits = FALSE)){
    return(list(enrichment_comparison=arrange(enrichment_comparison, delta_enrichment_empirical_pval), plot_overview=plot1, plot_comparePartitions=plot2, plot_key=plot_key))
  }else{
  return(list(enrichment_comparison=arrange(enrichment_comparison, delta_enrichment_empirical_pval), plot_overview=plot1, plot_comparePartitions=plot2))}
  
}
