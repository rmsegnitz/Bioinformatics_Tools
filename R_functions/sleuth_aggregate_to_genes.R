
#####################################################################
###                                                               ###
###     Aggregate Tx-level Sleuth results table to gene level     ###
###                                                               ###
#####################################################################

# Implement Sleuth's p-value aggreggation method for aggreggating tx level inference to gene level. 
#
# Details: Per Sleuth's method used by sleuth_results(), this function aggregates transcript results to the gene level by
# aggregating p-values using Lancaster weighted v-value aggregation. P values from transcript expression are weighted by the 
# mean tx expression in the tx level comparison.

# Input: transcript level results table returned by sleuth_results(), or a table in the same format.
# Output: gene-level results, equivalent to running sleuth_results(pval_aggregate = TRUE).

# take non aggregated results table and aggregate to gene-level
ag_results_to_gene<-function(tx_results_table, p_weight_func = identity){

  tx_results_table%>%
  group_by(ens_gene, ext_gene)%>%
  summarise(num_aggregated_transcripts = length(!is.na(pval)),
            sum_mean_obs_counts = sum(p_weight_func(mean_obs), na.rm = TRUE),
            pval = as.numeric(aggregation::lancaster(pval, p_weight_func(mean_obs))))%>%
  ungroup()%>%
  mutate(qval = p.adjust(pval, method = "BH"))%>%
  arrange(qval)

}




  