# plotting functions for Sleuth results
# Author: Max Segnitz, msegnitz@uw.edu
# Started January 2022
#
# © Richard M Segnitz 2022
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#

#'##################################
##  plot_sleuth_res_geneTx()   ##
#'##################################

#'# DESCRIPTION:
# A semi-generic plotting function for visualizing transcript and gene-level Sleuth results.
# Currently works to visualize RSTR effect in the RSTR/LTBI analysis. Could be better generalized 
# in the future

# plot_sleuth_res_geneTx() takes transcript and gene level results from Sleuth and returns
# a plot showing transcript level fold changes and corresponding gene-level p/FDR values.


# REQUIRED
# tx_results = (data.frame) A data frame of transcript level results returned by Sleuth. 
# gene_results = (data.frame) A data frame of gene-level results returned by Sleuth.
# gene = (character) HGNC formatted name of gene to by visualized.

# OPTIONAL
# palette = (vector, character) A user provided palette for plotting transcripts. Defaults to NULL & auto generates a palette.
# sig = (numeric) P value threshold at which to color transcripts. Defaults to NULL & automatically determines a pval threshold that maintains plotting colors <=6.
# gw_expression = (logical) indicate whether the expression scale should be set  across all results (as opposed to within the gene). default = F

########### DEFINE INPUTS ###############

plot_sleuth_res_geneTx<-
  function(tx_results, gene_results, gene, palette = NULL, sig = NULL, gw_expression=F){

require(ggthemes)

#---------------------------
# Format Data for plotting
#---------------------------

res_tx_sub<-
  tx_results%>%
  filter(ext_gene== gene, effect == "RSTR")%>%
  mutate(stim = ifelse(grepl("_TB", reference_level), "TB", "MEDIA"))%>%
  left_join(
    filter(gtf_df, gene_name==gene)%>%dplyr::select(transcript_id, transcript_name)%>%distinct(),
    by=c("target_id"="transcript_id"))



# Determine optimal p cutoff
if(is.null(sig)){
p_threshold_eval<-list()
for(i in c(0.05, 0.045,0.04, 0.035,  0.03, 0.025, 0.02, 0.015, 0.01, 0.001)){
  p_threshold_eval[[as.character(i)]]<-length(unique(filter(res_tx_sub, pval<i)$transcript_name))
}

sig<-as.numeric(names(p_threshold_eval[which(p_threshold_eval<=6)][1]))
}

if(is.na(sig)){stop("No sensible p threshold determined. Please provide p threshold (likely lower than 0.001)")}

non_sig_label<-paste0("no p<", sig)

res_tx_sub<-
  res_tx_sub%>%
  group_by(transcript_name)%>% # group by Tx
  mutate(tx_color_code = ifelse(min(pval)< sig, transcript_name, non_sig_label))%>%  # code as non-significant if no pval<0.05
  ungroup()%>% # ungroup before ordering color codes
  mutate(tx_color_code = factor(tx_color_code, levels = c(non_sig_label, gtools::mixedsort(unique(.$tx_color_code)[-grep(non_sig_label, unique(.$tx_color_code))], decreasing = T)))) # Order color coding smartly

res_gene_sub<-
  gene_results%>%
  filter(ext_gene== gene, effect == "RSTR")%>%
  mutate(stim = ifelse(grepl("_TB", reference_level), "TB", "MEDIA"))%>%
  left_join(
    filter(gtf_df, gene_name==gene)%>%dplyr::select(transcript_id, transcript_name)%>%distinct(),
    by=c("target_id"="transcript_id"))

#---------------------------
# Define Plotting palete
#---------------------------
if(is.null(palette)){
  if(length(unique(res_tx_sub$tx_color_code))<=7){
    
    plot_palette<-c("grey", colorblind_pal()(8)[-1]) # Generate color-blind friendly palette
    names(plot_palette) <- c(non_sig_label, as.character(unique(res_tx_sub$tx_color_code))[-grep(non_sig_label, unique(res_tx_sub$tx_color_code))]) # assign tx names to palette
    plot_palette<-plot_palette[-which(is.na(names(plot_palette)))] # drop unused colors
    
  } else stop(paste0("The number of transcripts with significant (p<0.05) results exceeds the default plotting palette (6 + grey/non-significant). Please provide a palette with at least ", length(unique(res_tx_sub$tx_color_code)), " colors, or provide a lower pvalue threshold for coloring transcripts."))
  
  
} else if(length(palette)>= length(unique(res_tx_sub$tx_color_code))){
  
  plot_palette=palette # assign palette
  names(plot_palette) <- c(non_sig_label, unique(res_tx_sub$tx_color_code)[-grep(non_sig_label, unique(res_tx_sub$tx_color_code))]) # assign tx names to palette
  plot_palette<-plot_palette[-which(is.na(names(plot_palette)))] # drop unused colors
  
  
  
  } else {stop("The number of transcripts with significant results exceeds the provided palette. Please provide palette with at least ", length(unique(res_tx_sub$tx_color_code)), "colors.")}


#---------------------------
#       Create Plots       
#---------------------------

# Tx logFC
tx_plot<-
  res_tx_sub%>%
  mutate(sign = ifelse(b>0, "pos", ifelse(b<0, "neg", "Z")))%>%
  arrange(tx_color_code)%>%
  
  ggplot(aes(x=-log10(pval), y=b, group=stim, color=tx_color_code))+
  geom_hline(yintercept = 0, linetype="dashed", size=0.25)+
  geom_vline(xintercept = 0, linetype="dashed", size=0.25)+
  geom_errorbar(aes(ymin = b-se_b, ymax = b+se_b), width=0, position=position_dodge(width=0.25))+
  geom_point(aes(size=mean_obs, shape=stim),alpha=0.5 ,position=position_dodge(width=0.25))+
  scale_shape_manual(values = c("MEDIA" = 1, "TB" = 16))+
  scale_color_manual(values = plot_palette)+
  scale_size_area()+
  ylab("Bias Adjusted Fold Change (b)")+
  xlab("-log10(Tx Pval)")+
  labs(shape = " ", color = "Transcript", size = "Mean Normalized \nExpression")+
  guides(shape = guide_legend(order = 3), 
         color = guide_legend(order = 2),
         size = guide_legend(order = 1))+
  ggtitle(gene)+
  theme_bw()+
  theme(aspect.ratio=1, 
        panel.grid = element_blank(), 
        strip.background = element_rect(fill="white"),
        plot.title = element_text(hjust=0.5))+
  facet_wrap(~stim)



geneP_FDR_plot<-
  res_gene_sub%>%
  dplyr::rename(p = pval, FDR = qval)%>%
  pivot_longer(c(p, FDR))%>%
  ggplot(aes(x=-log10(value), y=""))+
  geom_vline(xintercept = 0, linetype="dashed", size=0.25)+
  geom_point(aes(shape=name, fill=stim))+
  scale_shape_manual(values = c("p" = 21, "FDR" = 24))+
  scale_fill_manual(values = c("MEDIA" = "white", "TB" = "black"), guide="none")+
  # xlim(c(0,6))+
  ylab("")+
  xlab("-log10(Aggregate Gene P/FDR)")+
  labs(shape="")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank())

#---------------------------------------------------
# Define Plotting scale for expression if necessary
#---------------------------------------------------
if(gw_expression){
  
  # Determine expression limits & breaks
  exp_range<-
    range(filter(tx_results, ext_gene %in% plotting_genes_of_interest_clean)$mean_obs)
  
  exp_range_plot_breaks<-
    c(floor(floor(exp_range[1]*10)/5)*0.5, 
      floor(floor((exp_range[1] + diff(exp_range)/3)*10)/5)*0.5,
      floor(floor((exp_range[1] + 2*diff(exp_range)/3)*10)/5)*0.5,
      floor(ceiling(exp_range[2]*10)/5)*0.5)
  
  tx_plot<-
    tx_plot+
    scale_size_area(limits = range(exp_range_plot_breaks), breaks=exp_range_plot_breaks)
  
}



# Assemble multipanel

plot_layout <-"
AAAA#
AAAAC
AAAAC
AAAAC
AAAA#
AAAA#
BBBBD
"

tx_gene_results_plot<-
  tx_plot + theme(legend.position = "none")+
  geneP_FDR_plot+ theme(legend.position = "none")+
  get_legend(tx_plot)+
  get_legend(geneP_FDR_plot)+
  plot_layout(design=plot_layout)

return(tx_gene_results_plot)
  }
