

### WARNING: This script is in development and should not be used yet. 

sleuth_object<-so_int

samples<-sleuth_object$sample_to_covariates%>%filter(TB=="MEDIA")%>%dplyr::select(sample)%>%unlist()
samples<-NULL

targets<-SO_basic_PC$filter_df$target_id

subset_sleuth_object=function(sleuth_object, samples = NULL, targets = NULL){
  
  # Subset samples
  if(!is.null(samples)){
    
  
    samples_keep<-which(sleuth_object$kal[[1]]$abundance$sample %in% samples)
    
    so_sample_names<-vector(length = length(sleuth_object$kal))
    for(i in 1:length(sleuth_object$kal)){so_sample_names[i]<-sleuth_object$kal[[i]]$abundance$sample[1]}
    
    # Subset kalisto samples in kal object
    sleuth_object$kal<-sleuth_object$kal[c(which(so_sample_names %in% samples))]
    
    # Subset other internals
    sleuth_object$est_counts_sf<- sleuth_object$est_counts_sf[samples]
    
    sleuth_object$obs_raw<-filter(sleuth_object$obs_raw, sample %in% samples)
    sleuth_object$sample_to_covariates<-filter( sleuth_object$sample_to_covariates, sample %in% samples)
    sleuth_object$obs_norm<-filter(sleuth_object$obs_norm, sample %in% samples)
    sleuth_object$obs_norm_filt<-filter(sleuth_object$obs_norm_filt, sample %in% samples)
    
    # subset bootstrap quantification
    sleuth_object$bs_quants<- sleuth_object$bs_quants[samples]
    
    # subset bootstrap summary
    sleuth_object$bs_summary$obs_counts<-sleuth_object$bs_summary$obs_counts[, samples]
    sleuth_object$bs_summary$obs_tpm<-sleuth_object$bs_summary$obs_tpm[,samples]
    
    
    
  }
  if(!is.null(targets)){

    all_targets<-unique(sleuth_object$obs_norm$target_id)
    
    targets_keep<-which(sleuth_object$kal[[1]]$abundance$target_id %in% targets)
    
    for(i in length(sleuth_object$kal)){
      sleuth_object$kal[[i]]$abundance<-sleuth_object$kal[[i]]$abundance[targets_keep,]
      sleuth_object$kal[[i]]$bias_normalized<-sleuth_object$kal[[i]]$bias_normalized[targets_keep]
      sleuth_object$kal[[i]]$bias_observed<-sleuth_object$kal[[i]]$bias_observed[targets_keep]
    }
    
    # Subset other internals
    sleuth_object$obs_raw<-filter(sleuth_object$obs_raw, target_id %in% targets)
    sleuth_object$obs_norm<-filter(sleuth_object$obs_norm, target_id %in% targets)
    sleuth_object$obs_norm_filt<-filter(sleuth_object$obs_norm_filt, target_id %in% targets)
    
    # Set new filtered targets
    sleuth_object$filter_df<-data.frame(target_id = targets)
    sleuth_object$filter_bool<- all_targets %in% sleuth_object$filter_df$target_id
    names(sleuth_object$filter_bool)<-all_targets
    
    # subset bootstrap quantification
    for(i in 1:length(sleuth_object$bs_quants)){
      sleuth_object$bs_quants[[i]]$est_counts<-sleuth_object$bs_quants[[i]]$est_counts[targets,]
    }
    
    # subset bootstrap summary objects
    sleuth_object$bs_summary$obs_counts<-sleuth_object$bs_summary$obs_counts[targets,]
    sleuth_object$bs_summary$obs_tpm<-sleuth_object$bs_summary$obs_tpm[targets, ]
    sleuth_object$bs_summary$sigma_q_sq<-sleuth_object$bs_summary$sigma_q_sq[targets]
    sleuth_object$bs_summary$sigma_q_sq_tpm<-sleuth_object$bs_summary$sigma_q_sq_tpm[targets]
    
    
  }
  
  return(sleuth_object)
  
}





# Test subsetting approach

sleuth_object_test<-so_int

sleuth_object_test_filt<-subset_sleuth_object(sleuth_object_test, targets = targets)

### NOTE: CURRENTLY BROKEN & DOES NOT RESULT IN A USABLE OBJECT
# Error in me_model_by_row(obj, X, obj$bs_summary, which_var) : 
#   names(bs_summary[[sigma_var]]) and rownames(bs_summary[[which_var]]) are not equal:
#   Lengths (63801, 53015) differ (string compare on first 53015)
# 52428 string mismatches

# Check for consistency

sleuth_object_test_filt <- sleuth_fit(sleuth_object_test_filt, ~  M0_KCVSEX + M0_KCVAGE, 'reduced')


# Fit full model with covariates and condition | + M0_KCVSEX + M0_KCVAGE
sleuth_object_test_filt <- sleuth_fit(sleuth_object_test_filt, ~ RSTR + TB + RSTR:TB + M0_KCVSEX + M0_KCVAGE, 'full')

sleuth_object_test_filt <- sleuth_fit(sleuth_object_test_filt, ~ RSTR + TB + M0_KCVSEX + M0_KCVAGE, 'full_no_interaction')

# Perform DEG analysis via likelihood ratio test
sleuth_object_test_filt <- sleuth_lrt(sleuth_object_test_filt, 'reduced', 'full') # main model effect (tx sig for main effects)
sleuth_object_test_filt <- sleuth_lrt(sleuth_object_test_filt, 'full_no_interaction', 'full') # interaction effect (tx sig for interaction)


### Compare outputs

# compare pvals
data.frame(pval_og = SO_basic_PC$tests$lrt$`reduced:full`$pval, pval_test = sleuth_object_test_filt$tests$lrt$`reduced:full`$pval)%>%
  ggplot(aes(x=pval_og, y=pval_test))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="red")

# compare mean_obs
data.frame(mean_obs_og = SO_basic_PC$tests$lrt$`reduced:full`$mean_obs, mean_obs_test = sleuth_object_test_filt$tests$lrt$`reduced:full`$mean_obs)%>%
  ggplot(aes(x=mean_obs_og, y=mean_obs_test))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="red")

boxplot(x = (mean_obs_og = SO_basic_PC$tests$lrt$`reduced:full`$mean_obs-sleuth_object_test_filt$tests$lrt$`reduced:full`$mean_obs))

#NOTE: the means obs values are NEARLY identical, but not. Very interesting. 

