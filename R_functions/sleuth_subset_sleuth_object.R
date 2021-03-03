


subset_sleuth_samples=function(sleuth_object, samples = NULL){
  
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
  return(sleuth_object)
  }

subset_sleuth_targets=function(sleuth_object, targets = NULL){
  
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






