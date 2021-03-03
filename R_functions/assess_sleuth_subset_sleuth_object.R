
### WARNING: This script is in development and should not be used yet. 

sleuth_object<-SI_adv60_PC # assign sleuth import loaded elsewhere

samples<-sleuth_object$sample_to_covariates%>%filter(condition %in% c("RSTR_TB", "RSTR_MEDIA"))%>%dplyr::select(sample)%>%unlist()


# Fit models on non-subsetted sleuth object for comparison
sleuth_object_no_filt <- sleuth_fit(sleuth_object, ~  M0_KCVSEX + M0_KCVAGE, 'reduced')


# Fit full model with covariates and condition | + M0_KCVSEX + M0_KCVAGE
sleuth_object_no_filt <- sleuth_fit(sleuth_object_no_filt, ~ condition + M0_KCVSEX + M0_KCVAGE, 'full')

# Perform DEG analysis via likelihood ratio test
sleuth_object_no_filt <- sleuth_lrt(sleuth_object_no_filt, 'reduced', 'full') # main model effect (tx sig for main effects)




# Test subsetting approach for samples

sleuth_object_test<-sleuth_object

sleuth_object_sampleTest<- subset_sleuth_samples(sleuth_object_test, samples = samples)


# Fit models

sleuth_object_sampleTest_filt <- sleuth_fit(sleuth_object_sampleTest, ~  M0_KCVSEX + M0_KCVAGE, 'reduced')


# Fit full model with covariates and condition | + M0_KCVSEX + M0_KCVAGE
sleuth_object_sampleTest_filt <- sleuth_fit(sleuth_object_sampleTest_filt, ~ condition + M0_KCVSEX + M0_KCVAGE, 'full')

# Perform DEG analysis via likelihood ratio test
sleuth_object_sampleTest_filt <- sleuth_lrt(sleuth_object_sampleTest_filt, 'reduced', 'full') # main model effect (tx sig for main effects)

#------------

# Test subsetting for targets
targets<- sleuth_object_test_filt$tests$lrt$`reduced:full`$target_id[1:15000]

sleuth_object_targetTest_filt<-subset_sleuth_targets(sleuth_object_test, targets = targets)

# Fit models

sleuth_object_targetTest_filt <- sleuth_fit(sleuth_object_targetTest_filt, ~  M0_KCVSEX + M0_KCVAGE, 'reduced')

# Fit full model with covariates and condition | + M0_KCVSEX + M0_KCVAGE
sleuth_object_targetTest_filt <- sleuth_fit(sleuth_object_targetTest_filt, ~ condition + M0_KCVSEX + M0_KCVAGE, 'full')

# Perform DEG analysis via likelihood ratio test
sleuth_object_targetTest_filt <- sleuth_lrt(sleuth_object_targetTest_filt, 'reduced', 'full') # main model effect (tx sig for main effects)


### Compare outputs
results_no_filt<- sleuth_object_no_filt$tests$lrt$`reduced:full`
results_sample_filt<-sleuth_object_sampleTest_filt$tests$lrt$`reduced:full`
results_target_filt<-sleuth_object_targetTest_filt$tests$lrt$`reduced:full`

### Plots to compare estimates and p values
target_filtering_comp<-
  filter(results_no_filt, target_id %in% targets)%>% # filter to common targets
  mutate(analysis="No Filter")%>%
  bind_rows(mutate(results_target_filt,  analysis="Target Filter"))%>%
  dplyr::select(target_id,analysis, pval, qval, mean_obs, var_obs)%>%
  pivot_longer(c(pval, qval, mean_obs,var_obs))%>%
  pivot_wider(names_from = c(analysis))%>%
  mutate(name= factor(name, levels = c("mean_obs", "var_obs", "pval", "qval")))%>%
  ggplot(aes(x=`No Filter`, y= `Target Filter`))+
  geom_point(alpha=0.5)+
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color = "red")+
  facet_wrap(~name, scales="free")+
  labs(title = "Comparison With Target Filter")+
  theme(plot.title = element_text(hjust=0.5))

sample_filtering_comp<-
  mutate(results_no_filt, analysis="No Filter")%>%
  bind_rows(mutate(results_sample_filt,  analysis="Sample Filter"))%>%
  dplyr::select(target_id,analysis, pval, qval, mean_obs, var_obs)%>%
  pivot_longer(c(pval, qval, mean_obs,var_obs))%>%
  pivot_wider(names_from = c(analysis))%>%
  mutate(name= factor(name, levels = c("mean_obs", "var_obs", "pval", "qval")))%>%
  ggplot(aes(x=`No Filter`, y= `Sample Filter`))+
  geom_point(alpha=0.5)+
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color = "red")+
  facet_wrap(~name, scales="free")+
  labs(title = "Comparison With Sample Filter")+
  theme(plot.title = element_text(hjust=0.5))

# NOTE: Target Filters do reveal some odd behavior in that pvals change when the target list is filtered (I cannot think of why this would be)

