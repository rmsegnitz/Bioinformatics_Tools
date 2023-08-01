library(tidyverse)
library(kimma)
library(foreach)
source("R_functions/kmFit_resid.R")

# Unmodified kimma
km_test<-
       kmFit(dat = example.voom,
                         run_lm = TRUE, use_weights = FALSE,
                         subset_var = list("asthma"), subset_lvl = list(c("asthma")),
                         subset_genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
                         model = "~ virus + (1|ptID)", metrics = TRUE)

names(km_test)

# Modified, but not calling residuals
km_resid_test_noresid<-
  kmFit_resid(dat = example.voom,
              run_lm = TRUE, use_weights = FALSE,
              subset_var = list("asthma"), subset_lvl = list(c("asthma")),
              subset_genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
              model = "~ virus + (1|ptID)", metrics=T, residuals=F)

names(km_resid_test_noresid)

# Modified returning residuals
km_resid_test<-
  kmFit_resid(dat = example.voom,
        run_lm = TRUE, use_weights = FALSE,
        subset_var = list("asthma"), subset_lvl = list(c("asthma")),
        subset_genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
        model = "~ virus + (1|ptID)", metrics=TRUE, residuals=T)

names(km_resid_test)
view(km_resid_test$lm.residuals)
view(km_resid_test$lm.error)

# Compare to ensure accuracy

km_test$lm.fit
km_resid_test_noresid$lm.fit
km_resid_test$lm.fit

km_test$lm
km_resid_test_noresid$lm
km_resid_test$lm


# Manually fit a model to check accuracy of residuals

test_mod_simple<-
  lm(gene ~ virus, 
     data=filter(example.voom$targets, asthma=="asthma")%>%
       left_join(example.voom$E[km_resid_test$lm.fit$gene[1], ]%>%as.data.frame()%>%rownames_to_column("libID")%>%rename("gene"="."))
       )

summary(test_mod_simple)

residuals(test_mod_simple)

### TEST USE IN PRACTICE ###

km_residualized_1<-
  kmFit_resid(dat = example.voom,
              run_lm = TRUE, use_weights = TRUE,
              model = "~ virus + (1|ptID)", metrics=TRUE, residuals=TRUE)



example.voom_residualized_1<-example.voom
example.voom_residualized_1$E<-km_residualized_1$lm.residuals

km_residualized_1_post<-
  kmFit_resid(dat = example.voom_residualized_1,
              run_lm = TRUE, use_weights = TRUE,
              model = "~ asthma", metrics=TRUE)

km_residualized_1_post_comparison<-
  kmFit_resid(dat = example.voom,
              run_lm = TRUE, use_weights = TRUE,
              model = "~ virus + asthma + (1|ptID)", metrics=TRUE, residuals=TRUE)


# Test Use with lme

# Modified returning residuals
km_resid_test_lme<-
  kmFit_resid(dat = example.voom,
              run_lme = TRUE, use_weights = FALSE,
              subset_var = list("asthma"), subset_lvl = list(c("asthma")),
              subset_genes = c("ENSG00000250479","ENSG00000250510","ENSG00000255823"),
              model = "~ virus + (1|ptID)", metrics=TRUE, residuals=T)

names(km_resid_test_lme)
view(km_resid_test_lme$lme.residuals)

