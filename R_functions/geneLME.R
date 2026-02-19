
########################################################
# Scalable custom gene LMEs with contrast specification
########################################################

########################################################
# Starting point / framework
########################################################


geneLME_fit<-
  function(i, 
           dat=dat, 
           formula_obj=formula_obj, 
           model_weights=model_weights,
           run_contrast=run_contrast,
           contrast_vars = contrast_vars,
           contrast_var_2_levels=contrast_var_2_levels,
           contrasts_primary = contrasts_primary,
           contrasts_secondary = contrasts_secondary){
    
    gene_i <- rownames(dat$E)[i]
    
    # Initialize an empty status
    model_status <- "success"
    
    # Wrap the entire process in tryCatch
    result <- tryCatch({
      
      # --- FIT MODEL ---
      # We use a localized data frame to avoid passing the whole 'dat' object to emmeans later
      model_data <- dat$targets
      model_data$expression <- dat$E[i, ]
      
        if(is.null(model_weights)){
        lme_i <- lmer(
          formula_obj, 
          data = model_data,
          control = lmerControl(calc.derivs = FALSE, autoscale = TRUE))
        
        } else if(!is.null(model_weights) & model_weights){
          
          model_data$current_weights <- dat$weights[i, ]
          
          lme_i <- lmer(
            formula_obj, 
            weights = current_weights, 
            data = model_data,
            control = lmerControl(calc.derivs = FALSE, autoscale = TRUE)
          )
      }
      
      # Check for singularity manually to be explicit
      if(isSingular(lme_i)) stop("Boundary (singular) fit")
      
      # --- EXTRACT DATA (The only parts we keep) ---
      
      # 1. AIC
      aic_res <- data.frame(gene = gene_i, AIC = AIC(lme_i))
      
      # 2. ANOVA & Summary
      lme_i_anova <- car::Anova(lme_i) %>% broom.mixed::tidy()
      lme_i_summary <- summary(lme_i)$coefficients %>% 
        as.data.frame() %>% 
        rownames_to_column("variable")
      
      # 3. Process ANOVA Table
      lme_i_anova_tab <- lme_i_anova %>%
        rowwise()%>%
        mutate(
          gene = gene_i,
          predictor_class = case_when(
            grepl(":", term) ~ "interaction",
            #!term %in% colnames(model_data) ~ "other",
            is.numeric(model_data[[term]]) ~ "continuous",
            !(is.numeric(model_data[[term]])) & length(unique(model_data[[term]])) == 2 ~ "two-level-categorical",
            !(is.numeric(model_data[[term]])) & length(unique(model_data[[term]])) > 2 ~ "multi-level-categorical"
          ),
          Estimate_source = case_when(
            predictor_class %in% c("continuous", "two-level-categorical") ~ "lme_summary",
            predictor_class == "multi-level-categorical" ~ "seeContrasts",
            predictor_class == "interaction" &  
              length(grep(":", lme_i_summary$variable)) == 1 ~ "lme_summary",
            predictor_class == "interaction" &  
              length(grep(":", lme_i_summary$variable)) > 1 ~ "seeContrasts"
          ),
          Estimate = case_when(
            predictor_class == "continuous" ~ lme_i_summary$Estimate[match(term, lme_i_summary$variable)][[1]],
            predictor_class == "two-level-categorical" ~ lme_i_summary$Estimate[match(term, lme_i_summary$variable)][[1]],
            predictor_class == "interaction" &  Estimate_source=="lme_summary" ~ lme_i_summary$Estimate[grep(":", lme_i_summary$variable)][[1]],
            predictor_class == "interaction" & Estimate_source == "seeContrasts" ~ NA,
            .default = NA
          ),
          Estimate_SE = case_when(
            predictor_class == "continuous" ~ lme_i_summary$`Std. Error`[match(term, lme_i_summary$variable)][[1]],
            predictor_class == "two-level-categorical" ~ lme_i_summary$`Std. Error`[match(term, lme_i_summary$variable)][[1]],
            predictor_class == "interaction" &  Estimate_source=="lme_summary" ~ lme_i_summary$`Std. Error`[grep(":", lme_i_summary$variable)][[1]],
            predictor_class == "interaction" & Estimate_source == "seeContrasts" ~ NA,
            .default = NA
          )
        )
      
      # 4. Contrasts
      
      # --- CONTRASTS (Dynamic) ---
    
      if(!is.null(run_contrast)){
      
      spec_formula <- as.formula(paste("~", paste(contrast_vars, collapse = "|")))
      
      if(length(contrast_vars)==2 & !is.null(contrast_var_2_levels)){ # set levels of second term if multiple terms provided.
      by_list<-set_names(list(contrast_var_2_levels), paste(contrast_vars[2]))} else {
      by_list<-NULL
      }
      
      emm_1st <- emmeans(lme_i, spec = spec_formula, at=by_list, data=model_data) %>%
        contrast(method = contrasts_primary, adjust = "none")
      
      if(!is.null(contrasts_secondary) & length(contrasts_primary)>1){
      
      emm_2nd <- contrast(emm_1st, method = contrasts_secondary, adjust = "none")
      
      contrast_res <- bind_rows(
        as.data.frame(emm_1st) %>% mutate(contrast_order = "first_order"),
        as.data.frame(emm_2nd) %>% mutate(contrast_order = "second_order")
      ) %>% mutate(gene = gene_i)
      } else {
        contrast_res <-
          as.data.frame(emm_1st) %>% mutate(contrast_order = "first_order")
      }} else {contrast_res=NULL}
      
      # Return list of results. Once this function returns, 'lme_i' is destroyed.
      list(aic = aic_res, anova = lme_i_anova_tab, contrasts = contrast_res, model_status=setNames(c("success"), gene_i))
      
    }, error = function(e) {
      # This block triggers for actual errors AND our forced 'singular fit' stop
      err_msg <- as.character(e$message)
      
      list(
        aic = data.frame(gene = gene_i, AIC = NA, model_status = err_msg),
        # Return a 1-row ANOVA skeleton so the gene is not lost
        
        anova = data.frame(
          term = NA, statistic=NA, df=NA,  p.value =NA, gene=gene_i, 
          predictor_class=NA, Estimate_source=NA,  Estimate=NA,  Estimate_SE=NA, 
          model_status = err_msg),
        # Return a 1-row Contrast skeleton
        contrasts=data.frame(contrast =NA, estimate=NA, SE=NA, df=NA,  t.ratio=NA, p.value=NA,  
                             contrast_order=NA,   gene=gene_i, model_status = err_msg),
        model_status= setNames(paste(err_msg), gene_i)
      )
    })
    
    
    return(result)
  }



# Separate function for compiling results

geneLME_compiler<-function(fit){
  list(lme_anova_tab= map_dfr(fit, "anova"),   # gather anova tables
       lme_contrast=map_dfr(fit, "contrasts"), # gather contrast tables
       lme_fit=map_dfr(fit, "aic"), # gather model fit
       lme_err=map_chr(fit, "model_status") %>% # extract the status string for each gene & preserve names
         set_names(map_chr(fit, ~ names(.x$model_status))))
}



geneLME <- 
  function(dat, 
           formula_str, 
           model_weights=NULL,
           run_contrast=NULL,
           contrast_vars=NULL, 
           contrast_var_2_levels=NULL,
           contrasts_primary=NULL, 
           contrasts_secondary=NULL, 
           n_cores = NULL) {
    
#.     --- PRE-FLIGHT DATA VALIDATION ---
      
      # 1. Check Formula Variables
      # Extract all variables from the formula (excluding the LHS 'expression')
      required_vars <- all.vars(formula_obj)
      required_vars <- required_vars[required_vars != "expression"]
      
      missing_vars <- setdiff(required_vars, colnames(dat$targets))
      if (length(missing_vars) > 0) {
        stop(paste("The following variables in the formula are missing from dat$targets:", 
                   paste(missing_vars, collapse = ", ")))
      }
      
      # 2. Check Weights Alignment
      if (!is.null(model_weights)) {
        if (is.null(dat$weights)) {
          stop("model_weights is specified, but dat$weights is missing.")
        }
        if (nrow(dat$weights) != nrow(dat$E) || ncol(dat$weights) != ncol(dat$E)) {
          stop("Dimensions of dat$weights do not match dat$E. Weights must be a matrix of same size.")
        }
      }
      
      # 3. Check Contrast Variables & Levels
      if (!is.null(run_contrast)) {
        # Verify contrast_vars exist in targets
        missing_contrast_vars <- setdiff(contrast_vars, colnames(dat$targets))
        if (length(missing_contrast_vars) > 0) {
          stop(paste("contrast_vars missing from dat$targets:", 
                     paste(missing_contrast_vars, collapse = ", ")))
        }
        
        # Verify Level Filtering for the 'by' variable
        if (length(contrast_vars) == 2 && !is.null(contrast_var_2_levels)) {
          actual_levels <- unique(as.character(dat$targets[[contrast_vars[2]]]))
          invalid_levels <- setdiff(contrast_var_2_levels, actual_levels)
          
          if (length(invalid_levels) > 0) {
            stop(paste0("The specified levels for '", contrast_vars[2], 
                        "' do not exist in the data: ", paste(invalid_levels, collapse = ", ")))
          }
        }
      }
      
      # 4. Check for Character columns that should be Factors
      # emmeans and lmer handle strings, but it's safer to warn if a predictor has 100+ unique strings
      for (v in required_vars) {
        if (is.character(dat$targets[[v]]) && length(unique(dat$targets[[v]])) > (nrow(dat$targets)/2)) {
          warning(paste("Variable '", v, "' has a very high number of unique strings. Is it a unique ID instead of a factor?"))
        }
      }
      
      message("Input check complete: Data formatting looks good. Starting parallel LME fit...")
    
    
    
  
  # Set Parallel Plan
  workers <- if(is.null(n_cores)) parallel::detectCores() - 4 else n_cores
  plan(multisession, workers = workers)
  
  # Convert string formula to object
  
  formula_obj <- as.formula(paste0("expression", formula_str))
  
  # Run in parallel
  # We pass the contrast specs into the function
  fit_results <- future_lapply(
    1:nrow(dat$E),
    FUN =  function(x){
    geneLME_fit(i = x, 
    dat = dat, 
    formula_obj = formula_obj,
    model_weights=model_weights,
    run_contrast=run_contrast,
    contrast_vars = contrast_vars,
    contrast_var_2_levels=contrast_var_2_levels,
    contrasts_primary = contrasts_primary,
    contrasts_secondary = contrasts_secondary)},
    future.seed = TRUE
  )
  
  return(geneLME_compiler(fit_results))
}


########################################################
#     Testing
########################################################
# dat_sub<-RNAetc::subset_voom(dat, gene_keep = head(rownames(dat$E), 100))
# 
# test_mods<-
#   geneLME(
#     dat = dat_sub,
#     formula_str = "~ treatment*visit + rNANgUl + percent_duplication + median_cv_coverage + (1|ptID)",
#     model_weights = TRUE, 
#     run_contrast=TRUE,
#     contrast_vars = c("treatment", "visit"), 
#     contrast_var_2_levels=c("-1", "S9", "S13"),
#     contrasts_primary = list("Dual vs Placebo"=c(-1, 0, 1), "SLIT vs Placebo"=c(-1, 1, 0)),
#     contrasts_secondary = list("Second Order Test"= c(1, -1)),
#     n_cores = 10
#     
#   )


#PICKUP HERE:

# needs to debug scaling up and hadling of parallel tasks. We didn't sort this out perfectlt. 
# need to modify how contrasts are handled to be able to specify multiple levels of interaction contrasts. 

