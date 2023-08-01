# Modification of {kimma} differential expression model fitting to support estimation of model residuals 
# for downstream permutational analyses such as GSEA or network estimation.

# Author: Max Segnitz, msegnitz@uw.edu
# Started July 2023
#

#
# DESCRIPTION:
# Impliments minor augmentations to the core kimma model fitting components and internal functions to 
# extract and return model residuals if desired. When residuals=TRUE is specified, the kmFit results 
# include a GENE x SAMPLE matrix of model residuals.
#

#'##############################
### 1) kmFit_resid()   ####
#'##############################



## Modification of kmFit to return model residuals ##
kmFit_resid<-
function (dat = NULL, kin = NULL, patientID = "ptID", libraryID = "libID", 
          counts = NULL, meta = NULL, genes = NULL, weights = NULL, 
          subset_var = NULL, subset_lvl = NULL, subset_genes = NULL, 
          model, use_weights = FALSE, run_lm = FALSE, run_lme = FALSE, 
          run_lmerel = FALSE, metrics = FALSE, residuals=FALSE,  run_contrast = FALSE, 
          contrast_var = NULL, processors = NULL, p_method = "BH", 
          genotype_name = NULL, run.lmekin = NULL, subset.var = NULL, 
          subset.lvl = NULL, subset.genes = NULL, use.weights = FALSE, 
          run.lm = FALSE, run.lme = FALSE, run.lmerel = FALSE, run.contrast = FALSE, 
          contrast.var = NULL, p.method = NULL) 
{
  rowname <- libID <- variable <- statistic <- df <- pval <- group <- gene <- V1 <- V2 <- combo <- term <- p.value <- estimate <- contrast <- contrast.i <- weights.gene <- FDR <- contrast_ref <- contrast_lvl <- NULL
  if (any(!is.null(subset.var), !is.null(subset.lvl), !is.null(subset.genes), 
          use.weights, run.lm, run.lme, run.lmerel, run.contrast, 
          !is.null(contrast.var), !is.null(p.method))) {
    message("WARNING: Arguments with '.' have been deprecated. Please use '_' versions.")
  }
  if (!is.null(subset.var)) {
    subset_var <- subset.var
  }
  if (!is.null(subset.lvl)) {
    subset_lvl <- subset.lvl
  }
  if (!is.null(subset.genes)) {
    subset_genes <- subset.genes
  }
  if (use.weights) {
    use_weights <- use.weights
  }
  if (run.lm) {
    run_lm <- run.lm
  }
  if (run.lme) {
    run_lme <- run.lme
  }
  if (run.lmerel) {
    run_lmerel <- run.lmerel
  }
  if (run.contrast) {
    run_contrast <- run.contrast
  }
  if (!is.null(contrast.var)) {
    contrast_var <- contrast.var
  }
  if (!is.null(p.method)) {
    p_method <- p.method
  }
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    processors.to.use <- 2
  }else if (is.null(processors)) {
    processors.to.use <- parallel::detectCores() - 2
    if (processors.to.use == 0) {
      stop("Error processors: Default resulted in 0. Please correct.")
    }
  } else {
    processors.to.use <- processors
  }
  cl <- parallel::makeCluster(processors.to.use)
  if (!is.null(dat) & !(libraryID %in% colnames(dat$targets))) {
    stop("libraryID column not found in dat$targets.")
  }
  if (!is.null(meta) & !(libraryID %in% colnames(meta))) {
    stop("libraryID column not found in meta.")
  }
  if (!is.null(dat) & !(patientID %in% colnames(dat$targets))) {
    stop("patientID column not found in dat$targets.")
  }
  if (!is.null(meta) & !(patientID %in% colnames(meta))) {
    stop("patientID column not found in meta.")
  }
  if (is.null(subset_var) & !is.null(subset_lvl)) {
    stop("Sample subsetting has been selected. Please also provide subset_var")
  }
  if (!is.null(subset_var) & is.null(subset_lvl)) {
    stop("Sample subsetting has been selected. Please also provide subset_lvl")
  }
  if (run_lmerel & !grepl("\\|", model)) {
    stop("Kinship models require a random effect in the model as in (1 | ptID)")
  }
  if (is.null(kin) & run_lmerel) {
    stop("Kinship matrix is required to run lmerel")
  }
  if (!run_lm & !run_lme & !run_lmerel & !run_contrast) {
    stop("At least 1 model type must be selected. Please set one run parameter to TRUE.")
  }
  if (!run_lm & !run_lme & !run_lmerel & run_contrast) {
    stop("Contrast models must be run with an accompanying linear model.")
  }
  if (use_weights & is.null(weights) & is.null(dat$weights)) {
    stop("When use_weights is TRUE, must provide gene weights is dat object or separate data frame.")
  }
  if ("gene_weight" %in% c(colnames(meta), colnames(dat$targets))) {
    stop("Variable gene_weight is present in dat$targets or meta. This name is used for model weights. Please change variable name in your data.")
  }
  if (!is.null(run.lmekin)) {
    stop("run.lmekin no longer supported. Please use run_lmerel")
  }
  if (grepl("\n", model)) {
    stop("Model cannot contain hard returns \n. Please correct.")
  }
  if (!use_weights) {
    if (!is.null(weights) | !is.null(dat$weights)) {
      message("WARNING: To use weights provided in dat$weights or weights, set use_weights = TRUE\n")
    }
  }
  if (grepl("\\|", model)) {
    model.temp <- strsplit(gsub(" ", "", model), split = "\\+\\(1")[[1]][1]
    model.lm <- paste("expression", model.temp, sep = "")
  } else {
    model.lm <- paste("expression", model, sep = "")
    model.lm <- gsub(" ", "", model.lm)
  }
  model.lme <- paste("expression", gsub(" ", "", model), sep = "")
  if (run_lm) {
    message(paste("lm model:", model.lm))
  }
  if (run_lme | run_lmerel) {
    message(paste("lme/lmerel model:", model.lme))
  }
  if (run_contrast & is.null(contrast_var)) {
    contrast.temp <- strsplit(gsub(" |~", "", model), split = "\\+\\(1")[[1]][1]
    contrast.main <- strsplit(contrast.temp, split = "\\+|\\*")[[1]]
    if (grepl("\\*", model)) {
      contrast.interact <- strsplit(contrast.temp, split = "\\+")[[1]]
      contrast.interact <- contrast.interact[grepl("\\*", 
                                                   contrast.interact)]
      if (any(stringr::str_count(contrast.interact, "\\*") > 
              1)) {
        stop("Contrasts of triple interactions are not supported.")
      }
      contrast.interact <- gsub("\\*", ":", contrast.interact)
    }
    else {
      contrast.interact <- NULL
    }
    contrast_var <- c(contrast.main, contrast.interact)
    if (!is.null(genotype_name)) {
      contrast_var <- gsub("genotype", genotype_name, 
                           contrast_var)
    }
  }
  if (!is.null(contrast_var)) {
    run_contrast <- TRUE
  }
  if (run_contrast) {
    if (!is.null(dat)) {
      meta.temp <- dat$targets
    }
    else {
      meta.temp <- meta
    }
    for (v in contrast_var) {
      if (grepl("[*]|:", v)) {
        v.sep <- strsplit(v, split = "[*]|:")[[1]]
        v.class1 <- class(meta.temp[, v.sep[1]])
        v.class2 <- class(meta.temp[, v.sep[2]])
        if (all(c(v.class1, v.class2) %in% c("numeric", 
                                             "integer"))) {
          stop(paste("Contrast variable", v, "is numeric. Please specify only character/factor contrasts."))
        }
      }
      else {
        v.class <- class(meta.temp[, v])
        if (v.class %in% c("numeric", "integer")) {
          stop(paste("Contrast variable", v, "is numeric. Please specify only character/factor contrasts."))
        }
      }
    }
  }
  to.model.ls <- kimma:::kimma_cleaning(dat, kin, patientID, libraryID, 
                                counts, meta, genes, weights, subset_var, subset_lvl, 
                                subset_genes, model.lm, genotype_name)
  fit.results <- data.frame()
  doParallel::registerDoParallel(cl)
  
  
  fit.results <- data.table::rbindlist(fill = TRUE, foreach::foreach(gene = unique(to.model.ls[["to.model"]]$rowname), 
                                                                     .packages = c("dplyr", "tidyr",  "magrittr", "stats", "broom", 
                                                                                   "lme4", "car", "tibble", "lme4qtl", "utils", "emmeans", 
                                                                                   "data.table", "foreach", "doParallel"), .export = c("kimma_lm_resid", 
                                                                                                                                       "kimma_lme_resid", 
                                                                                                                                       "kimma_lmerel_resid", 
                                                                                                                                       "kmFit_contrast_resid")) %dopar% 
                                         {
                                           to.model.gene <- to.model.ls[["to.model"]] %>% dplyr::filter(rowname == 
                                                                                                          gene) %>% dplyr::arrange(patientID)
                                           results.lm.ls <- NULL
                                           if (run_lm) {
                                             results.lm.ls <- tryCatch({
                                               kimma_lm_resid(model.lm, to.model.gene, gene, use_weights,
                                                        metrics, residuals=residuals)
                                             }, error = function(e) {
                                               results.lm.ls[["error"]] <- data.frame(model = "lm", 
                                                                                      gene = gene, message = conditionMessage(e))
                                               return(results.lm.ls)
                                             })
                                           }
                                           results.lme.ls <- NULL
                                           if (run_lme) {
                                             results.lme.ls <- tryCatch({
                                               kimma_lme_resid(model.lme, to.model.gene, gene,
                                                         use_weights, metrics, residuals=residuals)
                                             }, error = function(e) {
                                               results.lme.ls[["error"]] <- data.frame(model = "lme", 
                                                                                       gene = gene, message = conditionMessage(e))
                                               return(results.lme.ls)
                                             
                                             })
                                           }
                                           results.kin.ls <- NULL
                                           if (run_lmerel) {
                                             results.kin.ls <- tryCatch({
                                               kimma_lmerel_resid(model.lme, to.model.gene, gene, 
                                                            to.model.ls[["kin.subset"]], use_weights, 
                                                            patientID, metrics, residuals=residuals,)
                                             }, error = function(e) {
                                               results.kin.ls[["error"]] <- data.frame(model = "lmerel", 
                                                                                       gene = gene, message = conditionMessage(e))
                                               return(results.kin.ls)
                                             
                                             })
                                           }
                                           contrast.lm <- NULL
                                           contrast.lme <- NULL
                                           contrast.kin <- NULL
                                           contrast.results <- NULL
                                           if (run_contrast) {
                                             if (!is.null(results.lm.ls[["results"]])) {
                                               contrast.lm <- tryCatch({
                                                 kimma:::kmFit_contrast(results.lm.ls[["fit"]], contrast_var, 
                                                                to.model.gene, genotype_name) %>% dplyr::mutate(model = "lm.contrast")
                                               }, error = function(e) {
                                                 contrast.lm.error <- data.frame(model = "lm.contrast", 
                                                                                 gene = gene, message = conditionMessage(e))
                                                 return(contrast.lm.error)
                                               })
                                             }
                                             if (!is.null(results.lme.ls[["results"]])) {
                                               contrast.lme <- tryCatch({
                                                 kimma:::kmFit_contrast(results.lme.ls[["fit"]], 
                                                                contrast_var, to.model.gene, genotype_name) %>% 
                                                   dplyr::mutate(model = "lme.contrast")
                                               }, error = function(e) {
                                                 contrast.lme.error <- data.frame(model = "lme.contrast", 
                                                                                  gene = gene, message = conditionMessage(e))
                                                 return(contrast.lme.error)
                                               })
                                             }
                                             if (!is.null(results.kin.ls[["results"]])) {
                                               contrast.kin <- tryCatch({
                                                 kimma:::kmFit_contrast(results.kin.ls[["fit"]], 
                                                                contrast_var, to.model.gene, genotype_name) %>% 
                                                   dplyr::mutate(model = "lmerel.contrast")
                                               }, error = function(e) {
                                                 contrast.kin.error <- data.frame(model = "lmerel.contrast", 
                                                                                  gene = gene, message = conditionMessage(e))
                                                 return(contrast.kin.error)
                                               })
                                             }
                                             contrast.results <- dplyr::bind_rows(contrast.lm, 
                                                                                  contrast.lme, contrast.kin) %>% dplyr::mutate(gene = gene)
                                           }
                                           if (any(is.character(results.lm.ls[["results"]]$estimate), 
                                                   is.character(results.lme.ls[["results"]]$estimate), 
                                                   is.character(results.kin.ls[["results"]]$estimate), 
                                                   is.character(contrast.results$estimate))) {
                                             results.lm.ls[["results"]]$estimate <- as.character(results.lm.ls[["results"]]$estimate)
                                             results.lme.ls[["results"]]$estimate <- as.character(results.lme.ls[["results"]]$estimate)
                                             results.kin.ls[["results"]]$estimate <- as.character(results.kin.ls[["results"]]$estimate)
                                             contrast.results$estimate <- as.character(contrast.results$estimate)
                                             results.lm.ls[["results"]]$statistic <- as.character(results.lm.ls[["results"]]$statistic)
                                             results.lme.ls[["results"]]$statistic <- as.character(results.lme.ls[["results"]]$statistic)
                                             results.kin.ls[["results"]]$statistic <- as.character(results.kin.ls[["results"]]$statistic)
                                             contrast.results$statistic <- as.character(contrast.results$statistic)
                                           }
                                           fit.results <- results.lm.ls[["results"]] %>% dplyr::bind_rows(results.lme.ls[["results"]]) %>% 
                                             dplyr::bind_rows(results.kin.ls[["results"]]) %>% 
                                             dplyr::bind_rows(contrast.results) %>% dplyr::bind_rows(results.lm.ls[["error"]]) %>% 
                                             dplyr::bind_rows(results.lme.ls[["error"]]) %>% 
                                             dplyr::bind_rows(results.kin.ls[["error"]]) %>% 
                                             dplyr::bind_rows(results.lm.ls[["metrics"]]) %>% 
                                             dplyr::bind_rows(results.lme.ls[["metrics"]]) %>% 
                                             dplyr::bind_rows(results.kin.ls[["metrics"]])%>%
                                             dplyr::bind_rows(tibble(
                                               model=case_when(run_lm ~"lm.residuals", run_lme ~ "lme.residuals", run_lmerel ~ "lmerel.residuals"),
                                               gene=gene,
                                               residuals=c(results.lm.ls[["residuals"]], results.lme.ls[["residuals"]], results.kin.ls[["residuals"]])[
                                                                         which(!is.null(c(results.lm.ls[["residuals"]], results.lme.ls[["residuals"]], results.kin.ls[["residuals"]])))]
                                                                       ))
                                           
                                              
                                           
                                         })
  
  parallel::stopCluster(cl)
  all <- length(unique(to.model.ls[["to.model"]]$rowname))
  
  message("Complete: ", all, " genes")
  if ("message" %in% colnames(fit.results)) {
    fail <- fit.results %>% tidyr::drop_na(message) %>% 
      dplyr::distinct(gene, message) %>% nrow()
    message("Failed: ", fail, " genes. See results[['model_error']]")
    error.results <- fit.results %>% dplyr::filter(!is.na(message)) %>% 
      dplyr::distinct(model, gene, message)
    fit.results <- fit.results %>% dplyr::filter(is.na(message)) %>% 
      dplyr::select(-message)
  }else {
    message("Failed: 0 genes")
    error.results <- NULL
  }
  
  
  kmFit.ls <- list()
  if (nrow(fit.results) > 0) {
    if (run_contrast) {
      kmFit.results <- fit.results %>% dplyr::group_by(model, 
                                                       variable, contrast_ref, contrast_lvl) %>% dplyr::mutate(FDR = stats::p.adjust(pval, 
                                                                                                                                     method = p_method)) %>% dplyr::ungroup() %>% 
        dplyr::select(model:statistic, contrast_ref, 
                      contrast_lvl, estimate, pval, FDR, dplyr::everything())
    }else {
      kmFit.results <- fit.results %>% dplyr::group_by(model, 
                                                       variable) %>% dplyr::mutate(FDR = stats::p.adjust(pval, 
                                                                                                         method = p_method)) %>% dplyr::ungroup() %>% 
        dplyr::select(model:statistic, estimate, pval, 
                      FDR, dplyr::everything())
      
    }
    kmFit.results.anno <- kmFit.results
    for (result.i in unique(kmFit.results.anno$model)) {
      result.temp <- dplyr::filter(kmFit.results.anno, 
                                   model == result.i)
      estimates <- unique(result.temp$estimate)
      estimates <- estimates[!is.na(estimates)]
      if (all(estimates != "seeContrasts") & !is.null(estimates)) {
        result.temp <- dplyr::filter(kmFit.results.anno, 
                                     model == result.i) %>% dplyr::mutate(estimate = as.numeric(estimate), 
                                                                          statistic = as.numeric(statistic))
      }
      kmFit.ls[[result.i]] <- result.temp %>% dplyr::select_if(function(x) any(!is.na(x)))%>%dplyr::select_if(function(x) !all(is.null(unlist(x))))
      
    }
  }
  
  if(residuals){
    temp_residuals<- do.call(rbind,  kmFit.ls[[grep("residuals", names(kmFit.ls))]]$residuals)
    rownames(temp_residuals)<-kmFit.ls[[grep("residuals", names(kmFit.ls))]]$gene
    colnames(temp_residuals)<-to.model.ls[["to.model"]]%>%arrange(patientID)%>%pull(libID)%>%unique()
    kmFit.ls[[grep("residuals", names(kmFit.ls))]]<-temp_residuals
  }
  
  
  if (!is.null(error.results)) {
    for (result.i in unique(error.results$model)) {
      kmFit.ls[[paste(result.i, "error", sep = ".")]] <- dplyr::filter(error.results, 
                                                                       model == result.i)
    }
  }
  
  if(!residuals){kmFit.ls<-kmFit.ls[which(!grepl("residuals", names(kmFit.ls)))]} # remove residuals matrix if not requested. This could/should be better implemented upstream
  return(kmFit.ls)
}

#'##############################
### 2) kimma_lm_resid()   ####
#'##############################

kimma_lm_resid <- function(model_lm, to_model_gene, gene, use_weights, metrics, residuals){
  #Place holder LM results
  p.lm <- NaN
  sigma.lm <- 0
  results.lm <- NULL
  .GlobalEnv$to_model_gene <- to_model_gene
  
  #Fit model
  if(use_weights){
    fit.lm <- stats::lm(model_lm, data=to_model_gene, weights=to_model_gene$gene_weight)
  } else{
    fit.lm <- stats::lm(model_lm, data=to_model_gene, weights=NULL)
  }
  
  p.lm <- broom::tidy(fit.lm)
  if(residuals){residuals.lm<-residuals(fit.lm)} else {residuals.lm=NULL}
  
  #Extract results
  results.lm <- data.frame(
    model = rep("lm", nrow(p.lm)),    #Label model as lm
    gene = rep(gene, nrow(p.lm)),     #gene name
    variable = p.lm$term,             #variables in model
    df = fit.lm$df.residual,          #degrees of freedom
    statistic = p.lm$statistic,       #test statistic
    estimate = p.lm$estimate,         #estimates in model
    pval = p.lm$p.value)              #P-value
  
  #Model fit metrics
  if(metrics){
    fit.metrics <- data.frame(
      model="lm.fit",
      gene=gene,
      sigma = stats::sigma(fit.lm),
      AIC = stats::AIC(fit.lm),
      BIC = stats::BIC(fit.lm),
      Rsq = summary(fit.lm)$r.squared,
      adj_Rsq = summary(fit.lm)$adj.r.squared)
    
    # estimate residuals
    if(residuals){
      fit.residuals<-
        data.frame(gene=residuals(fit.lm))} else {fit.residuals<-NULL}
    
  } else{
    fit.metrics <- NULL
    fit.residuals<-NULL
  }
  
  results.lm.ls <- list()
  results.lm.ls[["fit"]] <- fit.lm
  results.lm.ls[["metrics"]] <- fit.metrics
  results.lm.ls[["results"]] <- results.lm
  results.lm.ls[["residuals"]] <- fit.residuals
  
  return(results.lm.ls)
}



#'##############################
### 3) kimma_lme_resid()   ####
#'##############################

kimma_lme_resid <- function(model_lme, to_model_gene, gene, use_weights, metrics, residuals){
  rowname <- NULL
  #Place holder LME results
  p.lme <- NaN
  sigma.lme <- 0
  results.lme <- NULL
  .GlobalEnv$to_model_gene <- to_model_gene
  
  #Fit LME model
  if(use_weights){
    fit.lme <- lme4::lmer(model_lme, data=to_model_gene,
                          weights=to_model_gene$gene_weight)
  } else{
    fit.lme <- lme4::lmer(model_lme, data=to_model_gene, weights=NULL)
  }
  
  #Estimate P-value
  p.lme <- broom::tidy(car::Anova(fit.lme))
  p.rand <- as.data.frame(lmerTest::rand(fit.lme)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname != "<none>")
  #Get estimates
  est.lme <- as.data.frame(stats::coef(summary(fit.lme))) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname != "(Intercept)")
  

  #Extract results
  #If no 3+ level variables
  if(nrow(p.lme)==nrow(est.lme)){
    results.lme <- data.frame(
      model = "lme",                                  #Label model as lme
      gene = gene,                                    #gene name
      variable = c(p.lme$term, p.rand$rowname),       #variables in model
      statistic = c(p.lme$statistic, rep(NA, nrow(p.rand))),  #test statistic
      #df = NA,                                        #degrees of freedom
      estimate = c(est.lme$Estimate, p.rand$LRT),     #estimate in model
      pval = c(p.lme$p.value, p.rand$`Pr(>Chisq)`))   #P-value
  } else{
    #If 3+ variable
    results.lme <- data.frame(
      model = "lme",                                #Label model as lme
      gene = gene,                                  #gene name
      variable = c(p.lme$term, p.rand$rowname),     #variables in model
      statistic = "seeContrasts",                   #test statistic
      #df = "seeContrasts",                          #degrees of freedom
      estimate = "seeContrasts",                    #estimate in model
      pval = c(p.lme$p.value, p.rand$`Pr(>Chisq)`)) #P-value
  }
  
  if(metrics){
    #Calculate R-squared
    if(use_weights){
      null <- stats::glm(formula = expression ~ 1, data = to_model_gene,
                         weights = to_model_gene$gene_weight)
    } else{
      null <- stats::glm(formula = expression ~ 1, data = to_model_gene)
    }
    
    L0 <- as.vector(stats::logLik(null))
    L1 <- as.vector(stats::logLik(fit.lme))
    n <- stats::nobs(fit.lme)
    ret <- 1 - exp(-2 / n * (L1 - L0))
    max.r2 <- 1 - exp(2 / n * L0)
    
    #Model fit metrics
    fit.metrics <- data.frame(
      model="lme.fit",
      gene=gene,
      sigma = stats::sigma(fit.lme),
      AIC = stats::AIC(fit.lme),
      BIC = stats::BIC(fit.lme),
      Rsq = ret,
      adj_Rsq = ret / max.r2
    )
    
    # estimate residuals
    if(residuals){
      fit.residuals<-
        data.frame(gene=residuals(fit.lme))} else {fit.residuals<-NULL}
    
  } else{
    fit.metrics <- NULL
    fit.residuals<-NULL
  }
  
  results.lme.ls <- list()
  results.lme.ls[["fit"]] <- fit.lme
  results.lme.ls[["metrics"]] <- fit.metrics
  results.lme.ls[["results"]] <- results.lme
  results.lme.ls[["residuals"]] <- fit.residuals
  return(results.lme.ls)
  
}

#'##############################
### 4) kimma_lmerel_resid()   ##
#'##############################


kimma_lmerel_resid <- function(model_lme, to_model_gene, gene, kin_subset, use_weights,
                         patientID, metrics, residuals){
  rowname <- NULL
  #Place holder LME results
  p.kin <- NaN
  sigma.kin <- 0
  results.kin <- NULL
  .GlobalEnv$to_model_gene <- to_model_gene
  
  #Format kinship matrix to list
  kin.ls <- list()
  kin.ls[[patientID]] <- as.matrix(kin_subset)
  #Fit kin model
  if(use_weights){
    fit.kin <- lme4qtl::relmatLmer(model_lme, data=to_model_gene,
                                   relmat = kin.ls,
                                   weights=to_model_gene$gene_weight)
  } else{
    fit.kin <- lme4qtl::relmatLmer(model_lme, data=to_model_gene,
                                   relmat = kin.ls)
  }
  
  #Estimate P-value
  p.kin <- broom::tidy(car::Anova(fit.kin))
  p.rand <- as.data.frame(lmerTest::rand(fit.kin)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname != "<none>")
  #Get estimates
  est.kin <- as.data.frame(stats::coef(summary(fit.kin))) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname != "(Intercept)")
  

  #Extract results
  #If no 3+ level variables
  if(nrow(p.kin)==nrow(est.kin)){
    results.kin <- data.frame(
      model = "lmerel",                               #Label model as lme
      gene = gene,                                    #gene name
      variable = c(p.kin$term, p.rand$rowname),       #variables in model
      estimate = c(est.kin$Estimate, p.rand$LRT),     #estimate in model
      statistic = c(p.kin$statistic, rep(NA, nrow(p.rand))),  #test statistic
      #df = NA,                                        #degrees of freedom
      pval = c(p.kin$p.value, p.rand$`Pr(>Chisq)`))   #P-value
  } else{
    #If 3+ variable
    results.kin <- data.frame(
      model = "lmerel",                             #Label model as lme
      gene = gene,                                  #gene name
      variable = c(p.kin$term, p.rand$rowname),     #variables in model
      statistic = "seeContrasts",                   #test statistic
      #df = "seeContrasts",                          #degrees of freedom
      estimate = "seeContrasts",                    #estimate in model
      pval = c(p.kin$p.value, p.rand$`Pr(>Chisq)`)) #P-value
  }
  
  if(metrics){
    #Calculate R-squared
    if(use_weights){
      null <- stats::glm(formula = expression ~ 1, data = to_model_gene,
                         weights = to_model_gene$gene_weight)
    } else{
      null <- stats::glm(formula = expression ~ 1, data = to_model_gene)
    }
    
    L0 <- as.vector(stats::logLik(null))
    L1 <- as.vector(stats::logLik(fit.kin))
    n <- stats::nobs(fit.kin)
    ret <- 1 - exp(-2 / n * (L1 - L0))
    max.r2 <- 1 - exp(2 / n * L0)
    
    #Model fit metrics
    fit.metrics <- data.frame(
      model="lmerel.fit",
      gene=gene,
      sigma = stats::sigma(fit.kin),
      AIC = stats::AIC(fit.kin),
      BIC = stats::BIC(fit.kin),
      Rsq = ret,
      adj_Rsq = ret / max.r2)
    
    # estimate residuals
    if(residuals){
      fit.residuals<-
        data.frame(gene=residuals(fit.kin))} else {fit.residuals<-NULL}
    
  } else{
    fit.metrics <- NULL
    fit.residuals<-NULL
  }
  
  results.kin.ls <- list()
  results.kin.ls[["fit"]] <- fit.kin
  results.kin.ls[["metrics"]] <- fit.metrics
  results.kin.ls[["results"]] <- results.kin
  results.kin.ls[["residuals"]] <- fit.residuals
  return(results.kin.ls)
  
}

#'##############################
### 5) kmFit_contrast_resid()   ####
#'##############################


kmFit_contrast_resid<-
  function (fit, contrast_var, to_model_gene, genotype_name) 
  {
    contrast.i <- term <- p.value <- contrast_ref <- contrast_lvl <- contrast <- null.value <- estimate <- NULL
    contrast.result <- data.frame()
    for (contrast.i in contrast_var) {
      contrast.result.temp <- NULL
      i.split <- strsplit(contrast.i, split = ":")[[1]]
      contrast.is.numeric <- unlist(lapply(to_model_gene[, 
                                                         i.split], is.numeric))
      if (any(contrast.is.numeric)) {
        contrast.result.temp <- tryCatch({
          emmeans::emtrends(fit, adjust = "none", var = i.split[contrast.is.numeric], 
                            stats::as.formula(paste("pairwise~", i.split[!contrast.is.numeric], 
                                                    sep = "")))$contrasts %>% broom::tidy() %>% 
            dplyr::mutate(term = gsub(":", "*", contrast.i)) %>% 
            tidyr::separate(contrast, into = c("contrast_ref", 
                                               "contrast_lvl"), sep = " - ")
        }, error = function(e) {
          return(NULL)
        })
        if (is.data.frame(contrast.result.temp)) {
          contrast.result <- contrast.result.temp %>% 
            dplyr::bind_rows(contrast.result)
        }
      }
      else {
        contrast.result.temp <- tryCatch({
          emmeans::emmeans(fit, adjust = "none", stats::as.formula(paste("pairwise~", 
                                                                         contrast.i, sep = "")))$contrasts %>% broom::tidy() %>% 
            dplyr::mutate(term = gsub(":", "*", contrast.i)) %>% 
            tidyr::separate(contrast, into = c("contrast_ref", 
                                               "contrast_lvl"), sep = " - ")
        }, error = function(e) {
          return(NULL)
        })
        if (is.data.frame(contrast.result.temp)) {
          if (!is.null(genotype_name)) {
            if (grepl(genotype_name, contrast.i)) {
              contrast.result.temp <- contrast.result.temp %>% 
                dplyr::mutate(dplyr::across(c(contrast_ref, 
                                              contrast_lvl), ~gsub(genotype_name, 
                                                                   paste0(genotype_name, "_"), .)))
            }
          }
          contrast.result <- contrast.result.temp %>% 
            dplyr::bind_rows(contrast.result)
        }
      }
    }
    contrast.result.format <- contrast.result %>% dplyr::rename(variable = term, 
                                                                pval = p.value) %>% dplyr::select(-null.value) %>% dplyr::mutate(estimate = -estimate)
    return(contrast.result.format)
  }