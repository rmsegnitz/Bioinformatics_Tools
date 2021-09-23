## SCRIPT IN DEVELOPMENT ##

# Calculate and Visualize Johnson Neyman Intervals in Simple Slopes Analysis
# Script is adapted from the {interactions} package, written by Jacob Long (2019)
#   with changes faciliating variance estimation from categorical predictors.

# Author: Max Segnitz, msegnitz@uw.edu
# Started September 2021
#
# © Richard M Segnitz 2021
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#
# DESCRIPTION:
# Runs Johnson-Neyman interval analysis on model fit with interaction on main effect of interest. 
# Determines the interval of moderator x over which predictor has significant effect on outcome.

############
interactionJN<-
function (model, pred, modx, vmat = NULL, alpha = 0.05, plot = TRUE, 
          control.fdr = FALSE, line.thickness = 0.5, df = "residual", 
          digits = getOption("jtools-digits", 2), critical.t = NULL, 
          sig.color = "#00BFC4", insig.color = "#F8766D", mod.range = NULL, 
          title = "Johnson-Neyman plot") 
{
  pred <- quo_name(enexpr(pred))
  if (make.names(pred) != pred) 
    pred <- bt(pred)
  modx <- quo_name(enexpr(modx))
  if (make.names(modx) != modx) 
    modx <- bt(modx)
  if (modx == "NULL") {
    modx <- NULL
  }
  if (df == "residual") {
    df <- df.residual(model)
    if (is.null(df)) {
      warn_wrap("Tried to calculate residual degrees of freedom but result was\n                NULL. Using normal approximation instead.")
      df <- Inf
    }
  }
  else if (df == "normal") {
    df <- Inf
  }
  else if (!is.numeric(df)) {
    stop("df argument must be 'residual', 'normal', or a number.")
  }
  out <- list()
  out <- structure(out, pred = gsub("`", "", pred), modx = gsub("`", "", modx), 
                   alpha = alpha, plot = plot, digits = digits, control.fdr = control.fdr)
  get_coef <- function(mod) {
    if (inherits(mod, "merMod") | inherits(mod, "brmsfit")) {
      coef <- lme4::fixef(model)
      if (inherits(mod, "brmsfit")) {
        coefs <- coef[, 1, drop = TRUE]
        names(coefs) <- rownames(coef)
        coef <- coefs
      }
      return(coef)
    }
    else {
      coef(mod)
    }
  }
  
  # Check nature of predictor
  pred_class<-
    case_when(model.frame(model)[, pred]%>%is.numeric() ~ "numeric",
              model.frame(model)[, pred]%>%is.character() ~ "categorical", 
              model.frame(model)[, pred]%>%is.factor() ~ "categorical")
  
  # Define coef term if categorical
  if(pred_class=="categorical"){
    pred_term = names(as.list(get_coef(model)))[grep(pred, names(as.list(get_coef(model))))[1]]}
  
  # Check whether any categorical predictors are limited to 2 groups
  if(pred_class=="categorical"){
    if(model.frame(model)[, pred]%>%unique()%>%length()>2){
      stop("Categorical predictor with greater than 2 groups present. Predictors with groups >2 are not currently supported.")
    }
  }
  
  # Edit: find interaction term for categorical predictor

  intterm1 = case_when(pred_class=="numeric" ~ paste(pred, ":", modx, sep = ""),
                       pred_class=="categorical" ~ paste(names(as.list(get_coef(model)))[grep(pred, names(as.list(get_coef(model))))][1], 
                                                         modx, sep=":"))
  
  intterm1tf <- any(intterm1 %in% names(get_coef(model)))
  
  # Edit: find interaction term for categorical predictor
  
  intterm2 = case_when(pred_class == "numeric" ~ paste(modx, ":", pred, sep = ""),
                       pred_class == "categorical" ~ paste(modx, names(as.list(get_coef(model)))[grep(pred, names(as.list(get_coef(model))))][1] 
                                                           , sep=":"))

  intterm2tf <- any(intterm2 %in% names(get_coef(model)))
  
  
  coefs <- get_coef(model)
  inttermstf <- c(intterm1tf, intterm2tf)
  intterms <- c(intterm1, intterm2)
  intterm <- intterms[which(inttermstf)]
  modrange <- range(model.frame(model)[, modx])
  modrangeo <- range(model.frame(model)[, modx])
  modsd <- sd(model.frame(model)[, modx])
  if (is.null(mod.range)) {
    modrange[1] <- modrange[1] - modsd
    modrange[2] <- modrange[2] + modsd
  }
  else {
    modrange <- mod.range
  }
  if (modrange[1] >= modrangeo[1] & modrange[2] <= modrangeo[2]) {
    no_range_line <- TRUE
  }
  else {
    no_range_line <- FALSE
  }
  alpha <- alpha/2
  if (control.fdr == FALSE) {
    tcrit <- qt(alpha, df = df)
    tcrit <- abs(tcrit)
  }
  else if (control.fdr == TRUE) {
    predb <- case_when(pred_class == "numeric" ~ coefs[pred],
                       pred_class == "categorical" ~ coefs[pred_term])
                       
    intb <- coefs[intterm]
    vcovs <- vcov(model)
    vcov_pred <- ifelse(pred_class == "numeric", vcovs[pred, pred],
                   ifelse(pred_class == "categorical", vcovs[pred_term, pred_term], NA))
    
    vcov_int <- vcovs[intterm, intterm]
    
    vcov_pred_int <- ifelse(pred_class=="numeric", vcovs[pred, intterm],
                            ifelse(pred_class == "categorical", vcovs[pred_term, intterm], NA))
                            
    range_sequence <- seq(from = modrangeo[1], to = modrangeo[2], 
                          by = (modrangeo[2] - modrangeo[1])/1000)
    
    marginal_effects <- predb + intb * range_sequence
    me_ses <- sqrt(vcov_pred + (range_sequence^2) * vcov_int + 
                     2 * range_sequence * vcov_pred_int)
    ts <- marginal_effects/me_ses
    df <- df.residual(model)
    ps <- 2 * pmin(pt(ts, df = df), (1 - pt(ts, df = df)))
    multipliers <- seq_along(marginal_effects)/length(marginal_effects)
    ps_o <- order(ps)
    test <- 0
    i <- 1 + length(marginal_effects)
    while (test == 0 & i > 1) {
      i <- i - 1
      test <- min(ps[ps_o][1:i] <= multipliers[i] * (alpha *2))
    }
    tcrit <- abs(qt(multipliers[i] * alpha, df))
  }
  if (is.null(vmat)) {
    vmat <- vcov(model)
    covy3 <- vmat[intterm, intterm]
    covy1 <- ifelse(pred_class=="numeric" , vmat[pred, pred],
                    ifelse(pred_class=="categorical", vmat[pred_term, pred_term], NA))
    covy1y3 <- ifelse(pred_class=="numeric" , vmat[intterm, pred],
                      ifelse(pred_class=="categorical",vmat[intterm, pred_term], NA))
    y3 <- coefs[intterm]
    y1 <- ifelse(pred_class=="numeric", coefs[pred], ifelse(pred_class=="categorical", coefs[pred_term], NA))
  }
  else {
    covy3 <- vmat[intterm, intterm]
    covy1 <- ifelse(pred_class=="numeric" , vmat[pred, pred],
                    ifelse(pred_class=="categorical", vmat[pred_term, pred_term], NA))
    covy1y3 <- ifelse(pred_class=="numeric" , vmat[intterm, pred],
                      ifelse(pred_class=="categorical",vmat[intterm, pred_term], NA))
    
    y3 <- get_coef(model)[intterm]
    y1 <- ifelse(pred_class=="numeric", get_coef(model)[pred],
                 ifelse(pred_class=="categorical", get_coef(model)[pred_term], NA))
  }
  if (!is.null(critical.t)) {
    tcrit <- critical.t
  }
  a <- tcrit^2 * covy3 - y3^2
  b <- 2 * (tcrit^2 * covy1y3 - y1 * y3)
  c <- tcrit^2 * covy1 - y1^2
  discriminant <- function(a, b, c) {
    disc <- b^2 - 4 * a * c
    if (disc > 0) {
      out <- disc
    }
    else if (disc == 0) {
      return(NULL)
    }
    else {
      return(NULL)
    }
    return(out)
  }
  disc <- discriminant(a, b, c)
  if (is.null(disc)) {
    failed <- TRUE
  }
  else {
    failed <- FALSE
  }
  quadsolve <- function(a, b, c, disc) {
    x1 <- (-b + sqrt(disc))/(2 * a)
    x2 <- (-b - sqrt(disc))/(2 * a)
    result <- c(x1, x2)
    result <- sort(result, decreasing = FALSE)
    return(result)
  }
  if (!is.null(disc)) {
    bounds <- quadsolve(a, b, c, disc)
  }
  else {
    bounds <- c(-Inf, Inf)
  }
  names(bounds) <- c("Lower", "Higher")
  cbands <- function(x2, y1, y3, covy1, covy3, covy1y3, tcrit, 
                     predl, modx) {
    upper <- c()
    slopes <- c()
    lower <- c()
    slopesf <- function(i) {
      s <- y1 + y3 * i
      return(s)
    }
    upperf <- function(i, s) {
      u <- s + tcrit * sqrt((covy1 + 2 * i * covy1y3 + 
                               i^2 * covy3))
      return(u)
    }
    lowerf <- function(i, s) {
      l <- s - tcrit * sqrt((covy1 + 2 * i * covy1y3 + 
                               i^2 * covy3))
      return(l)
    }
    slopes <- sapply(x2, slopesf, simplify = "vector", USE.NAMES = FALSE)
    upper <- mapply(upperf, x2, slopes)
    lower <- mapply(lowerf, x2, slopes)
    out <- matrix(c(x2, slopes, lower, upper), ncol = 4)
    colnames(out) <- c(modx, predl, "Lower", "Upper")
    out <- as.data.frame(out)
    return(out)
  }
  x2 <- seq(from = modrange[1], to = modrange[2], length.out = 1000)
  predl <- paste("Slope of", pred)
  cbs <- cbands(x2, y1, y3, covy1, covy3, covy1y3, tcrit, 
                predl, modx)
  out$bounds <- bounds
  out <- structure(out, modrange = modrangeo)
  sigs <- which((cbs$Lower < 0 & cbs$Upper < 0) | (cbs$Lower > 
                                                     0 & cbs$Upper > 0))
  insigs <- setdiff(1:1000, sigs)
  cbs$Significance <- rep(NA, nrow(cbs))
  cbs$Significance <- factor(cbs$Significance, levels = c("Insignificant", 
                                                          "Significant"))
  index <- 1:1000 %in% insigs
  cbs$Significance[index] <- "Insignificant"
  index <- 1:1000 %in% sigs
  cbs$Significance[index] <- "Significant"
  out$cbands <- cbs
  index <- which(cbs$Significance == "Significant")[1]
  if (!is.na(index) & index != 0) {
    inside <- (cbs[index, modx] > bounds[1] && cbs[index, 
                                                   modx] < bounds[2])
    all_sig <- NULL
    if (is.na(which(cbs$Significance == "Insignificant")[1])) {
      all_sig <- TRUE
    }
    else {
      all_sig <- FALSE
    }
  }
  else {
    inside <- FALSE
    all_sig <- TRUE
  }
  out <- structure(out, inside = inside, failed = failed, 
                   all_sig = all_sig)
  cbso1 <- cbs[cbs[, modx] < bounds[1], ]
  cbso2 <- cbs[cbs[, modx] > bounds[2], ]
  cbsi <- cbs[(cbs[, modx] > bounds[1] & cbs[, modx] < bounds[2]), 
  ]
  alpha <- alpha * 2
  alpha <- gsub("0\\.", "\\.", as.character(alpha))
  pmsg <- paste("p <", alpha)
  plot <- ggplot2::ggplot() + ggplot2::geom_path(data = cbso1, 
                                                 ggplot2::aes(x = cbso1[, modx], y = cbso1[, predl], 
                                                              color = cbso1[, "Significance"]), size = line.thickness) + 
    ggplot2::geom_path(data = cbsi, ggplot2::aes(x = cbsi[, 
                                                          modx], y = cbsi[, predl], color = cbsi[, "Significance"]), 
                       size = line.thickness) + ggplot2::geom_path(data = cbso2, 
                                                                   ggplot2::aes(x = cbso2[, modx], y = cbso2[, predl], 
                                                                                color = cbso2[, "Significance"]), size = line.thickness) + 
    ggplot2::geom_ribbon(data = cbso1, ggplot2::aes(x = cbso1[, 
                                                              modx], ymin = cbso1[, "Lower"], ymax = cbso1[, "Upper"], 
                                                    fill = cbso1[, "Significance"]), alpha = 0.2) + 
    ggplot2::geom_ribbon(data = cbsi, ggplot2::aes(x = cbsi[, 
                                                            modx], ymin = cbsi[, "Lower"], ymax = cbsi[, "Upper"], 
                                                   fill = cbsi[, "Significance"]), alpha = 0.2) + ggplot2::geom_ribbon(data = cbso2, 
                                                                                                                       ggplot2::aes(x = cbso2[, modx], ymin = cbso2[, "Lower"], 
                                                                                                                                    ymax = cbso2[, "Upper"], fill = cbso2[, "Significance"]), 
                                                                                                                       alpha = 0.2) + ggplot2::scale_fill_manual(values = c(Significant = sig.color, 
                                                                                                                                                                            Insignificant = insig.color), labels = c("n.s.", pmsg), 
                                                                                                                                                                 breaks = c("Insignificant", "Significant"), drop = FALSE, 
                                                                                                                                                                 guide = ggplot2::guide_legend(order = 2)) + ggplot2::geom_hline(ggplot2::aes(yintercept = 0))
  if (no_range_line == FALSE) {
    plot <- plot + ggplot2::geom_segment(ggplot2::aes(x = modrangeo[1], 
                                                      xend = modrangeo[2], y = 0, yend = 0, linetype = "Range of\nobserved\ndata"), 
                                         lineend = "square", size = 1.25)
  }
  plot <- plot + ggplot2::scale_linetype_discrete(name = " ", 
                                                  guide = ggplot2::guide_legend(order = 1))
  if (out$bounds[1] < modrange[1]) {
  }
  else if (all_sig == FALSE) {
    plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = out$bounds[1]), 
                                       linetype = 2, color = sig.color)
  }
  if (out$bounds[2] > modrange[2]) {
  }
  else if (all_sig == FALSE) {
    plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = out$bounds[2]), 
                                       linetype = 2, color = sig.color)
  }
  plot <- plot + ggplot2::xlim(range(cbs[, modx])) + ggplot2::labs(title = title, 
                                                                   x = modx, y = predl) + ggplot2::scale_color_manual(name = "", 
                                                                                                                      values = c(Significant = sig.color, Insignificant = insig.color), 
                                                                                                                      guide = "none") + theme_apa(legend.pos = "right", legend.font.size = 10) + 
    ggplot2::theme(legend.key.size = ggplot2::unit(1, "lines"))
  out$plot <- plot
  if (control.fdr == TRUE) {
    out$t_value <- tcrit
  }
  
  # summarize results
  cbs_summary<-cbs
  cbs_summary$slope = ifelse(cbs[,2]>=0, "pos", "neg")
  cbs_summary %>% 
    mutate(predictor_significance = ifelse(Significance=="Significant", slope, 
                              "Non-significant"))%>%
    group_by(predictor_significance)%>%
    summarize(moderator_interval = paste(range(get(modx)), collapse = " - "))%>%
    ungroup()%>%
    mutate(test = "Johnson - Neyman", 
           predictor = pred, 
           moderator = modx, 
           alpha = alpha)%>%
    mutate(predictor_slope = ifelse(predictor_significance=="pos",  "+", 
                                       ifelse(predictor_significance=="neg", "-", "ns")),
           predictor_significance = recode(predictor_significance, pos = "Significant", neg = "Significant"))%>%
    dplyr::select(test, predictor, moderator, alpha, 
                  predictor_significance, predictor_slope, moderator_interval)-> cbs_summary
    
  out$summary<-cbs_summary
  
  class(out) <- "johnson_neyman"
  return(out)
  }
