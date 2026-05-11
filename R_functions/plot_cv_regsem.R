#######################################################################################
# visualize parameter-tuning / feature-selection for regularized SEM models in regsem
#######################################################################################

# Author: Max Segnitz, msegnitz@uw.edu
# Started February 2026
#
# © Richard M Segnitz 2026
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#
# DESCRIPTION:
# Takes as input an output obect from cv_regsem() or par_cv_regsem(). Visualizes parameter tuning effects on 
# feature shrinkage / exclusion.

################################################################################################################
# ============================================================
#  plot_cv_regsem()
#  Visualize lasso regularization path from a cv_regsem() fit
# ============================================================
#
#  Arguments:
#    fit_cv        : output object from regsem::cv_regsem()
#    treatment_var : character; name of the treatment variable as it
#                    appears in SEM parameter labels (e.g. "treatment")
#    outcome_var   : character; name of the outcome variable as it
#                    appears in SEM parameter labels (e.g. "TNSS_AUC")
#    mediator_suffix : character; suffix that identifies mediator
#                    variable names in the path labels (e.g. "modSEM")
#    z_suffix      : character; suffix appended to mediator names in
#                    the standardized parameterization (e.g. "_Z").
#                    Stripped for display. Set to "" if not applicable.
#    mediator_levels : optional character vector; factor level order
#                    for mediator names in plots. If NULL, order is
#                    derived from data.
#    lambda_max    : numeric; upper bound of lambda range to display.
#                    Default 0.15.
#    exclude_nonconv: logical; if TRUE, lambda values where conv == 1
#                    (non-convergent solutions) are excluded from all
#                    plots. Default TRUE.
#    lambda_min_sel: numeric; the selected/optimal lambda value,
#                    used to draw a reference band on plots. If NULL,
#                    no band is drawn.
#    band_halfwidth: numeric; half-width of the shaded lambda reference
#                    band. Default 0.0025.
#    interp_n      : integer; number of interpolation points for
#                    smoothing the BIC-coloured retention path.
#                    Default 100.
#    tag_subtitle  : character; subtitle for the assembled figure
#                    (e.g. "Dual Therapy"). Default "".
#
#  Returns:
#    A named list:
#      $reg_plot            : regularization path facet plot
#      $retention_grid      : per-mediator tile × lambda dot plot
#      $retention_path      : smoothed retention count × lambda plot
#      $combined            : patchwork-assembled multi-panel figure

plot_cv_regsem <- function(
    fit_cv,
    treatment_var    = "treatment",
    outcome_var      = "TNSS_AUC",
    mediator_suffix  = "modSEM",
    z_suffix         = "_Z",
    mediator_levels  = NULL,
    exclude_nonconv  = TRUE,
    lambda_max       = 0.15,
    lambda_min_sel   = NULL,
    band_halfwidth   = 0.0025,
    interp_n         = 100,
    tag_subtitle     = ""
) {
  
  # ---- 0. Dependencies -------------------------------------------------------
  requireNamespace("dplyr",    quietly = TRUE)
  requireNamespace("tidyr",    quietly = TRUE)
  requireNamespace("ggplot2",  quietly = TRUE)
  requireNamespace("patchwork",quietly = TRUE)
  requireNamespace("stringr",  quietly = TRUE)
  requireNamespace("tibble",   quietly = TRUE)
  requireNamespace("readr",    quietly = TRUE)
  
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(stringr); library(tibble); library(readr)
  
  # ---- 1. Identify non-convergent lambda values ------------------------------
  
  fits_df <- 
    as.data.frame(fit_cv$fits)%>%
    mutate(lambda = round(lambda, digits = 5))
  
  
  if (exclude_nonconv && "conv" %in% colnames(fits_df)) {
    nonconv_lambdas <- fits_df$lambda[fits_df$conv == 1]
    if (length(nonconv_lambdas) > 0) {
      message(sprintf(
        "plot_cv_regsem: excluding %d non-convergent solution(s) (conv == 1) at lambda = %s.",
        length(nonconv_lambdas),
        paste(round(nonconv_lambdas, 4), collapse = ", ")
      ))
    }
    fits_df <- filter(fits_df, conv != 1)
  }
  
  conv_lambdas <- fits_df$lambda
  
  # Function to count decimal places
  lambda_jump <- signif(fit_cv$call$jump, digits = 5)
  
  # ---- 2. Extract & tidy parameter estimates ---------------------------------
  
  # Build path strings for the two path types we care about:
  #   (a) treatment  -> mediator  (a-paths)
  #   (b) mediator_Z -> outcome   (b-paths)
  a_path_pattern <- paste0(treatment_var, " -> ", mediator_suffix)
  b_path_pattern <- paste0(z_suffix, " -> ", outcome_var)
  
  params_df <- as.data.frame(fit_cv$parameters) %>%
    dplyr::select(contains("->")) %>%
    rownames_to_column("lambda_round") %>%
    mutate(
      lambda = round(seq(from = min(fits_df$lambda), by = lambda_jump, length.out = n()),
                     digits = 5))%>%
    pivot_longer(
      cols      = contains("->"),
      names_to  = "parameter",
      values_to = "estimate"
    ) %>%
    select(lambda_round, lambda, parameter, estimate) %>%
    # Keep only a-paths and b-paths involving the mediator set
    filter(
      grepl(a_path_pattern, parameter, fixed = TRUE) |
        grepl(b_path_pattern, parameter, fixed = TRUE)
    ) %>%
    filter(grepl(mediator_suffix, parameter, fixed = TRUE))%>%
    filter(lambda %in% fits_df$lambda)
  
  # ---- 2. Derive mediator label from parameter string -----------------------
  # a-path label lives in the RHS  ("treatment -> modSEM_X_Z" → "modSEM_X_Z")
  # b-path label lives in the LHS  ("modSEM_X_Z -> TNSS_AUC" → "modSEM_X_Z")
  
  params_df <- params_df %>%
    rowwise() %>%
    mutate(
      mediator = case_when(
        grepl(outcome_var, parameter)  ~
          str_split(parameter, " -> ")[[1]][1],          # b-path: take LHS
        grepl(treatment_var, parameter) ~
          str_split(parameter, " -> ")[[1]][2]           # a-path: take RHS
      )
    ) %>%
    ungroup() %>%
    mutate(
      # Strip standardization suffix and convert to factor
      mediator = str_remove(mediator, fixed(z_suffix)),
      mediator = if (!is.null(mediator_levels)) {
        factor(mediator, levels = mediator_levels)
      } else {
        factor(mediator)
      },
      parameter_class = case_when(
        grepl(outcome_var,   parameter) ~ "mediator_outcome",   # b-path
        grepl(treatment_var, parameter) ~ "treatment_mediator"  # a-path
      )
    ) %>%
    select(-parameter)
  
  # ---- 3. Compute indirect effects and join model fit indices ----------------
  
  path_df <- params_df %>%
    group_by(lambda, mediator) %>%
    pivot_wider(
      names_from  = parameter_class,
      values_from = estimate
    ) %>%
    ungroup() %>%
    mutate(
      # Specific indirect effect = a-path × b-path
      ind_effect_coef = mediator_outcome * treatment_mediator
    ) %>%
    left_join(fits_df, by = "lambda") %>%
    mutate(dBIC = BIC - min(fits_df$BIC, na.rm = TRUE)) %>%
    filter(lambda <= lambda_max) %>%
    { if (exclude_nonconv && "conv" %in% colnames(fits_df))
      filter(., lambda %in% conv_lambdas) else . }
  
  # ---- 4. Regularization path plot ------------------------------------------
  # Three facets: b-path, a-path, indirect effect
  
  reg_plot <- path_df %>%
    pivot_longer(
      cols      = c(mediator_outcome, treatment_mediator, ind_effect_coef),
      names_to  = "parameter",
      values_to = "estimate"
    ) %>%
    mutate(
      parameter = recode(parameter,
                         "ind_effect_coef"    = "Specific Indirect Effect",
                         "mediator_outcome"   = paste0("Mediator Effect on ", outcome_var),
                         "treatment_mediator" = paste0(str_to_title(treatment_var), " Effect on Mediator")
      ),
      parameter = factor(parameter, levels = c(
        paste0("Mediator Effect on ", outcome_var),
        paste0(str_to_title(treatment_var), " Effect on Mediator"),
        "Specific Indirect Effect"
      ))
    ) %>%
    ggplot(aes(x = lambda, y = estimate)) +
    geom_hline(yintercept = 0, color = "grey50") +
    { if (!is.null(lambda_min_sel))
      geom_vline(xintercept = lambda_min_sel, linetype = "dashed") } +
    geom_path(aes(group = mediator, color = mediator), lwd = 1, alpha = 0.75) +
    labs(
      x     = "Lasso Penalty (\u03BB)",
      y     = "Parameter Estimate",
      color = "Mediator"
    ) +
    facet_wrap(~parameter, nrow = 1) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"))
  
  # ---- 5. Per-mediator feature retention tile + dot plot --------------------
  # Mediators ordered by how long they are retained across the lambda range
  
  retention_df <- path_df %>%
    mutate(
      retention = ifelse(abs(ind_effect_coef) > 0, "retained", "removed")
    ) %>%
    group_by(mediator) %>%
    mutate(retention_total = sum(retention == "retained")) %>%
    ungroup() %>%
    group_by(lambda) %>%
    mutate(n_i = sum(retention == "retained")) %>%
    ungroup()%>%
    arrange(retention_total, mediator) %>%
    mutate(mediator = factor(mediator, levels = unique(as.character(mediator))))
  
  retention_grid <- {
    p <- ggplot(retention_df, aes(y = mediator, x = lambda))
    
    if (!is.null(lambda_min_sel)) {
      p <- p + annotate("rect",
                        xmin = lambda_min_sel - band_halfwidth,
                        xmax = lambda_min_sel + band_halfwidth,
                        ymin = -Inf, ymax = Inf, fill = "lightgrey"
      )
    }
    
    p +
      geom_tile(color = "grey70", fill = NA) +
      geom_point(
        data = filter(retention_df, retention == "retained"),
        aes(size = abs(ind_effect_coef))
      ) +
      scale_size_continuous(range = c(0, 4)) +
      labs(
        x        = "\u03BB",
        y        = "",
        size     = "Indirect Effect\nEst."
      ) +
      theme_bw() +
      theme(
        aspect.ratio = nlevels(retention_df$mediator) / (lambda_max / lambda_jump),
        panel.grid   = element_blank()
      )
  }
  
  # ---- 6. Smoothed mediator retention count × lambda (interpolated BIC) -----
  
  retention_tally <- path_df %>%
    mutate(retention = ifelse(abs(ind_effect_coef) > 0, "retained", "removed")) %>%
    group_by(lambda, retention) %>%
    reframe(n = n(), BIC = unique(BIC)) %>%
    mutate(dBIC = BIC - min(fits_df$BIC, na.rm = TRUE)) %>%
    filter(retention == "retained")
  
  # Linear interpolation for smooth colour gradient along path
  interp_bic <- approx(x = retention_tally$lambda, y = retention_tally$dBIC, n = interp_n)
  interp_n_  <- approx(x = retention_tally$lambda, y = retention_tally$n,    n = interp_n)
  
  interp_df <- data.frame(
    lambda_i = interp_bic$x,
    dBIC_i   = interp_bic$y,
    n_i      = interp_n_$y
  )
  
  n_mediators <- nlevels(path_df$mediator)
  
  retention_path <- {
    p <- ggplot(interp_df, aes(x = lambda_i, y = n_i))
    
    if (!is.null(lambda_min_sel)) {
      p <- p + annotate("rect",
                        xmin = lambda_min_sel - band_halfwidth,
                        xmax = lambda_min_sel + band_halfwidth,
                        ymin = -Inf, ymax = Inf, fill = "lightgrey"
      )
    }
    
    p +
      geom_path(aes(color = dBIC_i), lwd = 2, lineend = "round") +
      scale_y_continuous(breaks = seq(0, n_mediators, 1)) +
      scale_x_continuous(
        minor_breaks = retention_tally$lambda,
        breaks       = seq(0, lambda_max, 0.05),
        labels       = seq(0, lambda_max, 0.05)
      ) +
      scale_color_viridis_c() +
      labs(
        x     = "\u03BB",
        y     = "Total Mediators Retained",
        color = "\u0394 BIC\nFrom Best-Fit"
      ) +
      theme_bw() +
      theme(panel.grid.minor.y = element_blank())
  }
  
  # ---- 7. Assemble multi-panel figure ----------------------------------------
  
  combined <- (
    retention_grid + labs(subtitle = "") + coord_cartesian(xlim = c(-0.001, lambda_max))
  ) / (
    retention_path + coord_cartesian(xlim = c(-0.001, lambda_max))
  ) +
    plot_layout(ncol = 1, nrow = 2, heights = c(1, 1)) +
    plot_annotation(
      tag_levels = "a",
      subtitle   = tag_subtitle,
      caption    = paste0(
        "a) Feature retention by lasso penalty.\n",
        "b) Model fit and total retained mediators by lasso penalty.\n",
        if (!is.null(lambda_min_sel)) "Shaded grey bar indicates selected lambda value." else ""
      ),
      theme = theme(plot.caption = element_text(hjust = 0))
    )
  
  # ---- 8. Return -------------------------------------------------------------
  
  invisible(list(
    reg_plot       = reg_plot,
    retention_grid = retention_grid,
    retention_path = retention_path,
    combined       = combined
  ))
}