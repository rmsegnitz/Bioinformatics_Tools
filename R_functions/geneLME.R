
########################################################
# Scalable custom gene LMEs with contrast specification
########################################################
# Merged from geneLME_dev.R (2026-02-20):
#   - Added geneLME_contrast_spec() helper
#   - Added Branch A: interaction contrast support via contrast_spec
#   - Added FDR adjustment (fdr_method argument)
#   - Added $contrast_spec in return value
#   - Soft-fail on wrong-length contrasts_secondary
#   - Added contrast_ref / contrast_lvl columns to lme_contrast output
#
# Merged from geneLME_dev2.R (2026-02-20):
#   - Singular fit → model_status flag ("singular_fit") instead of error.
#     Results returned for all genes; filter downstream on model_status.
#     model_status column now present in both lme_anova and lme_contrast.
#   - New geneLME_build_contrast_args() helper: pre-computes Branch A
#     contrasts_list and spec_lookup once before parallel dispatch,
#     eliminating nrow(contrast_spec) × n_genes per-gene rebuilds.
#   - geneLME_fit() signature: contrast_spec replaced by contrasts_list
#     and spec_lookup (Branch A). geneLME_dispatch() updated accordingly.
#   - Benchmarked (geneLME_benchmark2.html): dev2 estimates identical to
#     previous stable (r=1.0, MAD=0); ~1.1-1.2x faster; 100% direction
#     agreement with kimma; equal speed to kimma at 66 contrasts.
########################################################


########################################################
# geneLME_contrast_spec
# Helper: returns a reference template of available contrast levels,
# formatted for use with the contrast_spec argument of geneLME().
#
# Two modes depending on whether contrast_vars contains ":":
#
#   Interaction mode (e.g. contrast_vars = "treatment:visit"):
#     contrast_vars must be a single "var_a:var_b" string.
#     Returns a data.frame with columns contrast_ref and contrast_lvl,
#     one row per pairwise combination of interaction-level strings.
#     Filter this to the contrasts of interest, then pass as
#     contrast_spec to geneLME().
#
#   Single-variable mode (e.g. contrast_vars = c("treatment", "visit")):
#     contrast_vars is a character vector of one or more plain variable
#     names (no ":"). Returns a named list, one element per variable,
#     each a data.frame with a single column 'level' listing the sorted
#     unique values. The message printed for each variable explains how
#     it maps to geneLME() arguments:
#       - The first variable in contrast_vars → contrast_vars[1] in
#         geneLME(); its levels define the length and position order of
#         contrasts_primary vectors.
#       - Additional variables → contrast_vars[2], used as the 'by'
#         grouping variable; filter its levels via contrast_var_2_levels.
#     This list is a reference only — it is not passed to geneLME().
########################################################

geneLME_contrast_spec <- function(targets, contrast_vars) {

  has_interaction <- any(grepl(":", contrast_vars))
  has_plain       <- any(!grepl(":", contrast_vars))

  # Disallow mixing interaction and plain variable names in one call
  if (has_interaction && has_plain) {
    stop(
      "contrast_vars mixes interaction terms (containing ':') and plain variable names.\n",
      "Call geneLME_contrast_spec() separately for interaction and non-interaction variables."
    )
  }

  if (has_interaction) {

    # ---- Interaction mode: must be a single "var_a:var_b" string ----
    if (length(contrast_vars) != 1) {
      stop("Interaction mode requires a single 'var_a:var_b' string in contrast_vars.")
    }

    vars  <- strsplit(contrast_vars, ":")[[1]]
    var_a <- vars[1]
    var_b <- vars[2]

    if (!var_a %in% colnames(targets)) stop(paste0("Variable '", var_a, "' not found in targets."))
    if (!var_b %in% colnames(targets)) stop(paste0("Variable '", var_b, "' not found in targets."))

    # Coerce each variable to a factor, preserving any existing factor level order.
    # If the column is already a factor, levels() preserves the user-defined ordering.
    # If it is a plain character vector, we impose alphabetical order explicitly.
    # This ensures the interaction level string ordering — and therefore which member
    # of each pair lands in contrast_ref vs contrast_lvl — is fully deterministic and
    # consistent between geneLME_contrast_spec() and emmeans' internal ordering.
    fac_a <- if (is.factor(targets[[var_a]])) targets[[var_a]] else factor(targets[[var_a]], levels = sort(unique(targets[[var_a]])))
    fac_b <- if (is.factor(targets[[var_b]])) targets[[var_b]] else factor(targets[[var_b]], levels = sort(unique(targets[[var_b]])))

    # Build all interaction level strings in the same order emmeans will use.
    # interaction() with two ordered factors produces a factor whose levels are
    # var_a[1]:var_b[1], var_a[1]:var_b[2], ..., var_a[n]:var_b[m] — matching
    # emmeans' default grid ordering (vary the rightmost variable fastest).
    ixn_lvls <- levels(interaction(fac_a, fac_b, sep = " "))

    # combn() traverses ixn_lvls in order, so for every pair the first element
    # (contrast_ref) is always at a lower level-index than the second (contrast_lvl).
    # This guarantees: level(var_a in ref) <= level(var_a in lvl) AND
    #                  level(var_b in ref) <= level(var_b in lvl).
    # Users may freely swap ref/lvl in any row after filtering — geneLME() and
    # geneLME_fit() handle either direction correctly (the sign of the contrast
    # estimate will flip, but no error will occur).
    # Row-position indices for contrasts_secondary construction are added by geneLME()
    # to its $contrast_spec output element after the user's filtered spec is received.
    pairs <- combn(ixn_lvls, 2, simplify = FALSE)
    spec  <- data.frame(
      contrast_ref = sapply(pairs, `[[`, 1),
      contrast_lvl = sapply(pairs, `[[`, 2),
      stringsAsFactors = FALSE
    )

    message(
      nrow(spec), " pairwise combinations generated for '", contrast_vars, "'.\n",
      "Filter this data.frame to your contrasts of interest, then pass as contrast_spec to geneLME().\n",
      "geneLME() will attach an indexed copy of contrast_spec to its output ($contrast_spec)\n",
      "showing the row-position index used for contrasts_secondary vector construction."
    )

    return(spec)

  } else {

    # ---- Single/multi-variable mode ----
    # Validate all variables exist
    missing_vars <- setdiff(contrast_vars, colnames(targets))
    if (length(missing_vars) > 0) {
      stop(paste0("Variable(s) not found in targets: ", paste(missing_vars, collapse = ", ")))
    }

    result <- lapply(seq_along(contrast_vars), function(i) {
      v    <- contrast_vars[i]
      lvls <- sort(unique(as.character(targets[[v]])))
      spec <- data.frame(level = lvls, stringsAsFactors = FALSE)

      if (i == 1) {
        message(
          "'", v, "' — primary contrast variable (contrast_vars[1] in geneLME).\n",
          length(lvls), " levels (alphabetical order = position order for contrasts_primary vectors):\n",
          paste(seq_along(lvls), lvls, sep = ". ", collapse = "\n"), "\n",
          "→ contrast vectors passed to contrasts_primary must have length ", length(lvls), ",\n",
          "  with each element weighted by position (e.g. '", lvls[length(lvls)], " vs ", lvls[1],
          "' = c(", paste(c(-1, rep(0, length(lvls) - 2), 1), collapse = ", "), "))."
        )
      } else {
        message(
          "'", v, "' — secondary 'by' grouping variable (contrast_vars[", i, "] in geneLME).\n",
          length(lvls), " levels available:\n",
          paste(seq_along(lvls), lvls, sep = ". ", collapse = "\n"), "\n",
          "→ pass a subset of these levels to contrast_var_2_levels in geneLME() to restrict\n",
          "  which groups the primary contrasts are computed within."
        )
      }

      spec
    })

    names(result) <- contrast_vars
    return(result)
  }
}


########################################################
# geneLME_build_contrast_args
# Pre-computes Branch A contrast structures once, before parallel dispatch.
#
# Called by geneLME() when contrast_vars contains ":".
# Returns a list with two elements that are passed directly to workers:
#
#   $contrasts_list — named list of contrast vectors, one per row of
#     contrast_spec. Each vector has length == number of interaction cells,
#     named by emmeans level string ("TrtA V1", "TrtB V2", ...).
#     contrast_ref cell receives -1; contrast_lvl cell receives +1.
#     Name of each vector is "contrast_lvl - contrast_ref" — the string
#     emmeans will use as the 'contrast' column in its output.
#
#   $spec_lookup — data.frame(contrast, contrast_ref, contrast_lvl) for
#     joining ref/lvl back onto the emmeans first-order result by contrast
#     name. Avoids reconstructing this inside every geneLME_fit() call.
#
# The emmeans level ordering is derived from factor levels in targets (the
# same source emmeans uses internally), guaranteeing index alignment without
# requiring a fitted model object.
########################################################

geneLME_build_contrast_args <- function(targets, contrast_vars, contrast_spec) {

  vars  <- strsplit(contrast_vars, ":")[[1]]
  var_a <- vars[1]
  var_b <- vars[2]

  # Reproduce the exact factor-level ordering that emmeans will use.
  # emmeans builds its reference grid from the factor levels of the model
  # data (targets); using the same coercion logic as geneLME_contrast_spec()
  # guarantees the level strings here match those emmeans returns.
  fac_a <- if (is.factor(targets[[var_a]])) targets[[var_a]] else factor(targets[[var_a]], levels = sort(unique(targets[[var_a]])))
  fac_b <- if (is.factor(targets[[var_b]])) targets[[var_b]] else factor(targets[[var_b]], levels = sort(unique(targets[[var_b]])))

  # Build interaction level strings: "var_a_level var_b_level" in emmeans order.
  # emmeans varies the rightmost variable fastest, matching interaction() default.
  lvls_a    <- levels(fac_a)
  lvls_b    <- levels(fac_b)
  emm_lvls  <- as.vector(outer(lvls_a, lvls_b, paste))   # row = var_a, col = var_b

  # Template zero vector, one element per interaction cell
  default_vec        <- rep(0, length(emm_lvls))
  names(default_vec) <- emm_lvls

  # Build contrasts_list once — identical for all genes
  contrasts_list <- list()
  for (k in seq_len(nrow(contrast_spec))) {
    cv <- default_vec
    cv[contrast_spec$contrast_ref[k]] <- -1
    cv[contrast_spec$contrast_lvl[k]] <-  1
    contrast_name <- paste(contrast_spec$contrast_lvl[k],
                           contrast_spec$contrast_ref[k], sep = " - ")
    contrasts_list[[contrast_name]] <- cv
  }

  # Build spec_lookup once — identical for all genes
  spec_lookup <- data.frame(
    contrast     = paste(contrast_spec$contrast_lvl,
                         contrast_spec$contrast_ref, sep = " - "),
    contrast_ref = contrast_spec$contrast_ref,
    contrast_lvl = contrast_spec$contrast_lvl,
    stringsAsFactors = FALSE
  )

  list(contrasts_list = contrasts_list, spec_lookup = spec_lookup)
}


########################################################
# geneLME_fit
# Core per-gene fitting function. Called inside future_lapply.
# Receives only the minimal pre-extracted data needed for one gene,
# not the full EList object.
#
# Singular fit is no longer an error; model_status = "singular_fit"
# is recorded and results are returned for user-side filtering.
# Branch A receives pre-computed contrasts_list + spec_lookup instead
# of contrast_spec. The per-gene contrast vector loop is eliminated.
########################################################

geneLME_fit <-
  function(gene_name,
           expression_vec,         # named numeric vector: expression values for this gene
           weight_vec,             # named numeric vector or NULL: per-sample weights
           targets,                # data.frame: sample metadata (dat$targets)
           formula_str,            # character: formula RHS e.g. "~ treatment*visit + (1|ptID)"
           run_contrast,
           contrast_vars,
           contrast_var_2_levels,
           contrasts_list,         # Branch A: pre-computed named contrast vector list, or NULL
           spec_lookup,            # Branch A: pre-computed ref/lvl join table, or NULL
           contrasts_primary,      # Branch B: named list of contrast vectors, or NULL
           contrasts_secondary) {  # both branches: second-order contrast vectors, or NULL

    # Load worker-side packages quietly.
    # suppressPackageStartupMessages: silences "Loading required package: Matrix" etc.
    # suppressWarnings: silences "package 'lme4' was built under R version X.Y.Z".
    suppressPackageStartupMessages(suppressWarnings({
      library(lme4)
      library(emmeans)
      library(car)
      library(broom.mixed)
      library(dplyr)
      library(tibble)
    }))

    result <- tryCatch({

      # --- BUILD MODEL DATA ---
      model_data <- targets
      model_data$expression <- expression_vec

      # --- RECONSTRUCT FORMULA LOCALLY ---
      # Built from the raw string so its enclosing environment is this call frame.
      # lmer resolves 'weight_vec' and other names in local scope — no
      # environment stripping needed, no locked-environment errors.
      formula_obj <- as.formula(paste("expression", formula_str))

      # --- FIT MODEL ---
      # check.scaleX = "ignore": silences the "predictor variables on very
      # different scales" warning. This is expected with mixed RNA-seq covariates
      # (e.g. lib.size vs binary treatment indicator) and does not affect model fit.
      lme_ctrl <- lmerControl(calc.derivs = FALSE, check.scaleX = "ignore")
      if (is.null(weight_vec)) {
        lme_i <- lmer(formula_obj, data = model_data, control = lme_ctrl)
      } else {
        lme_i <- lmer(formula_obj, weights = weight_vec, data = model_data, control = lme_ctrl)
      }

      # --- SINGULAR FIT: flag, do not stop ---
      # isSingular() indicates the random effect variance hit its boundary (zero).
      # Fixed effect estimates and contrasts are still numerically valid; the user
      # should decide whether to retain or filter these genes downstream.
      # model_status records the flag; "success" indicates a clean fit.
      model_status <- if (isSingular(lme_i)) "singular_fit" else "success"

      # --- EXTRACT: AIC ---
      aic_res <- data.frame(gene = gene_name, AIC = AIC(lme_i))

      # --- EXTRACT: ANOVA & coefficient summary ---
      lme_i_anova   <- car::Anova(lme_i) %>% broom.mixed::tidy()
      lme_i_summary <- summary(lme_i)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("variable")

      # --- BUILD ANOVA TABLE ---
      lme_i_anova_tab <- lme_i_anova %>%
        rowwise() %>%
        mutate(
          gene         = gene_name,
          model_status = model_status,   # "success" or "singular_fit"
          predictor_class = case_when(
            grepl(":", term)                                                             ~ "interaction",
            is.numeric(model_data[[term]])                                              ~ "continuous",
            !is.numeric(model_data[[term]]) & length(unique(model_data[[term]])) == 2  ~ "two-level-categorical",
            !is.numeric(model_data[[term]]) & length(unique(model_data[[term]])) >  2  ~ "multi-level-categorical"
          ),
          Estimate_source = case_when(
            predictor_class %in% c("continuous", "two-level-categorical")                     ~ "lme_summary",
            predictor_class == "multi-level-categorical"                                      ~ "seeContrasts",
            predictor_class == "interaction" & length(grep(":", lme_i_summary$variable)) == 1 ~ "lme_summary",
            predictor_class == "interaction" & length(grep(":", lme_i_summary$variable)) >  1 ~ "seeContrasts"
          ),
          Estimate = case_when(
            predictor_class == "continuous"            ~ lme_i_summary$Estimate[match(term, lme_i_summary$variable)][[1]],
            predictor_class == "two-level-categorical" ~ lme_i_summary$Estimate[match(term, lme_i_summary$variable)][[1]],
            # Guard: grep(":", ...) returns integer(0) when the model has no interaction
            # coefficient (e.g. a main-effects-only formula). [[1]] on integer(0) errors,
            # so fall through to .default = NA_real_ via the length check.
            predictor_class == "interaction" & Estimate_source == "lme_summary" &
              length(grep(":", lme_i_summary$variable)) >= 1 ~
              lme_i_summary$Estimate[grep(":", lme_i_summary$variable)[1L]],
            .default = NA_real_
          ),
          Estimate_SE = case_when(
            predictor_class == "continuous"            ~ lme_i_summary$`Std. Error`[match(term, lme_i_summary$variable)][[1]],
            predictor_class == "two-level-categorical" ~ lme_i_summary$`Std. Error`[match(term, lme_i_summary$variable)][[1]],
            predictor_class == "interaction" & Estimate_source == "lme_summary" &
              length(grep(":", lme_i_summary$variable)) >= 1 ~
              lme_i_summary$`Std. Error`[grep(":", lme_i_summary$variable)[1L]],
            .default = NA_real_
          )
        )

      # --- CONTRASTS ---
      if (isTRUE(run_contrast)) {

        is_interaction <- any(grepl(":", contrast_vars))

        if (is_interaction) {

          # ---- BRANCH A: Interaction contrast via pre-computed contrasts_list ----
          # contrasts_list and spec_lookup were built once in geneLME() before
          # parallel dispatch — no per-gene reconstruction needed here.
          vars  <- strsplit(contrast_vars, ":")[[1]]
          var_a <- vars[1]
          var_b <- vars[2]

          emm_obj <- emmeans(
            lme_i,
            spec = as.formula(paste("~", var_a, "*", var_b)),
            data = model_data
          )

          emm_1st_A <- contrast(emm_obj, method = contrasts_list, adjust = "none")

          if (!is.null(contrasts_secondary)) {
            # suppressWarnings: emmeans internally calls lmer on the first-order contrast
            # estimates when computing second-order contrasts; this can emit a benign
            # lme4 scale warning when the estimate/SE values differ in magnitude from
            # the original predictors. The output is unaffected — the warning is cosmetic.
            emm_2nd_A <- suppressWarnings(
              contrast(emm_1st_A, method = contrasts_secondary, adjust = "none")
            )
            contrast_res <- bind_rows(
              # First-order: join ref/lvl from pre-computed spec_lookup by contrast name
              as.data.frame(emm_1st_A) %>%
                mutate(contrast_order = "first_order") %>%
                left_join(spec_lookup, by = "contrast"),
              # Second-order: contrasts-of-contrasts have no single ref/lvl pair
              as.data.frame(emm_2nd_A) %>%
                mutate(contrast_order  = "second_order",
                       contrast_ref    = NA_character_,
                       contrast_lvl    = NA_character_)
            ) %>% mutate(gene = gene_name, model_status = model_status)
          } else {
            contrast_res <- as.data.frame(emm_1st_A) %>%
              mutate(contrast_order = "first_order") %>%
              left_join(spec_lookup, by = "contrast") %>%
              mutate(gene = gene_name, model_status = model_status)
          }

        } else {

          # ---- BRANCH B: Non-interaction ----
          spec_formula <- as.formula(paste("~", paste(contrast_vars, collapse = "|")))

          by_list <- if (length(contrast_vars) == 2 && !is.null(contrast_var_2_levels)) {
            setNames(list(contrast_var_2_levels), contrast_vars[2])
          } else {
            NULL
          }

          emm_1st <- emmeans(lme_i, spec = spec_formula, at = by_list, data = model_data) %>%
            contrast(method = contrasts_primary, adjust = "none")

          if (!is.null(contrasts_secondary)) {
            # suppressWarnings: same benign lme4 scale warning as Branch A second-order step.
            emm_2nd <- suppressWarnings(
              contrast(emm_1st, method = contrasts_secondary, adjust = "none")
            )
            contrast_res <- bind_rows(
              # Branch B uses named contrast vectors, not a ref/lvl spec — set NA for both orders
              as.data.frame(emm_1st) %>%
                mutate(contrast_order = "first_order",
                       contrast_ref   = NA_character_,
                       contrast_lvl   = NA_character_),
              as.data.frame(emm_2nd) %>%
                mutate(contrast_order = "second_order",
                       contrast_ref   = NA_character_,
                       contrast_lvl   = NA_character_)
            ) %>% mutate(gene = gene_name, model_status = model_status)
          } else {
            contrast_res <- as.data.frame(emm_1st) %>%
              mutate(contrast_order = "first_order",
                     contrast_ref   = NA_character_,
                     contrast_lvl   = NA_character_,
                     gene           = gene_name,
                     model_status   = model_status)
          }
        }

      } else {
        contrast_res <- NULL
      }

      list(
        aic          = aic_res,
        anova        = lme_i_anova_tab,
        contrasts    = contrast_res,
        model_status = setNames(model_status, gene_name)
      )

    }, error = function(e) {
      err_msg <- conditionMessage(e)
      list(
        aic = data.frame(gene = gene_name, AIC = NA_real_),
        anova = data.frame(
          term            = NA_character_,
          statistic       = NA_real_,
          df              = NA_real_,
          p.value         = NA_real_,
          gene            = gene_name,
          model_status    = err_msg,
          predictor_class = NA_character_,
          Estimate_source = NA_character_,
          Estimate        = NA_real_,
          Estimate_SE     = NA_real_
        ),
        contrasts = data.frame(
          contrast       = NA_character_,
          estimate       = NA_real_,
          SE             = NA_real_,
          df             = NA_real_,
          t.ratio        = NA_real_,
          p.value        = NA_real_,
          contrast_order = NA_character_,
          contrast_ref   = NA_character_,
          contrast_lvl   = NA_character_,
          gene           = gene_name
        ),
        model_status = setNames(err_msg, gene_name)
      )
    })

    return(result)
  }


########################################################
# geneLME_compiler
# Aggregates list of per-gene results into named result tables,
# then appends FDR-adjusted p-values within each grouping unit.
#
# FDR grouping strategy:
#   lme_anova:    adjust within each model term (across all genes).
#                 Each term's p-values form one adjustment set.
#   lme_contrast: adjust within each contrast x contrast_order combination
#                 (across all genes). Branch B contrast labels already encode
#                 the 'by' variable level (e.g. "TrtC vs TrtA, visit = V2"),
#                 so grouping by contrast alone is sufficient.
#   NA p-values (failed gene models) are preserved as NA in p.value_adj.
#   singular_fit genes are included in the adjustment set (their p-values
#   are valid; users filter on model_status if they wish to exclude them).
########################################################

geneLME_compiler <- function(fit, fdr_method = "BH", contrast_spec = NULL) {

  lme_anova <- map_dfr(fit, "anova") %>%
    group_by(term) %>%
    mutate(p.value_adj = p.adjust(p.value, method = fdr_method)) %>%
    relocate(p.value_adj, .after = p.value)%>%
    ungroup()

  lme_contrast_raw <- map_dfr(fit, "contrasts")
  lme_contrast <- if (ncol(lme_contrast_raw) > 0 &&
                      all(c("contrast", "contrast_order") %in% colnames(lme_contrast_raw))) {
    lme_contrast_raw %>%
      group_by(contrast, contrast_order) %>%
      mutate(p.value_adj = p.adjust(p.value, method = fdr_method)) %>%
      relocate(p.value_adj, .after = p.value)%>%
      ungroup()
  } else {
    lme_contrast_raw   # no contrasts run; return as-is (empty or NULL-row stub)
  }

  list(
    lme_anova     = lme_anova,
    lme_contrast  = lme_contrast,
    lme_fit       = map_dfr(fit, "aic"),
    lme_err       = map_chr(fit, "model_status"),
    contrast_spec = contrast_spec   # indexed spec; NULL when no contrasts run
  )
}


########################################################
# geneLME_dispatch
# Runs future_lapply with explicit global declaration to prevent
# future's automatic environment scan from capturing large objects.
#
# Key design decisions:
#   1. Iterate over an integer index — a plain integer sequence carries
#      no environment baggage for future to scan.
#   2. All shared objects passed via future.globals, bypassing automatic
#      scanning entirely.
#
# contrast_spec replaced by contrasts_list + spec_lookup (Branch A)
# in both the function signature and future.globals.
########################################################

geneLME_dispatch <- function(gene_data_list,
                             targets_df,
                             formula_str,
                             run_contrast,
                             contrast_vars,
                             contrast_var_2_levels,
                             contrasts_list,
                             spec_lookup,
                             contrasts_primary,
                             contrasts_secondary) {
  n_genes <- length(gene_data_list)

  # Suppress "package 'X' was built under R version Y" warnings that
  # future re-raises on the main process after collecting worker results.
  # These are cosmetic version-mismatch notices — not actionable errors.
  withCallingHandlers(
    future.apply::future_lapply(
    seq_len(n_genes),
    FUN = function(i) {
      gene_data <- gene_data_list[[i]]
      geneLME_fit(
        gene_name             = gene_data$gene_name,
        expression_vec        = gene_data$expression_vec,
        weight_vec            = gene_data$weight_vec,
        targets               = targets_df,
        formula_str           = formula_str,
        run_contrast          = run_contrast,
        contrast_vars         = contrast_vars,
        contrast_var_2_levels = contrast_var_2_levels,
        contrasts_list        = contrasts_list,
        spec_lookup           = spec_lookup,
        contrasts_primary     = contrasts_primary,
        contrasts_secondary   = contrasts_secondary
      )
    },
    future.globals = list(
      gene_data_list        = gene_data_list,
      targets_df            = targets_df,
      formula_str           = formula_str,
      run_contrast          = run_contrast,
      contrast_vars         = contrast_vars,
      contrast_var_2_levels = contrast_var_2_levels,
      contrasts_list        = contrasts_list,
      spec_lookup           = spec_lookup,
      contrasts_primary     = contrasts_primary,
      contrasts_secondary   = contrasts_secondary,
      geneLME_fit           = geneLME_fit
    ),
    future.seed = TRUE
  ),
  warning = function(w) {
    if (grepl("was built under R version", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}


########################################################
# geneLME
# User-facing wrapper: validates inputs, sets up parallel plan,
# pre-extracts per-gene data, dispatches geneLME_fit in parallel.
#
# Calls geneLME_build_contrast_args() for Branch A to pre-compute
# contrasts_list and spec_lookup before the parallel stage.
########################################################

geneLME <-
  function(dat,
           formula_str,
           model_weights         = NULL,
           run_contrast          = NULL,
           contrast_vars         = NULL,
           contrast_var_2_levels = NULL,
           contrast_spec         = NULL,   # data.frame(contrast_ref, contrast_lvl) or NULL
                                           # required when contrast_vars contains ":"
           contrasts_primary     = NULL,
           contrasts_secondary   = NULL,
           fdr_method            = "BH",   # any method accepted by p.adjust()
           n_cores               = NULL) {

    # --- PRE-FLIGHT VALIDATION ---

    # Build a local formula object solely for variable checking — never
    # passed to workers. Workers reconstruct the formula from formula_str.
    formula_obj_local <- as.formula(paste("expression", formula_str))

    # Initialise objects populated during the contrast validation block below.
    indexed_contrast_spec <- NULL
    contrasts_list        <- NULL   # Branch A: pre-computed named contrast vectors
    spec_lookup           <- NULL   # Branch A: pre-computed ref/lvl join table

    # 1. Formula variables present in targets
    required_vars <- all.vars(formula_obj_local)
    required_vars <- required_vars[required_vars != "expression"]

    missing_vars <- setdiff(required_vars, colnames(dat$targets))
    if (length(missing_vars) > 0) {
      stop(paste(
        "The following variables in the formula are missing from dat$targets:",
        paste(missing_vars, collapse = ", ")
      ))
    }

    # 2. Weights alignment
    if (isTRUE(model_weights)) {
      if (is.null(dat$weights)) {
        stop("model_weights = TRUE but dat$weights is NULL.")
      }
      if (!identical(dim(dat$weights), dim(dat$E))) {
        stop("Dimensions of dat$weights do not match dat$E.")
      }
    }

    # 3. Contrast variables and spec
    if (isTRUE(run_contrast)) {

      if (is.null(contrast_vars)) {
        stop("run_contrast = TRUE but contrast_vars is NULL.")
      }

      is_interaction <- any(grepl(":", contrast_vars))

      # 3a. Formula-vs-contrast_vars consistency check.
      # Uses terms() so that both "a*b" and "a:b" formula syntax are handled correctly —
      # "a*b" expands to include "a:b" in term labels, so both forms are detected.
      formula_terms <- attr(terms(formula_obj_local), "term.labels")

      if (is_interaction) {
        # For interaction contrasts, verify the interaction term is actually in the model.
        # emmeans will silently run contrasts on additive margins if the interaction is absent,
        # which is statistically misleading.
        ixn_term <- contrast_vars[grep(":", contrast_vars)][1]  # e.g. "treatment:visit"
        if (!ixn_term %in% formula_terms) {
          stop(paste0(
            "contrast_vars specifies an interaction contrast for '", ixn_term, "', ",
            "but this interaction term is not present in formula_str.\n",
            "Either add the interaction to the formula (e.g. '~ treatment * visit + ...') ",
            "or change contrast_vars to a non-interaction term."
          ))
        }
      } else {
        # For non-interaction contrasts, verify each contrast variable is in the formula
        # as a main effect. Warn (not error) in case it is part of an interaction only.
        for (cv in contrast_vars) {
          if (!cv %in% formula_terms) {
            warning(paste0(
              "contrast_vars includes '", cv, "' but this term does not appear as a ",
              "main effect in formula_str. Contrasts may be computed from interaction ",
              "margins only — verify this is the intended model structure."
            ))
          }
        }
      }

      # 3b. contrast_spec validation (required for interaction contrasts)
      if (is_interaction && is.null(contrast_spec)) {
        stop(paste0(
          "contrast_vars specifies an interaction contrast ('", contrast_vars, "') ",
          "but contrast_spec is NULL.\n",
          "Use geneLME_contrast_spec(dat$targets, contrast_vars = '", contrast_vars,
          "') to generate a template, filter it to your contrasts of interest, ",
          "then pass it as contrast_spec."
        ))
      }

      if (!is.null(contrast_spec)) {
        if (!is.data.frame(contrast_spec)) {
          stop("contrast_spec must be a data.frame.")
        }
        if (!all(c("contrast_ref", "contrast_lvl") %in% colnames(contrast_spec))) {
          stop("contrast_spec must have columns 'contrast_ref' and 'contrast_lvl'.")
        }
        if (!is_interaction) {
          warning("contrast_spec is provided but contrast_vars contains no interaction term (':'). ",
                  "contrast_spec is designed for interaction contrasts — did you mean to use contrasts_primary?")
        }

        # Build the indexed contrast_spec that will be attached to the return value.
        # contrast_index here is simply 1:nrow(contrast_spec) — the actual row position
        # within the filtered spec passed by the user. This is what contrasts_secondary
        # vectors must index into (not any index from the full unfiltered template).
        indexed_contrast_spec <- contrast_spec %>%
          mutate(contrast_index = seq_len(n())) %>%
          select(contrast_index, everything())

        # Print the indexed spec as a reminder whenever contrasts_secondary is provided.
        if (is_interaction && !is.null(contrasts_secondary)) {
          n_first <- nrow(contrast_spec)
          message(
            "contrasts_secondary will be applied to the first-order interaction contrasts ",
            "in the order they appear in contrast_spec (", n_first, " contrasts):\n",
            paste(seq_len(n_first),
                  paste(contrast_spec$contrast_lvl, contrast_spec$contrast_ref, sep = " - "),
                  sep = ". ", collapse = "\n"),
            "\nEnsure contrasts_secondary vectors have length ", n_first,
            ", with each element corresponding to the contrast at that position.",
            "\nThe indexed contrast_spec is returned as $contrast_spec in the output."
          )

          # Soft-fail: wrong-length vectors produce silent NAs deep inside geneLME_fit().
          # Catch this here and return early with only $contrast_spec populated so the
          # user has the indexed reference they need to fix their vectors.
          bad_lens <- sapply(contrasts_secondary, length)
          bad_names <- names(bad_lens)[bad_lens != n_first]
          if (length(bad_names) > 0) {
            message(
              "\ncontrasts_secondary vector(s) have wrong length — returning early without running models.\n",
              "Expected length ", n_first, " (= nrow(contrast_spec) after filtering).\n",
              "Offending vector(s):\n",
              paste0("  '", bad_names, "': length ", bad_lens[bad_names], collapse = "\n"), "\n\n",
              "Common cause: vectors were built against nrow() of the full unfiltered\n",
              "spec_template rather than nrow() of the filtered contrast_spec passed here.\n",
              "Use $contrast_spec in the returned object to re-specify your secondary contrast vectors:\n",
              "  each vector must have length ", n_first, ", one element per row of $contrast_spec."
            )
            return(invisible(list(
              lme_anova    = NULL,
              lme_contrast = NULL,
              lme_fit      = NULL,
              lme_err      = NULL,
              contrast_spec = indexed_contrast_spec
            )))
          }
        }
      }

      # --- PRE-COMPUTE BRANCH A CONTRAST STRUCTURES ---
      # Build contrasts_list and spec_lookup once here, before parallel dispatch.
      # Workers receive these ready-made objects instead of rebuilding from
      # contrast_spec on every gene — eliminates nrow(contrast_spec) × n_genes
      # iterations of R-level vector construction.
      if (is_interaction) {
        contrast_args  <- geneLME_build_contrast_args(dat$targets, contrast_vars, contrast_spec)
        contrasts_list <- contrast_args$contrasts_list
        spec_lookup    <- contrast_args$spec_lookup
      }

      # 3c. contrast_vars present in targets (split on ":" for interaction)
      vars_to_check <- unique(unlist(strsplit(contrast_vars, ":")))
      missing_contrast_vars <- setdiff(vars_to_check, colnames(dat$targets))
      if (length(missing_contrast_vars) > 0) {
        stop(paste(
          "contrast_vars not found in dat$targets:",
          paste(missing_contrast_vars, collapse = ", ")
        ))
      }

      # 3d. contrast_var_2_levels validation (Branch B only)
      if (!is_interaction && length(contrast_vars) == 2 && !is.null(contrast_var_2_levels)) {
        actual_levels  <- unique(as.character(dat$targets[[contrast_vars[2]]]))
        invalid_levels <- setdiff(contrast_var_2_levels, actual_levels)
        if (length(invalid_levels) > 0) {
          stop(paste0(
            "Levels specified for '", contrast_vars[2],
            "' not found in data: ", paste(invalid_levels, collapse = ", ")
          ))
        }
      }

      # 3e. Branch B: build indexed contrast_spec from contrasts_primary names,
      #     and soft-fail on wrong-length contrasts_secondary vectors.
      if (!is_interaction && !is.null(contrasts_primary)) {
        indexed_contrast_spec <- data.frame(
          contrast_index = seq_along(contrasts_primary),
          contrast_name  = names(contrasts_primary),
          stringsAsFactors = FALSE
        )

        if (!is.null(contrasts_secondary)) {
          n_primary <- length(contrasts_primary)
          bad_lens  <- sapply(contrasts_secondary, length)
          bad_names <- names(bad_lens)[bad_lens != n_primary]
          if (length(bad_names) > 0) {
            message(
              "\ncontrasts_secondary vector(s) have wrong length — returning early without running models.\n",
              "Expected length ", n_primary, " (= number of contrasts_primary vectors).\n",
              "Offending vector(s):\n",
              paste0("  '", bad_names, "': length ", bad_lens[bad_names], collapse = "\n"), "\n\n",
              "Use $contrast_spec in the returned object to confirm the primary contrast ordering:\n",
              "  each contrasts_secondary vector must have length ", n_primary,
              ", one element per row of $contrast_spec."
            )
            return(invisible(list(
              lme_anova    = NULL,
              lme_contrast = NULL,
              lme_fit      = NULL,
              lme_err      = NULL,
              contrast_spec = indexed_contrast_spec
            )))
          }
        }
      }
    }

    # 4. fdr_method must be a recognised p.adjust method
    valid_fdr_methods <- p.adjust.methods   # exported character vector from stats
    if (!fdr_method %in% valid_fdr_methods) {
      stop(paste0(
        "'", fdr_method, "' is not a valid p.adjust method.\n",
        "Choose one of: ", paste(valid_fdr_methods, collapse = ", ")
      ))
    }

    # 5. Warn on likely ID columns used as predictors
    for (v in required_vars) {
      if (is.character(dat$targets[[v]]) &&
          length(unique(dat$targets[[v]])) > nrow(dat$targets) / 2) {
        warning(paste0(
          "Variable '", v, "' has many unique string values — is it an ID column rather than a predictor?"
        ))
      }
    }

    message("Input validation passed. Launching parallel LME fits...")

    # --- PRE-SLICE INTO PER-GENE LIST ---
    W_mat <- if (isTRUE(model_weights)) dat$weights else NULL

    gene_data_list <- lapply(rownames(dat$E), function(g) {
      list(
        gene_name      = g,
        expression_vec = dat$E[g, ],
        weight_vec     = if (!is.null(W_mat)) W_mat[g, ] else NULL
      )
    })

    targets_df <- dat$targets

    # --- PARALLEL PLAN ---
    workers <- if (is.null(n_cores)) max(1L, parallel::detectCores() - 4L) else n_cores
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    # --- DISPATCH via clean-scope helper ---
    fit_results <- geneLME_dispatch(
      gene_data_list        = gene_data_list,
      targets_df            = targets_df,
      formula_str           = formula_str,
      run_contrast          = run_contrast,
      contrast_vars         = contrast_vars,
      contrast_var_2_levels = contrast_var_2_levels,
      contrasts_list        = contrasts_list,
      spec_lookup           = spec_lookup,
      contrasts_primary     = contrasts_primary,
      contrasts_secondary   = contrasts_secondary
    )

    return(geneLME_compiler(fit_results, fdr_method = fdr_method,
                            contrast_spec = indexed_contrast_spec))
  }


########################################################
# Testing (see geneLME_test.R for full test suite)
########################################################

# --- Branch A: interaction contrast via contrast_spec ---
# # Step 1: generate level reference template
# spec_template <- geneLME_contrast_spec(
#   targets       = dat$targets,
#   contrast_vars = "treatment:visit"    # interaction string
# )
# # Step 2: filter to contrasts of interest
# my_spec <- spec_template %>%
#   dplyr::filter(...)   # e.g. same-visit cross-treatment, or within-treatment longitudinal
#
# test_mods_A <-
#   geneLME(
#     dat           = dat_sub,
#     formula_str   = "~ treatment * visit + age + sex + rNANgUl + (1|ptID)",
#     model_weights = TRUE,
#     run_contrast  = TRUE,
#     contrast_vars = "treatment:visit",  # must match an interaction term in formula_str
#     contrast_spec = my_spec,            # required for interaction contrasts
#     n_cores       = 10
#   )
#
# # Optional: second-order contrasts on top of Branch A first-order contrasts.
# # Vectors must have length == nrow(my_spec), ordered to match my_spec rows.
# # geneLME() will print the numbered list of first-order contrasts at validation time.
# test_mods_A2 <-
#   geneLME(
#     dat                 = dat_sub,
#     formula_str         = "~ treatment * visit + age + sex + rNANgUl + (1|ptID)",
#     model_weights       = TRUE,
#     run_contrast        = TRUE,
#     contrast_vars       = "treatment:visit",
#     contrast_spec       = my_spec,
#     contrasts_secondary = list(
#       "TrtA vs TrtB: V3 minus V2 effect" = c(1, 0, -1, 0, 0, 0),
#       "TrtA vs TrtC: V3 minus V2 effect" = c(0, 1, 0, -1, 0, 0)
#     ),
#     n_cores             = 10
#   )

# --- Branch B: non-interaction ---
# # Step 1 (optional): inspect available levels for each contrast variable
# geneLME_contrast_spec(dat$targets, contrast_vars = "treatment")  # reference only
# # Treatment levels (alphabetical): TrtA, TrtB, TrtC
# # contrasts_primary vectors have length 3: positions = [TrtA, TrtB, TrtC]
#
# test_mods_B <-
#   geneLME(
#     dat                   = dat_sub,
#     formula_str           = "~ treatment + visit + age + sex + rNANgUl + (1|ptID)",
#     model_weights         = TRUE,
#     run_contrast          = TRUE,
#     contrast_vars         = c("treatment", "visit"),
#     contrast_var_2_levels = c("V2", "V3"),
#     contrasts_primary     = list("TrtC vs TrtA" = c(-1, 0, 1),
#                                  "TrtB vs TrtA" = c(-1, 1, 0)),
#     contrasts_secondary   = list("TrtC vs TrtB" = c(1, -1)),
#     n_cores               = 10
#   )
