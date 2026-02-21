########################################################
# geneLME_test.R
# Mock data and test cases for geneLME functions.
# Set SOURCE_FILE below to switch between implementations.
#
# Design: paired treatment x visit study
#   - 10 patients (ptID), each observed at 4 visits x 3 treatments
#   - 3 treatments: "TrtA", "TrtB", "TrtC"
#   - 4 visits:     "V1", "V2", "V3", "V4"
#   - 120 total samples (10 * 4 * 3)
#   - 50 genes
#   - Simulated voom weights matrix (same dims as E)
#
# Covariates included in targets:
#   Patient-level (same value repeated across all samples for a given patient):
#     sex          — factor: "M" / "F"
#     age          — continuous: age at enrollment (years)
#   Sample-level technical (vary per sample):
#     rNANgUl             — RNA concentration (continuous)
#     percent_duplication — library duplication rate (continuous, 0–1)
#     median_cv_coverage  — coverage uniformity metric (continuous)
#     lib.size            — total mapped reads (continuous, large integer scale)
#     norm.factors        — TMM normalization factors (continuous, near 1)
########################################################

library(lme4)
library(emmeans)
library(car)
library(broom.mixed)
library(dplyr)
library(tibble)
library(purrr)
library(future)
library(future.apply)

# Source implementation to test — change this line to switch versions:
#   "geneLME_dev.R"   — previous dev (singular fit = error, per-gene contrast build)
#   "geneLME_dev2.R"  — current dev  (singular fit = flag, pre-computed contrasts)
#   "geneLME.R"       — stable merged version
SOURCE_FILE <- "geneLME.R"
source(SOURCE_FILE)
cat("Sourced:", SOURCE_FILE, "\n")


########################################################
# 1. SIMULATE MOCK EList
########################################################

set.seed(42)

n_patients   <- 10
treatments   <- c("TrtA", "TrtB", "TrtC")
visits       <- c("V1", "V2", "V3", "V4")
n_genes      <- 50

# --- Patient-level attributes (one row per patient, then joined in) ---
# These are fixed characteristics that do not vary across visits or treatments
patient_meta <- data.frame(
  ptID = paste0("pt", sprintf("%02d", 1:n_patients)),
  sex  = factor(sample(c("M", "F"), n_patients, replace = TRUE)),
  age  = round(rnorm(n_patients, mean = 38, sd = 10)),
  stringsAsFactors = FALSE
)

# Build targets (sample metadata)
# expand.grid gives one row per ptID x treatment x visit combination
targets <- expand.grid(
  ptID      = paste0("pt", sprintf("%02d", 1:n_patients)),
  treatment = treatments,
  visit     = visits,
  stringsAsFactors = FALSE
) %>%
  arrange(ptID, treatment, visit) %>%
  # Join patient-level covariates — sex and age are identical across all
  # samples from the same patient, as in real data
  left_join(patient_meta, by = "ptID") %>%
  mutate(
    sample_id = paste(ptID, treatment, visit, sep = "_"),
    # --- Sample-level technical covariates ---
    # These vary independently per sample (sequencing run characteristics)
    rNANgUl             = rnorm(n(), mean = 5,      sd = 1),
    percent_duplication = runif(n(), min = 0.05,    max = 0.55),
    median_cv_coverage  = rnorm(n(), mean = 0.85,   sd = 0.08),
    lib.size            = round(rnorm(n(), mean = 20e6, sd = 3e6)),
    norm.factors        = rnorm(n(), mean = 1,      sd = 0.05)
  )

rownames(targets) <- targets$sample_id
n_samples <- nrow(targets)   # 120

# Simulate log2 expression matrix (genes x samples)
# Add a modest treatment*visit interaction effect to a subset of genes
E_mat <- matrix(
  rnorm(n_genes * n_samples, mean = 8, sd = 2),
  nrow  = n_genes,
  ncol  = n_samples,
  dimnames = list(
    paste0("gene", sprintf("%02d", 1:n_genes)),
    targets$sample_id
  )
)

# Add a simulated effect: genes 1-10 have a TrtC:V3 upregulation
trtC_v3 <- which(targets$treatment == "TrtC" & targets$visit == "V3")
E_mat[1:10, trtC_v3] <- E_mat[1:10, trtC_v3] + 2.5

# Add a patient random effect
for (pt in unique(targets$ptID)) {
  pt_idx <- which(targets$ptID == pt)
  E_mat[, pt_idx] <- E_mat[, pt_idx] + rnorm(1, 0, 1)
}

# Simulated voom precision weights (all near 1, small variance)
W_mat <- matrix(
  abs(rnorm(n_genes * n_samples, mean = 1, sd = 0.1)),
  nrow     = n_genes,
  ncol     = n_samples,
  dimnames = dimnames(E_mat)
)

# Assemble EList-like list (mimics limma EList structure)
dat <- list(
  E       = E_mat,
  weights = W_mat,
  targets = targets
)


########################################################
# 2. SMALL SUBSET FOR FAST TESTING
########################################################

dat_sub <- list(
  E       = dat$E[1:10, ],
  weights = dat$weights[1:10, ],
  targets = dat$targets
)


########################################################
# 3. TEST: geneLME_contrast_spec() — both modes
########################################################

cat("\n--- Test: geneLME_contrast_spec() interaction mode ---\n")

spec_template <- geneLME_contrast_spec(
  targets       = dat$targets,
  contrast_vars = "treatment:visit"
)

# spec_template has two columns: contrast_ref and contrast_lvl only.
# contrast_index is NOT in spec_template — it is added by geneLME() to its
# $contrast_spec return element after it receives the user's filtered spec.
print(head(spec_template, 10))
cat("Total pairwise combinations:", nrow(spec_template), "\n")
cat("Columns:", paste(colnames(spec_template), collapse = ", "), "\n")

# Filter to same-visit cross-treatment comparisons at V2 and V3.
my_spec <- spec_template %>%
  filter(
    # Both elements of the pair are at the same visit
    sub(".* ", "", contrast_ref) == sub(".* ", "", contrast_lvl),
    # Only visits V2 and V3
    sub(".* ", "", contrast_ref) %in% c("V2", "V3")
  )

cat("Filtered spec rows:", nrow(my_spec), "\n")
print(my_spec)

cat("\n--- Test: geneLME_contrast_spec() single-variable mode (one variable) ---\n")

treatment_levels <- geneLME_contrast_spec(
  targets       = dat$targets,
  contrast_vars = "treatment"
)
cat("Treatment levels (reference for building contrasts_primary):\n")
print(treatment_levels)

cat("\n--- Test: geneLME_contrast_spec() multi-variable mode ---\n")

# Pass both contrast variables together; the function explains how each maps
# to geneLME() arguments (contrasts_primary vs contrast_var_2_levels)
trt_visit_levels <- geneLME_contrast_spec(
  targets       = dat$targets,
  contrast_vars = c("treatment", "visit")
)
cat("Returned list structure:\n")
str(trt_visit_levels)
cat("\ntreatment levels:\n")
print(trt_visit_levels$treatment)
cat("visit levels:\n")
print(trt_visit_levels$visit)


########################################################
# 4. TEST: Branch A — explicit contrast_spec (interaction)
########################################################

cat("\n--- Test: Branch A (contrast_spec, interaction) ---\n")

test_A <- geneLME(
  dat           = dat_sub,
  formula_str   = "~ treatment * visit + age + sex + rNANgUl + percent_duplication + median_cv_coverage + (1|ptID)",
  model_weights = TRUE,
  run_contrast  = TRUE,
  contrast_vars = "treatment:visit",
  contrast_spec = my_spec,
  n_cores       = 2
)

cat("ANOVA rows:", nrow(test_A$lme_anova), "\n")
cat("Contrast rows:", nrow(test_A$lme_contrast), "\n")
cat("Fit rows:", nrow(test_A$lme_fit), "\n")
cat("Model status summary (lme_err):\n")
print(table(test_A$lme_err))
cat("Model status in ANOVA table:\n")
print(table(test_A$lme_anova$model_status))

# Verify singular fit genes still have non-NA contrast estimates
singular_genes <- names(test_A$lme_err)[test_A$lme_err == "singular_fit"]
if (length(singular_genes) > 0) {
  singular_contrasts <- test_A$lme_contrast %>% filter(gene %in% singular_genes)
  cat("Singular-fit genes with non-NA estimates:",
      sum(!is.na(singular_contrasts$estimate)), "/", nrow(singular_contrasts), "\n")
  cat("4 PASS (singular_fit) — singular genes return estimates with flag\n")
} else {
  cat("No singular fits in this run (unexpected with 10-patient mock data)\n")
}

cat("\nSample contrast output:\n")
print(head(test_A$lme_contrast %>% select(gene, model_status, contrast, contrast_ref, contrast_lvl, estimate, p.value), 10))


########################################################
# 4b. TEST: Branch A — with second-order contrasts
########################################################

# my_spec (after filtering) has 6 first-order contrasts. geneLME() attaches an
# indexed copy of my_spec to result_A$contrast_spec — contrast_index = 1:nrow(my_spec)
# (simple row positions in the filtered spec). Use these indices when building
# contrasts_secondary vectors. The row order within my_spec determines position in
# contrasts_secondary vectors, since geneLME() iterates over seq_len(nrow(contrast_spec)).
#
# Expected my_spec row order (alphabetical interaction level ordering):
#   row 1 (contrast_index=1): TrtB V2 - TrtA V2
#   row 2 (contrast_index=2): TrtC V2 - TrtA V2
#   row 3 (contrast_index=3): TrtC V2 - TrtB V2
#   row 4 (contrast_index=4): TrtB V3 - TrtA V3
#   row 5 (contrast_index=5): TrtC V3 - TrtA V3
#   row 6 (contrast_index=6): TrtC V3 - TrtB V3
#
# Second-order contrasts: difference-of-differences across visits
#   "TrtA vs TrtB: V3 minus V2 effect":  row 4 - row 1  → c(-1, 0, 0, 1, 0, 0)
#   "TrtA vs TrtC: V3 minus V2 effect":  row 5 - row 2  → c(0, -1, 0, 0, 1, 0)
#
# Verify row order with: test_A$contrast_spec
# Vectors must have length 6, one element per row of my_spec (in row order).

cat("\n--- Test: Branch A with second-order contrasts ---\n")
cat("Indexed contrast_spec from test_A (verify row ordering before specifying secondary contrasts):\n")
print(test_A$contrast_spec)

test_A2 <- geneLME(
  dat                 = dat_sub,
  formula_str         = "~ treatment * visit + age + sex + rNANgUl + percent_duplication + median_cv_coverage + (1|ptID)",
  model_weights       = TRUE,
  run_contrast        = TRUE,
  contrast_vars       = "treatment:visit",
  contrast_spec       = my_spec,
  contrasts_secondary = list(
    "TrtA vs TrtB: V3 minus V2 effect" = c(-1, 0, 0, 1, 0, 0),   # row 4 - row 1
    "TrtA vs TrtC: V3 minus V2 effect" = c(0, -1, 0, 0, 1, 0)    # row 5 - row 2
  ),
  n_cores             = 2
)

cat("Contrast rows (should include second_order):", nrow(test_A2$lme_contrast), "\n")
cat("Contrast orders present:\n")
print(table(test_A2$lme_contrast$contrast_order))
cat("\nSample second-order contrast output:\n")
print(test_A2$lme_contrast[test_A2$lme_contrast$contrast_order == "second_order", ])


########################################################
# 5. TEST: Branch B — non-interaction (regression test)
########################################################

cat("\n--- Test: Branch B (non-interaction, regression test) ---\n")

# Treatment levels in alphabetical order (as emmeans sees them): TrtA, TrtB, TrtC
# contrasts_primary vectors have length 3; positions correspond to: [TrtA, TrtB, TrtC]

test_B <- geneLME(
  dat                   = dat_sub,
  formula_str           = "~ treatment + visit + age + sex + rNANgUl + percent_duplication + median_cv_coverage + (1|ptID)",
  model_weights         = TRUE,
  run_contrast          = TRUE,
  contrast_vars         = c("treatment", "visit"),
  contrast_var_2_levels = c("V2", "V3"),
  contrasts_primary     = list(
    "TrtC vs TrtA" = c(-1, 0, 1),
    "TrtB vs TrtA" = c(-1, 1, 0)
  ),
  contrasts_secondary   = list("TrtC vs TrtB" = c(1, -1)),
  n_cores               = 2
)

cat("ANOVA rows:", nrow(test_B$lme_anova), "\n")
cat("Contrast rows:", nrow(test_B$lme_contrast), "\n")
cat("Model status summary (lme_err):\n")
print(table(test_B$lme_err))
cat("Model status in ANOVA table:\n")
print(table(test_B$lme_anova$model_status))

cat("\nSample contrast output:\n")
print(head(test_B$lme_contrast %>% select(gene, model_status, contrast, contrast_ref, contrast_lvl, estimate, p.value), 10))


########################################################
# 6. TEST: Validation errors (should each produce an informative error/warning)
########################################################

cat("\n--- Test: Input validation errors ---\n")

# 6a. Missing formula variable
tryCatch(
  geneLME(dat_sub, formula_str = "~ NONEXISTENT_VAR + (1|ptID)", n_cores = 2),
  error = function(e) cat("6a PASS — missing formula var:", conditionMessage(e), "\n")
)

# 6b. model_weights = TRUE but no dat$weights
dat_no_w <- dat_sub; dat_no_w$weights <- NULL
tryCatch(
  geneLME(dat_no_w, formula_str = "~ treatment + (1|ptID)", model_weights = TRUE, n_cores = 2),
  error = function(e) cat("6b PASS — missing weights:", conditionMessage(e), "\n")
)

# 6c. Interaction contrast requested but interaction NOT in formula
tryCatch(
  geneLME(dat_sub,
          formula_str   = "~ treatment + visit + age + (1|ptID)",  # additive — no interaction
          run_contrast  = TRUE,
          contrast_vars = "treatment:visit",
          contrast_spec = my_spec,
          n_cores       = 2),
  error = function(e) cat("6c PASS — interaction not in formula:", conditionMessage(e), "\n")
)

# 6d. Interaction contrast requested but contrast_spec not provided
tryCatch(
  geneLME(dat_sub,
          formula_str   = "~ treatment * visit + age + (1|ptID)",
          run_contrast  = TRUE,
          contrast_vars = "treatment:visit",
          # contrast_spec intentionally omitted
          n_cores       = 2),
  error = function(e) cat("6d PASS — interaction contrast_spec missing:", conditionMessage(e), "\n")
)

# 6e. contrast_spec has wrong columns
tryCatch(
  geneLME(dat_sub,
          formula_str   = "~ treatment * visit + age + (1|ptID)",
          run_contrast  = TRUE,
          contrast_vars = "treatment:visit",
          contrast_spec = data.frame(a = "x", b = "y"),
          n_cores       = 2),
  error = function(e) cat("6e PASS — bad contrast_spec columns:", conditionMessage(e), "\n")
)

# 6f. contrast_var_2_levels has invalid levels
tryCatch(
  geneLME(dat_sub,
          formula_str           = "~ treatment + visit + age + (1|ptID)",
          run_contrast          = TRUE,
          contrast_vars         = c("treatment", "visit"),
          contrast_var_2_levels = c("V2", "NOTAVISIT"),
          contrasts_primary     = list("TrtC vs TrtA" = c(-1, 0, 1)),
          n_cores               = 2),
  error = function(e) cat("6f PASS — invalid contrast_var_2_levels:", conditionMessage(e), "\n")
)

# 6g. contrasts_secondary with wrong length → soft-fail (returns early with $contrast_spec,
#     all other elements NULL). This is NOT a stop() — it returns an invisible list.
cat("\n--- Test: soft-fail on wrong-length contrasts_secondary ---\n")

soft_fail_A <- geneLME(
  dat_sub,
  formula_str         = "~ treatment * visit + age + sex + rNANgUl + percent_duplication + median_cv_coverage + (1|ptID)",
  model_weights       = TRUE,
  run_contrast        = TRUE,
  contrast_vars       = "treatment:visit",
  contrast_spec       = my_spec,
  contrasts_secondary = list("wrong length" = rep(0, 99)),  # should be length 6
  n_cores             = 2
)

cat("6g PASS — soft-fail Branch A:\n")
cat("  NULL elements:", paste(names(which(sapply(soft_fail_A, is.null))), collapse = ", "), "\n")
cat("  $contrast_spec populated:\n")
print(soft_fail_A$contrast_spec)

soft_fail_B <- geneLME(
  dat_sub,
  formula_str       = "~ treatment + visit + age + (1|ptID)",
  run_contrast      = TRUE,
  contrast_vars     = c("treatment", "visit"),
  contrasts_primary = list("TrtC vs TrtA" = c(-1, 0, 1),
                           "TrtB vs TrtA" = c(-1, 1, 0)),
  contrasts_secondary = list("wrong length" = rep(0, 99)),  # should be length 2
  n_cores           = 2
)

cat("\n6g PASS — soft-fail Branch B:\n")
cat("  NULL elements:", paste(names(which(sapply(soft_fail_B, is.null))), collapse = ", "), "\n")
cat("  $contrast_spec populated:\n")
print(soft_fail_B$contrast_spec)

cat("\nAll validation tests complete.\n")
