# geneLME Project — Claude Session Notes

> **How to use this file:** At the start of each session, ask Claude to read this file.
> Claude should update the "Current Status" and "Session Log" sections at the end of any
> session where meaningful progress is made.

---

## Project Overview

**Goal:** A scalable, parallelized per-gene linear mixed effects (LME) modelling function for
RNA-seq expression data stored in a `limma` EList (`dat`) object. Key features:
- Fits one `lmer` model per gene in parallel via `future_lapply`
- Extracts ANOVA table, coefficient estimates, and emmeans-based contrasts
- Supports weighted models (voom weights from `dat$weights`)
- Supports flexible contrast specification including interaction terms

**Active worktree:** `.claude/worktrees/reverent-satoshi/`
**Key files:**
- `R_functions/geneLME.R` — current stable (merged from dev2 2026-02-20; fully featured)
- `R_functions/geneLME_dev.R` — former dev version (superseded; kept for reference)
- `R_functions/geneLME_dev2.R` — former dev2 (now merged into geneLME.R; kept for reference)
- `R_functions/geneLME_test.R` — mock data and full test suite (SOURCE_FILE variable selects version)
- `R_functions/geneLME_benchmark.Rmd` — v1 benchmark (geneLME vs kimma)
- `R_functions/geneLME_benchmark2.Rmd` — v2 benchmark (sign consistency + speed vs kimma; dev2 correctness section removed post-merge)
- `R_functions/CLAUDE_NOTES_geneLME.md` — this file

---

## Architecture Summary

The code is split into five functions:

| Function | Role |
|---|---|
| `geneLME_contrast_spec()` | Helper: generates a contrast level reference template. Two modes (see below). |
| `geneLME_build_contrast_args()` | Branch A helper: pre-computes `contrasts_list` + `spec_lookup` once before parallel dispatch. |
| `geneLME_fit()` | Per-gene worker: fits lmer, extracts ANOVA + contrasts. Runs inside `future_lapply`. |
| `geneLME_dispatch()` | Parallel launcher: calls `future_lapply` with explicit `future.globals` to avoid environment scanning overhead. |
| `geneLME_compiler()` | Aggregator: binds per-gene results into tidy data frames. |
| `geneLME()` | User-facing wrapper: validates inputs, pre-slices data, sets parallel plan, calls dispatch + compiler. |

### Key Design Decisions Made

1. **Per-gene pre-slicing** — `dat$E` and `dat$weights` are sliced into a `gene_data_list` of small bundles
   *before* entering the parallel stage. This avoids serializing the full EList to every worker.
2. **Formula passed as string** — `formula_str` (a plain character string) is passed to workers, which
   reconstruct the formula locally. Avoids locked-environment errors from passing formula objects with
   environment references across futures.
3. **Explicit `future.globals`** — `geneLME_dispatch()` explicitly lists all globals rather than relying
   on future's automatic environment scan, preventing accidental capture of large objects.
4. **`on.exit(plan(sequential))`** — restores the user's parallel plan after the run completes.
5. **Singular fit flagged, not errored** — `isSingular()` is checked and sets `model_status = "singular_fit"`.
   Results are returned for all genes; users filter downstream on `model_status`. Fixed effect estimates
   and contrasts are numerically valid even for singular fits — only the random effect structure is unreliable.

---

## Contrast System (`geneLME_dev.R`)

Two contrast branches in `geneLME_fit()`:

### Branch A: Interaction contrasts (`contrast_vars` contains `":"`)
- `contrast_vars` is a single string e.g. `"treatment:visit"`.
- **`contrast_spec`** (required): a data.frame with columns `contrast_ref` and `contrast_lvl`.
  Each row defines one pairwise contrast across interaction cells. Contrast vectors are built from the
  emmeans level ordering to guarantee index alignment.
- Use `geneLME_contrast_spec(dat$targets, contrast_vars = "treatment:visit")` to generate the full
  pairwise template, then filter to the rows of interest before passing as `contrast_spec`.
- The formula must contain the corresponding interaction term (`treatment:visit` or `treatment*visit`);
  `geneLME()` validation enforces this — emmeans would otherwise silently use additive margins.

### Branch B: Non-interaction contrasts
- `contrast_vars` is a character vector of 1–2 variable names.
- `contrasts_primary`: named list of contrast vectors. Vector length must match the number of levels
  of `contrast_vars[1]`, in alphabetical order (or the order emmeans sees them).
- `contrasts_secondary`: optional second-order contrast (contrast-of-contrasts).
- `contrast_var_2_levels`: optional filter on levels of the second `by` variable.

### Helper: `geneLME_contrast_spec(targets, contrast_vars)`
Two modes determined by whether `contrast_vars` contains `":"`:

- **Interaction mode** (`"treatment:visit"`): returns a data.frame with columns `contrast_ref` and
  `contrast_lvl` **only** (no `contrast_index`). One row per pairwise combination of interaction-level
  strings. User filters and passes as `contrast_spec` to `geneLME()`.
  - `contrast_index` is NOT part of the `geneLME_contrast_spec()` output. Instead, `geneLME()` attaches
    an indexed copy of the user's filtered spec to its return value as `$contrast_spec`, with
    `contrast_index = 1:nrow(contrast_spec)` (simple row positions within the filtered spec).
    This is the index users should reference when building `contrasts_secondary` vectors.
  - Ordering guarantee: both `var_a` and `var_b` components of `contrast_ref` are always at a
    level-index ≤ the corresponding component in `contrast_lvl` (factor level order for factors,
    alphabetical for plain characters). This is enforced by explicit factor coercion before
    `interaction(...) %>% levels()`.
  - Swap tolerance: users may manually swap ref/lvl in any row after filtering. `geneLME_fit()`
    handles either direction — the contrast estimate sign flips, but no error occurs.
- **Single-variable mode** (`"treatment"`): returns a data.frame with a single column `level`, listing
  unique sorted values of that variable. This is a reference only — used when building `contrasts_primary`
  vectors for Branch B. Not passed to `geneLME()` directly.

---

## Function Signatures (current `geneLME_dev.R`)

```r
geneLME_contrast_spec(targets, contrast_vars)
# Returns: data.frame(contrast_ref, contrast_lvl) [interaction mode]
#       or named list of data.frame(level) [single/multi-variable mode]

geneLME_fit(gene_name, expression_vec, weight_vec, targets,
            formula_str, run_contrast, contrast_vars,
            contrast_var_2_levels, contrast_spec,
            contrasts_primary, contrasts_secondary)

geneLME_dispatch(gene_data_list, targets_df, formula_str,
                 run_contrast, contrast_vars, contrast_var_2_levels,
                 contrast_spec, contrasts_primary, contrasts_secondary)

geneLME_compiler(fit, fdr_method = "BH", contrast_spec = NULL)
# Returns: list(lme_anova, lme_contrast, lme_fit, lme_err, contrast_spec)

geneLME(dat, formula_str,
        model_weights         = NULL,
        run_contrast          = NULL,
        contrast_vars         = NULL,
        contrast_var_2_levels = NULL,
        contrast_spec         = NULL,   # required when contrast_vars contains ":"
        contrasts_primary     = NULL,
        contrasts_secondary   = NULL,
        fdr_method            = "BH",   # any method accepted by p.adjust()
        n_cores               = NULL)
# Returns: list(lme_anova, lme_contrast, lme_fit, lme_err, contrast_spec)
#   $lme_contrast columns include contrast_ref and contrast_lvl:
#     - Branch A first-order: populated from contrast_spec (unambiguous sign reference)
#     - Branch A second-order, Branch B (all): NA (no single ref/lvl pair)
#   $contrast_spec: indexed copy of the filtered contrast_spec (Branch A) or
#                   data.frame(contrast_index, contrast_name) from contrasts_primary (Branch B);
#                   NULL when no contrasts run.
#   On soft-fail (wrong-length contrasts_secondary): all elements NULL except $contrast_spec.
```

Note: `contrast_by` parameter was removed in session 2026-02-19 (Branch A2 eliminated).

---

## Validation Checks in `geneLME()`

1. All formula variables present in `dat$targets`
2. `dat$weights` present and dimension-matched if `model_weights = TRUE`
3. `contrast_vars` not NULL when `run_contrast = TRUE`
4. **Formula-vs-`contrast_vars` consistency:**
   - Interaction contrasts: interaction term must exist in `formula_str` (uses `attr(terms(), "term.labels")`)
   - Non-interaction contrasts: each `contrast_vars` element warned if absent as main effect
5. `contrast_spec` required (not NULL) when `contrast_vars` contains `":"`
6. `contrast_spec` column names validated (`contrast_ref`, `contrast_lvl`)
7. All `contrast_vars` components present in `dat$targets`
8. `contrast_var_2_levels` values valid against actual data levels
9. `fdr_method` validated against `p.adjust.methods`
10. `contrasts_secondary` vector lengths validated (**soft-fail**, not hard stop):
    - Branch A: each vector must have length == `nrow(contrast_spec)` (the filtered spec, NOT the full template)
    - Branch B: each vector must have length == `length(contrasts_primary)`
    - On mismatch: emits a `message()` naming offending vectors (observed vs expected length),
      then returns `invisible(list(lme_anova=NULL, lme_contrast=NULL, lme_fit=NULL, lme_err=NULL, contrast_spec=indexed_contrast_spec))`
    - `$contrast_spec` in the returned list contains the indexed spec so the user can inspect
      row ordering and fix their vectors without needing to re-run everything
11. Warn on likely ID columns used as predictors

---

## Current Status

**Date of last update:** 2026-02-20 (Session 15)

**Stage:** `geneLME_dev2.R` merged into `geneLME.R` (now the single active stable file).
All tests pass. Warning suppression added. Documentation updated.

**Benchmark results (2,000 genes, 5 reps, 8 cores) — geneLME_benchmark2.html (Session 15):**

| Metric | Result |
|---|---|
| Sign consistency (geneLME vs kimma) | 3300/3300 (100%) same direction — no flipped pairs |
| Estimate r (direction-corrected) | 1.00000000 |
| Estimate MAD (direction-corrected) | 1.16e-15 (effectively 0 — floating point only) |
| geneLME  3 contrasts median | 17.4 s (1.85× faster than kimma) |
| geneLME  6 contrasts median | 18.1 s (1.78× faster than kimma) |
| geneLME 66 contrasts median | 32.6 s (0.99× — effectively equal to kimma) |
| kimma (66 contrasts) median | 32.1 s |

**Key findings:**
- geneLME and kimma always agree on contrast direction (100%) — no sign inversion issue
- kimma's `statistic` column is the **negated** t-statistic (opposite sign to estimate)
- geneLME is ~1.8× faster than kimma when running fewer contrasts (3–6); equal speed at 66
- Warning suppression confirmed working — no package version or rescaling warnings in output

**Possible future tasks:**
1. Build runtime scaling estimates: cores × genes × contrasts → expected runtime on this machine

---

## Known Issues / Open Questions

- [x] `geneLME.R` (stable) — interaction contrast support merged 2026-02-20
- [x] `geneLME_dev2.R` — tested (all tests pass) and benchmarked (geneLME_benchmark2.html generated)
- **kimma `statistic` column:** kimma's `$lme.contrast$statistic` is `–(estimate/SE)` — the **negated**
  t-statistic (opposite sign to estimate). To compare against geneLME `t.ratio`, use `-statistic`.
  This was discovered during benchmark2 sign-consistency investigation.
- [ ] `contrast_spec` level ordering with non-alphabetical factors: explicit factor coercion is now
      in place in `geneLME_contrast_spec()` (uses existing factor levels if present, otherwise
      alphabetical). However, `emmeans` must also see the same level ordering — emmeans inherits
      factor levels from the model data, so as long as `dat$targets` columns have matching factor
      levels, alignment is guaranteed. This has not been explicitly tested with a non-alphabetical
      user-defined factor ordering; worth a targeted test if this use case arises.
- [ ] With only 10 patients, `lmer` fits are singular (patient random effect collinear with
      patient-level covariates sex/age). All 10 genes show `isSingular = TRUE` in tests. This is
      expected with small mock N and does NOT affect correctness testing. Real data with larger N
      will not have this issue.
- [ ] Branch B `contrasts_primary` vector construction: users must supply vectors whose length matches
      the number of levels of `contrast_vars[1]` in the order emmeans sees them (alphabetical by
      default). This is implicit and could be a source of user error. Consider adding a runtime check
      that compares vector length to actual emmeans level count.
- [ ] Non-alphabetical factor level ordering: explicit factor coercion is now in place in
      `geneLME_contrast_spec()` (preserves existing factor levels, otherwise alphabetical).
      Alignment with `emmeans` internal ordering has not been explicitly tested with a
      user-defined non-alphabetical factor; worth a targeted test if this arises.

---

## Bugs Fixed

### Session 2026-02-19 — Bug 1: Branch A1 returned 0-row contrast result
**Root cause:** `vector("list", n)` pre-allocates n NULL slots indexed by integer. Subsequent
named assignment (`contrasts_list[["name"]] <- cv`) *appends* rather than filling existing slots,
leaving the first n slots as unnamed NULLs. `contrast()` silently skipped NULLs, leaving a
malformed result that `as.data.frame()` could not coerce.
**Fix:** Build `contrasts_list` as `list()` and assign each element by name from the start.

### Session 2026-02-19 — Bug 2: Branch A2 `"Nonconforming number of contrast coefficients"`
**Root cause:** The `"second"` direction (`var_b | var_a`) was incorrectly passing `contrast_var_2_levels`
keyed to `var_a`, injecting visit-level strings as treatment levels. Additionally, `contrasts_primary`
has length matching var_a's level count (3) but var_b has 4 levels — structurally incompatible.
**Resolution:** Branch A2 (`contrast_by`) was removed entirely per user decision (see Session Log).

---

## Session Log

### Session 15 — 2026-02-20
- **Completed merge of `geneLME_dev2.R` → `geneLME.R`** (picked up where Session 14 left off).
  Remaining edits applied to `geneLME()` function body:
  - Added `contrasts_list <- NULL` and `spec_lookup <- NULL` initialisations.
  - Added `geneLME_build_contrast_args()` pre-compute block (inside `run_contrast`, outside
    `!is.null(contrast_spec)` — called only when `is_interaction`).
  - Updated `plan()` calls: `plan(multisession)` → `future::plan(future::multisession, workers = workers)`;
    `on.exit(plan(sequential))` → `on.exit(future::plan(future::sequential), add = TRUE)`.
  - Updated `geneLME_dispatch()` call: `contrast_spec` → `contrasts_list` + `spec_lookup`.
- **Warning suppression** (from `Claud_next_steps.rtf`):
  - **lmer() rescaling warning** (`"Some predictor variables are on very different scales"`):
    Added `check.scaleX = "ignore"` to a shared `lme_ctrl <- lmerControl(...)` object used
    for both weighted and unweighted `lmer()` calls in `geneLME_fit()`.
  - **Package version warnings** (`"package 'lme4' was built under R version X.Y.Z"`):
    Two-pronged approach:
    1. Removed `future.packages = c(...)` from `future_lapply`; added explicit
       `suppressPackageStartupMessages(suppressWarnings({ library(lme4); ... }))` block at
       the top of `geneLME_fit()` — gives full control over how packages load on workers.
    2. Wrapped `future_lapply(...)` in `withCallingHandlers()` in `geneLME_dispatch()` with
       a handler that muffles any warning matching `"was built under R version"`, catching
       warnings that the future framework re-raises on the main process after collecting worker results.
- **All `geneLME_test.R` tests pass** against merged `geneLME.R` (SOURCE_FILE updated to `"geneLME.R"`):
  Branch A (10/10, 80 ANOVA + 60 contrast), Branch A2 (80 total rows, 60+20 orders),
  Branch B (10/10, 70+60), validation 6a–6f (all PASS), soft-fail 6g (both branches PASS).
  Worker-sourced package-version warnings no longer appear in test output.
- **Updated `geneLME_tutorial.Rmd`:**
  - `source("geneLME_dev.R")` → `source("geneLME.R")`
  - Key capabilities bullet: updated from "per-gene error capture" to `model_status` flagging description
  - Error Handling section: replaced old stop-based singular fit description with new `model_status`
    flagging explanation; updated code to check `!= "success"` for all non-success genes (not just failed)
    and added filter example (`result$lme_contrast %>% filter(model_status == "success")`).
- **Updated `geneLME_function_overview.md`:**
  - Companion file reference: `geneLME_dev.R` → `geneLME.R`; date updated to 2026-02-20
  - Function Map: added `geneLME_build_contrast_args()` row
  - `lme_err` output table entry: updated to describe `"singular_fit"` as a valid value
  - Error Handling section: rewrote to describe `model_status` flagging; removed old
    `"Boundary (singular) fit"` stop-based description
- **Updated `CLAUDE_NOTES_geneLME.md`:**
  - Key files: `geneLME_dev2.R` now "former dev2 (merged); kept for reference"
  - Architecture: `geneLME_build_contrast_args()` added to function table
  - Design decision 5: updated from "error" to "flagged, not errored"
  - Current Status: stage updated to reflect completed merge
- **Re-ran `geneLME_benchmark2.Rmd`** against merged `geneLME.R` only (dev2 correctness section
  removed; report now focuses solely on geneLME vs kimma). Fresh results (Session 15):
  - Sign consistency: 3300/3300 (100%) same direction, r=1.00000000, MAD=1.16e-15
  - Speed: geneLME 1.85× faster than kimma at 3 contrasts, 1.78× at 6, effectively equal (0.99×) at 66
  - Rendered `geneLME_function_overview.html` via pandoc for easy review

### Session 14 — 2026-02-20
- **Confirmed `geneLME.R` (stable) is correct** — does NOT contain dev2 changes (no
  `geneLME_build_contrast_args`, no `singular_fit` flagging); the concern about overwrite was unfounded.
- **Audited `geneLME_dev2.R` for unqualified package calls** — all non-base function calls
  are correctly handled: `future::plan()`, `future::multisession`, `future::sequential`, and
  `future.apply::future_lapply()` are fully qualified in `geneLME()`/`geneLME_dispatch()`;
  worker-side calls (`lmer`, `emmeans`, etc.) rely on `future.packages` which loads all required
  packages on each worker. No bare unqualified calls that would fail outside an interactive session.
- **All `geneLME_test.R` tests pass** against `geneLME_dev2.R`:
  Branch A (10/10 success, 80 ANOVA + 60 contrast rows), Branch A second-order (60+20 rows),
  Branch B (10/10 success, 70+60 rows), all 6 validation error tests (6a–6f), soft-fail (6g).
- **Debugged and fixed `geneLME_benchmark2.Rmd`:**
  1. `sign-inspect` chunk: kimma has no `contrast` string column → rewrote to compare
     `contrast_ref`/`contrast_lvl` pairs directly using canonical pair keys.
  2. Kimma t-statistic: discovered kimma's `statistic` column is the **negated** t-statistic
     (`–estimate/SE`), not the absolute value. Fixed reconstruction to `-statistic`.
  3. Inline `` `r ≈ 1.0` `` parsed as R code → converted to plain text.
  4. Multi-line `if/else` → wrapped in braces to prevent premature statement termination.
  5. `stable_success`: `lme_err` has no names → derived success genes from contrast table instead.
- **`geneLME_benchmark2.html` generated successfully** (44/44 chunks, no errors).
- **Benchmark results:** See Current Status table above. Summary:
  - 100% direction agreement with kimma (no flipped pairs)
  - dev2 estimates identical to stable (r=1.0, MAD=0)
  - dev2 ~1.1–1.2× faster than stable; equal speed to kimma at 66 contrasts

### Session 13 — 2026-02-20
- **Created `geneLME_benchmark2.Rmd`** — v2 benchmarking report covering:
  1. **Sign consistency (Section 1):** Direction-aware join between geneLME and kimma.
     Forward join (same ref/lvl) + flipped join (kimma ref/lvl swapped) combined.
     Scatter plots of estimates and t-statistics before and after direction correction.
     Summary table reports % pairs with direction agreement and corrected accuracy metrics.
  2. **dev2 correctness (Section 2):** geneLME.R vs geneLME_dev2.R on same 50-gene,
     6-contrast problem. Joined on (gene, contrast); reports r, MAD, max|Δ| for estimates.
     Also shows model_status distribution (singular_fit vs success) in dev2 output.
  3. **Speed comparison (Section 3):** Microbenchmark (5 reps each) sourcing stable and
     dev2 sequentially to avoid namespace collision. Bar chart + scaling line plot
     comparing all three methods (stable, dev2, kimma) at 3/6/66 contrasts.
     Summary table includes dev2 speedup ratios vs stable and vs kimma.
- Requires `gridExtra` package (added to libs chunk)

### Session 12 — 2026-02-20
- **Created `geneLME_dev2.R`** with two changes vs `geneLME.R`:
  1. **Singular fit → flag:** `isSingular()` check changed from `stop()` to setting
     `model_status = "singular_fit"`. Results returned for all genes; users filter
     downstream on `model_status`. `model_status` column now appears in both
     `lme_anova` and `lme_contrast` output.
  2. **Pre-computed Branch A contrast structures:** New `geneLME_build_contrast_args()`
     helper called once in `geneLME()` before parallel dispatch. Returns `contrasts_list`
     (named vector list) and `spec_lookup` (ref/lvl join table). Workers receive these
     ready-made; the per-gene `for` loop building contrast vectors is eliminated.
     `geneLME_dispatch()` and `geneLME_fit()` signatures updated accordingly
     (`contrast_spec` parameter replaced by `contrasts_list` + `spec_lookup`).
  - Level ordering in `geneLME_build_contrast_args()` uses the same factor-coercion
    logic as `geneLME_contrast_spec()` (preserves existing factor levels; alphabetical
    otherwise), ensuring alignment with emmeans' internal ordering without requiring
    a fitted model object.
- **Updated `geneLME_test.R`:**
  - Added `SOURCE_FILE` variable at top for easy switching between implementations
  - Test 4: added singular-fit verification block (checks non-NA estimates for flagged genes)
  - Tests 4 and 5: updated output display to show `model_status` column
  - Note: with 10-patient mock data all genes are expected to be `singular_fit`, not errors

### Session 11 — 2026-02-20
- Merged `geneLME_dev.R` → `geneLME.R` as the new stable version
- Added merge changelog to `geneLME.R` header
- Removed stale "PICKUP HERE" comment and dev-only header label
- Updated CLAUDE_NOTES key files and current status sections

### Session 1 — 2026-02-19 (first session)
- Identified correct active worktree (`reverent-satoshi`)
- Created `CLAUDE_NOTES_geneLME.md` (this file)
- Created `geneLME_test.R` with mock EList:
  - 10 patients × 3 treatments × 4 visits = 120 samples, 50 genes
  - Patient-level covariates: `sex` (factor), `age` (continuous)
  - Sample-level technical covariates: `rNANgUl`, `percent_duplication`,
    `median_cv_coverage`, `lib.size`, `norm.factors`
- Ran full test suite; identified and fixed two bugs in `geneLME_dev.R`
- All three contrast branches (A1, A2, B) passed initially after fixes

### Session 10 — 2026-02-19 (tenth context window)
- **Created `geneLME_benchmark.Rmd`** (→ `geneLME_benchmark.html`): kimma vs geneLME benchmarking report
  - **Mock data:** 10 patients × 3 treatments × 4 visits = 120 samples, 2000 genes; genes 1–100
    have +2.5 log2 TrtC:V3 effect; `libID = sample_id` added to targets for kimma compatibility
  - **Contrast subsets defined:**
    - `spec_3` (3): longitudinal V2→V3 within-treatment
    - `spec_6` (6): between-treatment within-visit at V2 and V3
    - `spec_full` (66): all pairwise — equivalent to kimma's default output
  - **Section 1 (Accuracy):** both methods on 50-gene subset with all 66 contrasts.
    Results joined on `(gene, contrast_ref, contrast_lvl)` — both methods now have these columns.
    Reports Pearson r and MAD for estimates and t-statistics; scatter plots included.
    Note: geneLME excludes singular-fit genes (treated as error); kimma returns estimates silently.
  - **Section 2 (Selective contrast efficiency):** `microbenchmark(times=5)` on 2000-gene full
    dataset comparing geneLME with 3/6/66 contrasts and kimma with 66 contrasts. Bar chart +
    scaling plot showing linear runtime scaling with n_contrasts for geneLME.
  - **Section 3 (Head-to-head):** reuses 66-contrast runs from Section 2 to compare
    geneLME vs kimma directly at equal contrast count. Reports speedup ratio.
  - **kimma interface notes:**
    - `libraryID = "libID"` required — added `libID` column to targets
    - `contrast_var = "treatment:visit"` runs all 66 pairwise automatically
    - Output: `$lme.contrast` with columns `gene, contrast_ref, contrast_lvl, estimate, statistic, pval, FDR, std.error, df`
    - `statistic` = t-ratio (same as geneLME `t.ratio`); `pval` = p-value; `std.error` = SE
  - **Key finding expected:** estimates identical (r ≈ 1.0, MAD ≈ 0) for shared non-singular genes;
    geneLME faster when running fewer contrasts; head-to-head speed depends on parallelism overhead

### Session 9 — 2026-02-19 (ninth context window)
- **Added `contrast_ref` and `contrast_lvl` to `lme_contrast` output** in `geneLME_fit()`:
  - Branch A first-order: built a `spec_lookup` data.frame keyed on the contrast name string
    (`paste(contrast_lvl, contrast_ref, sep = " - ")`), then `left_join`-ed onto the emmeans
    first-order result by `"contrast"`. This ensures `contrast_ref` = the −1 cell and
    `contrast_lvl` = the +1 cell, eliminating sign ambiguity.
  - Branch A second-order: `contrast_ref = NA_character_, contrast_lvl = NA_character_`
    (no single ref/lvl pair applies to a contrast-of-contrasts).
  - Branch B (all rows): `contrast_ref = NA_character_, contrast_lvl = NA_character_`
    (Branch B uses named coefficient vectors, not a ref/lvl spec).
  - Error stub data.frame in `tryCatch` handler: added both NA columns for schema consistency.
- **Tutorial Rmd updated:** added `contrast_ref`/`contrast_lvl` to `select()` calls in all
  contrast output display chunks; added explanatory prose in the `lme_contrast` description.
- **CLAUDE_NOTES updated:** added columns to return value documentation.
- Next up: kimma benchmarking (2000-gene mock, speed + estimate comparison).

### Session 8 — 2026-02-19 (eighth context window)
- **Redesigned `contrast_index` system** — root cause was that `contrast_index` in the
  `geneLME_contrast_spec()` output was being misused: users built `contrasts_secondary` vectors
  using the full template's row count as the vector length (234 elements) instead of the filtered
  spec's row count. `which(my_spec$contrast_index == ...)` correctly identified positions within
  `my_spec`, but the outer `rep(0, nrow(full_template))` produced vectors of the wrong length.
- **`geneLME_contrast_spec()` output simplified:** removed `contrast_index` column entirely.
  Output is now a two-column data.frame (`contrast_ref`, `contrast_lvl`) only.
- **`geneLME()` now appends `$contrast_spec` to its return value:**
  - Branch A: `contrast_spec %>% mutate(contrast_index = seq_len(n())) %>% select(contrast_index, everything())`
    — `contrast_index` is `1:nrow(contrast_spec)` (row positions within the *filtered* spec the user passed in)
  - Branch B: `data.frame(contrast_index = seq_along(contrasts_primary), contrast_name = names(contrasts_primary))`
  - NULL when no contrasts are run
- **Soft-fail on wrong-length `contrasts_secondary`** (changed from hard `stop()`):
  - Emits an informative `message()` naming offending vectors, observed length, expected length
  - Returns `invisible(list(lme_anova=NULL, lme_contrast=NULL, lme_fit=NULL, lme_err=NULL, contrast_spec=indexed_contrast_spec))`
  - User can inspect `result$contrast_spec` to fix their vectors without running any models
- **`geneLME_compiler()` updated:** now accepts `contrast_spec = NULL` argument; includes it in
  return list; added column-existence guard for no-contrast runs (prevented crash on empty contrast tibble)
- **Tutorial Rmd updated:**
  - `contrast-spec-interaction` chunk: added comment explaining `contrast_index` is NOT in
    `spec_template`, is added by `geneLME()` to `$contrast_spec`
  - Added `$contrast_spec` display to the Branch A output structure chunk
  - Branch A second-order chunk: updated comments to reference `result_A$contrast_spec` for
    row-order verification
  - Programmatic secondary contrast section: rewritten to use `result_A$contrast_spec` (the
    indexed spec from `geneLME()`) instead of `my_spec$contrast_index`
  - Added soft-fail validation test (Test 7) to Input Validation section
  - Quick-reference table: updated `contrast_spec` description; updated `contrasts_secondary`
    description to mention soft-fail
  - Key capabilities: added soft-fail bullet
- **`geneLME_test.R` updated:**
  - Test 3: removed `contrast_index` range check; updated comment to explain columns
  - Test 4b: updated comments to reference `test_A$contrast_spec` for indexed ordering;
    corrected `contrasts_secondary` vectors to match actual alphabetical row ordering
    (V2 rows 1–3 before V3 rows 4–6); added `print(test_A$contrast_spec)` before test_A2 call
  - Added test 6g: soft-fail validation for both Branch A and Branch B wrong-length vectors

### Session 7 — 2026-02-19 (seventh context window)
- **Root-caused the real-data NA/warning failure:**
  - User shared real-data `contrasts_secondary` vectors: length 234 (= nrow of full
    unfiltered spec_template). The filtered `contrast_spec` passed to `geneLME()` had far
    fewer rows. emmeans throws "Nonconforming number of contrast coefficients" which is
    silently caught by `tryCatch` inside `geneLME_fit()`, producing NA-filled rows.
  - The programmatic builder from Session 4 used `rep(0, nrow(defined_contrasts_spec))`
    where `defined_contrasts_spec` referred to the full template — should be
    `nrow(filtered_spec)`. The `which(my_spec$contrast_index == ...)` index lookup was
    correct; only the vector length was wrong.
- **Fix: added `contrasts_secondary` length validation to `geneLME()` pre-flight:**
  - Branch A (interaction): checked that every vector in `contrasts_secondary` has length
    == `nrow(contrast_spec)` (the filtered spec passed in, NOT the full template).
    Error message names the offending vectors, shows observed vs expected length, and
    explains the common cause (full template vs filtered spec confusion).
  - Branch B (non-interaction): checked that every vector has length ==
    `length(contrasts_primary)`.
  - Both checks fire before any parallel work launches — so the user gets a clear,
    immediate error instead of silent NAs in results.
- User's real-data fix: rebuild `contrasts_secondary` vectors using
  `rep(0, nrow(defined_contrasts_geneLME))` (filtered spec row count) and
  `which(defined_contrasts_geneLME$contrast_index == ...)` for index lookup.

### Session 6 — 2026-02-19 (sixth context window)
- **Diagnosed and fixed second-order contrast warning + spurious NA issue:**
  - Symptom: rescaling warnings ("Some predictor variables are on very different scales")
    and occasional NA-filled outputs when running with `contrasts_secondary`.
  - Root cause (warnings): `emmeans::contrast()` called on a first-order contrast object
    internally re-invokes lme4 machinery on the first-order estimate/SE values, which can be
    on different scales from the original predictors → benign lme4 scale warning, no effect
    on results. Confirmed output is correct by single-gene isolated test.
  - Root cause (NAs): not a code bug — was a stochastic singular-fit outcome in the small
    mock data (10 patients), not reproducible when re-run. `isSingular()` check is correct.
  - Fix 1: wrapped second-order `contrast()` calls in `suppressWarnings()` in **both** Branch A
    and Branch B of `geneLME_fit()`, with explanatory comment. Warning is cosmetic only.
  - Fix 2 (latent bug): both Branch A and Branch B were gating second-order contrasts with
    `length(contrasts_list) > 1` / `length(contrasts_primary) > 1`. This incorrectly
    prevented second-order contrasts when only one first-order contrast existed. Changed
    both to simply `!is.null(contrasts_secondary)` — the user's provision of
    `contrasts_secondary` is the correct and sufficient gate.
- Verified fix: all 10 genes succeed, 60 first-order + 20 second-order rows, all
  `p.value_adj` non-NA, zero warnings emitted.

### Session 5 — 2026-02-19 (fifth context window)
- Added FDR-adjusted p-values (`p.value_adj`) to both `lme_anova` and `lme_contrast` outputs:
  - **`geneLME_compiler()`**: now accepts `fdr_method = "BH"` argument; applies `p.adjust()`
    within grouped sets after binding all per-gene results:
    - `lme_anova`: grouped by `term` — each model term adjusted independently across genes
    - `lme_contrast`: grouped by `contrast × contrast_order` — each contrast (including
      first- and second-order separately) adjusted independently across genes
  - **`geneLME()`**: added `fdr_method = "BH"` argument (any `p.adjust.methods` value valid);
    validated against `p.adjust.methods` in pre-flight checks; passed through to compiler
  - NA p-values (failed genes) propagate to `p.value_adj` as NA — they are excluded from the
    effective adjustment set automatically by `p.adjust()`
- Updated `geneLME_tutorial.Rmd`:
  - Added `p.value_adj` to `select()` calls in ANOVA and contrast output display chunks
  - Updated descriptions of `lme_anova` and `lme_contrast` to explain grouping logic
  - Added new "FDR Adjustment" section (before Input Validation) with grouping table,
    method list, and BH vs Bonferroni comparison example
  - Added `fdr_method` row to quick-reference argument table
- CLAUDE_NOTES: updated function signatures (added `geneLME_compiler` signature, `fdr_method`
  arg to `geneLME`), validation check list (item 9), and current status

### Session 4 — 2026-02-19 (fourth context window)
- User provided a real-data pattern for programmatic secondary contrast construction using
  `contrast_index`. Reviewed approach and identified improvements:
  - Replaced `group_by` + `reframe` with `last()`/`first()` with a more robust `summarise`
    producing explicit `index_neg`/`index_pos` columns (two-column wide shape)
  - Removed unused `nest` + `rowwise` + dead `secondary_contrast` list column
  - Replaced opaque `[[i]][1,1][[1]]` nested indexing with `which(my_spec$contrast_index == ...)`
    for transparent index-to-row-position mapping
  - Replaced `for` loop with `setNames(lapply(...))` pattern consistent with how
    `contrasts_primary` is specified elsewhere
- Added new tutorial section "Programmatic second-order contrast construction via `contrast_index`"
  to `geneLME_tutorial.Rmd` (between Branch A second-order and Branch B sections)
- Updated quick-reference table to include `contrast_index` in `contrast_spec` column description

### Session 3 — 2026-02-19 (third context window)
- Three improvements to `geneLME_contrast_spec()` interaction mode:
  1. **Factor-level ordering enforced**: both `var_a` and `var_b` are explicitly coerced to
     factors before calling `interaction()`. Existing factor levels are preserved; plain character
     columns get alphabetical ordering imposed. This guarantees `contrast_ref` component levels ≤
     `contrast_lvl` component levels (for both variables), consistent with emmeans' internal ordering.
  2. **Swap tolerance documented**: `geneLME_fit()` Branch A already handled swapped ref/lvl
     gracefully (sign flip, no error). Added explicit comment in the contrast vector construction
     loop documenting this behaviour.
  3. **`contrast_index` column added**: `geneLME_contrast_spec()` interaction mode now returns
     a `contrast_index` column (integer, 1-based) as the first column. This index is stable — it
     reflects the row's position in the *full unfiltered template* and survives user `filter()` calls.
     Intended for downstream use (e.g. joining back to results, building `contrasts_secondary`
     with explicit index reference). Tutorial example to be added in a future session.
- Updated `geneLME_test.R`:
  - Added `contrast_index` range check in test 3
  - Updated test 4b comments to reflect corrected row ordering (alphabetical interaction levels:
    V2 rows before V3 rows) and corrected `contrasts_secondary` vectors accordingly
- Both `geneLME_dev.R` and `geneLME.R` preserved separately (user preference)

### Session 2 — 2026-02-19 (same day, second context window)
- Reviewed all three branches with user; demonstrated outputs
- **User decision:** eliminate Branch A2 (`contrast_by`) entirely
  - Rationale: specifying contrasts in both directions is error-prone; pairwise contrasts
    for the second variable can be computationally excessive with many levels
- **Design changes implemented:**
  1. Extended `geneLME_contrast_spec()` to two modes:
     - Interaction mode: returns `contrast_ref`/`contrast_lvl` pairs (unchanged)
     - Single-variable mode: returns `level` reference frame (new)
  2. Added formula-vs-`contrast_vars` consistency check using `attr(terms(), "term.labels")`:
     - Interaction contrast on additive model → hard error
     - Non-interaction contrast var absent as main effect → warning
  3. Made `contrast_spec` required (not optional) when `contrast_vars` contains `":"`
  4. Removed `contrast_by` parameter from all four functions
  5. Removed Branch A2 from `geneLME_fit()`; Branch A is now solely the `contrast_spec` approach
  6. Updated `geneLME_test.R`: removed A2 test; added single-variable `geneLME_contrast_spec()`
     tests; added 4 new validation error tests (6c–6f)
- Final test run: all 8 tests pass (2 positive + 6 validation errors)
  - Branch A: 10/10 genes `success`, 80 ANOVA rows, 60 contrast rows
  - Branch B: 10/10 genes `success`, 70 ANOVA rows, 60 contrast rows
  - Validation 6a–6f: all produce expected informative errors
