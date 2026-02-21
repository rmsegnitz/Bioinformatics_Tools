# geneLME — Function Overview

> **Companion to:** `geneLME.R`
> **Last updated:** 2026-02-20
> For session continuity notes see `CLAUDE_NOTES_geneLME.md`.
> For a worked tutorial with output see `geneLME_tutorial.html` (knitted from `geneLME_tutorial.Rmd`).

---

## Function Map

| Function | Role | Called by |
|---|---|---|
| `geneLME_contrast_spec()` | Pre-run helper: inspect available contrast levels and understand how to specify arguments | User |
| `geneLME_build_contrast_args()` | Branch A helper: pre-computes `contrasts_list` + `spec_lookup` once before parallel dispatch | `geneLME()` |
| `geneLME()` | User-facing entry point: validates inputs, slices data, sets parallel plan | User |
| `geneLME_dispatch()` | Launches `future_lapply` with explicit globals; clean parallel scope | `geneLME()` |
| `geneLME_fit()` | Per-gene worker: fits `lmer`, extracts ANOVA + contrasts | `geneLME_dispatch()` |
| `geneLME_compiler()` | Binds per-gene list results into four output data frames | `geneLME()` |

---

## `geneLME_contrast_spec(targets, contrast_vars)`

A pre-run helper to enumerate available levels and understand how to construct contrast arguments before calling `geneLME()`. Two modes:

### Interaction mode

`contrast_vars` is a single `"var_a:var_b"` string. Returns a `data.frame(contrast_ref, contrast_lvl)` with all pairwise combinations of interaction cells. The user filters this to their contrasts of interest and passes it directly as `contrast_spec` to `geneLME()`.

```r
spec_template <- geneLME_contrast_spec(dat$targets, contrast_vars = "treatment:visit")
# → data.frame with 66 rows (all pairs of 12 interaction cells)
#    contrast_ref contrast_lvl
# 1      TrtA V1      TrtB V1
# 2      TrtA V1      TrtC V1
# ...

my_spec <- spec_template %>%
  filter(sub(".* ", "", contrast_ref) == sub(".* ", "", contrast_lvl),  # same visit
         sub(".* ", "", contrast_ref) %in% c("V2", "V3"))               # only V2, V3
# → 6 rows: the 3 cross-treatment pairs at V2, and 3 at V3
```

### Single / multi-variable mode

`contrast_vars` is a character vector of one or more plain variable names (no `:`). Returns a **named list**, one `data.frame(level)` per variable. A message for each variable explains its role in subsequent `geneLME()` arguments. This list is a reference only — not passed to `geneLME()`.

```r
ref <- geneLME_contrast_spec(dat$targets, contrast_vars = c("treatment", "visit"))

# Message for treatment (position 1 — primary):
# 'treatment' — primary contrast variable (contrast_vars[1] in geneLME).
# 3 levels (alphabetical order = position order for contrasts_primary vectors):
# 1. TrtA   2. TrtB   3. TrtC
# → contrast vectors passed to contrasts_primary must have length 3
#   e.g. 'TrtC vs TrtA' = c(-1, 0, 1)

# Message for visit (position 2 — secondary 'by' variable):
# 'visit' — secondary 'by' grouping variable (contrast_vars[2] in geneLME).
# 4 levels: V1, V2, V3, V4
# → pass a subset to contrast_var_2_levels in geneLME() to restrict
#   which groups the primary contrasts are computed within.

ref$treatment   # data.frame: TrtA, TrtB, TrtC
ref$visit       # data.frame: V1, V2, V3, V4
```

**Note:** mixing interaction and plain variable names in one call raises an error — call separately.

---

## `geneLME()` — Main Function

### Signature

```r
geneLME(dat, formula_str,
        model_weights         = NULL,
        run_contrast          = NULL,
        contrast_vars         = NULL,
        contrast_var_2_levels = NULL,
        contrast_spec         = NULL,
        contrasts_primary     = NULL,
        contrasts_secondary   = NULL,
        fdr_method            = "BH",   # any method accepted by p.adjust()
        n_cores               = NULL)
```

### Input

`dat` — EList-like list with three elements:

| Element | Type | Description |
|---|---|---|
| `dat$E` | matrix (genes × samples) | log2 expression values |
| `dat$weights` | matrix (genes × samples), optional | voom precision weights |
| `dat$targets` | data.frame (samples × covariates) | sample metadata; all formula and contrast variables must be columns here |

### Output

Named list of five elements:

| Element | Contents |
|---|---|
| `lme_anova` | One row per model term per gene: `term`, `statistic`, `df`, `p.value`, `p.value_adj` (FDR within term), `gene`, `model_status`, `predictor_class`, `Estimate`, `Estimate_SE` |
| `lme_contrast` | One row per contrast per gene: `contrast`, `contrast_ref`, `contrast_lvl` (Branch A first-order only; NA otherwise), `estimate`, `SE`, `df`, `t.ratio`, `p.value`, `p.value_adj` (FDR within contrast × order), `contrast_order` (`"first_order"` / `"second_order"`), `gene`. Branch B also has a `visit` column. |
| `lme_fit` | One row per gene: `gene`, `AIC` |
| `lme_err` | Named character vector: gene → `"success"`, `"singular_fit"`, or unexpected error message. A `model_status` column carrying the same value also appears in `lme_anova` and `lme_contrast`. |
| `contrast_spec` | Indexed copy of the filtered `contrast_spec` (Branch A) or `data.frame(contrast_index, contrast_name)` from `contrasts_primary` names (Branch B); `NULL` when no contrasts run. On soft-fail (wrong-length `contrasts_secondary`), all other elements are NULL and only this is returned. |

### Validation Checks (pre-flight)

1. All formula variables present in `dat$targets`
2. `dat$weights` present and dimension-matched if `model_weights = TRUE`
3. `contrast_vars` not NULL when `run_contrast = TRUE`
4. **Interaction term present in formula** if `contrast_vars` contains `:` (hard error — `emmeans` would otherwise silently use additive margins, which is statistically misleading)
5. `contrast_spec` required (not NULL) when `contrast_vars` contains `:`
6. `contrast_spec` has correct columns (`contrast_ref`, `contrast_lvl`)
7. All `contrast_vars` component variables present in `dat$targets`
8. `contrast_var_2_levels` values valid against actual data levels
9. `fdr_method` validated against `p.adjust.methods`
10. **`contrasts_secondary` vector lengths** (soft-fail, not hard stop): each vector must have length == `nrow(contrast_spec)` (Branch A) or `length(contrasts_primary)` (Branch B). On mismatch returns early with only `$contrast_spec` populated.
11. Warn on likely ID columns used as predictors

---

## Branch A — Interaction Contrasts

**When to use:** model formula contains an interaction term (`treatment * visit` or `treatment:visit`) and you want specific pairwise contrasts across particular interaction cells.

### Step 1 — Generate and filter the contrast template

```r
spec_template <- geneLME_contrast_spec(dat$targets, "treatment:visit")
# Returns data.frame of all 66 pairwise interaction cell combinations

my_spec <- spec_template %>%
  filter(sub(".* ", "", contrast_ref) == sub(".* ", "", contrast_lvl),  # same visit only
         sub(".* ", "", contrast_ref) %in% c("V2", "V3"))               # visits V2 and V3 only
# 6 rows
```

### Step 2 — Run

```r
result_A <- geneLME(
  dat           = dat,
  formula_str   = "~ treatment * visit + age + sex + rNANgUl + (1|ptID)",
  model_weights = TRUE,
  run_contrast  = TRUE,
  contrast_vars = "treatment:visit",   # single ":" string → triggers Branch A
  contrast_spec = my_spec,             # required; each row defines one contrast
  n_cores       = 10
)
```

### Sample `lme_contrast` output (6 contrasts × n genes rows)

`contrast_ref` = the −1 cell; `contrast_lvl` = the +1 cell. Together they make the sign of
`estimate` unambiguous without parsing the contrast label string.

```
            contrast   contrast_ref  contrast_lvl  estimate    SE      df   t.ratio  p.value  contrast_order  gene
TrtB V2 - TrtA V2    TrtA V2        TrtB V2        -0.040    0.960  93.04  -0.042    0.967    first_order     gene01
TrtC V2 - TrtA V2    TrtA V2        TrtC V2         0.092    0.942  97.37   0.098    0.922    first_order     gene01
TrtC V2 - TrtB V2    TrtB V2        TrtC V2         0.132    0.993  89.10   0.133    0.894    first_order     gene01
TrtB V3 - TrtA V3    TrtA V3        TrtB V3         0.191    0.934  98.51   0.205    0.838    first_order     gene01
TrtC V3 - TrtA V3    TrtA V3        TrtC V3         1.625    0.945  90.82   1.721    0.089    first_order     gene01
TrtC V3 - TrtB V3    TrtB V3        TrtC V3         1.434    0.958  90.33   1.497    0.138    first_order     gene01
```

---

## Branch A with Second-Order Contrasts

At validation time, `geneLME()` prints the numbered first-order list (in `contrast_spec` row order) to guide vector construction. Secondary vectors must have length `nrow(contrast_spec)`.

```r
result_A2 <- geneLME(
  dat                 = dat,
  formula_str         = "~ treatment * visit + age + sex + rNANgUl + (1|ptID)",
  model_weights       = TRUE,
  run_contrast        = TRUE,
  contrast_vars       = "treatment:visit",
  contrast_spec       = my_spec,          # 6 rows, order defines vector positions
  contrasts_secondary = list(
    # my_spec rows: 1=TrtB V2-TrtA V2, 2=TrtC V2-TrtA V2, 3=TrtC V2-TrtB V2,
    #               4=TrtB V3-TrtA V3, 5=TrtC V3-TrtA V3, 6=TrtC V3-TrtB V3
    "TrtA vs TrtB: V3 vs V2" = c(1, 0, -1, 0, 0, 0),   # row1 − row3: V2 vs V3 delta
    "TrtA vs TrtC: V3 vs V2" = c(0, 1,  0,-1, 0, 0)    # row2 − row4
  ),
  n_cores             = 10
)
# lme_contrast: first_order (6 × n_genes) + second_order (2 × n_genes) rows
```

**Validation message printed:**
```
contrasts_secondary will be applied to the first-order interaction contrasts
in the order they appear in contrast_spec (6 contrasts):
1. TrtB V2 - TrtA V2
2. TrtC V2 - TrtA V2
3. TrtC V2 - TrtB V2
4. TrtB V3 - TrtA V3
5. TrtC V3 - TrtA V3
6. TrtC V3 - TrtB V3
Ensure contrasts_secondary vectors have length 6, with each element
corresponding to the contrast at that position.
```

---

## Branch B — Non-Interaction Contrasts

**When to use:** model is additive (no interaction), and you want named contrasts on the marginal means of one variable, optionally evaluated within specific levels of a second variable.

### Step 1 — Inspect levels

```r
ref <- geneLME_contrast_spec(dat$targets, contrast_vars = c("treatment", "visit"))
# treatment levels (alphabetical): TrtA[1], TrtB[2], TrtC[3]
# → contrasts_primary vectors must have length 3, positions = [TrtA, TrtB, TrtC]
# visit levels: V1, V2, V3, V4 → pass subset to contrast_var_2_levels
```

### Step 2 — Run

```r
result_B <- geneLME(
  dat                   = dat,
  formula_str           = "~ treatment + visit + age + sex + rNANgUl + (1|ptID)",
  model_weights         = TRUE,
  run_contrast          = TRUE,
  contrast_vars         = c("treatment", "visit"), # [1] = primary; [2] = by-variable
  contrast_var_2_levels = c("V2", "V3"),           # restrict output to these visit levels
  contrasts_primary     = list(
    "TrtC vs TrtA" = c(-1, 0, 1),  # positions: [TrtA=−1, TrtB=0, TrtC=+1]
    "TrtB vs TrtA" = c(-1, 1, 0)
  ),
  contrasts_secondary   = list(
    "TrtC vs TrtB" = c(1, -1)      # length = number of primary contrasts (2)
  ),
  n_cores               = 10
)
```

### Sample `lme_contrast` output (per gene: 2 primary × 2 visits + 1 secondary × 2 visits = 6 rows)

```
       contrast  visit   estimate    SE       df   t.ratio  p.value  contrast_order  gene
TrtC vs TrtA     V2      0.814    0.457  101.6   1.782     0.078    first_order     gene01
TrtB vs TrtA     V2      0.009    0.457  102.3   0.019     0.984    first_order     gene01
TrtC vs TrtA     V3      0.814    0.457  101.6   1.782     0.078    first_order     gene01
TrtB vs TrtA     V3      0.009    0.457  102.3   0.019     0.984    first_order     gene01
TrtC vs TrtB     V2      0.805    0.462   99.1   1.742     0.085    second_order    gene01
TrtC vs TrtB     V3      0.805    0.462   99.1   1.742     0.085    second_order    gene01
```

> **Note:** In an additive model, estimates are identical across V2 and V3 — `contrast_var_2_levels` restricts *which* visit cells `emmeans` evaluates at, but since there is no interaction the treatment effect is constant across visits. This filter is most useful to control the number of output rows and limit multiple testing burden.

---

## `lme_anova` Output Structure

One row per model term per gene. Key columns:

| Column | Description |
|---|---|
| `term` | Predictor name (e.g. `"treatment"`, `"treatment:visit"`, `"age"`) |
| `statistic` | Chi-square statistic from `car::Anova()` |
| `df` | Degrees of freedom |
| `p.value` | p-value |
| `gene` | Gene identifier |
| `model_status` | `"success"` or error message |
| `predictor_class` | `"continuous"`, `"two-level-categorical"`, `"multi-level-categorical"`, or `"interaction"` |
| `Estimate_source` | `"lme_summary"` (value in `Estimate`/`Estimate_SE`) or `"seeContrasts"` (NA — use contrast output for estimate) |
| `Estimate` | Coefficient estimate where available; NA for multi-level or complex interactions |
| `Estimate_SE` | Standard error of estimate where available |

---

## Error Handling and Singular Fits

Per-gene outcomes are tracked via a `model_status` column in both `lme_anova` and `lme_contrast`, and summarised in `lme_err` (named character vector, gene → status). Possible values:

- **`"success"`** — model converged cleanly; all results are reliable.
- **`"singular_fit"`** — `isSingular()` was `TRUE`; random effect variance hit its boundary (zero). Fixed effect estimates and contrasts are still returned and are numerically valid. Common with small N or highly collinear random/fixed effects. Filter downstream if desired.
- **Unexpected error string** — `tryCatch` caught an unexpected error; rows for that gene contain NAs.

Note: singular fits were previously treated as a hard `stop()` that excluded the gene entirely. They are now flagged and returned, giving users the choice to retain or filter.

```r
table(result$lme_err)
non_success <- names(result$lme_err)[result$lme_err != "success"]

# Filter to clean fits only downstream:
result$lme_contrast %>% filter(model_status == "success")
```
