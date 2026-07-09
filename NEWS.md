# coresynth (development version)

## New features

- **In-space placebo visualization** (Abadie, Diamond & Hainmueller 2010,
  Section 3.4): `mspe_ratio_pval()` now returns an `scm_placebo` object
  (still fully backward compatible with the previous plain-list fields) that
  additionally carries the placebo gap path for every donor unit. A new
  `plot.scm_placebo()` method renders the two companion figures from the
  paper: `type = "gaps"` overlays the treated unit's gap on the donor-pool
  placebo gaps (Figures 4-7, with `mspe_prune` implementing the paper's
  relative pre-treatment MSPE pruning at 20x/5x/2x), and `type = "ratios"`
  plots the post/pre-treatment MSPE ratio for every unit (Figure 8).
- **Plot style customization**: `plot.coresynth()` and `plot.scm_placebo()`
  gain `colors`, `vline`, `hline`, and (for `type = "weights"`) `fill`
  arguments. `vline`/`hline` accept a list of `geom_vline()`/`geom_hline()`
  aesthetic overrides merged onto the built-in defaults, or `NULL`/`FALSE` to
  suppress the reference line entirely — since a line already added to a
  returned `ggplot` object cannot be removed, only restyled or overplotted,
  suppression has to happen inside the plot method. `colors` accepts a named
  vector overriding individual series colors (unmentioned series keep their
  default). All defaults reproduce the previous appearance exactly, so
  existing calls are unaffected.

# coresynth 0.2.4

## Bug fixes

- `method = "mc"` treated missing `(id, time)` panel cells and `NA` outcomes as
  observed zeros instead of excluding them from the observation mask used by
  the Soft-Impute solver. Missing cells are now masked out the same way as
  treated post-adoption cells. **This changes numerical results for `mc` fits
  on unbalanced panels or panels with `NA` outcomes.**
- `method = "tasc"`'s EM loop already handles missing outcome cells (Kalman
  smoother plus per-unit M-steps), but its initialisation (`svd()`, the
  treated-unit loading OLS, and `var()`) failed on panels with missing cells.
  Initial values are now computed from a column-mean-imputed matrix; the EM
  loop itself still runs on the true, unimputed data.

## Input validation

Malformed panels previously produced confusing or misleading errors — for
example, a 0-row data frame (from an upstream filtering bug) was misclassified
as staggered adoption and failed with an unrelated "predictors not supported"
error, and missing `(id, time)` cells surfaced as a raw
`eig_sym(): decomposition failed` from the C++ layer. `panel_to_matrices()`
(shared by all six estimators) and `scm_design()` now validate and error
clearly on:

- 0-row input data
- `NA` unit or time identifiers
- `NA` or negative treatment indicator values
- no treated units, or (for sharp fits only) no control units — staggered
  fits may legitimately have no never-treated units when future adopters
  serve as clean controls
- missing `(id, time)` cells or non-finite outcomes, for the estimators that
  require a fully observed panel (`scm`, `sdid`, `gsc`, `si` — `mc` and
  `tasc` handle missing data by design)

`scm_fit()` also now rejects non-numeric (`factor`/`character`) outcome or
treatment columns and non-integer treatment values, instead of silently
coercing them (`as.integer(factor(...))` returns level codes, not the original
values). "All cohort-level fits failed" errors (SCM/SDID/GSC/SI staggered
paths) now point back to the preceding per-cohort warnings for diagnosis.

# coresynth 0.2.3

## Performance

- The inner simplex QP solvers (`solve_simplex_qp()` / `solve_simplex_qp_lr()`),
  used by every `scm_fit(method = "scm")` fit and by the `method = "sdid"` unit
  weights, now use FISTA with adaptive restart (O'Donoghue & Candès 2015,
  gradient scheme). The first-order solver's iteration count grows with the
  condition number of `Q = X0'VX0`, which is large for the collinear
  pre-treatment outcome panels typical of synthetic control; resetting the
  momentum term when it works against the gradient removes this slowdown. Each
  inner solve converges to the same optimum (verified against an exact QP solver
  to within `1e-10`), so the returned weights are unchanged.
- For the common outcomes-only case the fit is effectively identical
  (differences `~1e-5`). For `v_selection = "oos"` and penalised (`lambda_pen`)
  fits on ill-conditioned panels, the non-convex outer V search may now settle on
  a different local optimum, shifting results slightly; both the previous and the
  new solutions are valid SCM fits with the same objective. **This can change
  numerical results for `v_selection = "oos"` and penalised fits on poorly
  conditioned data.**

# coresynth 0.2.2

## Bug fixes

- `panel_to_matrices()` (and therefore every `scm_fit()` method) and
  `scm_design()` now error on duplicate `(id, time)` entries instead of silently
  keeping a single arbitrary row. A balanced panel requires each unit-time cell
  to be unique; duplicates were previously overwritten by the last matrix-index
  assignment, dropping data without warning. The error reports the number of
  offending rows and the first duplicated unit and time.
- `plot()` now renders `Date`/`POSIXct` time axes correctly. The time vector was
  coerced with `as.numeric()`, so dates appeared as days-since-epoch (e.g. 16000,
  20000) on the x-axis. `Date`/`POSIXct` values are now passed through unchanged
  so ggplot2 selects the appropriate date scale; only `character`/`factor` time
  values are coerced to numeric.

# coresynth 0.2.1

## Bug fixes

- `v_selection = "oos"` (outcomes-only case) previously fit candidate `W(V)` on
  the full pre-treatment outcome matrix and restricted only the MSPE evaluation
  to the validation window, allowing the V optimiser to fit the validation
  period indirectly (a data leak relative to Abadie (2021) S.3.2). The new
  `.scm_oos_outcomes()` implements the correct train/validation split:
  candidate `W(V)` are fitted on training-half outcomes only, `V*` minimises
  validation-half MSPE, and `W*` is refit with `V*` on the outcomes of the last
  `floor(T_pre/2)` pre-treatment periods. For OOS fits, `v_weights` now has
  `floor(T_pre/2)` entries and a new `v_rows` field records which periods they
  refer to. **This changes numerical results for `v_selection = "oos"`.**

## New features

- `scale_predictors` (default `TRUE`): predictor rows supplied via
  `predictors =` are now divided by their standard deviation across all units
  before optimisation, matching the Synth reference implementation (Abadie,
  Diamond & Hainmueller 2011, JSS). `predictor_table` continues to report
  values on the original scale. **This changes numerical results for SCM fits
  with user-supplied `predictors`** unless `scale_predictors = FALSE`.
- `placebo_in_time()`: in-time placebo (backdating) test for sharp SCM fits
  (Abadie, Diamond & Hainmueller 2015; Abadie & Vives-i-Bastida 2022).
- `loo_donors()`: leave-one-out donor robustness check with the predictor
  weights V held fixed (Abadie, Diamond & Hainmueller 2015, footnote 20).
- `build_predictor_matrices()` now errors with an informative message if a
  `pred()` time window produces missing or non-finite predictor values.

# coresynth 0.2.0

## New features

- **Conformal inference** (`conformal_inference()`): permutation-based p-values and
  confidence intervals following Chernozhukov, Wüthrich & Zhu (2021).
  Works with sharp fits across all supported estimation methods
  (`scm`, `sdid`, `gsc`, `mc`, `si`).  The counterfactual proxy is re-estimated
  under the null on all *T* periods (essential for finite-sample validity per
  CWZ S.2.2), and p-values are obtained via moving-block (cyclic-shift) permutation
  of the estimated residuals.  Confidence intervals are constructed by test
  inversion over a user-supplied or automatically chosen grid.
  Returns a `coresynth_inference` subclass compatible with `tidy()` and `glance()`.

## Minor improvements

- `panel_to_matrices()`: fill loop replaced by vectorised `match()` + matrix-index
  assignment; removes an O(n × (T + N)) bottleneck in the shared data-prep path.
- `tasc.cpp`: `safe_inv_sympd()` helper added so the Kalman filter degrades to
  `pinv` instead of aborting when the innovation covariance is not numerically PD.
- `%||%` null-coalescing helper centralised in `utils.R`; duplicate definitions in
  `broom.R` and `plot.R` removed.
- `check_sharp_adoption()` (unused internal function) removed.

# coresynth 0.1.0

First public release.

## Methods

- **SCM** (Abadie, Diamond & Hainmueller 2010): Synthetic Control Method with unified formula interface.
  Supports predictor variables via `pred()`, out-of-sample V selection (`v_selection = "oos"`),
  donor filtering (`donor_mspe_threshold`), penalised SCM (`lambda_pen`), and staggered adoption.
  Inference: MSPE ratio permutation test via `mspe_ratio_pval()`.
- **SDID** (Arkhangelsky et al. 2021): Synthetic Difference-in-Differences.
  Supports time-varying covariates (`covariates =`), sharp and staggered adoption.
  Inference: `sdid_inference()` with placebo / bootstrap / jackknife / jackknife_global.
- **GSC** (Xu 2017): Generalised Synthetic Control with interactive fixed effects.
  Supports time-varying covariates via the full EM algorithm, sharp and staggered adoption.
  Inference: parametric bootstrap (`gsc_boot()`) and non-parametric (`gsc_inference()`).
- **MC** (Athey et al. 2021): Matrix Completion via nuclear-norm regularisation (Soft-Impute).
  Supports sharp and staggered adoption.
- **TASC** (Rho et al. 2026): Time-Aware Synthetic Control via Kalman EM.
  Supports sharp and staggered adoption.
- **SI** (Agarwal et al. 2025): Synthetic Interventions via SI-PCR.
  Supports sharp, staggered, multi-arm (K > 1), and staggered × multi-arm.
  Inference: `si_inference()` with bootstrap / jackknife / jackknife_global.
- **SCM-Design** (Abadie & Zhao 2026): `scm_design()` with base / weakly_targeted / unit_level variants,
  blank-period permutation test, and split-conformal confidence intervals.

## Unified API

- Single `scm_fit(outcome ~ treatment | unit + time, data, method = ...)` entry point for all methods.
- `panel_to_tensor()` for multi-arm SI data preparation.
- `broom` integration: `tidy()`, `glance()`, `augment()` for all methods and inference objects.
- `plot.coresynth()`: trend, gap, and weights plots via ggplot2.
- `export_json()`: JSON export for reproducibility.

## Performance

All core optimisations implemented in C++ via RcppArmadillo:
50–70x faster than the **Synth** package for typical panel sizes (N_co ≤ 30).
`src/inference.cpp` placebo loops parallelised with OpenMP.
