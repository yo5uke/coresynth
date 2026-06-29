# Changelog

## coresynth 0.2.3

### Performance

- The inner simplex QP solvers (`solve_simplex_qp()` /
  `solve_simplex_qp_lr()`), used by every `scm_fit(method = "scm")` fit
  and by the `method = "sdid"` unit weights, now use FISTA with adaptive
  restart (O’Donoghue & Candès 2015, gradient scheme). The first-order
  solver’s iteration count grows with the condition number of
  `Q = X0'VX0`, which is large for the collinear pre-treatment outcome
  panels typical of synthetic control; resetting the momentum term when
  it works against the gradient removes this slowdown. Each inner solve
  converges to the same optimum (verified against an exact QP solver to
  within `1e-10`), so the returned weights are unchanged.
- For the common outcomes-only case the fit is effectively identical
  (differences `~1e-5`). For `v_selection = "oos"` and penalised
  (`lambda_pen`) fits on ill-conditioned panels, the non-convex outer V
  search may now settle on a different local optimum, shifting results
  slightly; both the previous and the new solutions are valid SCM fits
  with the same objective. **This can change numerical results for
  `v_selection = "oos"` and penalised fits on poorly conditioned data.**

## coresynth 0.2.2

CRAN release: 2026-06-26

### Bug fixes

- `panel_to_matrices()` (and therefore every
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md)
  method) and
  [`scm_design()`](https://yo5uke.com/coresynth/reference/scm_design.md)
  now error on duplicate `(id, time)` entries instead of silently
  keeping a single arbitrary row. A balanced panel requires each
  unit-time cell to be unique; duplicates were previously overwritten by
  the last matrix-index assignment, dropping data without warning. The
  error reports the number of offending rows and the first duplicated
  unit and time.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) now renders
  `Date`/`POSIXct` time axes correctly. The time vector was coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html), so dates
  appeared as days-since-epoch (e.g. 16000,
  20000. on the x-axis. `Date`/`POSIXct` values are now passed through
         unchanged so ggplot2 selects the appropriate date scale; only
         `character`/`factor` time values are coerced to numeric.

## coresynth 0.2.1

### Bug fixes

- `v_selection = "oos"` (outcomes-only case) previously fit candidate
  `W(V)` on the full pre-treatment outcome matrix and restricted only
  the MSPE evaluation to the validation window, allowing the V optimiser
  to fit the validation period indirectly (a data leak relative to
  Abadie (2021) S.3.2). The new `.scm_oos_outcomes()` implements the
  correct train/validation split: candidate `W(V)` are fitted on
  training-half outcomes only, `V*` minimises validation-half MSPE, and
  `W*` is refit with `V*` on the outcomes of the last `floor(T_pre/2)`
  pre-treatment periods. For OOS fits, `v_weights` now has
  `floor(T_pre/2)` entries and a new `v_rows` field records which
  periods they refer to. **This changes numerical results for
  `v_selection = "oos"`.**

### New features

- `scale_predictors` (default `TRUE`): predictor rows supplied via
  `predictors =` are now divided by their standard deviation across all
  units before optimisation, matching the Synth reference implementation
  (Abadie, Diamond & Hainmueller 2011, JSS). `predictor_table` continues
  to report values on the original scale. **This changes numerical
  results for SCM fits with user-supplied `predictors`** unless
  `scale_predictors = FALSE`.
- [`placebo_in_time()`](https://yo5uke.com/coresynth/reference/placebo_in_time.md):
  in-time placebo (backdating) test for sharp SCM fits (Abadie, Diamond
  & Hainmueller 2015; Abadie & Vives-i-Bastida 2022).
- [`loo_donors()`](https://yo5uke.com/coresynth/reference/loo_donors.md):
  leave-one-out donor robustness check with the predictor weights V held
  fixed (Abadie, Diamond & Hainmueller 2015, footnote 20).
- `build_predictor_matrices()` now errors with an informative message if
  a [`pred()`](https://yo5uke.com/coresynth/reference/pred.md) time
  window produces missing or non-finite predictor values.

## coresynth 0.2.0

CRAN release: 2026-06-12

### New features

- **Conformal inference**
  ([`conformal_inference()`](https://yo5uke.com/coresynth/reference/conformal_inference.md)):
  permutation-based p-values and confidence intervals following
  Chernozhukov, Wüthrich & Zhu (2021). Works with sharp fits across all
  supported estimation methods (`scm`, `sdid`, `gsc`, `mc`, `si`). The
  counterfactual proxy is re-estimated under the null on all *T* periods
  (essential for finite-sample validity per CWZ S.2.2), and p-values are
  obtained via moving-block (cyclic-shift) permutation of the estimated
  residuals. Confidence intervals are constructed by test inversion over
  a user-supplied or automatically chosen grid. Returns a
  `coresynth_inference` subclass compatible with
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html).

### Minor improvements

- `panel_to_matrices()`: fill loop replaced by vectorised
  [`match()`](https://rdrr.io/r/base/match.html) + matrix-index
  assignment; removes an O(n × (T + N)) bottleneck in the shared
  data-prep path.
- `tasc.cpp`: `safe_inv_sympd()` helper added so the Kalman filter
  degrades to `pinv` instead of aborting when the innovation covariance
  is not numerically PD.
- `%||%` null-coalescing helper centralised in `utils.R`; duplicate
  definitions in `broom.R` and `plot.R` removed.
- `check_sharp_adoption()` (unused internal function) removed.

## coresynth 0.1.0

First public release.

### Methods

- **SCM** (Abadie, Diamond & Hainmueller 2010): Synthetic Control Method
  with unified formula interface. Supports predictor variables via
  [`pred()`](https://yo5uke.com/coresynth/reference/pred.md),
  out-of-sample V selection (`v_selection = "oos"`), donor filtering
  (`donor_mspe_threshold`), penalised SCM (`lambda_pen`), and staggered
  adoption. Inference: MSPE ratio permutation test via
  [`mspe_ratio_pval()`](https://yo5uke.com/coresynth/reference/mspe_ratio_pval.md).
- **SDID** (Arkhangelsky et al. 2021): Synthetic
  Difference-in-Differences. Supports time-varying covariates
  (`covariates =`), sharp and staggered adoption. Inference:
  [`sdid_inference()`](https://yo5uke.com/coresynth/reference/sdid_inference.md)
  with placebo / bootstrap / jackknife / jackknife_global.
- **GSC** (Xu 2017): Generalised Synthetic Control with interactive
  fixed effects. Supports time-varying covariates via the full EM
  algorithm, sharp and staggered adoption. Inference: parametric
  bootstrap
  ([`gsc_boot()`](https://yo5uke.com/coresynth/reference/gsc_boot.md))
  and non-parametric
  ([`gsc_inference()`](https://yo5uke.com/coresynth/reference/gsc_inference.md)).
- **MC** (Athey et al. 2021): Matrix Completion via nuclear-norm
  regularisation (Soft-Impute). Supports sharp and staggered adoption.
- **TASC** (Rho et al. 2026): Time-Aware Synthetic Control via Kalman
  EM. Supports sharp and staggered adoption.
- **SI** (Agarwal et al. 2025): Synthetic Interventions via SI-PCR.
  Supports sharp, staggered, multi-arm (K \> 1), and staggered ×
  multi-arm. Inference:
  [`si_inference()`](https://yo5uke.com/coresynth/reference/si_inference.md)
  with bootstrap / jackknife / jackknife_global.
- **SCM-Design** (Abadie & Zhao 2026):
  [`scm_design()`](https://yo5uke.com/coresynth/reference/scm_design.md)
  with base / weakly_targeted / unit_level variants, blank-period
  permutation test, and split-conformal confidence intervals.

### Unified API

- Single
  `scm_fit(outcome ~ treatment | unit + time, data, method = ...)` entry
  point for all methods.
- `panel_to_tensor()` for multi-arm SI data preparation.
- `broom` integration:
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html),
  [`augment()`](https://generics.r-lib.org/reference/augment.html) for
  all methods and inference objects.
- [`plot.coresynth()`](https://yo5uke.com/coresynth/reference/plot.coresynth.md):
  trend, gap, and weights plots via ggplot2.
- [`export_json()`](https://yo5uke.com/coresynth/reference/export_json.md):
  JSON export for reproducibility.

### Performance

All core optimisations implemented in C++ via RcppArmadillo: 50–70x
faster than the **Synth** package for typical panel sizes (N_co ≤ 30).
`src/inference.cpp` placebo loops parallelised with OpenMP.
