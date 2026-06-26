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
