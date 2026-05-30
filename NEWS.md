# coresynth 0.2.0

## New features

- **Conformal inference** (`conformal_inference()`): permutation-based p-values and
  confidence intervals following Chernozhukov, Wüthrich & Zhu (2021).
  Works with sharp fits across all supported estimation methods
  (`scm`, `sdid`, `gsc`, `mc`, `si`).  The counterfactual proxy is re-estimated
  under the null on all *T* periods (essential for finite-sample validity per
  CWZ §2.2), and p-values are obtained via moving-block (cyclic-shift) permutation
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
