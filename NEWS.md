# coresynth (development version)

## Breaking changes

- **`colors`/`labels` in `plot.coresynth()` and `plot.scm_placebo()` are now
  keyed by one-word series identifiers**: `treated`, `synthetic`, and (with
  `show_donors > 0`) `donors` for trend plots; `treated` and `placebo` for
  the placebo plots. The previous keys were the displayed legend strings
  (`Treated`, `Synthetic Control`, `Placebo (donor pool)`), which required
  quoting for the multi-word names and conflated series identity with the
  legend text that `labels` itself relabels. Write
  `colors = c(treated = "black")` instead of `c(Treated = "black")`; an
  unknown key errors with the list of valid keys. The displayed legend text
  is unchanged.

- **`mspe_ratio_pval()` placebo refits now mirror the treated fit's
  predictor specification by default**: `use_covariates` defaults to `NULL`
  (auto) instead of `FALSE`. When the fit was estimated with a `predictors`
  specification, each placebo unit is now refit with that same specification
  -- the Abadie et al. (2010) / `Synth` convention, under which the treated
  and placebo test statistics are computed under one common spec. Previously
  the default silently switched the placebo refits to the outcomes-only
  spec, which (a) made the permutation compare statistics that were not
  exchangeable with the treated one, and (b) was extremely slow for long
  pre-periods because each placebo re-ran a `T_pre`-dimensional V
  optimisation. Outcomes-only fits are unaffected (auto resolves to the
  outcomes-only placebo path, as before). Pass `use_covariates = FALSE`
  explicitly to reproduce the old behaviour on covariate fits.

## Performance

- **The nested V/W coordinate descent is orders of magnitude faster.** The
  solver behind `scm_weights_cpp()` and the outcomes-only placebo loop in
  `scm_placebo_cpp()` now solves each inner simplex QP with a warm-started
  active-set method (KKT-verified exact solve, with FISTA fallback), applies
  V-coordinate changes to the Gram matrix as implicit rank-1 updates instead
  of rebuilding `X0' V X0` per grid point, and computes the Lipschitz
  constant once per coordinate-descent sweep instead of once per QP. On a
  monthly panel with `T_pre = 139`, `T_post = 19`, and 75 donors, the full
  75-unit outcomes-only placebo run drops from roughly an hour to about
  2 seconds, and a single outcomes-only `scm_fit()` from ~75 s to ~0.1 s.
  Per-unit results agree with the previous implementation within the
  coordinate-descent convergence tolerance: the estimator and its
  grid/accept/normalisation semantics are unchanged, and the inner QP
  solutions are now exact rather than first-order-approximate.

- **Covariate-spec placebo refits now run in parallel.** When the fit was
  estimated with a `predictors` specification, `mspe_ratio_pval()`
  previously refit each placebo unit in a sequential R loop; the loop now
  runs in C++ under OpenMP (new low-level routine `scm_placebo_x_cpp()`,
  the covariate counterpart of `scm_placebo_cpp()`). Each leave-one-out
  problem is solved by the same coordinate-descent core as before, so
  per-unit results are identical to machine precision -- only the wall
  time changes. The speedup is bounded by the slowest single placebo unit
  (iterations are distributed dynamically across cores), so it approaches
  the core count when per-unit costs are uniform and is smaller when one
  hard-to-fit donor dominates.

## New features

- **Treatment-line placement in plots**: `plot.coresynth()` (`type = "trend"`
  and `"gap"`) and `plot.scm_placebo()` (`type = "gaps"`) gain a
  `vline_offset` argument controlling where the vertical treatment line is
  drawn, in periods relative to the first post-treatment period (the previous
  fixed position). `vline_offset = -1` moves it to the last pre-treatment
  period, and fractional values such as `-0.5` interpolate between adjacent
  observed times, on numeric, `Date`, and `POSIXct` axes alike. The `vline`
  style list now also accepts an `xintercept` element for one or more
  absolute positions (previously an error), e.g.
  `vline = list(xintercept = "1989-01-01")` on a `Date` axis. Hiding the line
  still works via `vline = FALSE`.

- **Level-aligned trend and gap plots**: `plot.coresynth()` gains an `align`
  argument for `type = "trend"` and `"gap"`. `align = TRUE` shifts the
  synthetic series by its pre-treatment level gap to the treated series, so
  both are drawn on the same level. SDID's unit weights are estimated with
  the intercept concentrated out of the QP (Arkhangelsky et al. 2021), so its
  raw trend plot can show the synthetic control at a different level than the
  treated unit; `align = TRUE` uses the time-weight (lambda) weighted
  pre-period gap, which makes the average post-period gap in the plot equal
  the SDID estimate exactly.
- **Donor paths in trend plots**: `plot.coresynth(type = "trend")` gains a
  `show_donors` argument that draws the outcome paths of the `show_donors`
  donor units with the largest weights as thin background lines (`Inf` shows
  all). The new `"Donors"` series participates in `colors`/`labels`
  overrides.
- **SDID time weights in the weights plot**: for SDID fits,
  `plot.coresynth(type = "weights")` now shows two panels — donor unit
  weights (omega) and pre-period time weights (lambda) — in the same bar
  style. Other methods keep the single unit-weight panel.
- **Placebo SE and CI for SDID**: `sdid_inference(method = "placebo")` now
  also reports the placebo-distribution standard error (Clarke et al. 2023,
  Algorithm 4) and the corresponding normal-approximation confidence
  interval, for both sharp and staggered fits. The p-value is unchanged
  (permutation-based). `tidy()`/`glance()` pick the new columns up
  automatically.
- **Partially pooled staggered SCM** (Ben-Michael, Feller & Rothstein 2022,
  JRSS-B): `scm_fit(method = "scm")` on a staggered panel gains a `nu`
  argument. `nu = NULL` (default) keeps the previous behaviour (per-cohort
  V-optimised SCM). A numeric `nu` in `[0, 1]` minimises the convex
  combination `nu * (pooled pre-treatment imbalance)^2 +
  (1 - nu) * (per-cohort imbalance)^2` over all cohort weight vectors
  jointly, so that the aggregate ATT is anchored by the pooled fit
  (`nu = 0` reproduces separate per-cohort SCM with uniform lag weights,
  `nu = 1` is fully pooled). `nu = "auto"` selects the paper's heuristic
  value. The fit stores balance diagnostics in `fit$pooling`
  (`q_sep`, `q_pool`, and their separate-SCM baselines).
- **Intercept-shifted staggered SCM**: new `fixedeff` argument for staggered
  SCM fits. `fixedeff = TRUE` demeans every unit by its own pre-treatment
  mean within each cohort before fitting (Ben-Michael, Feller & Rothstein
  2022, Section 5.1; Doudchenko & Imbens 2017; Ferman & Pinto 2021), turning
  the estimator into a weighted difference-in-differences. Works with both
  the default path and the partially pooled path.
- **Wild bootstrap inference for staggered SCM**: new `scm_inference()`
  function implements the weighted multiplier (wild) bootstrap of
  Ben-Michael, Feller & Rothstein (2022, Section 5.3): donor weights are
  kept fixed and per-treated-unit effect contributions are perturbed with
  golden-ratio two-point multipliers. Returns a standard
  `coresynth_inference` object (works with `tidy()`/`glance()`). This is the
  first inference method available for staggered SCM fits.
- `solve_simplex_qp()` gains an optional `x0` warm-start argument, used by
  the partially pooled block coordinate descent to restart FISTA from the
  previous block solution. Together with an objective-based stopping rule
  this makes the pooled path roughly an order of magnitude faster on larger
  donor pools (N = 100: ~0.8 s to ~0.06 s per fit). Validated against the
  reference implementation `augsynth` (weights correlate at 1.0, identical
  heuristic `nu`, equal pooled imbalance at `nu = 1`).

## Bug fixes

- `plot.coresynth()` no longer embeds non-ASCII characters (Greek letters,
  the Unicode minus sign) in plot titles and subtitles. Some platforms
  cannot compute text metrics for these characters outside a UTF-8 locale,
  which crashed `R CMD check --run-donttest` on macOS.

# coresynth 0.3.0

## New features

- **Outcome-series accessors**: new exported generics `treated_outcomes()`,
  `synthetic_outcomes()`, and `donor_outcomes()` return the treated series,
  the estimated counterfactual, and the donor outcome matrix from any
  `coresynth` fit under a uniform interface, regardless of the estimation
  method. They replace the per-method field sniffing previously duplicated
  across `augment()`, `conformal_inference()`, and `plot()`; each returns
  `NULL` when the series is not stored in the fit (e.g. staggered fits,
  which keep their data per cohort in `fit$cohort_fits`).
- **Structural subclasses**: fits with staggered adoption now additionally
  inherit from `"coresynth_staggered"`, and multi-arm SI fits from
  `"coresynth_multiarm"`. S3 methods (`print()`, `summary()`, `tidy()`,
  `augment()`) dispatch on these subclasses instead of branching on internal
  flags. The `staggered`/`multi_arm` list fields are retained, so existing
  code that reads them keeps working, and all class checks remain
  `inherits()`-compatible.
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
  gain `colors`, `labels`, `vline`, `hline`, and (for `type = "weights"`)
  `fill` and `top_n` arguments. `top_n` restricts the weights bar chart to
  the `top_n` largest-weight donors (default `Inf` keeps every donor with a
  non-negligible weight), useful for large donor pools.
  `vline`/`hline` accept a list of `geom_vline()`/`geom_hline()`
  aesthetic overrides merged onto the built-in defaults, or `NULL`/`FALSE` to
  suppress the reference line entirely — since a line already added to a
  returned `ggplot` object cannot be removed, only restyled or overplotted,
  suppression has to happen inside the plot method. `colors` accepts a named
  vector overriding individual series colors (unmentioned series keep their
  default). `labels` accepts a named vector overriding the legend text of
  individual series (e.g. `labels = c(Treated = "California")`); keys are
  always the original series names, so `colors` and `labels` compose
  independently, and in `type = "ratios"` the treated unit's axis tick
  follows the relabeled legend entry. All defaults reproduce the previous
  appearance exactly, so existing calls are unaffected.

## Improvements

- `augment(fit, include_donors = TRUE)` now returns donor rows for
  `method = "si"` fits too (previously it warned that control outcomes were
  unavailable and returned treated rows only). Estimates are unaffected.
- `plot()` on a staggered fit now fails with a clear message explaining that
  staggered fits store their series per cohort, instead of an internal
  `data.frame` length error.
- The gap plot (`plot(fit, type = "gap")`) now labels its y-axis simply
  "Gap" and states the definition ("Treated − synthetic control") in a new
  subtitle, replacing the previous `Y_treated - Y_synthetic` axis label.
  The placebo gaps plot (`plot(<scm_placebo>, type = "gaps")`) uses the same
  "Gap" axis label — there each line is that unit's gap relative to its own
  synthetic control, so the old parenthetical was dropped rather than
  reworded. Cosmetic only; no computed values change.

## Bug fixes

- For `method = "tasc"`, `plot()`, `augment()`, and the `Y_synth` series of
  `export_json()` reported the average fitted value of **all** units as the
  counterfactual: TASC stores its fitted values as a full T x N matrix
  (`Y_hat`), unlike the other methods, and the shared extraction code
  averaged over every column. The treated unit's columns are now used
  (`synthetic_outcomes()` handles this per method), so the plotted/augmented
  counterfactual, gap, and residuals for `tasc` fits change; the ATT estimate
  itself was always computed from the correct per-unit gaps and is
  unaffected.
- `scm_design()` solved the `weakly_targeted` (eq. 9) and `unit_level`
  (eq. 10) designs by a sequential approximation: the treated weights `w`
  were always chosen to match the population average predictor vector alone,
  which made the `xi` penalty of eq. 10 effectively inert and left the
  eq. 9 objective jointly suboptimal. Both designs are now solved exactly
  for every candidate treated set — eq. 9 by two-block alternating
  minimisation of the jointly convex QP in `(w, v)`, eq. 10 by folding the
  per-unit synthetic-control losses into the treated-weight QP target —
  and the solutions have been verified against an independent exact QP
  solver. **This changes numerical results for
  `scm_design(design = "weakly_targeted")` and `design = "unit_level"` when
  `m >= 2`** (the default `m = 1` selects a single treated unit, so `w` is
  degenerate and both old and new solvers agree).

## Internal

- `conformal_inference()`'s counterfactual refit now dispatches on the fit's
  class (one S3 method per estimator) rather than an `if`-chain on the
  method string. Results are unchanged.

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
