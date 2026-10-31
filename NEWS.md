# coresynth (development version)

A review of the GSC, TASC, SI, SDID and MC (Matrix Completion) code against
their papers and against gsynth/fect/MCPanel. Estimates from `method = "gsc"`
with covariates, from `method = "tasc"`, and the intervals of
`sdid_inference()`, `gsc_inference()`, `si_inference()` and `gsc_boot()` change.

## Performance and computational improvements

- **Matrix Completion (`method = "mc"`, `soft_impute_fe_cpp()`) eliminates
  redundant matrix allocations.** Two-way fixed effects back-fitting sweeps
  now precompute residual row and column sums outside the sweep loop and update
  via matrix-vector products, eliminating 20 temporary $T \times N$ matrix
  allocations per SVD iteration. Low-rank reconstruction truncates to singular
  vectors with positive thresholded singular values rather than forming dense
  diagonal matrices and multiplying zeroed components. In-place matrix updates
  eliminate intermediate `repmat` matrices.

## Bug fixes

- **`conformal_inference()` supports Matrix Completion panels with missing values.**
  When a panel contains unobserved or `NA` cells outside treatment, the refit
  under the null now marks them as unobserved (`O = 0`) instead of assuming a
  dense observation mask.

- **`fit_mc_cpp()` adds comprehensive argument validation.** `fixedeff`,
  `lambda`, `cv_folds`, `cv_nlambda`, `cv_seed`, infinite outcome values, and
  degenerate all-missing masks are now validated with informative error messages.
  The CV spectral norm fallback masks post-treatment cells consistently.

- **Documentation cites the published JASA 2021 formulation of Matrix Completion.**
  Parameter documentation and internal rules reflect the unpenalised two-way fixed
  effects estimator in Section 4 (Equation 4.3) and Mazumder et al. (2010, JMLR).

- **GSC with covariates no longer stops before its estimates settle.** The
  loop that alternates the factors and the covariate coefficients updated
  `beta` by plain alternating least squares and stopped silently after 50
  iterations. When the covariates move with the factors -- the design of Xu
  (2017)'s own simulations -- that left `beta` at (0.79, 2.77) against a true
  (1, 3) and biased the effect by +0.06. `beta` is now updated with the
  loadings concentrated out (Bai 2009), which reaches the same fixed point in
  a few dozen iterations; the loop runs to a tolerance of `1e-8` within at
  most 1000 iterations and warns if it does not get there. `gsc_ife_cpp()`
  reports `n_iter` and `converged`.

- **`gsc_inference()` and `si_inference()` resample the treated units, not
  only the controls.** Holding the treated units fixed left their own noise
  out of the standard error: with one treated unit, a 95% interval covered
  the truth 30% of the time for GSC, and 43% with five treated units for SI.
  The bootstrap now draws treated and control units separately, and the
  jackknife drops each of them in turn. A fit with a single treated unit (in
  any arm) is an error that points to `gsc_boot()` or
  `conformal_inference()`, since no resampling of units can reach that
  unit's noise. Coverage in simulations is now 0.90-0.96 for GSC and
  0.85-0.93 for SI.

- **On staggered fits, the bootstrap draws the controls once for all
  cohorts.** It used to redraw them independently for each cohort, although
  under the default `control_group = "never_treated"` every cohort shares the
  same donors, which dropped the correlation between the cohort effects (a
  95% interval covered 42% of the time). Treated units are drawn from all
  cohorts together, so the cohorts and their weights are rebuilt from each
  draw, as in Clarke et al. (2024, Algorithm 2). Coverage is now 0.93-0.95.

- **`sdid_inference()` follows Arkhangelsky et al. (2021) and Clarke et al.
  (2024) for the bootstrap and the jackknife.** The bootstrap resampled the
  controls only, holding the treated mean and `zeta2` fixed, and on staggered
  fits redrew the controls separately for each cohort; the jackknife dropped
  only controls and re-estimated the weights. Their 95% intervals covered
  0.68-0.90 of the time. Now `method = "bootstrap"` (Algorithm 2) draws units
  with replacement from the whole panel, treated and control alike, redraws
  a sample without a treated unit or with fewer than two controls, and reruns
  the whole estimator -- the covariate partial-out, `zeta2` (unless fixed in
  `scm_fit()`) and both sets of weights -- with the cohorts rebuilt from the
  treated units drawn. `method = "jackknife"` (Algorithm 3) drops every unit
  in turn with the fitted unit and time weights held fixed. Coverage is now
  0.88-0.97. Both need more than one treated unit and point to
  `method = "placebo"` otherwise, which is unchanged. This is the same scheme
  `sdid_event_study()` already used, so the two now agree.

- **The jackknife of a staggered fit needs two treated units in every
  cohort** for `sdid_inference()`, `gsc_inference()` and `si_inference()`:
  dropping a cohort's only treated unit leaves its effect undefined (Clarke et
  al. 2024, Algorithm 3), and such a fit now raises an error that points to
  the bootstrap.

- **The per-cohort jackknife of a staggered multi-arm SI fit no longer
  inflates its standard error.** It took the variance of each cohort's
  weighted sum over arms rather than of its weighted mean, which made the
  standard error about 20 times too large.

- **`gsc_boot()` follows Xu (2017, Algorithm 2).** It drew independent
  normal errors with a single variance, which ignored serial correlation (an
  interval covered 77% of the time when the errors had an autocorrelation of
  0.5) and gave the treated units the in-sample error of the controls. It
  now resamples whole residual series: the controls' in-sample residuals,
  and for the treated units the prediction errors of control units left out
  as if they were treated. Coverage is 0.92-0.93 in both cases.

- **TASC estimates its state-space model by the exact EM algorithm.** The
  M-step used the smoothed state means only, leaving out their covariances.
  That understates the state noise, and the understatement compounded until
  the state noise collapsed to zero within the 20 iterations the fit ran.
  The M-step now follows Rho et al. (2026, Algorithm 7) -- smoothed
  covariances included, the initial state learned, diagonal noise
  covariances -- with the missing cells of Shumway & Stoffer (1982), and
  runs until the fitted values settle. `em_iter` is now the maximum number
  of iterations (default 2000), `em_tol` (default `1e-6`) sets the
  tolerance, and a fit that does not converge warns.

- **`conformal_inference()` keeps the covariates of a GSC fit.** It refitted
  the factor model without them, so an interval for a fit with covariates was
  computed on a different model.

- **`scm_design()` uses the order-statistic quantile of Abadie & Zhao (2026,
  Equation 18)** for its split-conformal intervals, instead of interpolating
  below it, which narrowed the intervals when there were few blank periods.
  When it samples permutations it no longer resets the random number
  generator of the calling session.

## New features

- **`method = "gsc"` takes additive fixed effects through `force`**: `"unit"`,
  `"time"` or `"two-way"` (default `"none"`), the extended model of Xu (2017,
  Remark 3). With additive effects `r = 0` is allowed, and the
  cross-validation of `r` starts from it. Estimates match gsynth to machine
  precision for each of these choices.

- **TASC gives each unit its own intercept** in the observation equation
  (`intercept = TRUE`, the default). Without one, the states have to carry
  each unit's level, and on panels far from zero the EM drifted for tens of
  thousands of iterations, moving the effect by up to 0.4 depending on when
  it stopped. `intercept = FALSE` fits the paper's model as written.

## Internal

- Staggered GSC fits the never-treated donors once rather than once per
  cohort, and `sdid_inference()`, `gsc_inference()` and `si_inference()`
  share one resampling routine. GSC and SI staggered fits record
  `control_group`; SDID staggered fits record `zeta2_user` and `zeta_t`, and
  the raw outcomes and covariates when they have covariates.

# coresynth 1.1.0

## New features

- **`sdid_event_study()` computes period-by-period SDID treatment effects with
  cluster or unit bootstrap inference.** For a sharp `method = "sdid"` fit, it
  returns the event-study series aligned on the `lambda`-weighted pre-treatment
  gap -- so the post-treatment mean of the period effects equals `fit$estimate`
  (the ATT) -- alongside bootstrap standard errors and both normal-approximation
  and percentile confidence intervals for every period and for the overall ATT
  from the same bootstrap draws:

  - `cluster` accepts a column name in `data` (or one-sided formula `~ col`) or
    a named vector keyed by unit ID (positional unnamed vectors are rejected so
    unit-to-cluster alignment never depends on row order). When `cluster = NULL`
    (default), individual units are resampled.
  - Each bootstrap draw resamples **both treated and control clusters/units**
    with replacement, treats duplicated units as distinct columns, and
    re-estimates `sigma_hat`, `zeta2`, unit weights `omega`, and time weights
    `lambda`. Draws in which either the treated or control group is empty (or
    has fewer than two controls) are discarded, and `n_boot_valid` and
    `n_boot_dropped` report the valid and discarded counts.
  - `boot_gaps` (`n_boot_valid x T` matrix) and `boot_att` (length
    `n_boot_valid` vector) store the replication draws for full reproducibility
    with `seed`.
  - `plot()`, `plot_data()`, `tidy()`, and `glance()` methods format
    and plot the event-study series with point estimates, confidence intervals,
    zero and treatment reference lines, and observed-period axis breaks that
    handle non-consecutive time indices cleanly.

- **`sdid_inference()` documents and enforces its control-only bootstrap
  scope.** Its `method = "bootstrap"` resamples control units only with the
  treated mean series and `zeta2` held fixed; passing `cluster` to
  `sdid_inference()` now raises an informative error pointing to
  `sdid_event_study()` (and notes when a non-bootstrap `method` was selected).

# coresynth 1.0.1

## Bug fixes

- **`gsc_boot()` reports a confidence interval for the ATT.** Its
  replicates are generated with no treatment effect, so they are draws of the
  estimation error (Xu 2017, Algorithm 2), and the interval is the estimate
  minus their quantiles. `ci_lower`/`ci_upper` were the raw quantiles
  instead -- an acceptance region around zero that never contained a
  non-zero estimate. **`ci_lower`/`ci_upper` change**; `p_value`, `se` and
  `boot_dist` do not.

- **`scm_fit()` rejects a transformed outcome or treatment in the formula.**
  The columns are read from `data` by name, so `log(y) ~ d | id + time` was
  fitted as `y` and `y ~ I(t > 1990) | id + time` used `t` as the treatment,
  without a word. Both are now errors that ask for the variable to be created
  in `data` first.

- **The parallel placebo loops cannot take the R session down.** An
  exception thrown inside an OpenMP region terminates the process; the
  outcomes-only SCM and the SDID placebo loops now catch it and mark that
  placebo `NaN`, as the predictor-spec SCM loop already did.

- **OpenMP honours a thread limit.** The placebo loops and the multi-start
  search used every core, which the CRAN Repository Policy forbids during
  checks. `options(coresynth.threads = n)` now caps the thread count, and
  `R CMD check --as-cran` runs on at most two threads (see
  `?coresynth-package`). Results do not depend on the thread count.

- **`R CMD check` runs the test suite.** `tests/testthat.R` was missing, so
  the checks (locally and on CI) never ran the tests.

## Performance

- **SDID weights are solved without forming the donor Gram matrix.** The unit
  weight QP has one row per pre-treatment period and a column per donor, so
  with more donors than periods the gradient and the step size now come from
  the factor instead of an N-by-N matrix and its eigendecomposition. Same
  iterates up to round-off (estimates, standard errors and p-values agree to
  ten digits); on a 300-donor panel the placebo, jackknife and bootstrap
  methods of `sdid_inference()` run about 6 times faster.

- **The TASC Kalman filter avoids an N-by-N inverse per period.** The
  observation noise is diagonal, so the gain reduces to an r-by-r solve by the
  push-through identity. Same moments up to round-off; a 200-unit fit runs
  about 10 times faster.

# coresynth 1.0.0

## Breaking changes

- **Staggered fits of `"scm"`, `"sdid"`, `"gsc"` and `"si"` now take their
  donors from the never-treated units only (`control_group =
  "never_treated"`).** The old default, `"clean"`, also admitted units that
  adopt after the cohort. Each cohort's effect is averaged from its adoption
  date to the end of the panel, so those units were treated themselves during
  part of that window, and their post-adoption outcomes entered the cohort's
  counterfactual. With a constant effect of 3 in a three-cohort simulation
  the aggregate ATT came out 0.24 to 0.59 too low depending on the method,
  against a bias of at most 0.04 with never-treated donors. Neither source
  design does this: Clarke et al. (2024, Algorithm 1) subset to pure
  controls, and Ben-Michael, Feller & Rothstein (2022, S.2.2) admit a donor
  only if it is untreated through the last event time estimated. **Staggered
  estimates from these four methods change.** `"clean"` still runs but warns
  and will be removed in the next major release. A staggered panel with fewer
  than two never-treated units is now an error that points to `"mc"` and
  `"tasc"`, which impute from untreated cells and need no donor pool.

- **Covariate partial-out now includes unit and time fixed effects.** The
  `covariates` argument of `"sdid"` (sharp and staggered) and of staggered
  `"scm"` with `cov_method = "partial_out"` implements the projection of
  Kranz (2022) that Clarke et al. (2024, S.2.2) adopt: regress the outcome
  on the covariates *and two-way fixed effects* over the untreated cells,
  then subtract `X * beta`. The fixed effects were missing, so `beta` also
  absorbed how the covariates' unit levels and trends line up with the
  outcome's, and subtracting it wrote a spurious differential trend back into
  the series. In a simulation with a true coefficient of 1, `beta_hat` was
  3.75 and the staggered SCM and SDID estimates of an effect of 2 were 1.65
  and 1.18; they are now 1.00, 2.01 and 2.01. Sharp SDID also regressed on
  the control units' pre-treatment cells only; it now uses every untreated
  cell, the default of Kranz's own implementation, as staggered fits did. **Estimates with `covariates`
  change.** A covariate that varies only across units or only over time is
  absorbed by the fixed effects; it now gets a zero coefficient and a warning
  (for staggered SCM, `cov_method = "balance"` balances such covariates
  instead).

- **`scm_inference()` defaults to the leave-one-unit-out jackknife** of
  Ben-Michael, Feller & Rothstein (2022, Appendix A.1), which drops each unit
  of the panel, donors included, refits with `nu` held fixed, and uses the
  spread of the estimates. It refits once per unit -- about 0.7 s at 100
  units and 3 s at 200 -- skipping donors that carry no weight in any block,
  whose removal leaves the solution unchanged, and warm-starting the rest
  from the full fit. It warns below 3 treated units, where it over-rejected
  mildly (8.5% at a nominal 5%). The result carries the leave-one-out
  estimates in `jack_ests`.

- **`method = "wild_bootstrap"` now perturbs every unit, as Otsu & Rai
  (2017) prescribe.** The bootstrap BFR (2022, S.5.3) use adapts Otsu &
  Rai's wild bootstrap for matching estimators, whose key step is to
  resample "not only the regression residuals of treated and untreated
  observations but also the ones to estimate the average treatment effects"
  (their Eq. 9). The previous implementation perturbed only the treated
  units' deviations from the ATT, so it missed the error of the synthetic
  control that every treated unit of a cohort shares, and under-covered: a
  nominal 95% interval covered 80% with 8 donors and 8 treated units per
  cohort, and 45%, 63%, 82% and 91% with 2, 3, 4 and 8 treated units. Each
  donor now carries its own multiplier on its residual, weighted by how much
  the synthetic controls use it, and only the treated units' heterogeneity
  is centred at the ATT -- as in Otsu & Rai, and unlike BFR's Eq. (13) as
  printed and `augsynth`, which centre every unit and so widen the interval
  in proportion to the effect. Residuals are against the two-way model of
  the untreated outcome, so a common period shock does not widen it either.
  On the same designs it covers 93% to 99.5% (93% with 2 treated units); it
  is conservative where outcomes carry unit-specific trends the two-way
  model does not absorb. It needs no refitting, which makes it the option
  for panels too large for the jackknife, and now warns below 3 treated
  units rather than 5. **Its intervals change.**

- Staggered SCM fits store every unit's adoption date (`T_adopt`), and
  `scm_balance_frontier()` and the new jackknife rebuild the panel from it.
  Rebuilding it from `cohort_fits`, as before, turned the treated units of a
  cohort the fit had skipped (for too short a pre-treatment window) into
  never-treated donors for the refits.

- **`"gsc"`, `"mc"` and `"si"` now choose their tuning parameter the way their
  source papers say to.** All three picked it by a fixed rule that appears
  nowhere in the literature: `r = 2` factors for GSC, a penalty of 1% of the
  spectral norm for MC, and `k = floor(sqrt(min(T_pre, N_co)))` components for
  SI. Each paper specifies a data-driven selection instead, and each is now
  the default, so **estimates from all three change** unless the parameter is
  supplied explicitly.

  - `"gsc"`: `r = NULL` runs the leave-one-period-out cross-validation of
    Xu (2017, Algorithm 1) — for each candidate `r`, hold out one
    pre-treatment period, refit the treated units' factor loadings on the
    rest, and score the prediction. Xu presents this as central to the method
    ("it avoids specification searches"). On panels built with a known number
    of factors it recovered `r` exactly for 1 through 4, where the old fixed
    `r = 2` was right only by coincidence. The scoring table is kept in
    `fit$r_cv`. Candidates start at 1: `r = 0` is not representable by a
    rank-`r` truncated SVD.
  - `"mc"`: `lambda = NULL` runs the K-fold cross-validation of Athey et al.
    (2021, S.4), drawing folds so the observed fraction inside a fold matches
    the panel's and scoring each penalty on the held-out cells. It uses its
    own RNG stream, so the penalty is reproducible and the caller's random
    state is left untouched; `cv_folds`, `cv_nlambda` and `cv_seed` control
    it, and `fit$lambda_cv` reports whether the search ran.
  - `"si"`: `k = NULL` applies the optimal hard threshold of Gavish & Donoho
    (2014) to the pre-treatment donor block, one of the three data-driven
    rules Agarwal et al. name. It costs a single SVD, where cross-validating
    `k` would cost a PCR fit per candidate. On the same known-rank panels it
    recovered the rank exactly for 1 through 4; the old heuristic returned 4
    in every case. A staggered fit applies the threshold to each cohort's own
    donor block, and clips a supplied `k` to that cohort's dimensions, so
    `fit$k` is one value per cohort, named by adoption date.

- **Staggered `"mc"` and `"tasc"` fits now report the same estimand as every
  other method.** Both estimators recover the counterfactual by imputation:
  their observation mask already blanks each treated unit from its own
  adoption date, so a staggered panel needs no cohort-by-cohort refitting and
  the single fit was, and remains, a valid staggered fit. What was wrong was
  the aggregation. They averaged the per-unit post-treatment gaps with one
  weight per treated unit, while `"scm"`, `"sdid"`, `"gsc"`, and `"si"` use
  the `N_treated,g * T_post,g` weights of Clarke et al. (2024). The two agree
  only when every cohort has the same size and the same number of
  post-treatment periods; otherwise the same panel yielded a materially
  different number from `"mc"` and `"tasc"` than from the other four, with
  nothing to indicate it. They now use the cohort weights as well, so
  **staggered `"mc"` and `"tasc"` estimates change**. Sharp fits are
  unaffected and return bit-identical results.

  The fits also declare themselves properly. They gain `staggered = TRUE`,
  `cohort_estimates`, and `cohort_fits`, and inherit from
  `"coresynth_staggered"`. Consequently `glance()` reports `staggered = TRUE`
  rather than `FALSE`, `tidy()` returns the per-cohort table instead of an
  empty data frame, `augment()` returns one series per cohort with `.period`
  keyed to that cohort's own adoption date, and `plot()` errors and points at
  `augment()`, as it already did for the other staggered fits. Most
  consequentially, `conformal_inference()` rejects them: a staggered `"mc"`
  fit previously passed its sharp-only guard unnoticed and returned a p-value
  computed against a single global treatment date.

- **`method = "scm"` now fits panels where several units adopt at the same
  date, instead of erroring.** `scm_fit()` refused any sharp panel with more
  than one treated unit, while the *same* units with adoption dates one period
  apart went down the staggered path and fitted without complaint. Simultaneous
  adoption is one cohort, and S.4.1 of Ben-Michael, Feller & Rothstein (2022)
  shows that the across-unit variation putting unit-level fits into the ATT
  error bound vanishes when the adoption date is shared, so the bound is
  minimised by pooling those units: one synthetic control for their average.
  That is the step their Appendix A.2 hierarchy already takes within each
  cohort, and the one the staggered path here already took.

  The fit reports the treated average in `Y_treat` and keeps the raw per-unit
  series in `Y_treat_mat`, alongside `n_treated`, which `print()`, `summary()`,
  `glance()` and `export_json()` all report. Donor weights stay a single
  simplex vector, so `predictors`, `v_window`, `lambda_pen`,
  `v_selection = "oos"`, `plot()`, `placebo_in_time()`, `loo_donors()` and
  `augment_scm()` all work as they do for one treated unit. Each `pred()` row
  of the treated predictor vector is the average over the pooled units, under
  the same weights the outcome rows are pooled with, so the predictor and
  outcome sides describe one target series rather than two.

  `mspe_ratio_pval()` does not, and now says so, pointing at `scm_rank_pval()`
  below. It compares the treated unit against placebo runs on single donors,
  and an average of several treated units is less noisy than any one donor,
  which inflates its pre/post MSPE ratio. Simulated under no treatment effect,
  the test rejected at 28% with three pooled units and 53% with six, against a
  nominal 10%.

- **Staggered `"scm"` defaults to the pooling heuristic and the intercept
  shift.** `nu` defaulted to separate per-cohort SCM and `fixedeff` to
  `FALSE`. Ben-Michael, Feller & Rothstein argue against exactly that
  combination: separate SCM is "unlikely to yield good results" (abstract),
  S.4.2 advises against values of `nu` near 0 or 1, and S.6 reports that the
  intercept-shifted estimator "dominates either of the alternatives". Their
  own implementation, `augsynth`, defaults to both. So **staggered SCM
  estimates change**; `nu = "separate"` and `fixedeff = FALSE` restore the old
  behaviour, and `nu` gained `"separate"` as a name for the older per-cohort
  path so that `0` can mean what it means in the paper.

  Reproducing their S.6 design -- 49 units, 39 periods, selection into
  treatment through the unit effects, about 20 treated units, a sharp null --
  bears this out. Under the two-way fixed effects model the heuristic
  improved on separate SCM in bias, MAD and RMSE (MAD 0.0395 against 0.0413);
  under the factor model full pooling was much worse than the heuristic (MAD
  0.1492 against 0.0607), as the paper reports. The gain depends on selection
  into treatment: in a panel where treatment is assigned at random the
  heuristic is the worse choice, which is worth knowing before reaching for
  it.

  The default falls back to `"separate"` where the pooled solver cannot go
  (`donor_mspe_threshold`, `lambda_pen`, `v_selection = "oos"`); naming a
  numeric `nu` there is still an error. The fallback reports itself in a
  message naming the argument that caused it, rather than changing the
  estimator in silence, and `fit$nu` records which of the two ran: the number
  the pooled solver settled on, or `"separate"`. The heuristic itself was
  checked against `augsynth` and agrees to four or five decimal places,
  across cohorts of equal and unequal size alike.

- **`qp_solver` now defaults to `"auto"`, which is the Wolfe solver wherever
  the inner QP can be degenerate.** The inner SCM problem projects the treated unit onto the
  convex hull of the donor columns in predictor space, so by Caratheodory's
  theorem an optimum supported on at most `k + 1` donors always exists, and
  Wolfe (1976) returns such a solution. `"active_set"` returns an exact
  optimum too, but when the QP is degenerate — fewer predictors than donors
  in the support, where the optimal set is a face rather than a point — it
  breaks the tie on round-off, so its weights can move with the arithmetic.
  This was opt-in since 0.4.0 precisely because it changes those weights;
  1.0.0 is where that change belongs.

  On a degenerate predictor spec the difference is not only which optimum is
  reported. In a 30-donor test panel with one or two aggregated predictors,
  Wolfe reached a *lower* pre-treatment loss than the active-set path (2.75
  against 3.25, and 2.79 against 3.31) on 3 and 2 non-zero donors rather than
  13 and 17 — the sparser, more stable inner solutions give the non-convex
  outer V search a cleaner objective to descend.

  Degeneracy needs few rows next to the donor pool, so `"auto"` is Wolfe for
  every predictor fit and for outcomes-only fits with fewer than half as many
  pre-treatment periods as donors. The other outcomes-only fits have a unique
  optimum in practice, which both solvers return, and there `"auto"` keeps
  the active set: Wolfe cold-starts each of the coordinate sweep's inner
  problems and was three times slower for identical weights.

  `"auto"` also falls back to `"active_set"` on the paths that run their own QP
  outside the sharp inner solver (`v_selection = "oos"`, `lambda_pen`,
  `v_optim = "bfgs"`, and staggered adoption). Naming `"wolfe"` explicitly on
  one of those is still an error rather than a silent downgrade, and
  `fit$qp_solver` records the value actually used.

- **The default `vline_offset` in `plot.coresynth()` and
  `plot.scm_placebo()` is now `-1` instead of `0`**, moving the treatment
  line from the first post-treatment period to the last pre-treatment one.
  At offset 0 the treatment is already in effect, so the treated path and its
  counterfactual have usually separated by the time the line is reached and
  the figure reads as though the gap opened before treatment. With `-1`
  everything left of the line is the period the fit was asked to match. Pass
  `vline_offset = 0` for the previous position.

- **`scm_fit()` now takes `...` ahead of its estimator-specific arguments**,
  which must therefore be named in full; partial matches such as `v_opt =`
  and any use of the fourth or later argument by position no longer work.
  This is what makes `lambda` reach `method = "mc"`: while `lambda_pen` (an
  SCM-only argument) sat before `...`, R matched `lambda` to it partially, so
  MC's nuclear-norm penalty could not be set through `scm_fit()` at all and
  the value silently went to an argument MC ignores.

- **`scm_fit()` rejects arguments the chosen estimator does not take.**
  One formula interface serves six estimators, so an argument meant for
  another one — `predictors` or `nu` with `method = "mc"`, `covariates` or
  `fixedeff` with `"tasc"`, `control_group` with either — or simply
  misspelled, was absorbed by `...` and changed nothing about the fit that
  came back. These now error and name the methods the argument applies to.

- **Sharp `method = "scm"` rejects `covariates` instead of ignoring it.** The
  partial-out ran inside the staggered branch only, so on a single-cohort
  panel the argument reached the fit function, was never read, and returned an
  estimate identical to the one without it -- the same silent no-op the
  argument check above was added to stop, one level deeper. It is not an
  oversight that SCM has no partial-out: SCM balances covariates rather than
  residualising them. Abadie, Diamond & Hainmueller (2010, S.2.3) put them in
  the V-weighted predictor distance, and Ben-Michael, Feller & Rothstein
  (2022, S.5.2) add a separate balance term with its own weight, noting in
  S.4 that "where Abadie et al. (2010, 2015) balance auxiliary covariates, we
  focus exclusively on lagged outcomes". Partialling out is the device of
  Clarke et al. (2024, S.2.2), which SDID needs because it has nothing to
  balance with.

  A sharp fit does have something: `predictors = list(pred("x", <pre-treatment
  times>, "mean"))` balances a covariate on its pre-treatment aggregate, which
  is how BFR carry a time-varying covariate into the balance problem too. The
  error names that route. Staggered SCM keeps the partial-out, because
  `predictors` is unavailable there.

- **`pool_cohorts` now takes effect on the per-cohort staggered path.** It was
  read only by the partially pooled solver, so `pool_cohorts = FALSE` was
  silently discarded whenever the fit went down the `nu = "separate"` path --
  including the automatic fallback that `donor_mspe_threshold`, `lambda_pen`
  and `v_selection = "oos"` trigger. The fit then carried cohort-level blocks
  the caller had asked not to have, and `scm_rank_pval()`, which needs one
  synthetic control per treated unit, rejected it with a message telling the
  caller to pass the argument they had just passed. Both staggered paths now
  build their blocks the same way, and `fit$pool_cohorts` records it.
  Sharp fits reject the argument, as they already did `nu` and `fixedeff`.

- **Partially pooled staggered SCM weights each cohort by its size.** The two
  imbalance measures `nu` trades off, and the heuristic that sets it, treated
  every cohort as equally important however many units it held. Appendix A.2
  of Ben-Michael, Feller & Rothstein (2022) -- the hierarchy that pools units
  sharing an adoption date, which is what `pool_cohorts = TRUE` implements --
  writes the cohort problem on the *sum* over the cohort's units: its target
  is `sum_i 1{T_i = T(g)} Y_i` and its weights satisfy `sum_i gamma_ig = n_g`.
  A cohort's residual is therefore `n_g` times its mean gap, so a cohort of
  `n_g` units counts `n_g` times towards the pooled imbalance and `n_g^2`
  towards the separate one. Weighting them all equally is neither, and it is
  also not what the fit aggregates with, since the cohort ATTs combine with
  the `N_treated,g * T_post,g` weights of Clarke et al. (2024).

  **Staggered `"scm"` estimates change where the cohorts differ in size.**
  Panels whose cohorts are all the same size, and panels with one treated unit
  per adoption date, are unaffected -- the weighting is a common factor there,
  and the objective is normalised. `augsynth`, the authors' implementation,
  already followed the appendix; on a panel with cohorts of 1 and 4 its
  heuristic returned 0.94114 where this one returned 0.93563, and it now
  returns 0.94113. Across three further panels with cohorts of 1, 2 and 8 the
  two now agree to four decimal places, where they previously differed by as
  much as 0.025.

- **`weights` now reaches the objective, not just the aggregation.** A
  staggered fit multiplied the supplied mass into each cohort's aggregation
  weight while still fitting its weights to the *unweighted* average of that
  cohort's treated units, so with unequal weights inside one cohort the
  reported ATT was neither the equally weighted quantity nor the weighted one.
  The mass now sets the block's target series, its weight in the partially
  pooled objective, and its `T_post`-multiplied aggregation weight alike, and
  `scm_inference()` resamples each treated unit at its own mass rather than
  splitting the block's evenly. `scm_balance_frontier()` refits on the same
  weighted objective, where it previously ignored `weights` entirely.
  **Staggered `"scm"` estimates change when `weights` is supplied**; equal
  weights reproduce the unweighted fit exactly.

- **`"mc"` estimates unit and time fixed effects outside the nuclear norm.**
  The S.4 model these fits implemented, `Y = L + e`, leaves the panel's level
  inside the penalised matrix, so the nuclear norm shrinks it along with the
  factor structure and the imputed counterfactual is pulled towards zero. On a
  panel shifted to a level of 20 the ATT came out at +0.230 against a true 2
  (z = 4.5 over 40 draws). S.8.1 of the same paper keeps unit and time effects
  out of the penalty -- "we regularize the estimates of L*, but do not wish to
  regularize the estimates of the fixed effects" -- and MCPanel, the authors'
  own implementation, estimates them by default. Doing so removes the bias
  (-0.020, z = -0.4) and cuts RMSE by about 30%. The previous behaviour
  matched MCPanel with its fixed effects switched off, and the new behaviour
  matches its default, to within simulation error in both cases.

  `fixedeff = FALSE` recovers the S.4 baseline. The fit gains `fe_unit`,
  `fe_time` and `L_lowrank`; `L_hat` remains the completed matrix and now
  includes the fixed effects. `conformal_inference()` refits under whichever
  model the fit used.

  This also corrects the note published with the cross-validated penalty
  above: the small-panel upward bias recorded as a known limitation was this
  misspecification, not a property of nuclear-norm regularisation. The
  debiasing of Choi & Yuan (2024, Eq. 3.1) was tested against it and only
  halves it (+0.123, z = 2.5), adding nothing once the fixed effects are in.

## New features

- **`scm_rank_pval()` tests the aggregate effect when several units are
  treated** (Abadie & L'Hour 2021, S.3.2). Each treated unit contributes its
  own statistic, the identities of the treated are reassigned at random, and
  the sum of the ranks of the real treated units is compared with its
  permutation distribution. Ranking is what makes this work where
  `mspe_ratio_pval()` cannot: every statistic being compared comes from a
  single unit, so the averaged treated series never enters. Simulated with no
  treatment effect it rejected at 7.0% and 10.0% against a nominal 10% (three
  and six treated units), where `mspe_ratio_pval()` on the same fits reached
  28% and 53%; with an effect of 1.5 its power was 77% and 85%.

- **`mspe_ratio_pval()` can discard ill-fitting placebo runs.** Abadie,
  Diamond & Hainmueller (2010) report their permutation distribution at
  cutoffs of 20, 5 and 2 times the treated unit's pre-treatment MSPE, on the
  grounds that a placebo run the synthetic control could not track before
  treatment says nothing about how rare a large post-treatment gap is. That
  cutoff existed only in `plot.scm_placebo()`, where it changed the picture
  but not the p-value. `mspe_prune` now applies it to the test, and
  `n_placebo_pruned` reports how many runs it dropped. This is a different
  control from `mspe_threshold`, which is a floor guarding the ratio's
  denominator against a near-zero pre-treatment MSPE.

  `statistic = "post_mspe"` adds the alternative of Abadie (2021, S.3.5):
  rank the post-treatment MSPE itself, rather than its ratio to the
  pre-treatment MSPE, once the ill-fitting runs are gone. On a panel carrying
  five donors no synthetic control can track, that statistic gives p = 0.150
  unpruned -- the unfittable donors sit at the top of the distribution on size
  alone -- and p = 0.029 at `mspe_prune = 5`. The ratio statistic reaches
  p = 0.025 without any cutoff, which is the advantage ADH claim for it in
  their Figure 8.

  `plot()` and `plot_data()` on the result now default to whichever cutoff the
  test used, so the figure and the p-value describe the same set of placebo
  runs. This holds for `type = "ratios"` as well as `type = "gaps"`: the
  subtitle reports the number of runs the p-value rests on, so a figure
  drawing the discarded ones contradicts its own caption. Pass
  `mspe_prune = Inf` to draw them all regardless.

  `plot(type = "ratios")` also draws whichever statistic the test ranked, so
  a `statistic = "post_mspe"` result gives the post-treatment MSPE rather than
  its ratio, with the axis and title to match. `plot_data(type = "ratios")`
  returns both statistics, in `ratio` and `post_mspe`.

- **`pool_cohorts = FALSE` fits one synthetic control per treated unit.**
  Staggered SCM pools units sharing an adoption date into a cohort, which is
  the Appendix A.2 hierarchy of Ben-Michael, Feller & Rothstein (2022). Their
  S.4 formulation instead gives every treated unit its own weight vector, and
  that is what `pool_cohorts = FALSE` fits; `cohort_fits` then names the unit
  behind each block, which is what `scm_rank_pval()` needs. The default stays
  with the appendix, where S.4.1 shows there is no across-unit variation to
  exploit inside a cohort. (`augsynth` defaults the other way.)

  Those blocks are labelled by unit rather than by adoption date, because
  several of them share one: `cohort_estimates` carries a `unit` column,
  `tidy()` returns `unit_<id>` terms, and `augment()` keys `.unit` on it.
  Labelling them by cohort would have given two treated units the same key,
  which merges their series wherever it is grouped or drawn on. A fit whose
  blocks are cohorts keeps the labels it had, and the `unit` column names the
  unit of a one-unit cohort and is `NA` for a pooled one. Estimators without
  the notion carry no such column.

- **`weights` sets how much each treated unit counts towards the ATT.**
  Abadie (2021, S.8) notes that with regions or countries a population- or
  GDP-weighted average effect is often the quantity of interest. A staggered
  fit multiplies the supplied mass by each block's `T_post`; a sharp fit makes
  the pooled target a weighted average of the treated series. Left `NULL` the
  units count equally, as before.

- **`scm_balance_frontier()` and `plot()` trace the pooling trade-off.**
  Ben-Michael, Feller & Rothstein call this plot "an important tool for using
  partially pooled SCM in practice" (S.4.2). It refits across a grid of `nu`
  and reports the imbalance of the treated average, which the aggregate ATT
  rests on, against the root mean square of the per-block imbalances, which
  the block-level effects rest on. The curve is strongly convex in practice,
  so it shows where more pooling stops paying for itself: on a test panel,
  moving `nu` from 0 to 0.5 cut the pooled imbalance by 21% for a 6% rise in
  per-block imbalance, while the last quarter cost 39% for a further 8%.

  This is worth looking at because the `nu` heuristic does not always pick
  well. Where treatment is assigned without selection on the unit effects,
  every block already fits to the noise floor and the block-level gaps are
  near-independent, so their average shrinks by about `1/sqrt(J)` and the
  heuristic returns roughly that -- a number driven by noise averaging rather
  than by imbalance pooling could remove. Separating the two would need an
  estimate of the noise floor, which is neither in the paper nor in
  `augsynth`, so the default is unchanged and the frontier is offered instead.

- **`print()` and `summary()` describe staggered fits as staggered.** Both
  previously used the sharp output, which reports `T_pre` and `T_post` from
  the panel-level split — the earliest adoption date, correct for one cohort
  and wrong for the rest — with nothing to say the fit was staggered at all.
  `print()` now labels the method `(staggered)` and gives the cohort count in
  place of that split, and `summary()` prints the per-cohort table of
  estimates, weights, and treatment windows. Multi-arm SI fits show the
  cohort-arm table on the same terms.

- **`cov_method = "balance"` fits the auxiliary-covariate form of staggered
  SCM** (Ben-Michael, Feller & Rothstein 2022, S.5.2, Equation 11). Every SCM
  paper handles covariates by balancing them -- Abadie, Diamond & Hainmueller
  put them in the V-weighted predictor distance, and BFR note in S.4 that
  "where Abadie et al. (2010, 2015) balance auxiliary covariates, we focus
  exclusively on lagged outcomes" before reintroducing them this way.
  Partialling out, which is what `covariates` has always done here, is the
  device of Clarke et al. (2024, S.2.2) that SDID needs because it has
  nothing to balance with; staggered SCM borrowed it only because
  `predictors` is unavailable on that path.

  Each covariate's pre-treatment average now enters the partially pooled
  objective as a second imbalance measure, pooled and per-block like the
  lagged outcomes and carrying the same cohort masses, with `xi` weighting
  the two sides. `xi` defaults to the paper's rule, the sample variance of
  the donor units' outcomes before the last adoption date; both normalisation
  constants are the combined measure at the `nu = 0` solution, as their
  footnote 7 defines. `fit$pooling` gains `q_cov_sep`, `q_cov_pool` and
  `n_covariates`, and `scm_balance_frontier()` traces the objective including
  them and returns them as columns of its own. Read its relative columns
  accordingly: they divide by the combined constant, so they are not 1 at
  `nu = 0`, and the two outcome imbalances can both fall as `nu` rises
  because the covariate side is the third quantity being traded off.

  It needs the partially pooled solver, so it is refused alongside
  `nu = "separate"`, `donor_mspe_threshold`, `lambda_pen` and
  `v_selection = "oos"` rather than falling back silently. `xi` is refused on
  the same grounds wherever the balance term is not in the objective -- a
  sharp fit, or `cov_method = "partial_out"` -- instead of being accepted and
  ignored. `"partial_out"` remains the default, so no existing fit changes.

  As `xi` falls to zero the fit collapses onto the outcomes-only one at a
  rate linear in `xi`, which is the check that the covariate side reduces
  exactly to the existing path; raising it trades outcome balance for
  covariate balance monotonically, which is the trade-off the equation
  describes.

## Bug fixes

- **Outcomes-only SCM is fast again.** Since 0.4.1 the default fit and its
  in-space placebo test had run about 3 times slower on panels with many
  pre-treatment periods, and 1.0.0's Wolfe default had made it 10 to 14
  times slower. Two causes, both fixed without changing any result (weights,
  estimates and placebo statistics agree with the previous build to 1e-13):
  the exact inner solve that seeds the coordinate sweep started the active
  set from all donors at once, paying a face solve per donor it dropped, and
  now reaches the same optimum from the support a gradient presolve finds;
  and the sweep over a panel with at least half as many periods as donors
  now carries the donor Gram matrix with rank-1 updates, where the
  predictor-oriented factor form rebuilt each face's Gram from scratch. On a
  100-donor, 96-period panel a fit takes 0.013 s and its placebo test
  0.4 s, against 0.19 s and 5.6 s before.

- **The partially pooled staggered SCM solver no longer decomposes an
  N-by-N matrix per block update.** Each block's QP matrix is the Gram of a
  factor with one row per lag, whose small Gram gives the same Lipschitz
  constant; `solve_simplex_qp()` now accepts it. Estimates are unchanged (to
  1e-14); a 200-unit fit takes 0.08 s instead of 0.13 s, and its jackknife 3 s
  instead of 27 s.

- **A staggered fit with `covariates` reports its treated series on the scale
  it was fitted on.** The covariates are partialled out first, and every other
  series the fit carries -- the per-cohort targets, `Y_all`/`Y_mat`, and the
  sharp path's own `Y_treat` -- is on that residualised scale, but the
  staggered `Y_treat` was the raw outcome. `treated_outcomes()` therefore
  returned one series while `augment()` returned another. Affects `"scm"` and
  `"sdid"`; without `covariates` the two scales coincide and nothing changes.

- **`mspe_ratio_pval()` names `mspe_post_placebo` by donor.** Its companion
  `mspe_pre_placebo` was named and it was not, so the two could not be lined
  up by name although the documentation presents them as a pair.

- **`loo_donors()` re-solves with the inner solver the fit used.** It always
  called the active-set method, so on a degenerate predictor spec -- where the
  optimal set is a face rather than a point and the two solvers pick different
  points of it -- dropping a donor also changed which optimum was reported,
  mixing the two effects the exercise is meant to separate.

- **`conformal_inference()` no longer substitutes a retired rule for a tuning
  parameter a fit does not carry.** The refit has to impose the model the fit
  used, but `r`, `k` and `lambda` each fell back to the fixed heuristic that
  was this release's default before it was replaced -- `r = 2`, `k =
  floor(sqrt(min(T_pre, N_co)))`, 1% of the spectral norm -- which would score
  the null against a differently specified counterfactual. Every current fit
  records its own, so the fallback was unreachable; it is now an error naming
  the parameter. The fields are also read exactly rather than with `$`, which
  would partial-match `r` to `r_cv`, or `lambda` to `lambda_cv`, if the exact
  name were ever absent.

# coresynth 0.5.0

## New features

- **`linetypes` in `plot.coresynth()` and `plot.scm_placebo()`** makes the
  line type of each series settable, the third member of the
  `colors`/`labels` family and keyed the same way. Line types were previously
  hard-coded, so the treated series was always solid and the synthetic
  control always dashed; drawing both solid -- the usual choice when the two
  paths are already distinguished by color, or when the figure is going into
  a document that separates series by color alone -- meant rebuilding the
  plot from `plot_data()`.

  For `type = "trend"` it takes a named vector over `"treated"`,
  `"synthetic"`, and (with `show_donors > 0`) `"donors"`, e.g.
  `linetypes = c(synthetic = "solid")`; series not mentioned keep their
  default. For `type = "gap"` it takes a single line type for the gap line,
  mirroring how `colors` takes a single color there. For
  `plot.scm_placebo(type = "gaps")` the keys are `"treated"` and
  `"placebo"`. Values are `ggplot2` line type names (`"solid"`, `"dashed"`,
  `"dotted"`, ...) or the equivalent integer codes `0:6`; an unknown key
  errors with the list of valid keys, as for `colors` and `labels`. Bar
  charts (`"weights"`, `"pred_weights"`) and the point chart
  (`"ratios"`) ignore the argument.

  All defaults reproduce the previous appearance exactly. The argument sits
  after `labels` in both signatures, so the arguments that follow it
  (`vline`, `vline_offset`, `hline`, ...) shift by one position for anyone
  passing them positionally rather than by name.

## Changed

- **A near miss on the 0.4.1 outcomes-only routing is now reported.**
  Whether a per-period outcome spec is the outcomes-only fit depends on
  the pre-treatment periods `data` carries, not on the span of the
  `pred()` list, so the same `predictors` argument switches paths when
  `data` reaches further back -- widening the panel for a plot is enough.
  Sharp SCM fits whose predictors are all single-period outcome lags but
  do not cover the pre-treatment window exactly now say so, report the
  coverage, and point at `v_window` and at subsetting `data` as the two
  ways to move the window the fit is judged on. Both fits were already
  correct for the model each one states; what was missing was any sign
  that the model had changed.

- **`pred()` time windows are checked against the panel.** A window
  naming times the panel does not carry silently aggregated the overlap
  alone while keeping its full label, and a window reaching past the
  treatment date silently aggregated values the treatment had already
  moved. Both now warn, naming the offending times. A window with no
  overlap at all keeps the existing hard error. `scm_design()` gained
  the same non-finite predictor guard `scm_fit()` already had -- such a
  window previously surfaced as "no feasible treatment set found", which
  named neither the predictor nor the window.

# coresynth 0.4.1

## Changed

- **Sharp SCM fits with `predictors = NULL` now emit a one-line message
  naming the default model** ("using the outcome in each of the k
  pre-treatment periods as predictors"), so the implied specification is
  visible rather than only documented. The message is suppressed where it
  would be noise: staggered fits (always outcomes-only, no alternative) and
  `v_selection = "oos"` (its own split is described under that argument).
  Wrap the call in `suppressMessages()` to silence it.

- **A predictor specification consisting solely of the outcome at each
  single pre-treatment period is now fitted through the outcomes-only
  path**, returning a fit identical to `predictors = NULL` (with a message).
  Such a spec builds the same X0/X1 as the outcomes-only fit -- it is the
  standard way reference implementations express an outcomes-only SCM
  (`Synth` requires it: with no predictors at all, `dataprep()` errors) --
  but previously it entered the predictor path (per-row SD scaling +
  multi-start V search), which could return slightly different weights for
  the same model and ran orders of magnitude slower. Verified on
  Proposition 99: `predictors = NULL`, the per-period spec, and `Synth`'s
  per-period `special.predictors` run now agree (weights within 3e-3 of
  `Synth`, identical between the two coresynth forms). Pass
  `v_optim = "multistart"` explicitly to force the predictor-path optimiser
  for such a spec.

## New features

- **`scm_fit(qp_solver = "wolfe")` selects a Caratheodory-sparse inner
  solver.** `X0 w` traces the convex hull of the donor columns in the
  k-dimensional predictor space as `w` ranges over the simplex, so the inner
  QP is a projection onto a polytope in `R^k`, and by Caratheodory's theorem
  an optimum supported on at most `k + 1` donors always exists. The Wolfe
  (1976) min-norm-point method returns such a solution: it grows a small
  affinely independent corral one donor at a time instead of shedding
  coordinates from the full pool, terminates finitely, and is KKT-verified,
  so it is an exact solver and not an approximation.

  Use it when the predictor specification has fewer predictors than donors
  in the support. There the inner QP's optimal set is a face rather than a
  point, the default active-set solver breaks the tie on round-off (so its
  weights are not stable across arithmetic), and the weights it returns can
  spread over most of the donor pool -- 221 of 240 donors in one measured
  case -- which is neither interpretable nor reproducible. Where the QP is
  *not* degenerate the two solvers agree: on Proposition 99 an outcomes-only
  fit returns the same six donors and the same pre-treatment loss to
  1e-15 under either.

  Opt-in for now, since it changes reported weights for degenerate specs;
  it is intended to become the default in a future major release. Not
  available with `v_selection = "oos"`, `lambda_pen`, `v_optim = "bfgs"`, or
  staggered adoption, which solve their own QP outside the sharp inner
  solver; those combinations error rather than silently ignoring the
  argument. `mspe_ratio_pval()` mirrors the fit's choice so the permutation
  test stays symmetric.

## Bug fixes

- **The inner SCM QP now returns a KKT-exact optimum.**
  `scm_inner_weights_cpp()` -- which supplies the reported weights for
  `v_selection = "oos"`, `v_optim = "bfgs"`, `loo_donors()`,
  `conformal_inference()`, `scm_design()`, and the coordinate sweep's
  starting candidate -- ran accelerated projected gradient to a 1e-6 step
  tolerance with no optimality gate. On rank-deficient problems that left it
  visibly short of the optimum: measured KKT residuals ran to 1e-2 relative,
  meaning the returned weights did not solve the problem they were
  documented to solve. It now runs the same warm-started active-set solver
  the main path uses, with the gradient method only as a fallback; residuals
  across the regression suite dropped from 1e-5..1e-2 to below 1e-11.

  This moves reported weights for specifications whose inner QP is
  degenerate. Outcomes-only fits improve slightly (Proposition 99
  pre-treatment loss 7.220088 -> 7.220081). Specifications where the treated
  unit's predictor vector lies strictly inside the donor hull -- so the
  predictor distance is zero for a continuum of weights, e.g. two or three
  aggregated predictor windows over a large donor pool -- can move
  substantially, because *which* zero-distance solution is returned is
  arbitrary and the outcome path attached to it is arbitrary with it. That
  is a property of such specifications, not of the solver; `qp_solver =
  "wolfe"` above makes it visible and stable.

## Performance

- **Predictor-based SCM and its in-space placebo test are 100-10,000x
  faster on large donor pools.** The inner simplex QP no longer forms the
  N_co x N_co Gram `Q = X0' V X0`: it carries the k x N_co factor
  `B = diag(sqrt(V)) X0` instead. `Q` has rank at most k, so with few
  predictors and many donors it is maximally rank deficient -- the regime
  where forming it costs the most and says the least. The active-set face
  solve now runs in the null space of the sum constraint *inside the
  factor*, at O(k^2 m) instead of O(m^3) (same minimiser, same rank
  truncation: `pinv(B_A P) = Z pinv(B_A Z)` for any orthonormal basis Z of
  `{z : 1'z = 0}`), the dual-feasibility gradient costs O(k N_co) instead of
  O(N_co^2), and the Lipschitz bound comes from the k x k Gram instead of an
  N_co x N_co eigendecomposition. The NNLS rescue no longer needs an
  eigendecomposition at all, since the factor *is* the least-squares form it
  wants. At N_co = 240, T_pre = 20: a single fit with one predictor drops
  from 3.78 s to under 0.01 s and its 240-donor placebo test from 315 s to
  0.03 s; with three predictors, from 5.79 s to 0.05 s and 463 s to 2.1 s.
  Outcomes-only placebo runs improve more modestly (1.9 s to 1.1 s), as that
  path was never the rank-deficient one.

- **`v_optim = "multistart"` and the coordinate sweep short-circuit when
  there is a single predictor.** With k = 1 the normalised V is the single
  point V = 1: every start in the multi-start set is the same point and
  every grid candidate rescales the inner QP by a positive constant, leaving
  its argmin unchanged. The search is therefore vacuous and now reduces to
  one exact inner solve, which is also strictly more accurate than the old
  path's cold FISTA solution.

- **The outcomes-only bordered-KKT face solve now falls back to the exact
  null-space solve instead of to FISTA.** That saddle system is singular
  once the face exceeds the rank of its Gram, and the previous code then
  gave up and let the caller fall back to FISTA, which stops at its 1e-6
  tolerance -- so the outer V search read an inexact objective. It now falls
  through to the exact routes. This makes `v_selection = "oos"` measurably
  better fitting (3.5% lower validation RMSPE on Proposition 99).

  Note on reproducibility: when the inner QP is rank deficient (fewer
  predictors than donors in the support) its solution set is a face, not a
  point, and the active set's entering-variable choice breaks ties on
  round-off. Weights from predictor-based fits can therefore differ from
  0.4.0 while remaining exact optima of the same problem; the pre-treatment
  fit itself is unchanged to round-off for outcomes-only specs and moves by
  at most a few times 1e-5 for in-sample predictor specs under the default
  `v_optim`. Outcomes-only fits, staggered SCM, SDID, GSC, MC, TASC, SI,
  `scm_design()`, and `conformal_inference()` are bit-identical.

- **Predictor matrix construction is no longer a bottleneck for
  predictor-based SCM fits.** `build_predictor_matrices()` aggregated each
  predictor row with one full scan of the input data frame per unit
  (O(rows x units x data size)); with many `pred()` entries -- e.g. one
  outcome lag per pre-treatment year -- this R-level loop dominated the fit
  time (1.75 s of a 3.9 s fit at N_co = 150, T_pre = 80, 15k rows). Each
  predictor row is now computed with a single grouped aggregation
  (0.06 s for the same spec, ~29x). Output is bit-identical to the previous
  implementation across ops (`mean`/`median`/`sum`), id types, NA handling,
  and the empty-window error path, so every fit is unchanged. The same
  per-unit-scan pattern was removed from `build_covariate_array()` (the
  covariate paths of GSC, SDID, and staggered SCM) and from
  `scm_design()`'s predictor builder, likewise with bit-identical output.

- **Outcomes-only SCM fits are ~2-2.5x faster again.** Hardening the inner
  active-set QP for the predictor multi-start search (0.4.0) replaced the
  cheap bordered-KKT face solve with a scale-robust null-space projection on
  every pivot. That robustness is only needed when `predictors` are supplied
  (rank-deficient V, scale invariance under predictor rescaling); the
  outcomes-only path -- dense V, well-conditioned faces, no rescaling -- paid
  the cost without the benefit, running roughly twice as slow as 0.3.x
  despite returning the identical fit. The outcomes-only inner QP now uses
  the cheap bordered-KKT solve again (`scm_weights_cpp(cheap_face = TRUE)`,
  wired automatically for outcomes-only fits, staggered SCM, the outcomes-only
  in-space placebo, `placebo_in_time()`, and OOS V selection). Predictor-based
  fits and the multi-start search keep the null-space solver, so their results
  and scale invariance are untouched. Outcomes-only results are bit-identical
  to the 0.2/0.3 line and differ from 0.4.0 only at round-off
  (<= 5e-14 in weights), i.e. not observably. Verified against `Synth` on
  Proposition 99 (ATT, synthetic series, pre-RMSPE, donor weights all agree)
  and across a degeneracy/noise sweep.

# coresynth 0.4.0

## Breaking changes

- **Predictor-based SCM fits now default to a multi-start outer V
  optimisation** (`v_optim = "auto"`, the new default, resolves to
  `"multistart"` when `predictors` are supplied and to the single-start
  `"coord_descent"` for outcomes-only fits). The outer problem -- choose V
  so that the implied donor weights minimise pre-treatment outcome MSPE --
  is non-convex, and the previous single start from uniform V could settle
  in a poor basin: on the Proposition 99 predictor specification it
  reported a pre-treatment SSR of 137.8 where 58.5 is attainable (below the
  reference `Synth`/`tidysynth` solution's 78.1), halving the treated
  unit's MSPE ratio. The multi-start search screens a fixed,
  fully deterministic start set (uniform, one-hot per predictor, 100
  fixed-seed draws) with warm-started exact inner QPs, then runs the
  leaders through a Nelder-Mead / coordinate-descent refinement pipeline;
  its solution is never worse (in pre-treatment loss) than the single-start
  path, and no R RNG state is consumed. `mspe_ratio_pval()` applies the
  same search to every placebo refit so the permutation test stays
  symmetric. Numeric results for predictor-based fits therefore change
  (toward better pre-treatment fits); pass `v_optim = "coord_descent"` to
  reproduce the previous behaviour. Outcomes-only fits are unchanged.
  Note that the old `"auto"` (choose `"bfgs"` when `k <= 15`) is gone;
  `"bfgs"` itself remains available.

## New features

- **`plot(fit, type = "pred_weights")`** draws the predictor/variable weight
  matrix `V` of a sharp SCM fit as a horizontal bar chart, the companion
  to the existing donor `type = "weights"` chart. Bars are labelled with the
  predictor names (or `V1..V_Tpre` for outcomes-only fits) and honour the
  `fill` and `top_n` arguments. Staggered SCM and the other estimators do not
  estimate a `V` matrix and raise an informative error.
- **`plot_data()`** returns the tidy `data.frame` behind any `plot()` view,
  taking the same `type` argument (`"trend"`, `"gap"`, `"weights"`,
  `"pred_weights"` for a fit; `"gaps"`, `"ratios"` for an `scm_placebo`
  object). It is the data-side companion to `plot()`: extract the plotted
  series, weights, or placebo paths to relabel the simplified series names,
  feed them into a table, or build a bespoke figure. Columns use plain names
  (`time`, `value`, `series`, ...), and the arguments that change which rows
  or values appear (`align`, `top_n`, `show_donors`, `mspe_prune`) are
  honoured. The cosmetic near-zero-weight filter that `type = "weights"`
  applies to the chart is not used here, so every donor is returned (use
  `top_n` to subset).
- **`v_window` argument in `scm_fit()`** (sharp SCM fits): a vector of
  pre-treatment time values over which the outer V optimisation evaluates
  the pre-treatment fit, e.g. `v_window = 1975:1988`. The default (`NULL`)
  keeps the current behaviour of evaluating on all pre-treatment periods.
  The window restricts only the outer evaluation loss -- predictor matrices
  and the reported full-window `loss` are untouched -- and
  `mspe_ratio_pval()` mirrors it in every placebo refit. Cannot be combined
  with `v_selection = "oos"`, which manages its own train/validation split.
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
  also reports the placebo-distribution standard error (Clarke et al. 2024,
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

## Deprecated

- **`v_optim = "bfgs"` is deprecated** and will be removed in a future
  release. It is a single-start L-BFGS-B outer optimiser with no advantage
  over the alternatives now that the default is multi-start: use
  `v_optim = "multistart"` (or the `"auto"` default) for predictor-based
  fits, and `"coord_descent"` for outcomes-only fits. Passing `"bfgs"` still
  works but now emits a deprecation warning. `"coord_descent"` is **not**
  deprecated: besides being a selectable optimiser, it is the engine for
  outcomes-only and staggered fits and the reference against which the
  multi-start never-worse guarantee is defined.

## Performance

- **The multi-start outer search runs at interactive speed.** Outer-loss
  evaluations are KKT-gated warm-started active-set solves; the active set
  gained Bland's-rule anti-cycling (engaged only beyond 30 pivots, so
  previously terminating pivot paths are untouched) and now solves each
  face subproblem in the null space of the sum constraint via a
  Cholesky/eigendecomposition hybrid that is robust to the rank-deficient
  metrics (k < |active set|) the multi-start screen hits constantly -- and
  whose branch decision cannot be flipped by last-bit input perturbations,
  preserving the scale-invariance of `scale_predictors`. On a degenerate
  face the warm active-set solve can churn to its pivot cap, and the old
  cold-FISTA fallback (thousands of iterations on ill-conditioned metrics)
  made a single slow donor dominate placebo wall time; a finite-termination
  Lawson-Hanson NNLS now seeds the exact face solve on those cases in
  microseconds, gated by the same KKT check so the returned solution is
  unchanged. The rescue is confined to the multi-start internals -- every
  other path (outcomes-only, `coord_descent`, `oos`, staggered,
  `loo_donors()`, conformal) keeps its historical solver bit-for-bit, and
  the multi-start uniform-start leg stays on the plain path so the
  never-worse guarantee holds. Candidate refinement pipelines run in
  parallel, and the placebo battery flattens the (donor x candidate) task
  list into one schedule -- donor-level parallelism alone lets the slowest
  donor pin a thread while the others go idle. On the Proposition 99
  8-predictor spec a full multi-start fit takes ~0.06s and the complete
  38-donor `mspe_ratio_pval()` battery ~1.7s (every donor's fit verified
  never worse than the single-start path). A side effect of the exact face solves:
  inner QPs on rank-deficient faces are now solved exactly where the
  previous code fell back to a loosely-converged FISTA iterate.
  Outcomes-only sharp, staggered, and placebo results are bit-identical to
  the previous release; `v_selection = "oos"` outcomes fits (where the
  training half has fewer rows than donors) can shift slightly, toward
  better fits.
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
