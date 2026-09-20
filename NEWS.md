# coresynth 1.0.0

## Breaking changes

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
  - `"si"`: `k = NULL` applies the optimal hard threshold of Donoho & Gavish
    (2014) to the pre-treatment donor block, one of the three data-driven
    rules Agarwal et al. name. It costs a single SVD, where cross-validating
    `k` would cost a PCR fit per candidate. On the same known-rank panels it
    recovered the rank exactly for 1 through 4; the old heuristic returned 4
    in every case.

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
  series in `Y_treat_mat`, alongside `n_treated`. Donor weights stay a single
  simplex vector, so `predictors`, `v_window`, `lambda_pen`,
  `v_selection = "oos"`, `plot()`, `placebo_in_time()`, `loo_donors()` and
  `augment_scm()` all work as they do for one treated unit.

  `mspe_ratio_pval()` does not, and now says so. It compares the treated unit
  against placebo runs on single donors, and an average of several treated
  units is less noisy than any one donor, which inflates its pre/post MSPE
  ratio. Simulated under no treatment effect, the test rejected at 28% with
  three pooled units and 53% with six, against a nominal 10%.

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
  numeric `nu` there is still an error. The heuristic itself was checked
  against `augsynth` and agrees to three or four decimal places.

- **`qp_solver` now defaults to `"auto"`, which is the Wolfe solver wherever
  it is wired.** The inner SCM problem projects the treated unit onto the
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
  outer V search a cleaner objective to descend. Outcomes-only fits, where
  the QP is not degenerate, are unchanged to machine precision.

  `"auto"` falls back to `"active_set"` on the paths that run their own QP
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
  runs. Pass `mspe_prune = Inf` to draw them all regardless.

- **`pool_cohorts = FALSE` fits one synthetic control per treated unit.**
  Staggered SCM pools units sharing an adoption date into a cohort, which is
  the Appendix A.2 hierarchy of Ben-Michael, Feller & Rothstein (2022). Their
  S.4 formulation instead gives every treated unit its own weight vector, and
  that is what `pool_cohorts = FALSE` fits; `cohort_fits` then names the unit
  behind each block, which is what `scm_rank_pval()` needs. The default stays
  with the appendix, where S.4.1 shows there is no across-unit variation to
  exploit inside a cohort. (`augsynth` defaults the other way.)

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
