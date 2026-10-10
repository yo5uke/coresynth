#' Fit a Synthetic Control Method Model
#'
#' Unified formula interface for Synthetic Control and related causal
#' inference methods.  The formula syntax is:
#'
#'   `outcome ~ treatment | unit_id + time_id`
#'
#' @param formula A `Formula` object, e.g. `y ~ D | unit + time`.
#' @param data    A `data.frame` in **long** format (one row per unit-time).
#' @param method  One of `"scm"`, `"sdid"`, `"gsc"`, `"mc"`, `"tasc"`, `"si"`.
#' @param predictors A `list()` of [pred()] specifications that define the
#'   predictor matrix for SCM (see Abadie et al. 2010, S.2.3). Each [pred()]
#'   entry aggregates one or more variables over a time window. Pass `NULL`
#'   (default) to use all pre-treatment outcome periods as predictors; a
#'   one-line message states this default when it applies.
#'   A specification consisting solely of the outcome variable at each
#'   single pre-treatment period (one `pred()` per period, jointly covering
#'   the full pre-treatment window) defines the same predictor matrix and is
#'   therefore fitted through the same outcomes-only path as `NULL`,
#'   returning an identical fit; supply `v_optim = "multistart"` explicitly
#'   to force the predictor-path optimiser instead. That window is whatever
#'   `data` carries before the treatment date, not the span of the `pred()`
#'   list: adding earlier periods to `data` (to widen a plot, say) leaves the
#'   same list short of full coverage and moves the fit onto the predictor
#'   path, so a message reports the mismatch when it happens.
#'   Restricting the periods the fit is judged on is a separate choice from
#'   `predictors` -- subset `data`, or keep `data` whole and pass
#'   `v_window`.
#'   `pred()` windows are expected to lie inside the pre-treatment periods
#'   present in `data`: times the panel does not carry, and times at or after
#'   the treatment date, each raise a warning because the row then aggregates
#'   something other than what the window reads as.
#'   Applies to `method = "scm"` only. Predictor rows are scaled by their
#'   standard deviation across all units before optimisation, matching the
#'   Synth reference implementation (ADH 2011, JSS); pass
#'   `scale_predictors = FALSE` to disable.
#'   When several units adopt at the same date they are pooled into one
#'   cohort, and each predictor row of the treated vector is their average
#'   under the same `weights` the outcome rows are pooled with, so both sides
#'   of the fit describe one target series.
#' @param covariates An optional named `list` of additional time-varying
#'   covariates to partial out before estimation. Each element is a character
#'   string naming a column in `data`. Supported for `method = "sdid"`,
#'   `"gsc"`, and **staggered** `method = "scm"`.
#'   Partialling out is the device of Clarke et al. (2024, S.2.2), which SDID
#'   needs because it has no way to balance a covariate. SCM does: Abadie,
#'   Diamond & Hainmueller (2010, S.2.3) put covariates into the V-weighted
#'   predictor distance, and Ben-Michael, Feller & Rothstein (2022, S.5.2) add
#'   a separate balance term for them. A sharp SCM fit therefore rejects
#'   `covariates` and asks for `predictors = list(pred(...))` instead; the
#'   staggered path borrows the partial-out only because `predictors` is
#'   unavailable there.
#' @param v_selection V matrix selection method for `method = "scm"`.
#'   `"insample"` (default) follows Abadie et al. (2010): V is chosen by
#'   minimising in-sample pre-treatment MSPE. `"oos"` follows Abadie (2021)
#'   S.3.2 / ADH (2015): the pre-treatment window is split into a training
#'   half and a validation half. In the default outcomes-as-predictors case,
#'   candidate W(V) are fitted on training-half outcomes only, V* minimises
#'   the validation-half MSPE, and the final W* is refit with V* on the
#'   outcomes of the last `floor(T_pre/2)` pre-treatment periods (so
#'   `v_weights` has `floor(T_pre/2)` entries). With user-supplied
#'   `predictors`, the predictor matrix is fixed and only the MSPE evaluation
#'   window is restricted to the validation half; lag your [pred()] windows
#'   to the training period for a fully out-of-sample exercise.
#' @param donor_mspe_threshold Donor pool filtering threshold (Abadie 2021 S.4).
#'   For `method = "scm"` only. Each donor's individual pre-treatment MSPE
#'   (using that donor alone as the counterfactual) is divided by the minimum
#'   such MSPE across all donors. Donors whose ratio exceeds this threshold are
#'   excluded from estimation. `Inf` (default) disables filtering.
#' @param lambda_pen Penalised SCM parameter (Abadie & L'Hour 2021, JASA).
#'   For `method = "scm"` only. `NULL` (default) runs standard unpenalised SCM.
#'   `"auto"` selects the penalty via out-of-sample pre-treatment MSPE on the
#'   same validation window as `v_selection = "oos"`. A non-negative number
#'   uses that value directly.
#' @param v_optim Outer V-optimisation method for `method = "scm"`.
#'   The outer problem (choose V so that the implied donor weights W(V)
#'   minimise pre-treatment outcome MSPE) is non-convex, and a single local
#'   search can settle in a poor basin when `predictors` are supplied.
#'   `"auto"` (default) therefore selects `"multistart"` for predictor-based
#'   fits and `"coord_descent"` for outcomes-only fits (where the single
#'   start is empirically reliable). `"multistart"` runs a deterministic
#'   multi-start search: a fixed start set (uniform V, one-hot V per
#'   predictor, and 100 fixed-seed random draws) is screened at one inner QP
#'   each, the leaders are polished by coordinate descent, and the winner is
#'   refined by a Nelder-Mead pass. Its solution is never worse (in
#'   pre-treatment loss) than `"coord_descent"`, at roughly the cost of a
#'   handful of single-start fits, and is fully reproducible (no RNG state
#'   is consumed). `"coord_descent"` is the classic single-start coordinate
#'   descent with an 11-point grid, and is also the outcomes-only and
#'   staggered engine and the multi-start never-worse reference.
#'   `"bfgs"` (a single-start L-BFGS-B) is **deprecated** and will be removed
#'   in a future release: it has no advantage over `"multistart"` (which
#'   dominates it on predictor fits) or `"coord_descent"`.
#'   [mspe_ratio_pval()] mirrors a multi-start fit in its placebo refits so
#'   the permutation test stays symmetric.
#' @param qp_solver Inner-QP solver for `method = "scm"`.
#'   `"wolfe"` is the Wolfe (1976) min-norm-point method, which exploits the
#'   fact that the inner QP is a projection onto the convex hull of the donor
#'   columns in the k-dimensional predictor space: by Caratheodory's theorem
#'   an optimum supported on at most `k + 1` donors always exists, and Wolfe
#'   returns such a solution. `"active_set"` is the warm-started active-set
#'   method. Both return an exact optimum; they differ in *which* optimum when
#'   the QP is degenerate (fewer predictors than donors in the support), where
#'   the optimal set is a face rather than a point. `"wolfe"` is preferred
#'   there because its solution is sparse and reproducible, while
#'   `"active_set"` breaks those ties on round-off and can return different
#'   weights for the same fit under a different build.
#'   `"auto"` (default) is `"wolfe"` where that degeneracy can arise -- sharp
#'   fits with `predictors`, and outcomes-only fits with fewer than half as
#'   many pre-treatment periods as donors -- and `"active_set"` for the other
#'   outcomes-only fits, whose QP has a unique optimum in practice, so both
#'   solvers return the same weights and the active set gets there several
#'   times faster. It is also `"active_set"` on the paths that run their own
#'   solver: `v_selection = "oos"`, `lambda_pen`, `v_optim = "bfgs"`, and
#'   staggered adoption. Naming `"wolfe"` explicitly
#'   on one of those paths is an error rather than a silent downgrade. The
#'   resolved value is recorded in `fit$qp_solver`.
#' @param v_window Optional vector of pre-treatment time values (matching the
#'   time index in `data`) over which the outer V optimisation evaluates the
#'   pre-treatment fit, for `method = "scm"` (sharp fits only). `NULL`
#'   (default) evaluates on all pre-treatment periods. The window restricts
#'   only the outer evaluation loss: predictor matrices (or the full
#'   pre-treatment outcome rows in the outcomes-only case) still enter the
#'   inner QP unchanged, and the reported `loss` and [mspe_ratio_pval()]
#'   MSPE components always cover the full pre-treatment window. Cannot be
#'   combined with `v_selection = "oos"`, which manages its own
#'   train/validation split. Use it to judge the fit on part of the
#'   pre-treatment period while keeping the rest of `data` available for
#'   plotting and post-estimation.
#' @param nu How strongly the adoption cohorts are pooled in **staggered**
#'   SCM fits (Ben-Michael, Feller & Rothstein 2022, JRSS-B). Cohort weight
#'   vectors are chosen jointly to minimise
#'   `nu * (normalised pooled pre-treatment imbalance)^2 +
#'   (1 - nu) * (normalised per-cohort imbalance)^2`, trading the fit of each
#'   cohort against the fit of their average, which is what the aggregate ATT
#'   rests on.
#'   `"auto"` (the default) is the paper's heuristic, the ratio of the pooled
#'   to the average per-cohort imbalance of the separate solution; it is also
#'   what `augsynth`, the authors' implementation, uses when `nu` is unset.
#'   A number in `[0, 1]` sets it directly, with `0` separate and `1` fully
#'   pooled -- the paper advises against both extremes. `"separate"` is
#'   `0` reached through the older per-cohort path, which keeps the
#'   V optimisation and donor filtering that the pooled solver has no room
#'   for; it is also what the default falls back to when
#'   `donor_mspe_threshold`, `lambda_pen`, or `v_selection = "oos"` rules the
#'   pooled solver out. Balance diagnostics are stored in `fit$pooling`.
#'
#'   Staggered panels only. Under a single adoption date there are no cohorts
#'   to pool across, and several units sharing one date are already fully
#'   pooled into that cohort: S.4.1 of the paper shows the across-unit
#'   variation that would justify doing otherwise is zero when the adoption
#'   date is shared.
#' @param pool_cohorts For **staggered** `method = "scm"`: whether treated
#'   units sharing an adoption date are averaged into one block before the
#'   weights are fitted. `TRUE` (default) is the Appendix A.2 hierarchy
#'   of Ben-Michael, Feller & Rothstein (2022), whose S.4.1 shows there is no
#'   across-unit variation to exploit within a cohort. `FALSE` is their S.4
#'   formulation, one weight vector per treated unit, which is what the
#'   rank-based permutation test in [scm_rank_pval()] needs. `augsynth` makes
#'   the unit-level form its default; the paper's own appendix does not.
#'   It applies to both staggered paths -- the partially pooled solver and
#'   `nu = "separate"` -- so it still takes effect where `donor_mspe_threshold`,
#'   `lambda_pen` or `v_selection = "oos"` have sent the fit down the latter.
#'   Sharp fits have a single cohort and reject it.
#' @param weights Optional named numeric vector over the treated units giving
#'   the mass each carries in the aggregate ATT, for `method = "scm"`. `NULL`
#'   (default) weights them equally, which for a shared adoption date is the
#'   `N_treated * T_post` weighting used throughout. Abadie (2021, S.8) notes
#'   that with regions or countries a population- or GDP-weighted average
#'   effect is often the quantity of interest; supply those figures here.
#'   The mass is used consistently: a block's target series is the
#'   correspondingly weighted average of its treated units, that same mass is
#'   what the block counts for in the partially pooled objective (in place of
#'   the cohort size of Ben-Michael, Feller & Rothstein 2022, Appendix A.2)
#'   and in its `T_post`-multiplied aggregation weight, and
#'   [scm_inference()] resamples each unit at its own mass. Only the ratios
#'   matter -- the vector is rescaled internally -- so equal weights reproduce
#'   the unweighted fit exactly.
#' @param fixedeff Unpenalised additive level terms. For `method = "mc"` it
#'   estimates unit and time fixed effects outside the nuclear norm, as in
#'   Athey et al. (2021, S.8.1), and defaults to `TRUE`; `FALSE` recovers the
#'   S.4 baseline where the level is absorbed into the penalised matrix.
#'   For `method = "scm"` it defaults to `TRUE` on staggered panels, matching
#'   `augsynth` and the S.6 finding that the intercept-shifted estimator
#'   "dominates either of the alternatives"; it does not apply to sharp fits.
#'   If `TRUE`, staggered SCM demeans every unit by its own
#'   pre-treatment mean within each cohort before fitting (intercept shift;
#'   Ben-Michael, Feller & Rothstein 2022, Section 5.1), which turns the
#'   estimator into a weighted difference-in-differences and typically
#'   improves fit when outcome levels differ across units. The reported
#'   `Y_synth` is shifted back to the raw outcome scale. Works with both
#'   the default and the partially pooled path. For `method = "scm"` on
#'   staggered panels only.
#' @param ...     Estimator-specific arguments forwarded to the chosen method:
#'   `lambda` plus `cv_folds`, `cv_nlambda` and `cv_seed` (`"mc"`), `r`
#'   (`"gsc"`, `"tasc"`), `em_iter` and `fix_A` (`"tasc"`), `zeta2` and
#'   `zeta_t` (`"sdid"`), `k` (`"si"`), `control_group` (`"scm"`, `"sdid"`,
#'   `"gsc"`, `"si"`), and `scale_predictors`, `cov_method` and `xi`
#'   (`"scm"`). A name the chosen
#'   estimator does not take is an error rather than a silently ignored
#'   argument. Every argument after `...` must be named in full: this is what
#'   lets `lambda` reach `method = "mc"` instead of partially matching
#'   `lambda_pen`.
#'
#'   `control_group` chooses each cohort's donors in a **staggered** fit of
#'   `"scm"`, `"sdid"`, `"gsc"` or `"si"`. `"never_treated"` (default) uses
#'   the units never treated in the panel, the pure controls of Clarke et al.
#'   (2024, Algorithm 1); at least two are required. `"clean"` is deprecated:
#'   it adds the units adopting after the cohort, but each cohort's effect runs
#'   to the end of the panel, so those units' post-adoption outcomes enter its
#'   counterfactual and bias the effect towards zero. `"mc"` and `"tasc"`
#'   impute from every untreated cell and need no never-treated units.
#'
#'   `cov_method` chooses how `covariates` enter a **staggered** `"scm"` fit.
#'   `"partial_out"` (default) residualises the outcome on them first, the
#'   device of Clarke et al. (2024, S.2.2). `"balance"` is the formulation of
#'   Ben-Michael, Feller & Rothstein (2022, S.5.2, Equation 11): each
#'   covariate's pre-treatment average enters the partially pooled objective
#'   as a second imbalance measure, weighted against the lagged outcomes by
#'   `xi`. It requires the partially pooled solver, so it cannot be combined
#'   with `nu = "separate"`, `donor_mspe_threshold`, `lambda_pen` or
#'   `v_selection = "oos"`. `xi` defaults to the paper's rule, the sample
#'   variance of the donor units' outcomes before the last adoption date.
#'
#'   Three tuning parameters that govern model complexity are selected from
#'   the data when left at their `NULL` default, each by the rule its source
#'   paper specifies: `r` for `"gsc"` by the leave-one-period-out
#'   cross-validation of Xu (2017, Algorithm 1), `lambda` for `"mc"` by the
#'   K-fold cross-validation of Athey et al. (2021), and `k` for `"si"` by the
#'   optimal hard threshold of Gavish & Donoho (2014). Supplying a value skips
#'   the search. The `r` of `"tasc"` is a different quantity -- the dimension
#'   of the state vector its Kalman smoother carries, not a factor count
#'   estimated from a donor block -- and keeps its fixed default of `2`.
#'
#' @return An object of classes `c("coresynth_<method>", "coresynth")`.
#'   Fits with staggered adoption additionally inherit from
#'   `"coresynth_staggered"`, and multi-arm SI fits from
#'   `"coresynth_multiarm"`; S3 methods such as [tidy()] and [augment()]
#'   dispatch on these subclasses.
#'   All methods return at minimum:
#'   * `method`: estimator name
#'   * `estimate`: average treatment effect (ATT)
#'   * `times`: time index vector
#'   * `T_pre`: number of pre-treatment periods
#'   * `Y_treat`: treated unit outcome series
#'   * `gap`: treatment effect series (Y_treat - counterfactual)
#'
#'   Staggered `method = "scm"` fits add `cohort_estimates`, `cohort_fits`,
#'   and the settings the fit resolved: `nu` (the number the partially pooled
#'   solver used, or `"separate"` for the per-cohort path), `fixedeff`,
#'   `pool_cohorts` and `qp_solver`.
#'
#' @export
#'
#' @examples
#' # Synthetic balanced panel: 10 units over 20 periods, unit 1 treated
#' # after period 15.
#' set.seed(1)
#' panel <- expand.grid(unit = 1:10, year = 1:20)
#' panel$treated <- as.integer(panel$unit == 1 & panel$year > 15)
#' panel$gdp <- panel$unit + 0.5 * panel$year +
#'   rnorm(nrow(panel)) + 3 * panel$treated
#'
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "sdid")
#' summary(fit)
#'
#' \donttest{
#' # Visualise the estimated gap (requires ggplot2)
#' plot(fit, type = "gap")
#' }
scm_fit <- function(
  formula,
  data,
  method = c("scm", "sdid", "gsc", "mc", "tasc", "si"),
  ...,
  predictors = NULL,
  covariates = NULL,
  v_selection = c("insample", "oos"),
  donor_mspe_threshold = Inf,
  lambda_pen = NULL,
  v_optim = c("auto", "coord_descent", "bfgs", "multistart"),
  qp_solver = c("auto", "active_set", "wolfe"),
  v_window = NULL,
  nu = NULL,
  fixedeff = NULL,
  pool_cohorts = TRUE,
  weights = NULL
) {
  v_selection <- match.arg(v_selection)
  v_optim     <- match.arg(v_optim)
  qp_solver   <- match.arg(qp_solver)
  method <- match.arg(method)

  if (identical(v_optim, "bfgs")) {
    warning(
      "v_optim = \"bfgs\" is deprecated and will be removed in a future ",
      "release. It is a single-start outer optimiser with no advantage over ",
      "the alternatives: use v_optim = \"multistart\" (or the \"auto\" ",
      "default) for predictor-based fits, and \"coord_descent\" for ",
      "outcomes-only fits.",
      call. = FALSE
    )
  }

  if (!is.null(v_window) && method != "scm") {
    stop("'v_window' applies to method = \"scm\" only.", call. = FALSE)
  }

  .check_method_args(
    method    = method,
    supplied  = setdiff(names(match.call(expand.dots = FALSE)), c("", "...")),
    dot_names = names(list(...))
  )

  # Parse Formula
  f_parts <- Formula::Formula(formula)

  if (length(f_parts)[2] < 2) {
    stop(
      "Formula must specify unit and time after '|', e.g. y ~ D | unit + time"
    )
  }

  # The outcome and treatment are read from `data` by column name, so a term
  # such as log(y) or I(t > 1990) would otherwise be reduced to its variable
  # and fitted untransformed.
  for (part in list(list(lhs = 1, rhs = 0, what = "outcome"),
                    list(lhs = 0, rhs = 1, what = "treatment"))) {
    f_part <- formula(f_parts, lhs = part$lhs, rhs = part$rhs)
    term   <- f_part[[2L]]   # lhs of y ~ 0, rhs of ~ D
    if (!is.name(term))
      stop(sprintf(paste0(
        "The %s must be a column name, not the expression '%s'. Create the ",
        "transformed variable in `data` first, e.g. ",
        "data$log_y <- log(data$y)."),
        part$what, deparse1(term)), call. = FALSE)
  }

  y_var <- all.vars(formula(f_parts, lhs = 1, rhs = 0))
  d_var <- all.vars(formula(f_parts, lhs = 0, rhs = 1))
  idx_vars <- all.vars(formula(f_parts, lhs = 0, rhs = 2))

  if (length(y_var) != 1) {
    stop("Exactly one outcome variable required.")
  }
  if (length(d_var) != 1) {
    stop("Exactly one treatment variable required.")
  }
  if (length(idx_vars) != 2) {
    stop("Must specify exactly two index variables: unit_id and time_id.")
  }

  id_var <- idx_vars[1]
  time_var <- idx_vars[2]

  # Validate columns exist
  for (v in c(y_var, d_var, id_var, time_var)) {
    if (!v %in% names(data)) {
      stop(paste0("Variable '", v, "' not found in data."))
    }
  }

  y_raw <- data[[y_var]]
  d_raw <- data[[d_var]]

  # as.numeric()/as.integer() on a factor silently return level codes, so
  # factor/character columns must be rejected rather than coerced.
  if (is.factor(y_raw) || is.character(y_raw)) {
    stop(sprintf(paste0(
      "Outcome variable '%s' is of type %s; it must be numeric. Convert it ",
      "explicitly, e.g. as.numeric(as.character(x))."),
      y_var, class(y_raw)[1L]), call. = FALSE)
  }
  if (is.factor(d_raw) || is.character(d_raw)) {
    stop(sprintf(paste0(
      "Treatment variable '%s' is of type %s; it must be numeric, integer, ",
      "or logical (0 = control, 1 = treated). Convert it explicitly, e.g. ",
      "as.integer(as.character(x))."),
      d_var, class(d_raw)[1L]), call. = FALSE)
  }
  if (is.numeric(d_raw) &&
      any(abs(d_raw - round(d_raw)) > 1e-8, na.rm = TRUE)) {
    stop(sprintf(paste0(
      "Treatment variable '%s' contains non-integer values; it must be 0/1 ",
      "(or 0..K integer arms for method = 'si')."),
      d_var), call. = FALSE)
  }

  y <- as.numeric(y_raw)
  d <- as.integer(round(as.numeric(d_raw)))
  id <- data[[id_var]]
  time <- data[[time_var]]

  # Dispatch
  res <- switch(
    method,
    "scm" = fit_scm_cpp(
      y,
      d,
      id,
      time,
      data = data,
      id_var = id_var,
      time_var = time_var,
      outcome_var = y_var,
      predictors = predictors,
      covariates = covariates,
      v_selection = v_selection,
      donor_mspe_threshold = donor_mspe_threshold,
      lambda_pen = lambda_pen,
      v_optim = v_optim,
      qp_solver = qp_solver,
      v_window = v_window,
      nu = nu,
      fixedeff = fixedeff,
      pool_cohorts = pool_cohorts,
      weights = weights,
      ...
    ),
    "sdid" = fit_sdid_cpp(
      y, d, id, time,
      covariates = covariates,
      data       = data,
      id_var     = id_var,
      time_var   = time_var,
      ...
    ),
    "gsc" = fit_gsc_cpp(
      y,
      d,
      id,
      time,
      data = data,
      id_var = id_var,
      time_var = time_var,
      covariates = covariates,
      ...
    ),
    "mc" = fit_mc_cpp(y, d, id, time,
                      fixedeff = fixedeff %||% TRUE, ...),
    "tasc" = fit_tasc_cpp(y, d, id, time, ...),
    "si" = fit_si_cpp(y, d, id, time, ...),
    stop(paste0("Unknown method: '", method, "'"))
  )

  new_coresynth(res, method)
}

# scm_fit() is one interface over six estimators, so an argument meant for a
# different one -- or simply misspelled -- reaches a fit function that has no
# formal for it and is absorbed by its `...`, leaving the fit silently
# unchanged. These two tables turn that into an error.
#
# Arguments that scm_fit() takes as its own formals but only forwards to some
# estimators:
.SCM_FIT_METHOD_ONLY <- list(
  predictors           = "scm",
  v_selection          = "scm",
  donor_mspe_threshold = "scm",
  lambda_pen           = "scm",
  v_optim              = "scm",
  qp_solver            = "scm",
  nu                   = "scm",
  pool_cohorts         = "scm",
  weights              = "scm",
  fixedeff             = c("scm", "mc"),
  covariates           = c("scm", "sdid", "gsc")
)

# Estimator-specific arguments travel through `...`; anything not a formal of
# the chosen fit function is rejected. The plumbing that scm_fit() fills in
# itself, and its own formals, are excluded so they never show up as choices.
.method_dot_args <- function(method) {
  fit_fun <- switch(
    method,
    "scm"  = fit_scm_cpp,  "sdid" = fit_sdid_cpp, "gsc" = fit_gsc_cpp,
    "mc"   = fit_mc_cpp,   "tasc" = fit_tasc_cpp, "si"  = fit_si_cpp
  )
  setdiff(
    names(formals(fit_fun)),
    c("...", "y", "d", "id", "time", "data", "id_var", "time_var",
      "outcome_var", names(formals(scm_fit)))
  )
}

.check_method_args <- function(method, supplied, dot_names) {
  for (nm in intersect(supplied, names(.SCM_FIT_METHOD_ONLY))) {
    ok <- .SCM_FIT_METHOD_ONLY[[nm]]
    if (!method %in% ok)
      stop(sprintf(
        "'%s' applies to method = %s only, not \"%s\", where it would be ignored.",
        nm, paste0("\"", ok, "\"", collapse = " / "), method), call. = FALSE)
  }

  if (length(dot_names) == 0L) return(invisible(NULL))
  if (any(!nzchar(dot_names)))
    stop("arguments passed through `...` to scm_fit() must be named.",
         call. = FALSE)

  known <- .method_dot_args(method)
  bad   <- setdiff(dot_names, known)
  if (length(bad))
    stop(sprintf(
      "unused argument(s) for method = \"%s\": %s. %s",
      method, paste(bad, collapse = ", "),
      if (length(known))
        paste0("This method takes ", paste(known, collapse = ", "),
               " beyond the shared arguments.")
      else "This method takes no estimator-specific arguments."),
      call. = FALSE)
  invisible(NULL)
}

# Internal constructor: attach the class tags that drive S3 dispatch.
# Structural variants get their own subclass so methods (print, tidy,
# augment, ...) can dispatch instead of branching on isTRUE(res$staggered) /
# isTRUE(res$multi_arm). "coresynth_multiarm" precedes "coresynth_staggered"
# so multi-arm methods win and can delegate via NextMethod().
new_coresynth <- function(res, method) {
  cls <- character(0L)
  if (isTRUE(res$multi_arm)) cls <- c(cls, "coresynth_multiarm")
  if (isTRUE(res$staggered)) cls <- c(cls, "coresynth_staggered")
  class(res) <- c(cls, paste0("coresynth_", method), "coresynth")
  res
}

#' @export
print.coresynth <- function(x, ...) {
  cat("=== coresynth fit ===\n")
  cat("Method :", toupper(x$method), "\n")
  cat("Estimate (ATT):", round(x$estimate, 4), "\n")
  .cat_pooled_treated(x)
  cat("Pre-treatment periods:", x$T_pre, "\n")
  invisible(x)
}

# Several units adopting at one date are pooled into a single synthetic
# control, so say so: otherwise nothing in the output distinguishes that fit
# from a one-treated-unit one.
.cat_pooled_treated <- function(x) {
  n <- x$n_treated %||% 1L
  if (n > 1L)
    cat("Treated units :", n, "(pooled into one synthetic control",
        if (is.null(x$treated_weights)) "as a simple average)"
        else "as a weighted average)", "\n")
  invisible(NULL)
}

#' @export
print.coresynth_staggered <- function(x, ...) {
  cat("=== coresynth fit ===\n")
  cat("Method :", toupper(x$method), "(staggered)\n")
  cat("Estimate (ATT):", round(x$estimate, 4), "\n")
  # T_pre is the panel-level split (the earliest adoption date), which holds
  # for one cohort only, so report the cohort count instead and leave the
  # per-cohort windows to summary() and tidy().
  cat("Cohorts:", nrow(x$cohort_estimates) %||% NA_integer_, "\n")
  invisible(x)
}

#' @export
print.coresynth_multiarm <- function(x, ...) {
  cat("=== coresynth fit ===\n")
  cat("Method :", toupper(x$method), "\n")
  stag_label <- if (inherits(x, "coresynth_staggered")) " (staggered)" else ""
  cat(sprintf("Multi-arm SI%s (K = %d arms)\n", stag_label, length(x$arm_levels)))
  cat("Per-arm ATT:",
      paste(names(x$arm_estimates), round(x$arm_estimates, 4),
            sep = "=", collapse = "  "), "\n")
  cat("Estimate (ATT):", round(x$estimate, 4), "\n")
  cat("Pre-treatment periods:", x$T_pre, "\n")
  invisible(x)
}

#' @export
summary.coresynth <- function(object, ...) {
  cat("=== coresynth summary ===\n")
  cat("Method :", toupper(object$method), "\n")
  .summary_coresynth_body(object)
  invisible(object)
}

#' @export
summary.coresynth_staggered <- function(object, ...) {
  cat("=== coresynth summary ===\n")
  cat("Method :", toupper(object$method), "(staggered)\n")
  .summary_cohort_table(object)
  cat("Aggregate ATT:", round(object$estimate, 6), "\n")
  invisible(object)
}

# Per-cohort estimates, weights, and treatment windows. Shared by
# summary.coresynth_staggered and the multi-arm summary, whose fits are
# staggered too whenever they carry a cohort table.
.summary_cohort_table <- function(object) {
  ce <- object$cohort_estimates %||% object$cohort_arm_estimates
  if (is.null(ce) || nrow(ce) == 0L) return(invisible(NULL))
  # Under pool_cohorts = FALSE a row is a treated unit, not a cohort, and
  # several rows share an adoption date.
  if (isFALSE(object$pool_cohorts))
    cat("Blocks  :", nrow(ce), "(one per treated unit) across",
        length(unique(ce$cohort)), "cohorts\n")
  else
    cat("Cohorts :", nrow(ce), "\n")
  ce$estimate <- round(ce$estimate, 6)
  ce$weight   <- round(ce$weight, 4)
  print(ce, row.names = FALSE)
  invisible(NULL)
}

#' @export
summary.coresynth_multiarm <- function(object, ...) {
  cat("=== coresynth summary ===\n")
  cat("Method :", toupper(object$method), "\n")
  stag_label <- if (inherits(object, "coresynth_staggered")) " (staggered)" else ""
  cat(sprintf("Multi-arm SI%s (K = %d treatment arms)\n",
              stag_label, length(object$arm_levels)))
  cat("Per-arm ATT:\n")
  print(round(object$arm_estimates, 6))
  if (inherits(object, "coresynth_staggered")) {
    .summary_cohort_table(object)
    cat("Aggregate ATT:", round(object$estimate, 6), "\n")
  } else {
    .summary_coresynth_body(object)
  }
  invisible(object)
}

# Shared tail of summary.coresynth / summary.coresynth_multiarm.
.summary_coresynth_body <- function(object) {
  .cat_pooled_treated(object)
  cat(
    "Periods : T_pre =",
    object$T_pre,
    "| T_post =",
    length(object$times) - object$T_pre,
    "\n"
  )
  cat("ATT estimate:", round(object$estimate, 6), "\n")
  if (!is.null(object$unit_weights)) {
    cat("Unit weights (non-zero donors):\n")
    w <- object$unit_weights
    print(round(w[w > 1e-4], 4))
  }
  if (!is.null(object$predictor_table)) {
    cat("\nPredictor balance:\n")
    pt <- object$predictor_table
    pt$treated <- round(pt$treated, 4)
    pt$synthetic <- round(pt$synthetic, 4)
    print(pt, row.names = FALSE)
  }
  invisible(object)
}
