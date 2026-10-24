# ── Internal: control-unit IFE fit and the treated-unit step ─────────────────
# Xu (2017) estimates the factors and beta from the control units alone
# (Step 1) and only then each treated unit's loadings from its own
# pre-treatment periods (Step 2). Splitting the two lets one control-unit fit
# serve every cohort of a staggered panel and every refit of the inference
# routines. `X` arguments are T x N x p covariate cubes, or NULL. `force`
# names the additive fixed effects of Xu (2017, Remark 3).
.gsc_fit_controls <- function(Y_co, X_co, r, force = "none") {
  TT <- nrow(Y_co)
  p  <- if (is.null(X_co)) 0L else dim(X_co)[3L]
  if (p == 0L) X_co <- array(0.0, c(TT, ncol(Y_co), 0L))
  gsc_ife_cpp(Y_co, matrix(0.0, TT, 0L), r, X_co, array(0.0, c(TT, 0L, p)),
              force = force)
}

# Whether `force` gives each unit its own effect, which the treated units
# then estimate from their pre-treatment periods alongside their loadings.
.gsc_unit_fe <- function(force) force %in% c("unit", "two-way")

# Counterfactual over all T periods for treated columns Y_tr whose first
# T_pre_g periods are untreated, covariate contribution included.
.gsc_counterfactual <- function(ctrl, Y_tr, X_tr, T_pre_g) {
  TT  <- nrow(Y_tr)
  n   <- ncol(Y_tr)
  p   <- length(ctrl$beta)
  pre <- seq_len(T_pre_g)
  X_pre <- if (p > 0L) X_tr[pre, , , drop = FALSE] else
    array(0.0, c(T_pre_g, n, 0L))
  force <- ctrl[["force"]] %||% "none"
  tr <- gsc_treated_cpp(ctrl$F, ctrl$beta, Y_tr[pre, , drop = FALSE], X_pre,
                        ctrl[["mu"]] %||% 0, ctrl[["xi"]] %||% numeric(0),
                        .gsc_unit_fe(force))
  Y_hat <- tr$Y_tr_hat
  for (j in seq_len(p))
    Y_hat <- Y_hat + matrix(X_tr[, , j], TT, n) * ctrl$beta[j]
  list(Y_hat = Y_hat, L_tr = tr$L_tr, alpha_tr = drop(tr$alpha_tr))
}

# ── Internal: cohort-by-cohort GSC for staggered adoption ────────────────────
# Clarke et al. (2024). Each cohort's treated loadings come from its own
# pre-treatment window; the factors come from the donors over all T periods.
# Under the default "never_treated" every cohort has the same donors, so the
# control-unit fit is done once. The deprecated "clean" adds the units
# adopting later, which differ by cohort, so it fits per cohort.
.fit_gsc_staggered <- function(pan, Y, r, covariate_array = NULL,
                                control_group = "never_treated",
                                force = "none") {
  TT      <- nrow(Y)
  T_adopt <- pan$T_adopt
  idx_tr  <- pan$idx_treat
  idx_co  <- pan$idx_control
  slice_X <- function(idx)
    if (is.null(covariate_array)) NULL else covariate_array[, idx, , drop = FALSE]

  shared <- NULL
  if (control_group == "never_treated") {
    shared <- tryCatch(
      .gsc_fit_controls(Y[, idx_co, drop = FALSE], slice_X(idx_co), r, force),
      error = function(e) stop("GSC staggered: the control-unit fit failed: ",
                               conditionMessage(e), call. = FALSE))
    .gsc_warn_unconverged(!isTRUE(shared$converged), "the control-unit fit")
  }

  cohorts <- sort(unique(T_adopt[idx_tr]))

  cohort_list <- lapply(cohorts, function(g) {
    idx_tr_g <- idx_tr[T_adopt[idx_tr] == g]
    T_pre_g  <- g - 1L
    T_post_g <- TT - T_pre_g

    if (control_group == "never_treated") {
      idx_co_g <- idx_co
    } else {
      future_tr <- idx_tr[!is.na(T_adopt[idx_tr]) & T_adopt[idx_tr] > g]
      idx_co_g  <- c(idx_co, future_tr)
    }

    N_co_g <- length(idx_co_g)
    if (N_co_g < r || T_pre_g < r + .gsc_unit_fe(force)) {
      warning(sprintf(
        "GSC staggered: cohort g=%d skipped (T_pre=%d, N_co=%d, r=%d).",
        g, T_pre_g, N_co_g, r), call. = FALSE)
      return(NULL)
    }

    ctrl <- shared
    if (is.null(ctrl)) {
      ctrl <- tryCatch(
        .gsc_fit_controls(Y[, idx_co_g, drop = FALSE], slice_X(idx_co_g), r,
                          force),
        error = function(e) {
          warning(sprintf("GSC staggered: cohort g=%d failed: %s",
                          g, conditionMessage(e)), call. = FALSE)
          NULL
        })
      if (is.null(ctrl)) return(NULL)
      .gsc_warn_unconverged(!isTRUE(ctrl$converged), sprintf("cohort g=%d", g))
    }

    Y_treat_g <- Y[, idx_tr_g, drop = FALSE]
    cf        <- .gsc_counterfactual(ctrl, Y_treat_g, slice_X(idx_tr_g), T_pre_g)
    post_rows <- (T_pre_g + 1L):TT
    att_g     <- mean(Y_treat_g[post_rows, , drop = FALSE] -
                        cf$Y_hat[post_rows, , drop = FALSE])

    list(
      cohort    = g,
      n_treated = length(idx_tr_g),
      T_pre     = T_pre_g,
      T_post    = T_post_g,
      estimate  = att_g,
      weight    = as.numeric(length(idx_tr_g)) * T_post_g,
      idx_tr    = idx_tr_g,
      idx_co    = idx_co_g,
      Y_tr_hat  = cf$Y_hat,
      Y_treat   = Y_treat_g,
      F         = ctrl$F,
      L_co      = ctrl$L_co,
      L_tr      = cf$L_tr,
      alpha_tr  = cf$alpha_tr,
      beta      = ctrl$beta
    )
  })

  valid <- !vapply(cohort_list, is.null, logical(1L))
  if (!any(valid)) {
    stop("All cohort-level GSC fits failed; see the preceding warnings for ",
         "per-cohort reasons.", call. = FALSE)
  }
  .aggregate_cohorts(cohort_list[valid])
}

# gsc_ife_cpp() caps the covariate loop at max_iter. A fit that stops at the
# cap carries a beta that has not settled, which moves the counterfactual, so
# say so rather than return it silently. `n` counts the fits that hit the cap.
.gsc_warn_unconverged <- function(n, where) {
  if (n > 0L)
    warning(sprintf(paste0(
      "GSC: the covariate loop did not converge (%s); beta and the effect ",
      "estimate may not have settled."), where), call. = FALSE)
  invisible(NULL)
}

# ── Internal: Xu (2017) Algorithm 1 -- cross-validate the number of factors ──
# Leave-one-period-out over the pre-treatment window. For each candidate r the
# IFE model is estimated once on the control units; then each pre-treatment
# period s is held out in turn, the treated units' factor loadings are refit by
# OLS on the remaining pre-treatment periods, and period s is predicted from
# them. The r with the smallest mean squared prediction error wins.
#
# With additive fixed effects the candidates start at r = 0, the effects
# alone; the treated units' unit effects are refit alongside their loadings.
# Without them r = 0 would leave nothing but the covariates, so the
# candidates start at 1.
#
# Dividing by T_pre alone, rather than by T_pre * N_treated, follows the
# paper's Step 3. The treated count is constant across candidates, so it does
# not affect which r wins.
.gsc_cv_r <- function(Y_co, Y_tr_pre, X_co, X_tr_pre, r_max = 5L,
                      force = "none") {
  T_pre   <- nrow(Y_tr_pre)
  N_co    <- ncol(Y_co)
  N_tr    <- ncol(Y_tr_pre)
  p       <- dim(X_tr_pre)[3L]
  unit_fe <- .gsc_unit_fe(force)

  # Holding out one period leaves T_pre - 1 to identify the loadings, and the
  # unit effect when there is one, by OLS
  lo <- if (identical(force, "none")) 1L else 0L
  hi <- min(as.integer(r_max), T_pre - 1L - unit_fe, N_co)
  if (hi < lo)
    stop("GSC cannot cross-validate the number of factors with T_pre = ",
         T_pre, " pre-treatment period(s) and ", N_co, " control unit(s). ",
         "Supply r explicitly.", call. = FALSE)

  cands <- lo:hi
  n_unconv <- 0L
  mspe  <- vapply(cands, function(r) {
    res <- tryCatch(gsc_ife_cpp(Y_co, Y_tr_pre, r, X_co, X_tr_pre,
                                force = force),
                    error = function(e) NULL)
    if (is.null(res)) return(NA_real_)
    if (!isTRUE(res$converged)) n_unconv <<- n_unconv + 1L
    Z <- res$F[seq_len(T_pre), , drop = FALSE]           # T_pre x r
    if (unit_fe) Z <- cbind(1, Z)

    # Take the covariate part and the grand mean and time effects out of the
    # treated pre-treatment outcomes, so what the held-out loadings have to
    # explain is the unit-specific part.
    R_tr <- Y_tr_pre - res$mu - drop(res$xi)[seq_len(T_pre)]
    if (p > 0L)
      for (j in seq_len(p))
        R_tr <- R_tr - matrix(X_tr_pre[, , j], T_pre, N_tr) * res$beta[j]

    err <- vapply(seq_len(T_pre), function(s) {
      if (ncol(Z) == 0L) return(sum(R_tr[s, ]^2))
      Z_ms <- Z[-s, , drop = FALSE]
      lam  <- tryCatch(
        solve(crossprod(Z_ms), crossprod(Z_ms, R_tr[-s, , drop = FALSE])),
        error = function(e) NULL)                         # ncol(Z) x N_tr
      if (is.null(lam)) return(NA_real_)
      sum((R_tr[s, ] - drop(Z[s, , drop = FALSE] %*% lam))^2)
    }, numeric(1L))

    if (anyNA(err)) NA_real_ else sum(err) / T_pre
  }, numeric(1L))

  if (all(is.na(mspe)))
    stop("GSC cross-validation failed for every candidate number of factors. ",
         "Supply r explicitly.", call. = FALSE)
  .gsc_warn_unconverged(n_unconv, sprintf(
    "%d of the %d candidate fits of the factor-count cross-validation",
    n_unconv, length(cands)))

  list(r     = cands[which.min(mspe)],
       table = data.frame(r = cands, mspe = mspe))
}

#' Fit Generalized Synthetic Control (Xu, 2017)
#'
#' Uses the Interactive Fixed Effects (IFE) model:
#'   Y_it(0) = x_it'beta + lambda_i' f_t + (mu + alpha_i + xi_t) + eps_it
#' Factors f_t and control loadings Lambda_co are extracted from control units
#' by the iteration of Bai (2009), alternating a truncated SVD with the
#' regression for beta. The additive effects in parentheses are those of the
#' extended model in Xu (2017, Remark 3), chosen by `force`.
#'
#' @param y    Outcome vector (long format)
#' @param d    Treatment indicator (long format)
#' @param id   Unit identifier (long format)
#' @param time Time identifier (long format)
#' @param r    Number of latent factors. `NULL` (default) selects it by the
#'   leave-one-period-out cross-validation of Xu (2017, Algorithm 1).
#' @param data Long-format data frame. Required when `covariates` is non-NULL.
#' @param id_var   Name of the unit identifier column in `data`.
#' @param time_var Name of the time identifier column in `data`.
#' @param covariates Character vector of time-varying covariate column names in
#'   `data`. When non-NULL, `data`, `id_var`, and `time_var` must be supplied.
#' @param force Additive fixed effects: `"none"` (default), `"unit"`, `"time"`
#'   or `"two-way"`.
#' @return A list of class `coresynth`.
#' @noRd
fit_gsc_cpp <- function(
  y,
  d,
  id,
  time,
  r = NULL,
  data = NULL,
  id_var = NULL,
  time_var = NULL,
  covariates = NULL,
  control_group = c("never_treated", "clean"),
  force = c("none", "unit", "time", "two-way"),
  ...
) {
  control_group <- match.arg(control_group)
  force         <- match.arg(force)
  pan <- panel_to_matrices(y, d, id, time)
  .check_panel_complete(pan$Y, "GSC")
  .check_control_group(control_group, pan, "GSC")

  r_min <- if (force == "none") 1L else 0L
  if (!is.null(r) && (!is.numeric(r) || length(r) != 1L || r < r_min ||
                      r != as.integer(r)))
    stop(if (r_min == 1L)
           "r must be NULL (cross-validated) or a single positive integer."
         else
           "r must be NULL (cross-validated) or a single non-negative integer.",
         call. = FALSE)
  unit_fe <- .gsc_unit_fe(force)

  # Xu (2017) Algorithm 1 picks the number of factors before the effect is
  # estimated. Every treated unit is untreated throughout the global
  # pre-treatment window, so one cross-validation serves sharp and staggered
  # panels alike.
  r_cv <- NULL
  if (is.null(r)) {
    T_pre_cv <- pan$T_pre
    Y_co_cv  <- pan$Y[, pan$idx_control, drop = FALSE]
    Y_tr_cv  <- pan$Y[seq_len(T_pre_cv), pan$idx_treat, drop = FALSE]
    if (!is.null(covariates) && !is.null(data)) {
      if (is.null(id_var) || is.null(time_var))
        stop("'id_var' and 'time_var' must be provided with covariates.",
             call. = FALSE)
      X_co_cv <- build_covariate_array(data, id_var, time_var, covariates,
                                       pan$units[pan$idx_control], pan$times)
      X_tr_cv <- build_covariate_array(
        data, id_var, time_var, covariates,
        pan$units[pan$idx_treat], pan$times
      )[seq_len(T_pre_cv), , , drop = FALSE]
    } else {
      X_co_cv <- array(0.0, c(nrow(pan$Y), length(pan$idx_control), 0L))
      X_tr_cv <- array(0.0, c(T_pre_cv,    length(pan$idx_treat),   0L))
    }
    cv   <- .gsc_cv_r(Y_co_cv, Y_tr_cv, X_co_cv, X_tr_cv, force = force)
    r    <- cv$r
    r_cv <- cv$table
  }

  # ── Staggered adoption path ──────────────────────────────────────────────
  if (!pan$is_sharp) {
    use_cov_st <- !is.null(covariates) && !is.null(data)
    if (use_cov_st) {
      if (is.null(id_var) || is.null(time_var))
        stop("'id_var' and 'time_var' must be provided with covariates.",
             call. = FALSE)
      X_all_arr <- build_covariate_array(data, id_var, time_var, covariates,
                                          pan$units, pan$times)
    } else {
      X_all_arr <- NULL
    }
    res_st <- .fit_gsc_staggered(pan, pan$Y, r,
                                  covariate_array = X_all_arr,
                                  control_group   = control_group,
                                  force           = force)
    res_st$method              <- "gsc"
    res_st$staggered           <- TRUE
    res_st$r                   <- r
    res_st$r_cv                <- r_cv
    res_st$force               <- force
    res_st$Y_treat             <- pan$Y[, pan$idx_treat, drop = FALSE]
    res_st$Y_synth             <- NULL
    res_st$gap                 <- NULL
    res_st$times               <- pan$times
    res_st$T_pre               <- pan$T_pre
    res_st$unit_weights        <- NULL
    res_st$covariates          <- covariates
    res_st$Y_all               <- pan$Y            # T × N — needed by gsc_inference()
    res_st$idx_co              <- pan$idx_control
    res_st$idx_tr              <- pan$idx_treat
    res_st$covariate_array_all <- if (use_cov_st) X_all_arr else NULL
    res_st$control_group       <- control_group
    return(res_st)
  }

  Y <- pan$Y
  T_pre <- pan$T_pre
  idx_tr <- pan$idx_treat
  idx_co <- pan$idx_control

  if (length(idx_co) < r) {
    stop("Need at least r control units for GSC.")
  }
  if (T_pre < r + unit_fe) {
    stop("Need at least r pre-treatment periods for GSC, plus one for the ",
         "unit effect when force includes it.")
  }

  Y_co_all <- Y[, idx_co, drop = FALSE] # T x N_co
  Y_tr_pre <- Y[seq_len(T_pre), idx_tr, drop = FALSE] # T_pre x N_tr

  # ── Build covariate arrays ────────────────────────────────────────────────
  use_cov <- !is.null(covariates) && !is.null(data)
  if (use_cov) {
    if (is.null(id_var) || is.null(time_var)) {
      stop(
        "'id_var' and 'time_var' must be provided with covariates.",
        call. = FALSE
      )
    }
    co_units <- pan$units[idx_co]
    tr_units <- pan$units[idx_tr]
    X_co_arr <- build_covariate_array(
      data,
      id_var,
      time_var,
      covariates,
      co_units,
      pan$times
    ) # T x N_co x p
    X_tr_arr <- build_covariate_array(
      data,
      id_var,
      time_var,
      covariates,
      tr_units,
      pan$times
    ) # T x N_tr x p
    # Pre-treatment slice of X_tr for Xu (2017) Step 2 loading estimation
    X_tr_pre_arr <- X_tr_arr[seq_len(T_pre), , , drop = FALSE] # T_pre x N_tr x p
    p <- length(covariates)
  } else {
    T_full <- nrow(Y_co_all)
    N_co_n <- ncol(Y_co_all)
    X_co_arr <- array(0.0, dim = c(T_full, N_co_n, 0L)) # empty cube (p=0)
    X_tr_arr <- NULL
    X_tr_pre_arr <- NULL
    p <- 0L
  }

  # Build empty X_tr_pre cube for the no-covariate case
  if (is.null(X_tr_pre_arr)) {
    N_tr_n <- length(idx_tr)
    X_tr_pre_arr <- array(0.0, dim = c(T_pre, N_tr_n, 0L))
  }

  res <- gsc_ife_cpp(Y_co_all, Y_tr_pre, r, X_co_arr, X_tr_pre_arr,
                     force = force)
  .gsc_warn_unconverged(!isTRUE(res$converged), "the fit")
  beta <- res$beta # numeric(0) when p=0

  # ── Add beta contribution to Y_tr_hat and Y_hat_co ───────────────────────
  Y_tr_hat <- res$Y_tr_hat # T x N_tr, factors and additive effects
  Y_hat_co <- res$Y_co_hat # T x N_co, on the same basis

  if (use_cov && p > 0L) {
    T_all <- nrow(Y_tr_hat)
    N_tr_n <- ncol(Y_tr_hat)
    for (j in seq_len(p)) {
      cov_tr <- matrix(X_tr_arr[,, j], T_all, N_tr_n)
      Y_tr_hat <- Y_tr_hat + cov_tr * beta[j]
      Y_hat_co <- Y_hat_co + X_co_arr[,, j] * beta[j]
    }
  }

  Y_treat <- Y[, idx_tr, drop = FALSE] # T x N_tr

  # ATT: average post-treatment residual across treated units and post periods
  post_rows <- (T_pre + 1L):nrow(Y)
  gap_post <- Y_treat[post_rows, , drop = FALSE] -
    Y_tr_hat[post_rows, , drop = FALSE]
  att <- mean(gap_post)

  list(
    method = "gsc",
    r = r,
    r_cv = r_cv,
    force = force,
    F = res$F,
    L_co = res$L_co,
    L_tr = res$L_tr,
    mu = res$mu,
    alpha_co = drop(res$alpha_co),
    alpha_tr = drop(res$alpha_tr),
    xi = drop(res$xi),
    singular_values = res$singular_values,
    beta = beta, # numeric(0) when no covariates
    Y_tr_hat = Y_tr_hat, # T x N_tr (beta included if covariates)
    Y_treat = Y_treat,
    gap = Y_treat - Y_tr_hat,
    times = pan$times,
    T_pre = T_pre,
    unit_weights = NULL, # GSC uses loadings, not simplex weights
    estimate = att,
    Y_co_all = Y_co_all, # T x N_co — needed by gsc_boot()
    Y_hat_co = Y_hat_co, # T x N_co, full prediction (beta included)
    covariate_array_co = if (use_cov) X_co_arr else NULL,
    covariate_array_tr = if (use_cov) X_tr_arr else NULL,
    covariates = covariates # variable names or NULL
  )
}

#' Parametric Bootstrap Inference for GSC (Xu 2017, Algorithm 2)
#'
#' Uncertainty for the ATT of a sharp GSC fit, conditional on the covariates,
#' factors and loadings, by resampling the residuals of the estimated IFE
#' model (Xu 2017, S.3.2 and Algorithm 2). It is the procedure the paper
#' proposes when few units are treated, where resampling the treated units
#' themselves cannot reach their noise.
#'
#' Each replicate rebuilds the panel without a treatment effect from the
#' fitted values and re-estimates it. Residuals are resampled as whole time
#' series, which keeps their serial correlation. The control units take the
#' in-sample residuals of the IFE fit. The treated units take prediction
#' errors instead, which are larger: before the bootstrap, a control unit is
#' set aside as if it were treated, the model is refitted on a resample of
#' the remaining controls, and the left-out unit's observed series minus its
#' predicted counterfactual is one prediction error. Covariates enter through
#' the fitted values and are re-estimated in each replicate.
#'
#' @param fit   A `coresynth` object from [scm_fit()] with `method = "gsc"`.
#' @param B     Number of bootstrap replications, and of prediction errors
#'   drawn for the treated units (default 499L).
#' @param alpha Significance level for the confidence interval (default 0.05).
#' @param seed  RNG seed for reproducibility (default NULL).
#' @return A list with:
#'   * `p_value`:   Two-sided p-value: mean(|ATT*| >= |ATT_obs|)
#'   * `ci_lower`:  Lower bound of the (1-alpha)*100% CI for the ATT,
#'     `ATT_obs - q(1 - alpha/2)` where `q` are quantiles of `boot_dist`
#'   * `ci_upper`:  Upper bound, `ATT_obs - q(alpha/2)`
#'   * `se`:        Bootstrap standard error
#'   * `boot_dist`: Numeric vector of length B (bootstrap ATT* values)
#'   * `att_obs`:   Observed ATT from the original fit
#' @export
gsc_boot <- function(fit, B = 499L, alpha = 0.05, seed = NULL) {
  if (!inherits(fit, "coresynth_gsc")) {
    stop("gsc_boot() requires a coresynth object with method = 'gsc'.")
  }
  if (inherits(fit, "coresynth_staggered")) {
    stop("gsc_boot() does not support staggered adoption fits. ",
         "Use gsc_inference() instead.", call. = FALSE)
  }
  if (is.null(fit$Y_co_all)) {
    stop("fit is missing Y_co_all. Re-run scm_fit() with the current version.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  r       <- fit$r
  force   <- fit[["force"]] %||% "none"
  T_pre   <- fit$T_pre
  att_obs <- fit$estimate
  Y_co    <- fit$Y_co_all
  X_co    <- fit$covariate_array_co      # T x N_co x p, or NULL
  X_tr    <- fit$covariate_array_tr      # T x N_tr x p, or NULL
  TT      <- nrow(Y_co)
  N_co    <- ncol(Y_co)
  N_tr    <- ncol(fit$Y_treat)
  post    <- (T_pre + 1L):TT
  slice_X <- function(X, idx) if (is.null(X)) NULL else X[, idx, , drop = FALSE]

  # Fitted values without a treatment effect, covariate part included
  Y_hat_co <- fit$Y_hat_co %||% (fit$F %*% t(fit$L_co))
  Y_hat_tr <- fit$Y_tr_hat
  res_co   <- Y_co - Y_hat_co

  n_unconverged <- 0L
  fit_ctrl <- function(Y, X) {
    ctrl <- tryCatch(.gsc_fit_controls(Y, X, r, force), error = function(e) NULL)
    if (!is.null(ctrl) && !isTRUE(ctrl$converged))
      n_unconverged <<- n_unconverged + 1L
    ctrl
  }

  # Step 1: prediction errors, each from a control unit left out as if treated
  pred_err <- matrix(NA_real_, TT, B)
  for (m in seq_len(B)) {
    i    <- sample.int(N_co, 1L)
    rest <- seq_len(N_co)[-i]
    co_m <- rest[sample.int(length(rest), N_co, replace = TRUE)]
    ctrl <- fit_ctrl(Y_co[, co_m, drop = FALSE], slice_X(X_co, co_m))
    if (is.null(ctrl)) next
    cf <- tryCatch(
      .gsc_counterfactual(ctrl, Y_co[, i, drop = FALSE], slice_X(X_co, i), T_pre),
      error = function(e) NULL)
    if (!is.null(cf)) pred_err[, m] <- Y_co[, i] - cf$Y_hat[, 1L]
  }
  pred_err <- pred_err[, colSums(is.na(pred_err)) == 0L, drop = FALSE]
  if (ncol(pred_err) == 0L)
    stop("gsc_boot(): every leave-one-control-out fit failed, so there are ",
         "no prediction errors to resample. Consider reducing r.",
         call. = FALSE)

  # Steps 2-3: resample residual series onto the fitted values and re-estimate
  boot_dist <- numeric(B)
  for (b in seq_len(B)) {
    Y_star_co <- Y_hat_co + res_co[, sample.int(N_co, N_co, replace = TRUE),
                                   drop = FALSE]
    Y_star_tr <- Y_hat_tr + pred_err[, sample.int(ncol(pred_err), N_tr,
                                                  replace = TRUE), drop = FALSE]
    ctrl <- fit_ctrl(Y_star_co, X_co)
    cf   <- if (is.null(ctrl)) NULL else tryCatch(
      .gsc_counterfactual(ctrl, Y_star_tr, X_tr, T_pre),
      error = function(e) NULL)
    boot_dist[b] <- if (is.null(cf)) NA_real_ else
      mean(Y_star_tr[post, , drop = FALSE] - cf$Y_hat[post, , drop = FALSE])
  }

  n_fail <- sum(is.na(boot_dist))
  if (n_fail > 0.1 * B) {
    warning(
      sprintf(
        "%.0f%% of bootstrap replications failed (NA). ",
        100 * n_fail / B
      ),
      "Consider reducing r or checking model fit."
    )
  }
  .gsc_warn_unconverged(n_unconverged, sprintf(
    "%d of the bootstrap refits", n_unconverged))

  valid <- boot_dist[!is.na(boot_dist)]
  # Each replicate is generated with no effect, so ATT* is a draw of the
  # estimation error ATT_hat - ATT (Xu 2017, Algorithm 2). The interval for
  # the ATT therefore subtracts the error quantiles from the estimate; the raw
  # quantiles of ATT* are an acceptance region around zero, not a CI.
  err_q <- quantile(valid, c(alpha / 2, 1 - alpha / 2), names = FALSE)

  list(
    p_value = mean(abs(valid) >= abs(att_obs)),
    ci_lower = att_obs - err_q[2],
    ci_upper = att_obs - err_q[1],
    se = sd(valid),
    boot_dist = boot_dist,
    att_obs = att_obs
  )
}

# ── Internal: GSC effect for a draw of units ─────────────────────────────────
# The estimator of fit_gsc_cpp() applied to never-treated columns `co` and
# treated columns `tr` of ctx$Y; see .unit_inference(). Cohorts are rebuilt
# from the treated units drawn and combined with the Clarke et al. (2024)
# weights, and the treated units of each cohort get their loadings from that
# cohort's pre-treatment window.
.gsc_unit_estimate <- function(ctx, co, tr) {
  TT      <- nrow(ctx$Y)
  adopt   <- ctx$T_adopt[tr]
  cohorts <- sort(unique(adopt))
  slice_X <- function(idx)
    if (is.null(ctx$X)) NULL else ctx$X[, idx, , drop = FALSE]
  unconverged <- 0L
  fit_ctrl <- function(idx) {
    if (length(idx) < ctx$r) return(NULL)
    ctrl <- tryCatch(
      .gsc_fit_controls(ctx$Y[, idx, drop = FALSE], slice_X(idx), ctx$r,
                        ctx$force),
      error = function(e) NULL)
    if (!is.null(ctrl) && !isTRUE(ctrl$converged))
      unconverged <<- unconverged + 1L
    ctrl
  }

  shared <- NULL
  if (identical(ctx$control_group, "never_treated")) {
    shared <- fit_ctrl(co)
    if (is.null(shared)) return(NULL)
  }
  est <- w <- numeric(length(cohorts))
  for (k in seq_along(cohorts)) {
    g    <- cohorts[k]
    tr_g <- tr[adopt == g]
    ctrl <- shared %||% fit_ctrl(c(co, tr[adopt > g]))
    if (is.null(ctrl)) return(NULL)
    cf <- tryCatch(
      .gsc_counterfactual(ctrl, ctx$Y[, tr_g, drop = FALSE], slice_X(tr_g),
                          g - 1L),
      error = function(e) NULL)
    if (is.null(cf)) return(NULL)
    post   <- g:TT
    est[k] <- mean(ctx$Y[post, tr_g, drop = FALSE] -
                     cf$Y_hat[post, , drop = FALSE])
    w[k]   <- length(tr_g) * length(post)
  }
  list(estimate = sum(w * est) / sum(w), cohort = cohorts,
       cohort_estimate = est, cohort_weight = w, unconverged = unconverged)
}

# The panel a GSC fit was estimated on, laid out for .gsc_unit_estimate().
.gsc_inference_context <- function(fit) {
  if (inherits(fit, "coresynth_staggered")) {
    if (is.null(fit$Y_all))
      stop("fit does not contain Y_all. Re-estimate with the current version.",
           call. = FALSE)
    tr <- unlist(lapply(fit$cohort_fits, `[[`, "idx_tr"), use.names = FALSE)
    return(list(
      Y = fit$Y_all, X = fit$covariate_array_all, r = fit$r,
      T_adopt = .adopt_from_cohorts(ncol(fit$Y_all), fit$cohort_fits),
      co = fit$idx_co, tr = tr, control_group = .fit_control_group(fit),
      force = fit[["force"]] %||% "none"))
  }
  if (is.null(fit$Y_co_all))
    stop("fit does not contain Y_co_all. Re-estimate with the current version.",
         call. = FALSE)
  N_co <- ncol(fit$Y_co_all)
  N_tr <- ncol(fit$Y_treat)
  X <- NULL
  if (!is.null(fit$covariate_array_co)) {
    p <- dim(fit$covariate_array_co)[3L]
    X <- array(0.0, c(nrow(fit$Y_co_all), N_co + N_tr, p))
    X[, seq_len(N_co), ]        <- fit$covariate_array_co
    X[, N_co + seq_len(N_tr), ] <- fit$covariate_array_tr
  }
  list(
    Y = cbind(fit$Y_co_all, fit$Y_treat), X = X, r = fit$r,
    T_adopt = c(rep(NA_integer_, N_co), rep(fit$T_pre + 1L, N_tr)),
    co = seq_len(N_co), tr = N_co + seq_len(N_tr),
    control_group = "never_treated", force = fit[["force"]] %||% "none")
}

#' Non-parametric Inference for GSC (Xu 2017)
#'
#' Estimates SE and confidence intervals for the ATT by re-estimating the
#' fit on resampled units: a bootstrap of the treated and the control units,
#' drawn separately, or a jackknife that drops one unit at a time. Works for
#' both sharp and staggered GSC fits.
#'
#' Xu (2017, S.3.2) notes that resampling units is valid when many units are
#' treated and proposes the parametric bootstrap of [gsc_boot()] for when few
#' are: the treated units' own noise varies across replicates only through
#' which treated units are drawn. With a single treated unit there is no such
#' variation, so this function stops and points to [gsc_boot()] and
#' [conformal_inference()].
#'
#' On a staggered fit the bootstrap draws the treated units from all cohorts
#' together, so cohort sizes vary across replicates, and draws the controls
#' once per replicate for every cohort; cohorts are rebuilt from the draw and
#' combined with their `N_treated * T_post` weights. `"jackknife"` drops each
#' cohort's own donors and treated units and combines the per-cohort
#' variances by the delta method, which ignores the covariance that shared
#' donors induce between cohorts; `"jackknife_global"` drops every unit from
#' the whole panel in turn and keeps it. Both jackknives drop treated units,
#' so every cohort needs at least two (Clarke et al. 2024, Algorithm 3).
#'
#' Note: `gsc_boot()` performs a *parametric* bootstrap under H0 (hypothesis
#' testing). `gsc_inference()` provides *non-parametric* SE and CIs suitable
#' for inference about the ATT magnitude.
#'
#' @param fit   A `coresynth` object from [scm_fit()] with `method = "gsc"`.
#' @param method `"bootstrap"` (default), `"jackknife"`, or
#'   `"jackknife_global"` (staggered fits only).
#' @param n_boot Number of bootstrap replications (default 499L; ignored for jackknife).
#' @param level  Confidence level (default 0.95).
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param seed  RNG seed for reproducibility (default NULL).
#' @return A list of class `coresynth_inference`.
#' @export
gsc_inference <- function(
  fit,
  method      = c("bootstrap", "jackknife", "jackknife_global"),
  n_boot      = 499L,
  level       = 0.95,
  alternative = c("two.sided", "greater", "less"),
  seed        = NULL
) {
  method      <- match.arg(method)
  alternative <- match.arg(alternative)

  if (!inherits(fit, "coresynth_gsc"))
    stop("gsc_inference() requires a coresynth fit with method = 'gsc'.",
         call. = FALSE)
  if (method == "jackknife_global" && !inherits(fit, "coresynth_staggered"))
    stop("gsc_inference() method='jackknife_global' requires a staggered fit.",
         call. = FALSE)

  ctx <- .gsc_inference_context(fit)
  if (!is.null(seed)) set.seed(seed)
  .unit_inference(
    fit, ctx, .gsc_unit_estimate, method, n_boot, level, alternative,
    label = "gsc_inference()",
    one_treated_hint = paste0(
      "Use gsc_boot(), the parametric bootstrap of Xu (2017, Algorithm 2) ",
      "for few treated units, or conformal_inference()."))
}
