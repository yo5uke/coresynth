# ── Internal: cohort-by-cohort GSC for staggered adoption ────────────────────
# Clarke et al. (2023); Arkhangelsky et al. (2021) Appendix S.8.
# Factors estimated independently per cohort (gsc_ife_cpp per cohort).
# Y_co_g uses all T periods (pragmatic: future-treated contamination treated as
# noise, same rationale as SDID staggered). Covariates handled via per-cohort
# EM (subset of global T×N×p array).
.fit_gsc_staggered <- function(pan, Y, r, covariate_array = NULL,
                                control_group = "clean") {
  TT      <- nrow(Y)
  T_adopt <- pan$T_adopt
  idx_tr  <- pan$idx_treat
  idx_co  <- pan$idx_control
  p       <- if (!is.null(covariate_array)) dim(covariate_array)[3L] else 0L

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
    if (N_co_g < r || T_pre_g < r) {
      warning(sprintf(
        "GSC staggered: cohort g=%d skipped (T_pre=%d, N_co=%d, r=%d).",
        g, T_pre_g, N_co_g, r), call. = FALSE)
      return(NULL)
    }

    pre_rows  <- seq_len(T_pre_g)
    post_rows <- (T_pre_g + 1L):TT

    Y_co_g     <- Y[, idx_co_g, drop = FALSE]
    Y_tr_pre_g <- Y[pre_rows, idx_tr_g, drop = FALSE]
    Y_treat_g  <- Y[, idx_tr_g, drop = FALSE]

    if (p > 0L) {
      X_co_g     <- covariate_array[, idx_co_g, , drop = FALSE]
      X_tr_pre_g <- covariate_array[pre_rows, idx_tr_g, , drop = FALSE]
    } else {
      X_co_g     <- array(0.0, dim = c(TT, N_co_g, 0L))
      X_tr_pre_g <- array(0.0, dim = c(T_pre_g, length(idx_tr_g), 0L))
    }

    res <- tryCatch(
      gsc_ife_cpp(Y_co_g, Y_tr_pre_g, r, X_co_g, X_tr_pre_g),
      error = function(e) {
        warning(sprintf("GSC staggered: cohort g=%d failed: %s",
                        g, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(res)) return(NULL)

    Y_tr_hat_g <- res$Y_tr_hat
    if (p > 0L) {
      X_tr_full_g <- covariate_array[, idx_tr_g, , drop = FALSE]
      for (j in seq_len(p)) {
        Y_tr_hat_g <- Y_tr_hat_g +
          matrix(X_tr_full_g[,, j], TT, length(idx_tr_g)) * res$beta[j]
      }
    }

    att_g <- mean(Y_treat_g[post_rows, , drop = FALSE] -
                    Y_tr_hat_g[post_rows, , drop = FALSE])

    list(
      cohort    = g,
      n_treated = length(idx_tr_g),
      T_pre     = T_pre_g,
      T_post    = T_post_g,
      estimate  = att_g,
      weight    = as.numeric(length(idx_tr_g)) * T_post_g,
      idx_tr    = idx_tr_g,
      idx_co    = idx_co_g,
      Y_tr_hat  = Y_tr_hat_g,
      Y_treat   = Y_treat_g,
      F         = res$F,
      L_co      = res$L_co,
      L_tr      = res$L_tr,
      beta      = res$beta
    )
  })

  valid <- !vapply(cohort_list, is.null, logical(1L))
  if (!any(valid)) stop("All cohort-level GSC fits failed.", call. = FALSE)
  r_list <- cohort_list[valid]

  w   <- vapply(r_list, `[[`, numeric(1L), "weight")
  tau <- vapply(r_list, `[[`, numeric(1L), "estimate")
  att <- sum(w * tau) / sum(w)

  cohort_df <- data.frame(
    cohort    = vapply(r_list, `[[`, integer(1L), "cohort"),
    n_treated = vapply(r_list, `[[`, integer(1L), "n_treated"),
    T_pre     = vapply(r_list, `[[`, integer(1L), "T_pre"),
    T_post    = vapply(r_list, `[[`, integer(1L), "T_post"),
    estimate  = tau,
    weight    = w / sum(w),
    stringsAsFactors = FALSE
  )

  list(estimate = att, cohort_estimates = cohort_df, cohort_fits = r_list)
}

#' Fit Generalized Synthetic Control (Xu, 2017)
#'
#' Uses the Interactive Fixed Effects (IFE) model:
#'   Y_it(0) = x_it'beta + lambda_i' f_t + eps_it
#' Factors f_t and control loadings Lambda_co are extracted from control units
#' via an EM loop that alternates truncated SVD (E-step) with panel OLS for
#' beta (M-step). When no covariates are supplied (covariates = NULL), the
#' plain 3-step SVD estimator (beta = 0) is used.
#'
#' @param y    Outcome vector (long format)
#' @param d    Treatment indicator (long format)
#' @param id   Unit identifier (long format)
#' @param time Time identifier (long format)
#' @param r    Number of latent factors (default 2)
#' @param data Long-format data frame. Required when `covariates` is non-NULL.
#' @param id_var   Name of the unit identifier column in `data`.
#' @param time_var Name of the time identifier column in `data`.
#' @param covariates Character vector of time-varying covariate column names in
#'   `data`. When non-NULL, `data`, `id_var`, and `time_var` must be supplied.
#' @return A list of class `coresynth`.
#' @noRd
fit_gsc_cpp <- function(
  y,
  d,
  id,
  time,
  r = 2,
  data = NULL,
  id_var = NULL,
  time_var = NULL,
  covariates = NULL,
  control_group = c("clean", "never_treated"),
  ...
) {
  control_group <- match.arg(control_group)
  pan <- panel_to_matrices(y, d, id, time)

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
                                  control_group   = control_group)
    res_st$method              <- "gsc"
    res_st$staggered           <- TRUE
    res_st$r                   <- r
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
    return(res_st)
  }

  Y <- pan$Y
  T_pre <- pan$T_pre
  idx_tr <- pan$idx_treat
  idx_co <- pan$idx_control

  if (length(idx_co) < r) {
    stop("Need at least r control units for GSC.")
  }
  if (T_pre < r) {
    stop("Need at least r pre-treatment periods for GSC.")
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

  res <- gsc_ife_cpp(Y_co_all, Y_tr_pre, r, X_co_arr, X_tr_pre_arr)
  beta <- res$beta # numeric(0) when p=0

  # ── Add beta contribution to Y_tr_hat and compute Y_hat_co ───────────────
  Y_tr_hat <- res$Y_tr_hat # T x N_tr, factor part
  Y_hat_co <- res$F %*% t(res$L_co) # T x N_co, factor part

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
    F = res$F,
    L_co = res$L_co,
    L_tr = res$L_tr,
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

#' Parametric Bootstrap Inference for GSC (Xu 2017 S.3)
#'
#' Generates the null distribution of the ATT under H0 (no treatment effect)
#' by parametric resampling from the estimated IFE factor model. Under H0,
#' both the control panel and treated unit are generated from the fitted
#' factor model with homoskedastic noise. When the fit includes covariate
#' adjustment (beta), the covariate contribution is included in the simulated
#' DGP and re-estimated in each bootstrap replicate.
#'
#' @param fit   A `coresynth` object from [scm_fit()] with `method = "gsc"`.
#' @param B     Bootstrap replications (default 499L).
#' @param alpha Significance level for the confidence interval (default 0.05).
#' @param seed  RNG seed for reproducibility (default NULL).
#' @return A list with:
#'   * `p_value`:   Two-sided p-value: mean(|ATT*| >= |ATT_obs|)
#'   * `ci_lower`:  Lower bound of (1-alpha)*100% bootstrap CI
#'   * `ci_upper`:  Upper bound of (1-alpha)*100% bootstrap CI
#'   * `se`:        Bootstrap standard error
#'   * `boot_dist`: Numeric vector of length B (bootstrap ATT* values)
#'   * `att_obs`:   Observed ATT from the original fit
#' @export
gsc_boot <- function(fit, B = 499L, alpha = 0.05, seed = NULL) {
  if (!inherits(fit, "coresynth") || fit$method != "gsc") {
    stop("gsc_boot() requires a coresynth object with method = 'gsc'.")
  }
  if (isTRUE(fit$staggered)) {
    stop("gsc_boot() does not support staggered adoption fits. ",
         "Use gsc_inference() instead.", call. = FALSE)
  }
  if (is.null(fit$Y_co_all)) {
    stop("fit is missing Y_co_all. Re-run scm_fit() with the current version.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  F_hat <- fit$F
  L_co_hat <- fit$L_co
  L_tr_hat <- fit$L_tr
  r <- fit$r
  T_pre <- fit$T_pre
  att_obs <- fit$estimate
  Y_co_all <- fit$Y_co_all
  beta <- fit$beta # numeric(0) if no covariates
  X_tr_arr <- fit$covariate_array_tr # T x N_tr x p  or NULL
  X_co_arr <- fit$covariate_array_co # T x N_co x p  or NULL

  T <- nrow(F_hat)
  N_co <- ncol(Y_co_all)
  N_tr <- nrow(L_tr_hat)
  T_post <- T - T_pre
  pre <- seq_len(T_pre)
  post <- (T_pre + 1L):T

  # Full prediction for controls (includes beta contribution if covariates used)
  Y_hat_co <- if (!is.null(fit$Y_hat_co)) {
    fit$Y_hat_co
  } else {
    F_hat %*% t(L_co_hat)
  }

  sigma2_hat <- mean((Y_co_all - Y_hat_co)^2)
  sd_hat <- sqrt(sigma2_hat)

  # Full prediction for treated pre/post (includes beta contribution)
  Y_hat_tr_pre <- F_hat[pre, , drop = FALSE] %*% t(L_tr_hat) # T_pre x N_tr
  Y_hat_tr_post <- F_hat[post, , drop = FALSE] %*% t(L_tr_hat) # T_post x N_tr
  if (length(beta) > 0L && !is.null(X_tr_arr)) {
    T_pre_n <- length(pre)
    T_post_n <- length(post)
    for (j in seq_along(beta)) {
      Y_hat_tr_pre <- Y_hat_tr_pre +
        matrix(X_tr_arr[pre, , j], T_pre_n, N_tr) * beta[j]
      Y_hat_tr_post <- Y_hat_tr_post +
        matrix(X_tr_arr[post, , j], T_post_n, N_tr) * beta[j]
    }
  }

  # Covariate cubes to pass to gsc_ife_cpp in each replicate
  X_co_arr_boot <- if (!is.null(X_co_arr)) {
    X_co_arr
  } else {
    array(0.0, dim = c(T, N_co, 0L))
  }
  # Pre-treatment slice of treated covariates (for Xu 2017 Step 2 demeaning)
  X_tr_pre_boot <- if (!is.null(X_tr_arr)) {
    # T_pre x N_tr x p
    X_tr_arr[seq_len(T_pre), , , drop = FALSE]
  } else {
    array(0.0, dim = c(T_pre, N_tr, 0L))
  }

  boot_dist <- numeric(B)

  for (b in seq_len(B)) {
    eps_co_star <- matrix(rnorm(T * N_co, 0, sd_hat), T, N_co)
    eps_tr_pre_star <- matrix(rnorm(T_pre * N_tr, 0, sd_hat), T_pre, N_tr)
    eps_tr_post_star <- matrix(rnorm(T_post * N_tr, 0, sd_hat), T_post, N_tr)

    Y_star_co <- Y_hat_co + eps_co_star
    Y_star_tr_pre <- Y_hat_tr_pre + eps_tr_pre_star

    res_star <- tryCatch(
      gsc_ife_cpp(Y_star_co, Y_star_tr_pre, r, X_co_arr_boot, X_tr_pre_boot),
      error = function(e) NULL
    )
    if (is.null(res_star)) {
      boot_dist[b] <- NA_real_
      next
    }

    Y_star_tr_post <- Y_hat_tr_post + eps_tr_post_star
    Y_hat_star_post <- res_star$Y_tr_hat[post, , drop = FALSE] # factor part

    # Add beta_star covariate contribution to Y_hat_star_post
    beta_star <- res_star$beta
    if (length(beta_star) > 0L && !is.null(X_tr_arr)) {
      T_post_n <- length(post)
      for (j in seq_along(beta_star)) {
        Y_hat_star_post <- Y_hat_star_post +
          matrix(X_tr_arr[post, , j], T_post_n, N_tr) * beta_star[j]
      }
    }

    boot_dist[b] <- mean(Y_star_tr_post - Y_hat_star_post)
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

  valid <- boot_dist[!is.na(boot_dist)]
  ci <- quantile(valid, c(alpha / 2, 1 - alpha / 2), names = FALSE)

  list(
    p_value = mean(abs(valid) >= abs(att_obs)),
    ci_lower = ci[1],
    ci_upper = ci[2],
    se = sd(valid),
    boot_dist = boot_dist,
    att_obs = att_obs
  )
}

# ── Internal helper: re-estimate one GSC cohort with given control indices ─────
.refit_gsc_cohort <- function(Y_all, cf, r, cov_array_all) {
  TT        <- nrow(Y_all)
  T_pre_g   <- cf$T_pre
  pre_rows  <- seq_len(T_pre_g)
  post_rows <- (T_pre_g + 1L):TT
  N_co_g    <- length(cf$idx_co)

  Y_co_g   <- Y_all[, cf$idx_co, drop = FALSE]
  Y_tr_pre <- Y_all[pre_rows, cf$idx_tr, drop = FALSE]

  p <- if (!is.null(cov_array_all)) dim(cov_array_all)[3L] else 0L
  if (p > 0L) {
    X_co_g     <- cov_array_all[, cf$idx_co, , drop = FALSE]
    X_tr_pre_g <- cov_array_all[pre_rows, cf$idx_tr, , drop = FALSE]
    X_tr_full  <- cov_array_all[, cf$idx_tr, , drop = FALSE]
  } else {
    X_co_g     <- array(0.0, dim = c(TT, N_co_g, 0L))
    X_tr_pre_g <- array(0.0, dim = c(T_pre_g, length(cf$idx_tr), 0L))
    X_tr_full  <- NULL
  }

  res <- tryCatch(
    gsc_ife_cpp(Y_co_g, Y_tr_pre, r, X_co_g, X_tr_pre_g),
    error = function(e) NULL
  )
  if (is.null(res)) return(NA_real_)

  Y_tr_hat <- res$Y_tr_hat  # T × N_tr_g (factor part)
  if (p > 0L) {
    for (j in seq_len(p))
      Y_tr_hat <- Y_tr_hat +
        matrix(X_tr_full[,, j], TT, length(cf$idx_tr)) * res$beta[j]
  }

  Y_treat_post <- Y_all[post_rows, cf$idx_tr, drop = FALSE]
  mean(Y_treat_post - Y_tr_hat[post_rows, , drop = FALSE])
}

#' Non-parametric Inference for GSC (Xu 2017)
#'
#' Estimates SE and confidence intervals for the ATT via non-parametric cluster
#' bootstrap or jackknife over control units. Works for both sharp and staggered
#' GSC fits. For staggered fits, bootstrap resamples each cohort's control pool
#' independently, and jackknife uses a per-cohort LOO with delta-method variance
#' aggregation.
#'
#' Note: `gsc_boot()` performs a *parametric* bootstrap under H0 (hypothesis
#' testing). `gsc_inference()` provides *non-parametric* SE and CIs suitable
#' for inference about the ATT magnitude.
#'
#' @param fit   A `coresynth` object from [scm_fit()] with `method = "gsc"`.
#' @param method `"bootstrap"` (default) or `"jackknife"`.
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

  if (!inherits(fit, "coresynth") || !identical(fit$method, "gsc"))
    stop("gsc_inference() requires a coresynth fit with method = 'gsc'.",
         call. = FALSE)

  tau_hat   <- fit$estimate
  alpha     <- 1 - level
  staggered <- isTRUE(fit$staggered)

  if (!is.null(seed)) set.seed(seed)

  # ── Sharp path ──────────────────────────────────────────────────────────────
  if (!staggered) {
    if (method == "jackknife_global")
      stop("gsc_inference() method='jackknife_global' requires a staggered fit.",
           call. = FALSE)
    if (is.null(fit$Y_co_all))
      stop("fit does not contain Y_co_all. Re-estimate with the current version.",
           call. = FALSE)

    Y_co  <- fit$Y_co_all
    r     <- fit$r
    T_pre <- fit$T_pre
    N_co  <- ncol(Y_co)
    TT    <- nrow(Y_co)
    post_rows  <- (T_pre + 1L):TT
    Y_tr_pre   <- fit$Y_treat[seq_len(T_pre), , drop = FALSE]

    X_co_arr   <- if (!is.null(fit$covariate_array_co)) fit$covariate_array_co else
                    array(0.0, dim = c(TT, N_co, 0L))
    X_tr_pre   <- if (!is.null(fit$covariate_array_tr))
                    fit$covariate_array_tr[seq_len(T_pre), , , drop = FALSE] else
                    array(0.0, dim = c(T_pre, ncol(fit$Y_treat), 0L))
    X_tr_full  <- fit$covariate_array_tr
    p          <- dim(X_co_arr)[3L]

    .refit_sharp <- function(co_idx) {
      N_b <- length(co_idx)
      X_co_b <- if (p > 0L) X_co_arr[, co_idx, , drop = FALSE] else
                  array(0.0, dim = c(TT, N_b, 0L))
      res <- tryCatch(
        gsc_ife_cpp(Y_co[, co_idx, drop = FALSE], Y_tr_pre, r, X_co_b, X_tr_pre),
        error = function(e) NULL
      )
      if (is.null(res)) return(NA_real_)
      Y_hat_post <- res$Y_tr_hat[post_rows, , drop = FALSE]
      if (p > 0L && !is.null(X_tr_full)) {
        N_tr <- ncol(Y_tr_pre)
        T_post <- length(post_rows)
        for (j in seq_len(p))
          Y_hat_post <- Y_hat_post +
            matrix(X_tr_full[post_rows,, j], T_post, N_tr) * res$beta[j]
      }
      mean(fit$Y_treat[post_rows, , drop = FALSE] - Y_hat_post)
    }

    if (method == "bootstrap") {
      boot_ests <- replicate(n_boot, .refit_sharp(sample(N_co, replace = TRUE)))
      n_fail <- sum(is.na(boot_ests))
      if (n_fail > 0.1 * n_boot)
        warning(sprintf("%.0f%% of bootstrap replications failed.", 100 * n_fail / n_boot))
      valid    <- boot_ests[!is.na(boot_ests)]
      se       <- sd(valid)
      ci_lower <- unname(quantile(valid, alpha / 2))
      ci_upper <- unname(quantile(valid, 1 - alpha / 2))
      z        <- tau_hat / max(se, .Machine$double.eps)
      p_value  <- switch(alternative,
        two.sided = 2 * pnorm(-abs(z)),
        greater   = pnorm(-z),
        less      = pnorm(z)
      )
      return(structure(list(
        estimate    = tau_hat, se = se, p_value = p_value,
        ci_lower    = ci_lower, ci_upper = ci_upper,
        method      = method, staggered = FALSE, n_controls = N_co,
        alternative = alternative, boot_ests = boot_ests
      ), class = "coresynth_inference"))
    }

    # jackknife
    jack_ests <- vapply(seq_len(N_co),
                        function(i) .refit_sharp(seq_len(N_co)[-i]),
                        numeric(1))
    jack_var <- (N_co - 1) / N_co * sum((jack_ests - mean(jack_ests))^2)
    se       <- sqrt(jack_var)
    z        <- tau_hat / max(se, .Machine$double.eps)
    ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
    ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se
    p_value  <- switch(alternative,
      two.sided = 2 * pnorm(-abs(z)),
      greater   = pnorm(-z),
      less      = pnorm(z)
    )
    return(structure(list(
      estimate    = tau_hat, se = se, p_value = p_value,
      ci_lower    = ci_lower, ci_upper = ci_upper,
      method      = method, staggered = FALSE, n_controls = N_co,
      alternative = alternative, boot_ests = NULL
    ), class = "coresynth_inference"))
  }

  # ── Staggered path ──────────────────────────────────────────────────────────
  if (is.null(fit$Y_all))
    stop("fit does not contain Y_all. Re-estimate with the current version.",
         call. = FALSE)

  Y_all       <- fit$Y_all
  r           <- fit$r
  cov_all     <- fit$covariate_array_all  # NULL if no covariates
  cohort_fits <- fit$cohort_fits
  w_g         <- vapply(cohort_fits, `[[`, numeric(1L), "weight")  # unnorm
  W           <- sum(w_g)

  .boot_one_gsc <- function() {
    att_total <- 0
    for (cf in cohort_fits) {
      N_co_g  <- length(cf$idx_co)
      idx_b   <- cf$idx_co[sample(N_co_g, replace = TRUE)]
      cf_b    <- cf
      cf_b$idx_co <- idx_b
      att_g_b <- .refit_gsc_cohort(Y_all, cf_b, r, cov_all)
      att_total <- att_total + cf$weight * att_g_b
    }
    att_total / W
  }

  if (method == "bootstrap") {
    boot_ests <- replicate(n_boot, .boot_one_gsc())
    n_fail <- sum(is.na(boot_ests))
    if (n_fail > 0.1 * n_boot)
      warning(sprintf("%.0f%% of staggered bootstrap replications failed.",
                      100 * n_fail / n_boot))
    valid    <- boot_ests[!is.na(boot_ests)]
    se       <- sd(valid)
    ci_lower <- unname(quantile(valid, alpha / 2))
    ci_upper <- unname(quantile(valid, 1 - alpha / 2))
    z        <- tau_hat / max(se, .Machine$double.eps)
    p_value  <- switch(alternative,
      two.sided = 2 * pnorm(-abs(z)),
      greater   = pnorm(-z),
      less      = pnorm(z)
    )
    n_co_mean <- round(mean(vapply(cohort_fits,
                                   function(cf) length(cf$idx_co), numeric(1L))))
    return(structure(list(
      estimate    = tau_hat, se = se, p_value = p_value,
      ci_lower    = ci_lower, ci_upper = ci_upper,
      method      = method, staggered = TRUE, n_controls = n_co_mean,
      alternative = alternative, boot_ests = boot_ests
    ), class = "coresynth_inference"))
  }

  if (method == "jackknife_global") {
    # Global jackknife: drop one unique control across ALL cohorts simultaneously.
    # Captures cross-cohort correlation ignored by per-cohort LOO.
    all_co    <- sort(unique(unlist(lapply(cohort_fits, `[[`, "idx_co"))))
    orig_ests <- vapply(cohort_fits, `[[`, numeric(1L), "estimate")
    jack_ests <- vapply(all_co, function(i) {
      att_total <- 0
      for (k in seq_along(cohort_fits)) {
        cf <- cohort_fits[[k]]
        if (i %in% cf$idx_co) {
          keep <- cf$idx_co[cf$idx_co != i]
          if (length(keep) < 2L) return(NA_real_)
          cf_loo        <- cf
          cf_loo$idx_co <- keep
          att_g <- .refit_gsc_cohort(Y_all, cf_loo, r, cov_all)
          if (is.na(att_g)) return(NA_real_)
        } else {
          att_g <- orig_ests[k]
        }
        att_total <- att_total + cf$weight * att_g
      }
      att_total / W
    }, numeric(1L))
    valid <- jack_ests[!is.na(jack_ests)]
    N_v   <- length(valid)
    se    <- sqrt((N_v - 1L) / N_v * sum((valid - mean(valid))^2))
    z        <- tau_hat / max(se, .Machine$double.eps)
    ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
    ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se
    p_value  <- switch(alternative,
      two.sided = 2 * pnorm(-abs(z)),
      greater   = pnorm(-z),
      less      = pnorm(z)
    )
    return(structure(list(
      estimate    = tau_hat, se = se, p_value = p_value,
      ci_lower    = ci_lower, ci_upper = ci_upper,
      method      = method, staggered = TRUE, n_controls = length(all_co),
      alternative = alternative, boot_ests = NULL
    ), class = "coresynth_inference"))
  }

  # staggered jackknife: per-cohort LOO + delta-method
  var_components <- vapply(cohort_fits, function(cf) {
    N_co_g <- length(cf$idx_co)
    if (N_co_g < 2L) return(0)
    jack_g <- vapply(seq_len(N_co_g), function(i) {
      cf_loo <- cf
      cf_loo$idx_co <- cf$idx_co[-i]
      .refit_gsc_cohort(Y_all, cf_loo, r, cov_all)
    }, numeric(1))
    valid_g <- jack_g[!is.na(jack_g)]
    if (length(valid_g) < 2L) return(0)
    (N_co_g - 1) / N_co_g * sum((valid_g - mean(valid_g))^2)
  }, numeric(1))

  var_st   <- sum((w_g / W)^2 * var_components)
  se       <- sqrt(var_st)
  z        <- tau_hat / max(se, .Machine$double.eps)
  ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
  ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se
  p_value  <- switch(alternative,
    two.sided = 2 * pnorm(-abs(z)),
    greater   = pnorm(-z),
    less      = pnorm(z)
  )
  n_co_mean <- round(mean(vapply(cohort_fits,
                                 function(cf) length(cf$idx_co), numeric(1L))))
  structure(list(
    estimate    = tau_hat, se = se, p_value = p_value,
    ci_lower    = ci_lower, ci_upper = ci_upper,
    method      = method, staggered = TRUE, n_controls = n_co_mean,
    alternative = alternative, boot_ests = NULL
  ), class = "coresynth_inference")
}
