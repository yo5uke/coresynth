# ── Relabeling Permutation Test (Liu & Xu 2026, S.5 & Algorithm 6) ──────────
#
# Finite-sample exact randomization/relabeling test of the sharp null of no
# effect in any treated cell (H0: tau_it = 0 for all treated cells). Reassigns
# the observed treatment pattern (including adoption dates under staggered
# adoption) uniformly at random across all units and compares the observed
# aggregate ATT (or |ATT|) against the relabeled distribution.

# Build a method-independent context for relabeling_test() holding the full
# T x N panel, observed adoption-row vector T_adopt (NA = never treated),
# original treated/control column indices, and fixed hyperparameters.
.relabel_context <- function(fit) {
  method    <- fit$method
  staggered <- inherits(fit, "coresynth_staggered")

  if (identical(method, "scm")) {
    if (staggered) {
      pan <- .pan_from_staggered_fit(fit)
      sep <- identical(fit$nu, "separate")
      if (sep && is.null(fit$separate_args))
        stop("this staggered SCM fit lacks 'separate_args'; re-run scm_fit().",
             call. = FALSE)
      return(list(
        method        = "scm",
        staggered     = TRUE,
        pan           = pan,
        Y             = pan$Y,
        T_adopt       = pan$T_adopt,
        idx_tr        = pan$idx_treat,
        idx_co        = pan$idx_control,
        sep           = sep,
        nu            = fit$nu,
        fixedeff      = isTRUE(fit$fixedeff),
        pool_cohorts  = fit$pool_cohorts %||% TRUE,
        control_group = fit$control_group %||% "never_treated",
        separate_args = fit$separate_args,
        Z_balance     = fit$Z_balance,
        treated_w     = fit$treated_weights
      ))
    }

    Y_tr <- fit[["Y_treat_mat"]] %||% as.matrix(fit[["Y_treat"]])
    Y_co <- donor_outcomes(fit)
    if (is.null(Y_tr) || is.null(Y_co))
      stop("SCM fit lacks treated/donor matrices; re-run scm_fit().",
           call. = FALSE)
    n_tr <- ncol(Y_tr)
    if (n_tr > 1L && !is.null(fit[["X0_mat"]]))
      stop(paste0(
        "relabeling_test() is not supported for sharp predictor-based SCM ",
        "with multiple treated units because only the pooled treated ",
        "predictor vector is stored."), call. = FALSE)
    X_all <- if (!is.null(fit[["X0_mat"]])) {
      cbind(fit[["X1_vec"]], fit[["X0_mat"]])
    } else NULL

    T_pre   <- as.integer(fit[["T_pre"]])
    T_adopt <- c(rep(T_pre + 1L, n_tr), rep(NA_integer_, ncol(Y_co)))
    return(list(
      method    = "scm",
      staggered = FALSE,
      Y         = cbind(Y_tr, Y_co),
      T_adopt   = T_adopt,
      idx_tr    = seq_len(n_tr),
      idx_co    = n_tr + seq_len(ncol(Y_co)),
      T_pre     = T_pre,
      X_all     = X_all,
      fit       = fit,
      treated_w = fit[["treated_weights"]]
    ))
  }

  if (identical(method, "sdid")) {
    ctx <- .sdid_inference_context(fit)
    # Keep user-supplied zeta2 fixed; when data-driven (NULL), recompute
    # equivariantly from each relabeled control pool (Proposition 9).
    return(list(
      method    = "sdid",
      staggered = staggered,
      sdid_ctx  = ctx,
      Y         = ctx$Y,
      T_adopt   = ctx$T_adopt,
      idx_tr    = ctx$tr,
      idx_co    = ctx$co
    ))
  }

  if (identical(method, "gsc")) {
    ctx <- .gsc_inference_context(fit)
    ctx$cache <- NULL
    return(list(
      method    = "gsc",
      staggered = staggered,
      gsc_ctx   = ctx,
      Y         = ctx$Y,
      T_adopt   = ctx$T_adopt,
      idx_tr    = ctx$tr,
      idx_co    = ctx$co
    ))
  }

  if (identical(method, "si")) {
    ctx <- .si_inference_context(fit)
    return(list(
      method    = "si",
      staggered = staggered,
      si_ctx    = ctx,
      Y         = ctx$Y,
      T_adopt   = ctx$T_adopt,
      idx_tr    = ctx$tr,
      idx_co    = ctx$co
    ))
  }

  if (method %in% c("mc", "tasc")) {
    Y_tr <- as.matrix(fit[["Y_treat"]])
    Y_co <- donor_outcomes(fit)
    if (is.null(Y_tr) || is.null(Y_co))
      stop(sprintf(
        "%s fit lacks treated/donor matrices; re-run scm_fit().",
        toupper(method)), call. = FALSE)
    n_tr <- ncol(Y_tr)
    n_co <- ncol(Y_co)
    Y    <- cbind(Y_tr, Y_co)

    if (staggered) {
      T_adopt <- rep(NA_integer_, n_tr + n_co)
      for (cf in fit$cohort_fits) {
        pos <- match(cf$idx_tr, fit$idx_tr)
        if (anyNA(pos))
          stop("Could not map staggered cohort units back to Y_treat.",
               call. = FALSE)
        T_adopt[pos] <- as.integer(cf$cohort)
      }
    } else {
      T_pre   <- as.integer(fit[["T_pre"]])
      T_adopt <- c(rep(T_pre + 1L, n_tr), rep(NA_integer_, n_co))
    }

    return(list(
      method    = method,
      staggered = staggered,
      Y         = Y,
      T_adopt   = T_adopt,
      idx_tr    = seq_len(n_tr),
      idx_co    = if (n_co > 0L) n_tr + seq_len(n_co) else integer(0L),
      fit       = fit
    ))
  }

  stop(sprintf("relabeling_test() does not support method = '%s'.", method),
       call. = FALSE)
}

# Re-estimate aggregate ATT for a reassigned treatment pattern (tr_b, adopt_b).
# Returns NA_real_ when a reassigned fit fails.
.relabel_fit_att <- function(ctx, tr_b, co_b, adopt_b) {
  method <- ctx$method
  TT     <- nrow(ctx$Y)

  if (identical(method, "scm")) {
    if (isTRUE(ctx$staggered)) {
      pan_b <- ctx$pan
      pan_b$T_adopt     <- adopt_b
      pan_b$idx_treat   <- tr_b
      pan_b$idx_control <- co_b
      pan_b$T_pre       <- min(adopt_b[tr_b], na.rm = TRUE) - 1L
      D <- matrix(0L, TT, ncol(ctx$Y))
      for (j in tr_b) D[adopt_b[j]:TT, j] <- 1L
      pan_b$D <- D

      tw_b <- if (is.null(ctx$treated_w)) NULL else {
        stats::setNames(as.numeric(ctx$treated_w),
                        as.character(pan_b$units[tr_b]))
      }

      f <- tryCatch(suppressWarnings(suppressMessages(
        if (isTRUE(ctx$sep)) {
          a <- ctx$separate_args
          .fit_scm_staggered(
            pan_b, unit_weights = tw_b, Y_matrix = pan_b$Y,
            v_selection = a$v_selection, v_optim = a$v_optim,
            donor_mspe_threshold = a$donor_mspe_threshold,
            lambda_pen = a$lambda_pen, control_group = ctx$control_group,
            pool_cohorts = ctx$pool_cohorts, fixedeff = ctx$fixedeff)
        } else {
          .fit_scm_staggered_pooled(
            pan_b, Y_matrix = pan_b$Y, unit_weights = tw_b, nu = ctx$nu,
            fixedeff = ctx$fixedeff, control_group = ctx$control_group,
            pool_cohorts = ctx$pool_cohorts, Z = ctx$Z_balance)
        }
      )), error = function(e) NULL)
      return(if (is.null(f) || !is.finite(f$estimate)) NA_real_ else f$estimate)
    }

    # Sharp SCM
    pre  <- seq_len(ctx$T_pre)
    post <- (ctx$T_pre + 1L):TT
    Y_co_pre  <- ctx$Y[pre, co_b, drop = FALSE]
    Y_co_post <- ctx$Y[post, co_b, drop = FALSE]
    Y_tr_mat  <- ctx$Y[, tr_b, drop = FALSE]
    y_tr_all  <- if (length(tr_b) == 1L) {
      as.numeric(Y_tr_mat)
    } else if (is.null(ctx$treated_w)) {
      rowMeans(Y_tr_mat)
    } else {
      w_tr <- as.numeric(ctx$treated_w)
      w_tr <- w_tr / sum(w_tr)
      drop(Y_tr_mat %*% w_tr)
    }
    X_B <- if (!is.null(ctx$X_all)) ctx$X_all[, co_b, drop = FALSE] else NULL
    x_i <- if (!is.null(ctx$X_all)) ctx$X_all[, tr_b[1L]] else NULL

    w <- tryCatch(
      .uc_scm_fit_weights(Y_co_pre, y_tr_all[pre], X_B = X_B, x_i = x_i,
                          fit = ctx$fit),
      error = function(e) NULL
    )
    if (is.null(w) || any(!is.finite(w))) return(NA_real_)
    return(mean(y_tr_all[post] - drop(Y_co_post %*% w)))
  }

  if (identical(method, "sdid")) {
    c_b <- ctx$sdid_ctx
    c_b$T_adopt <- adopt_b
    res <- .sdid_unit_estimate(c_b, co_b, tr_b)
    return(if (is.null(res) || !is.finite(res$estimate)) NA_real_ else res$estimate)
  }

  if (identical(method, "gsc")) {
    c_b <- ctx$gsc_ctx
    c_b$T_adopt <- adopt_b
    c_b$cache   <- NULL
    res <- .gsc_unit_estimate(c_b, co_b, tr_b)
    return(if (is.null(res) || !is.finite(res$estimate)) NA_real_ else res$estimate)
  }

  if (identical(method, "si")) {
    c_b <- ctx$si_ctx
    c_b$T_adopt <- adopt_b
    res <- .si_unit_estimate(c_b, co_b, tr_b)
    return(if (is.null(res) || !is.finite(res$estimate)) NA_real_ else res$estimate)
  }

  # MC and TASC: single full-panel imputation with reassigned mask
  N <- ncol(ctx$Y)
  if (identical(method, "mc")) {
    lam      <- .conformal_tuning(ctx$fit[["lambda"]], "lambda", "mc")
    fixedeff <- isTRUE(ctx$fit[["fixedeff"]] %||% FALSE)
    Y_obs <- ctx$Y
    O     <- matrix(1.0, nrow = TT, ncol = N)
    for (j in tr_b) {
      t0 <- adopt_b[j]
      if (!is.na(t0) && t0 <= TT) O[t0:TT, j] <- 0.0
    }
    na_mask <- !is.finite(Y_obs)
    if (any(na_mask)) {
      O[na_mask]     <- 0.0
      Y_obs[na_mask] <- 0.0
    }
    if (any(rowSums(O) == 0) || any(colSums(O) == 0)) return(NA_real_)
    M_hat <- tryCatch(
      soft_impute_fe_cpp(Y_obs, O, lam, estimate_fe = fixedeff)$M,
      error = function(e) NULL
    )
    if (is.null(M_hat)) return(NA_real_)
  } else {
    # TASC
    r_tasc  <- as.integer(ctx$fit[["r"]] %||% 2L)
    em_iter <- as.integer(ctx$fit[["em_iter"]] %||% 2000L)
    em_tol  <- ctx$fit[["em_tol"]] %||% 1e-6
    fix_A   <- isTRUE(ctx$fit[["fix_A"]])
    interc  <- if (is.null(ctx$fit[["intercept"]])) TRUE else isTRUE(ctx$fit[["intercept"]])
    T_pre_b <- min(adopt_b[tr_b], na.rm = TRUE) - 1L
    if (T_pre_b < 1L) return(NA_real_)

    Y_obs <- ctx$Y
    for (j in tr_b) {
      t0 <- adopt_b[j]
      if (!is.na(t0) && t0 <= TT) Y_obs[t0:TT, j] <- NA_real_
    }
    res <- tryCatch(
      .tasc_em_fit(Y_obs, T_pre_b, idx_tr = tr_b, idx_co = co_b,
                   r = r_tasc, em_iter = em_iter, em_tol = em_tol,
                   fix_A = fix_A, intercept = interc,
                   warn_unconverged = FALSE),
      error = function(e) NULL
    )
    if (is.null(res)) return(NA_real_)
    M_hat <- res$Y_hat_all
  }

  # Clarke et al. (2024) aggregation across reassigned cohorts
  cohorts <- sort(unique(adopt_b[tr_b]))
  est_g   <- numeric(length(cohorts))
  w_g     <- numeric(length(cohorts))
  for (k in seq_along(cohorts)) {
    g    <- cohorts[k]
    tr_g <- tr_b[adopt_b[tr_b] == g]
    post <- g:TT
    est_g[k] <- mean(ctx$Y[post, tr_g, drop = FALSE] -
                       M_hat[post, tr_g, drop = FALSE], na.rm = TRUE)
    w_g[k]   <- length(tr_g) * length(post)
  }
  if (any(!is.finite(est_g)) || sum(w_g) <= 0) return(NA_real_)
  sum(w_g * est_g) / sum(w_g)
}

#' Relabeling Permutation Test for Sharp Nulls
#'
#' Implements the Monte Carlo relabeling test of Liu & Xu (2026, Section 5
#' and SM Algorithm 6) for the sharp null hypothesis of no treatment effect in
#' any treated cell (\eqn{H_0: \tau_{it} = 0}). Works for sharp fits
#' (including multiple treated units) and staggered adoption fits across all six
#' estimators (`"scm"`, `"sdid"`, `"gsc"`, `"mc"`, `"tasc"`, and `"si"`),
#' including `"mc"` and `"tasc"` staggered panels where no resampling-based
#' inference routine is otherwise available.
#'
#' In each of `B` draws, the observed adoption-date pattern across the `K`
#' treated units is reassigned uniformly at random without replacement among all
#' `N` units in the panel, and the estimator is refitted with its tuning
#' parameters (`r`, `lambda`, `k`, `nu`, user-supplied `zeta2`, etc.) held
#' fixed at the original fit's values (for SDID with data-driven `zeta2`,
#' the ridge penalty is recomputed equivariantly from each relabeled control
#' pool per Proposition 9).
#'
#' @param fit A sharp or staggered `coresynth` object from [scm_fit()]
#'   (excluding multi-arm SI fits).
#' @param B Number of random treatment reassignments (default `200L`).
#' @param statistic Test statistic: `"abs_att"` (default, absolute aggregate
#'   ATT for two-sided testing) or `"att"` (signed aggregate ATT). One-sided
#'   alternatives (`"greater"` or `"less"`) require `statistic = "att"`.
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param seed Optional RNG seed for reproducibility (default `NULL`; caller's
#'   `.Random.seed` is preserved).
#'
#' @return A list of class `c("relabeling_test", "coresynth_inference")`
#'   compatible with [tidy()] and [glance()], containing `estimate`, `se`
#'   (`NA_real_`), `p_value`, `ci_lower` (`NA_real_`), `ci_upper` (`NA_real_`),
#'   `method` (`"relabeling"`), `n_controls`, `n_treated`, `alternative`,
#'   `staggered`, `boot_ests` (vector of valid relabeled ATT estimates), `B`,
#'   `n_valid`, `statistic`, and `base_method`.
#'
#' @references Liu, Z., & Xu, Y. (2026). Unit-Level Conformal Inference for
#'   Panel Counterfactuals. *arXiv:2611.19209*.
#'
#' @seealso [unit_conformal()] for single-treated-unit confidence sets,
#'   [scm_rank_pval()] for the rank-based SCM permutation test.
#' @export
relabeling_test <- function(
  fit,
  B           = 200L,
  statistic   = c("abs_att", "att"),
  alternative = c("two.sided", "greater", "less"),
  seed        = NULL
) {
  statistic   <- match.arg(statistic)
  alternative <- match.arg(alternative)

  if (!inherits(fit, "coresynth"))
    stop("relabeling_test() requires a coresynth object from scm_fit().",
         call. = FALSE)
  if (inherits(fit, "coresynth_multiarm"))
    stop("relabeling_test() does not support multi-arm SI fits; use si_inference().",
         call. = FALSE)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1 ||
      B != trunc(B))
    stop("`B` must be a single positive integer.", call. = FALSE)
  B <- as.integer(B)
  if (alternative != "two.sided" && identical(statistic, "abs_att"))
    stop("One-sided tests ('greater'/'less') require statistic = 'att', not 'abs_att'.",
         call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed)))
    stop("`seed` must be NULL or a single finite number.", call. = FALSE)

  ctx <- .relabel_context(fit)
  N   <- ncol(ctx$Y)
  K   <- length(ctx$idx_tr)
  if (N - K < 1L && !fit$method %in% c("mc", "tasc"))
    stop("relabeling_test() requires at least 1 control unit.", call. = FALSE)
  if (N < 2L)
    stop("relabeling_test() requires at least 2 units.", call. = FALSE)

  adopt_slots <- ctx$T_adopt[ctx$idx_tr]
  obs_att     <- fit$estimate

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed)
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  ests <- numeric(B)
  for (b in seq_len(B)) {
    tr_b <- sample.int(N, K)
    co_b <- setdiff(seq_len(N), tr_b)
    adopt_b       <- rep(NA_integer_, N)
    adopt_b[tr_b] <- adopt_slots
    ests[b] <- .relabel_fit_att(ctx, tr_b, co_b, adopt_b)
  }

  valid  <- ests[is.finite(ests)]
  n_fail <- B - length(valid)
  if (length(valid) == 0L)
    stop("relabeling_test(): all ", B, " relabeled refits failed.",
         call. = FALSE)
  if (n_fail > 0L)
    warning(sprintf(paste0(
      "relabeling_test(): %d of %d relabeled refits failed; counted in ",
      "favor of the null per Proposition 6."),
      n_fail, B), call. = FALSE)

  # Per Proposition 6, failed draws receive an extreme score in favor of H0
  # while keeping the denominator fixed at B + 1.
  n_extreme <- switch(alternative,
    two.sided = sum(abs(valid) >= abs(obs_att) * (1 - 1e-12)),
    greater   = sum(valid >= obs_att - 1e-12 * max(1, abs(obs_att))),
    less      = sum(valid <= obs_att + 1e-12 * max(1, abs(obs_att)))
  )
  p_value <- (1 + n_fail + n_extreme) / (B + 1)

  structure(list(
    estimate    = obs_att,
    se          = NA_real_,
    p_value     = p_value,
    ci_lower    = NA_real_,
    ci_upper    = NA_real_,
    method      = "relabeling",
    n_controls  = length(ctx$idx_co),
    n_treated   = K,
    alternative = alternative,
    staggered   = isTRUE(ctx$staggered),
    boot_ests   = valid,
    B           = B,
    n_valid     = length(valid),
    statistic   = statistic,
    base_method = ctx$method
  ), class = c("relabeling_test", "coresynth_inference"))
}

#' @export
print.relabeling_test <- function(x, digits = 4, ...) {
  cat(sprintf("Relabeling Permutation Test (Liu & Xu 2026, base = %s%s, alternative = %s)\n",
              x$base_method, if (isTRUE(x$staggered)) ", staggered" else "",
              x$alternative))
  cat("  Estimate (ATT):", round(x$estimate, digits), "\n")
  cat("  p-value       :", round(x$p_value, digits), "\n")
  cat(sprintf("  Draws         : %d valid / %d requested\n",
              x$n_valid, x$B))
  cat("  Treated units :", x$n_treated, "| Controls:", x$n_controls, "\n")
  invisible(x)
}
