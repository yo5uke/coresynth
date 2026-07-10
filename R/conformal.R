# ── Conformal Inference (Chernozhukov, Wuthrich & Zhu 2021, JASA) ─────────────
# Permutation inference based on permuting blocks of estimated residuals.
# The counterfactual proxy is re-estimated under the null (post-treatment
# treated outcomes imputed as Y1_t - tau0) using ALL T periods, which is
# essential for the finite-sample validity of the procedure (CWZ 2021 S.2.2).

# Internal: extract the donor matrix (T x N_co) and treated series (length T)
# from a sharp coresynth fit. Returns NULL if the fit is unsupported.
.conformal_extract <- function(fit) {
  Yco <- donor_outcomes(fit)
  if (is.null(Yco)) return(NULL)
  list(Yco = Yco, y1 = treated_outcomes(fit))
}

# Internal generic: re-estimate the treated counterfactual P_hat (length T)
# under the null, using all T periods of the (null-imputed) treated series y1
# and the observed donor matrix Yco. Dispatches on the estimator class.
.conformal_refit <- function(fit, Yco, y1) UseMethod(".conformal_refit")

.conformal_refit.default <- function(fit, Yco, y1) {
  stop(sprintf("conformal_inference(): method '%s' is not supported.",
               fit$method), call. = FALSE)
}

.conformal_refit.coresynth_scm <- function(fit, Yco, y1) {
  # Canonical SC over all T: simplex weights minimising ||y1 - Yco w||^2.
  TT <- nrow(Yco)
  V  <- rep(1 / TT, TT)
  w  <- drop(scm_inner_weights_cpp(Yco, y1, V))
  drop(Yco %*% w)
}

.conformal_refit.coresynth_si <- function(fit, Yco, y1) {
  TT   <- nrow(Yco)
  N_co <- ncol(Yco)
  k <- fit$k %||% max(1L, floor(sqrt(min(TT, N_co))))
  k <- as.integer(min(k, min(TT, N_co)))
  res <- si_pcr_cpp(Yco, Yco, matrix(y1, ncol = 1L), k)
  drop(res$Y_hat)
}

.conformal_refit.coresynth_sdid <- function(fit, Yco, y1) {
  # SDID synthetic series (omega-weighted donors + concentrated intercept),
  # unit weights estimated over all T under the null.
  TT    <- nrow(Yco)
  zeta2 <- fit$zeta2 %||% (sqrt(TT) * var(diff(y1)))
  mu_co <- colMeans(Yco)
  mu_y1 <- mean(y1)
  Yco_c <- sweep(Yco, 2L, mu_co, "-")
  y1_c  <- y1 - mu_y1
  omega  <- drop(sdid_unit_weights_cpp(Yco_c, y1_c, zeta2))
  omega0 <- mu_y1 - sum(omega * mu_co)
  omega0 + drop(Yco %*% omega)
}

.conformal_refit.coresynth_gsc <- function(fit, Yco, y1) {
  TT   <- nrow(Yco)
  N_co <- ncol(Yco)
  r <- fit$r %||% 2L
  r <- as.integer(min(r, min(TT, N_co)))
  empty_co <- array(0.0, dim = c(TT, N_co, 0L))
  empty_tr <- array(0.0, dim = c(TT, 1L, 0L))
  res <- gsc_ife_cpp(Yco, matrix(y1, ncol = 1L), r, empty_co, empty_tr)
  drop(res$Y_tr_hat)
}

.conformal_refit.coresynth_mc <- function(fit, Yco, y1) {
  TT   <- nrow(Yco)
  N_co <- ncol(Yco)
  M <- cbind(Yco, y1)                       # T x (N_co + 1), treated last
  O <- matrix(1.0, nrow = TT, ncol = N_co + 1L)  # fully observed under null
  lambda <- fit$lambda %||% (0.01 * svd(M, nu = 0, nv = 0)$d[1])
  L <- soft_impute_cpp(M, O, lambda)
  drop(L[, N_co + 1L])
}

# Internal: residual vector (length T) under the null H0: tau = tau0.
.conformal_residuals <- function(fit, Yco, y1_obs, T_pre, tau0) {
  TT   <- length(y1_obs)
  post <- (T_pre + 1L):TT
  y1_null <- y1_obs
  y1_null[post] <- y1_null[post] - tau0          # impute counterfactual outcome
  P_hat <- .conformal_refit(fit, Yco, y1_null)
  y1_null - P_hat
}

# Internal: moving-block (cyclic) permutation p-value for a residual vector.
.conformal_pval <- function(u, T_pre, q, alternative) {
  TT   <- length(u)
  post <- (T_pre + 1L):TT
  Sstat <- function(v) {
    up <- v[post]
    if (alternative == "two.sided") {
      mean(abs(up)^q)^(1 / q)
    } else {
      mean(up)                                   # signed: average post residual
    }
  }
  S0 <- Sstat(u)
  base <- seq_len(TT) - 1L
  S_perm <- vapply(0:(TT - 1L), function(j) {
    Sstat(u[((base + j) %% TT) + 1L])
  }, numeric(1L))
  switch(alternative,
    two.sided = mean(S_perm >= S0),
    greater   = mean(S_perm >= S0),
    less      = mean(S_perm <= S0)
  )
}

#' Conformal Inference for Synthetic Control Estimators
#'
#' Implements the permutation-based conformal inference procedure of
#' Chernozhukov, Wuthrich & Zhu (2021, JASA). The test inverts a sharp null
#' \eqn{H_0: \tau = \tau_0} by imputing the treated post-treatment
#' counterfactual as \eqn{Y_{1t} - \tau_0}, re-estimating the counterfactual
#' proxy on **all** \eqn{T} periods (imposing the null), and computing a
#' moving-block permutation p-value from the estimated residuals. A confidence
#' interval is obtained by test inversion over a grid of candidate \eqn{\tau_0}.
#'
#' Supported for **sharp** (single-cohort) fits with `method` in
#' `c("scm", "sdid", "gsc", "mc", "si")`. Staggered, multi-arm, and `tasc`
#' fits are not supported (use `sdid_inference()`, `gsc_inference()`, or
#' `si_inference()` instead).
#'
#' @param fit A `coresynth` object from [scm_fit()].
#' @param tau0 Null value of the ATT for the reported p-value (default 0).
#' @param q Exponent of the \eqn{S_q} test statistic
#'   (`S_q = (T_post^{-1} \sum |u_t|^q)^{1/q}`). Default 1, robust to
#'   heavy-tailed data (CWZ 2021). Used only for `alternative = "two.sided"`;
#'   one-sided tests use the signed mean post-treatment residual.
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param ci Logical; construct a confidence interval by test inversion
#'   (default `TRUE`).
#' @param level Confidence level for the interval (default 0.95).
#' @param grid Optional numeric vector of candidate \eqn{\tau_0} values for test
#'   inversion. When `NULL` (default), a grid of `n_grid` points is centred on
#'   the point estimate with half-width `grid_mult` times the pre-treatment
#'   residual standard deviation.
#' @param n_grid Number of grid points when `grid = NULL` (default 200).
#' @param grid_mult Half-width multiplier when `grid = NULL` (default 4).
#' @param ... Unused.
#'
#' @return A list of class `c("conformal_inference", "coresynth_inference")`
#'   with `estimate`, `se` (`NA`; conformal has no SE), `p_value` (at `tau0`),
#'   `ci_lower`, `ci_upper`, `method` (`"conformal"`), `n_controls`,
#'   `alternative`, `staggered` (`FALSE`), plus `tau0`, `q`, `grid`, and
#'   `p_grid` (p-values along the grid). Compatible with [tidy()] / [glance()].
#'
#' @references Chernozhukov, V., Wuthrich, K., & Zhu, Y. (2021). An Exact and
#'   Robust Conformal Inference Method for Counterfactual and Synthetic
#'   Controls. *Journal of the American Statistical Association*, 116(536),
#'   1849-1864.
#'
#' @export
conformal_inference <- function(
  fit,
  tau0        = 0,
  q           = 1,
  alternative = c("two.sided", "greater", "less"),
  ci          = TRUE,
  level       = 0.95,
  grid        = NULL,
  n_grid      = 200L,
  grid_mult   = 4,
  ...
) {
  alternative <- match.arg(alternative)
  if (!inherits(fit, "coresynth"))
    stop("conformal_inference() requires a coresynth object.", call. = FALSE)
  if (inherits(fit, "coresynth_staggered") || inherits(fit, "coresynth_multiarm"))
    stop("conformal_inference() supports sharp (single-cohort) fits only.\n",
         "  For staggered/multi-arm fits use sdid_inference(), gsc_inference(), ",
         "or si_inference().", call. = FALSE)

  method <- fit$method
  if (!method %in% c("scm", "sdid", "gsc", "mc", "si"))
    stop(sprintf("conformal_inference() does not support method = '%s'.\n", method),
         "  Supported: 'scm', 'sdid', 'gsc', 'mc', 'si'.", call. = FALSE)

  ex <- .conformal_extract(fit)
  if (is.null(ex))
    stop("fit does not contain the donor matrix needed for conformal inference. ",
         "Re-estimate with the current version of coresynth.", call. = FALSE)

  Yco    <- ex$Yco
  y1_obs <- ex$y1
  T_pre  <- fit$T_pre
  TT     <- length(y1_obs)
  if (TT - T_pre < 1L)
    stop("conformal_inference() requires at least one post-treatment period.",
         call. = FALSE)

  p_at <- function(t0) {
    u <- .conformal_residuals(fit, Yco, y1_obs, T_pre, t0)
    .conformal_pval(u, T_pre, q, alternative)
  }

  p_value <- p_at(tau0)
  alpha   <- 1 - level

  ci_lower <- NA_real_
  ci_upper <- NA_real_
  grid_used <- NULL
  p_grid    <- NULL

  if (isTRUE(ci)) {
    if (is.null(grid)) {
      gap_pre <- (y1_obs - .conformal_refit(fit, Yco, y1_obs))[seq_len(T_pre)]
      s <- stats::sd(gap_pre)
      if (!is.finite(s) || s < 1e-8) s <- max(abs(fit$estimate), 1)
      half <- grid_mult * s
      grid_used <- seq(fit$estimate - half, fit$estimate + half, length.out = n_grid)
    } else {
      grid_used <- sort(as.numeric(grid))
    }
    p_grid <- vapply(grid_used, p_at, numeric(1L))
    in_ci  <- p_grid >= alpha
    if (any(in_ci)) {
      ci_lower <- min(grid_used[in_ci])
      ci_upper <- max(grid_used[in_ci])
    }
  }

  structure(list(
    estimate    = fit$estimate,
    se          = NA_real_,
    p_value     = p_value,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    method      = "conformal",
    n_controls  = ncol(Yco),
    alternative = alternative,
    staggered   = FALSE,
    boot_ests   = NULL,
    tau0        = tau0,
    q           = q,
    grid        = grid_used,
    p_grid      = p_grid,
    base_method = method
  ), class = c("conformal_inference", "coresynth_inference"))
}

#' @export
print.conformal_inference <- function(x, digits = 4, ...) {
  cat(sprintf("Conformal Inference (CWZ 2021, base = %s, alternative = %s)\n",
              x$base_method, x$alternative))
  cat("  Estimate (ATT):", round(x$estimate, digits), "\n")
  cat(sprintf("  p-value (H0: tau = %s): %s\n",
              format(x$tau0), round(x$p_value, digits)))
  if (!is.na(x$ci_lower))
    cat("  CI            : [", round(x$ci_lower, digits), ", ",
        round(x$ci_upper, digits), "]\n", sep = "")
  else if (!is.null(x$grid))
    cat("  CI            : empty at this level (widen 'grid' or raise level)\n")
  cat("  Controls      :", x$n_controls, "\n")
  invisible(x)
}
