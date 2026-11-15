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
# `prep` is what .conformal_prepare() returned for the same fit and donors.
.conformal_refit <- function(fit, Yco, y1, prep = NULL)
  UseMethod(".conformal_refit")

# Internal generic: the part of the refit that the donors alone determine.
# The test inversion refits once per grid point with only y1 changing, so an
# estimator whose donor-side step does not involve y1 computes it here once.
# NULL (the default) means the refit has no such step.
.conformal_prepare <- function(fit, Yco) UseMethod(".conformal_prepare")

.conformal_prepare.default <- function(fit, Yco) NULL

# The refit has to impose the model the fit used, so a tuning parameter the
# fit does not carry cannot be guessed: substituting some other rule would
# score the null against a differently specified counterfactual. Every
# current fit stores its own, so this only catches objects saved by an
# earlier version. Callers read the field with [[ ]], not $: a GSC fit also
# carries r_cv and an MC fit lambda_cv, which $ would partial-match to if the
# exact name were ever absent.
.conformal_tuning <- function(value, name, method) {
  if (is.null(value))
    stop(sprintf(paste0(
      "this %s fit does not record '%s', so conformal_inference() cannot ",
      "refit under the same model. Re-run scm_fit() with the current ",
      "version of coresynth."), method, name), call. = FALSE)
  value
}

.conformal_refit.default <- function(fit, Yco, y1, prep = NULL) {
  stop(sprintf("conformal_inference(): method '%s' is not supported.",
               fit$method), call. = FALSE)
}

.conformal_refit.coresynth_scm <- function(fit, Yco, y1, prep = NULL) {
  # Canonical SC over all T: simplex weights minimising ||y1 - Yco w||^2.
  TT <- nrow(Yco)
  V  <- rep(1 / TT, TT)
  w  <- drop(scm_inner_weights_cpp(Yco, y1, V))
  drop(Yco %*% w)
}

# PCR is linear in the treated series: the counterfactual is H y1 with
# H = Yco pinv_k(Yco), which the donors and k fix.
.conformal_prepare.coresynth_si <- function(fit, Yco) {
  TT   <- nrow(Yco)
  N_co <- ncol(Yco)
  k <- as.integer(min(.conformal_tuning(fit[["k"]], "k", "si"), TT, N_co))
  si_pcr_cpp(Yco, Yco, diag(TT), k)$Y_hat
}

.conformal_refit.coresynth_si <- function(fit, Yco, y1, prep = NULL) {
  H <- prep %||% .conformal_prepare(fit, Yco)
  drop(H %*% y1)
}

.conformal_refit.coresynth_sdid <- function(fit, Yco, y1, prep = NULL) {
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

# Step 1 of Xu (2017) -- factors, beta and the additive effects -- uses the
# control units only, so it is the same fit at every null value.
.conformal_prepare.coresynth_gsc <- function(fit, Yco) {
  TT   <- nrow(Yco)
  N_co <- ncol(Yco)
  r <- as.integer(min(.conformal_tuning(fit[["r"]], "r", "gsc"), TT, N_co))
  # Covariates stay in the model: y1 is the treated units' average, whose
  # covariate part under the linear term is the average covariate path.
  X_tr <- NULL
  if (!is.null(fit$covariate_array_tr)) {
    p    <- dim(fit$covariate_array_tr)[3L]
    X_tr <- array(apply(fit$covariate_array_tr, c(1L, 3L), mean), c(TT, 1L, p))
  }
  list(ctrl = .gsc_fit_controls(Yco, fit$covariate_array_co, r,
                                fit[["force"]] %||% "none"),
       X_tr = X_tr)
}

.conformal_refit.coresynth_gsc <- function(fit, Yco, y1, prep = NULL) {
  prep <- prep %||% .conformal_prepare(fit, Yco)
  drop(.gsc_counterfactual(prep$ctrl, matrix(y1, ncol = 1L), prep$X_tr,
                           nrow(Yco))$Y_hat)
}

.conformal_refit.coresynth_mc <- function(fit, Yco, y1, prep = NULL) {
  TT   <- nrow(Yco)
  N_co <- ncol(Yco)
  M <- cbind(Yco, y1)                       # T x (N_co + 1), treated last
  O <- matrix(1.0, nrow = TT, ncol = N_co + 1L)  # fully observed under null
  na_mask <- is.na(M)
  if (any(na_mask)) {
    O[na_mask] <- 0.0
    M[na_mask] <- 0.0
  }
  lambda <- .conformal_tuning(fit[["lambda"]], "lambda", "mc")
  # Refit under the same model the fit used: leaving the fixed effects out
  # here would score the null against a differently specified counterfactual.
  L <- soft_impute_fe_cpp(M, O, lambda,
                          estimate_fe = isTRUE(fit[["fixedeff"]] %||% FALSE))$M
  drop(L[, N_co + 1L])
}

# Internal: residual vector (length T) under the null H0: tau = tau0.
.conformal_residuals <- function(fit, Yco, y1_obs, T_pre, tau0, prep = NULL) {
  TT   <- length(y1_obs)
  post <- (T_pre + 1L):TT
  y1_null <- y1_obs
  y1_null[post] <- y1_null[post] - tau0          # impute counterfactual outcome
  P_hat <- .conformal_refit(fit, Yco, y1_null, prep)
  y1_null - P_hat
}

# Internal: moving-block (cyclic) permutation p-value for a residual vector.
.conformal_pval <- function(u, T_pre, q, alternative) {
  TT   <- length(u)
  post <- (T_pre + 1L):TT
  Sstat <- function(v) {
    up <- v[post]
    if (alternative == "two.sided") {
      mean(abs(up)^q, na.rm = TRUE)^(1 / q)
    } else {
      mean(up, na.rm = TRUE)                     # signed: average post residual
    }
  }
  S0 <- Sstat(u)
  base <- seq_len(TT) - 1L
  S_perm <- vapply(0:(TT - 1L), function(j) {
    Sstat(u[((base + j) %% TT) + 1L])
  }, numeric(1L))
  switch(alternative,
    two.sided = mean(S_perm >= S0, na.rm = TRUE),
    greater   = mean(S_perm >= S0, na.rm = TRUE),
    less      = mean(S_perm <= S0, na.rm = TRUE)
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
#' The interval collects the candidates whose p-value exceeds `1 - level`,
#' as in Algorithm 1 of the paper, and is reported as the range they span.
#' The default grid is only a first guess at where that range lies, so while
#' a side of it still ends on an accepted candidate the search moves outwards,
#' in steps that double, until the test rejects; the bound is then resolved to
#' the spacing of the initial grid. A warning is raised if the search stops
#' with a bound still on its edge. That is what happens when the interval is
#' unbounded: the permutation p-value is a multiple of \eqn{1/T}, so a panel
#' with fewer than \eqn{1/(1 - level)} periods rejects no value at all. A grid
#' supplied through `grid` is used as given. A one-sided test leaves the side
#' of its alternative open by construction, so that bound is the end of the
#' grid and is not flagged.
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
#'   inversion, used as given. When `NULL` (default), a grid of `n_grid`
#'   points is centred on the point estimate with half-width `grid_mult` times
#'   the pre-treatment residual standard deviation, and is extended outwards
#'   while the interval reaches its edge.
#' @param n_grid Number of points in the initial grid when `grid = NULL`
#'   (default 200); the grid returned is longer when the search was extended.
#' @param grid_mult Half-width multiplier of the initial grid when
#'   `grid = NULL` (default 4).
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

  prep <- .conformal_prepare(fit, Yco)
  p_at <- function(t0) {
    u <- .conformal_residuals(fit, Yco, y1_obs, T_pre, t0, prep)
    .conformal_pval(u, T_pre, q, alternative)
  }

  p_value <- p_at(tau0)
  alpha   <- 1 - level

  ci_lower <- NA_real_
  ci_upper <- NA_real_
  grid_used <- NULL
  p_grid    <- NULL

  if (isTRUE(ci)) {
    # CWZ (2021, Algorithm 1) keep the candidates whose p-value exceeds alpha.
    # The p-value is a multiple of 1/T, so whenever alpha * T is an integer
    # some candidates sit exactly on alpha; the margin keeps the rounding in
    # 1 - level from deciding which side they fall on.
    accepted <- function(p) p > alpha + 1e-9
    # The identity shift always counts, so no p-value is below 1/T.
    can_reject <- !accepted(1 / TT)
    # A one-sided test never rejects in the direction of its alternative, so
    # that side of the interval is open by construction.
    bounded <- c(lower = alternative != "less", upper = alternative != "greater")

    auto <- is.null(grid)
    if (auto) {
      gap_pre <- (y1_obs - .conformal_refit(fit, Yco, y1_obs, prep))[seq_len(T_pre)]
      s <- stats::sd(gap_pre, na.rm = TRUE)
      if (!is.finite(s) || s < 1e-8) s <- max(abs(fit$estimate), 1)
      half <- grid_mult * s
      grid_used <- seq(fit$estimate - half, fit$estimate + half, length.out = n_grid)
    } else {
      grid_used <- sort(as.numeric(grid))
    }
    p_grid <- vapply(grid_used, p_at, numeric(1L))

    # The default grid is a guess at where the interval lies: a few
    # pre-treatment residual standard deviations around the estimate, which a
    # tight pre-treatment fit makes too narrow. While a bounded side still
    # ends on an accepted candidate, the search moves outwards in steps that
    # double, and the bound it then brackets is resolved back to the spacing
    # of the initial grid by bisection.
    if (auto && can_reject) {
      step  <- 2 * half / max(as.integer(n_grid) - 1L, 1L)
      n_add <- max(4L, as.integer(n_grid) %/% 4L)
      extra_x <- numeric(0L)
      extra_p <- numeric(0L)
      for (side in c(-1, 1)) {
        if (!bounded[[if (side < 0) "lower" else "upper"]]) next
        widened <- FALSE
        for (r in seq_len(6L)) {
          edge_i <- if (side < 0) 1L else length(grid_used)
          if (!accepted(p_grid[edge_i])) break
          new   <- grid_used[edge_i] + side * step * 2^r * seq_len(n_add)
          p_new <- vapply(new, p_at, numeric(1L))
          if (side < 0) {
            grid_used <- c(rev(new), grid_used)
            p_grid    <- c(rev(p_new), p_grid)
          } else {
            grid_used <- c(grid_used, new)
            p_grid    <- c(p_grid, p_new)
          }
          widened <- TRUE
        }
        acc <- which(accepted(p_grid))
        if (!widened || length(acc) == 0L) next
        i_in  <- if (side < 0) min(acc) else max(acc)
        i_out <- i_in + as.integer(side)
        if (i_out < 1L || i_out > length(grid_used)) next
        a <- grid_used[i_in]
        b <- grid_used[i_out]
        while (abs(b - a) > step * (1 + 1e-8)) {
          m  <- (a + b) / 2
          pm <- p_at(m)
          extra_x <- c(extra_x, m)
          extra_p <- c(extra_p, pm)
          if (accepted(pm)) a <- m else b <- m
        }
      }
      if (length(extra_x) > 0L) {
        ord       <- order(c(grid_used, extra_x))
        p_grid    <- c(p_grid, extra_p)[ord]
        grid_used <- c(grid_used, extra_x)[ord]
      }
    }

    in_ci <- accepted(p_grid)
    if (any(in_ci)) {
      ci_lower <- min(grid_used[in_ci])
      ci_upper <- max(grid_used[in_ci])
      # A bound still on the first or last candidate is where the search
      # stopped, not where the test started rejecting.
      n_g  <- length(grid_used)
      edge <- c(if (in_ci[1L]  && bounded[["lower"]]) "lower",
                if (in_ci[n_g] && bounded[["upper"]]) "upper")
      if (length(edge) > 0L) {
        warning(
          "conformal_inference(): ",
          if (length(edge) == 2L) "both confidence bounds sit"
          else paste("the", edge, "confidence bound sits"),
          " on the edge of the search grid, so the interval is cut off ",
          "there and extends beyond it. ",
          if (!can_reject)
            sprintf(paste0(
              "With %d periods the permutation p-value cannot fall below ",
              "1/%d = %.3f, so no value is rejected at level %s and the ",
              "interval is unbounded; lower `level` or use a longer panel."),
              TT, TT, 1 / TT, format(level))
          else if (auto)
            paste0("The search already reached ",
                   formatC(grid_used[1L], format = "g", digits = 4), " to ",
                   formatC(grid_used[n_g], format = "g", digits = 4),
                   " without the test rejecting, so the interval may be ",
                   "unbounded on that side.")
          else "Widen `grid`.",
          call. = FALSE)
      }
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
