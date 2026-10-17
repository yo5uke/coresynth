#' Fit Synthetic Difference-in-Differences (Arkhangelsky et al. 2021)
#'
#' Implements Algorithm 1 from Arkhangelsky et al. (2021) exactly:
#'
#' 1. **Unit weights** omega (with implicit intercept omega_0 concentrated out):
#'    zeta = (N_tr * T_post)^{1/4} * sigma_hat, where sigma_hat^2 = Var of
#'    first differences of control outcomes in pre-period (Eq. 2.2).
#'    Column-demeaning Y_co_pre/Y_tr_pre eliminates omega_0 from the QP.
#'
#' 2. **Time weights** lambda (with implicit intercept lambda_0 concentrated out):
#'    Row-demeaning Y_co_pre and demeaning Ybar_co_post eliminates lambda_0.
#'    Tiny zeta_t = 1e-6 * sigma_hat for uniqueness (Eq. 2.3, footnote 3).
#'
#' 3. **SDID estimate** (closed-form DiD; omega_0/lambda_0 cancel algebraically):
#'    tau = (Ybar_tr_post - lambda'Y_tr_pre) - (omega'Ybar_co_post - omega'Y_co_pre'lambda)
#'
#' For staggered adoption, applies the cohort-based approach from Arkhangelsky
#' et al. (2021) Appendix Section 8: SDID is run separately for each adoption
#' cohort and the resulting estimates are averaged with weights proportional to
#' N_treated * T_post for that cohort.
#'
#' @param y     Outcome vector (long format)
#' @param d     Treatment indicator (long format)
#' @param id    Unit identifier (long format)
#' @param time  Time identifier (long format)
#' @param zeta2 Override zeta^2 (penalty BEFORE T_pre multiplication). Default: auto.
#' @param zeta_t Base tiny ridge scale for time weight uniqueness (default 1e-6).
#' @param covariates Character vector of time-varying covariate column names in
#'   `data`. When provided, covariates are partialled out of Y before
#'   computing SDID weights and the ATT, by the "projected" method of Clarke
#'   et al. (2024, S.2.2) and Kranz (2022): the coefficients come from a
#'   regression of Y on the covariates with unit and time fixed effects over
#'   every untreated (Wit = 0) cell, sharp and staggered alike. Requires
#'   `data`, `id_var`, and `time_var` to be supplied (done automatically when
#'   called via [scm_fit()]).
#' @param data        Full long-format data frame (passed by [scm_fit()]).
#' @param id_var      Name of the unit identifier column (passed by [scm_fit()]).
#' @param time_var    Name of the time identifier column (passed by [scm_fit()]).
#' @param control_group For staggered adoption: which units serve as controls
#'   for each cohort. `"never_treated"` (default) uses the pure controls of
#'   Clarke et al. (2024, Algorithm 1). `"clean"` (deprecated) adds units
#'   whose adoption date is later than cohort g, whose post-adoption outcomes
#'   then contaminate the cohort's counterfactual.
#' @return A list of class `coresynth`.
#' @noRd

# ── Internal: covariate partial-out ──────────────────────────────────────────
# Kranz (2022), as adopted by Clarke et al. (2024, S.2.2, "projected"): beta
# comes from the two-way fixed-effects regression
#   Y_it = X_it beta + mu_i + lambda_t + e_it
# on untreated cells only, and Y - X beta is what the estimator then sees.
# The fixed effects are not optional: without them beta also absorbs how the
# covariate's unit levels and time trend line up with the outcome's, and
# subtracting X beta then writes that spurious slope back into every unit's
# series -- including a differential trend wherever treated units' covariates
# trend differently, which is the usual reason to adjust for them at all.
#
# The cells form an unbalanced set, so the within transformation is done by
# alternating projections (sweep unit means, then time means, until both
# vanish); on a balanced block it is exact after one pass. By Frisch-Waugh-
# Lovell, OLS on the demeaned columns is the fixed-effects slope.
.twfe_partial_out <- function(Y, X_arr, cells, label) {
  p   <- dim(X_arr)[3]
  K   <- nrow(cells)
  y_v <- Y[cells]
  X_m <- matrix(NA_real_, K, p)
  for (k in seq_len(p))
    X_m[, k] <- X_arr[cbind(cells, rep(k, K))]

  ok <- stats::complete.cases(X_m) & is.finite(y_v)
  if (sum(ok) < p + 1L)
    stop("Too few complete cases for covariate partial-out in ", label, ".",
         call. = FALSE)
  X_m  <- X_m[ok, , drop = FALSE]
  y_v  <- y_v[ok]
  f_t  <- factor(cells[ok, 1L])
  f_i  <- factor(cells[ok, 2L])

  demean <- function(v) {
    scale <- max(1, max(abs(v)))
    for (it in seq_len(10000L)) {
      v   <- v - (rowsum(v, f_i, reorder = TRUE) / tabulate(f_i))[f_i]
      m_t <- (rowsum(v, f_t, reorder = TRUE) / tabulate(f_t))[f_t]
      v   <- v - m_t
      if (max(abs(m_t)) < 1e-12 * scale) break
    }
    v
  }
  y_w <- demean(y_v)
  X_w <- apply(X_m, 2L, demean)
  if (!is.matrix(X_w)) X_w <- matrix(X_w, ncol = p)

  # A covariate with no variation left after the fixed effects (time-
  # invariant, or a common time path) is absorbed by them and has no slope to
  # estimate. The ridge below would return ~0 for it anyway; say so.
  var_raw <- apply(X_m, 2L, function(v) sum((v - mean(v))^2))
  var_w   <- colSums(X_w^2)
  absorbed <- var_w <= 1e-10 * pmax(var_raw, .Machine$double.eps)
  if (any(absorbed)) {
    cov_names <- dimnames(X_arr)[[3L]] %||% paste0("covariate ", seq_len(p))
    warning("Covariate(s) ", paste(cov_names[absorbed], collapse = ", "),
            " vary only across units or only over time, so the unit and ",
            "time fixed effects of the partial-out regression absorb them ",
            "and they are not adjusted for.", call. = FALSE)
  }

  beta_hat <- numeric(p)
  keep <- !absorbed
  if (any(keep))
    beta_hat[keep] <- drop(solve(
      crossprod(X_w[, keep, drop = FALSE]) + diag(1e-8, sum(keep)),
      crossprod(X_w[, keep, drop = FALSE], y_w)))

  Y_tilde <- Y
  for (k in seq_len(p))
    Y_tilde <- Y_tilde - X_arr[, , k] * beta_hat[k]

  list(Y_tilde = Y_tilde, beta_hat = beta_hat)
}

# The sample is every untreated cell, as in Kranz's reference implementation
# (xsynthdid::adjust.outcome.for.x, whose default rows are treatment == 0):
#   - control units x all T periods
#   - treated unit j x periods 1..(T_adopt[j]-1)
# The same rule serves sharp and staggered panels.
# Returns Y_tilde = Y - X * beta_hat (all units, all periods).
.sdid_partial_out_staggered <- function(Y, X_arr, pan) {
  TT      <- nrow(Y)
  T_adopt <- pan$T_adopt
  co_cells <- as.matrix(expand.grid(row = seq_len(TT), col = pan$idx_control))
  tr_cells <- do.call(rbind, lapply(pan$idx_treat, function(j) {
    t0 <- T_adopt[j]
    if (is.na(t0) || t0 <= 1L) return(NULL)
    cbind(row = seq_len(t0 - 1L), col = j)
  }))
  .twfe_partial_out(Y, X_arr, rbind(co_cells, tr_cells),
                    "staggered adoption")
}

# ── Internal: SDID from pre-built matrices ────────────────────────────────────
# Core estimation logic, callable from both sharp and staggered paths.
# Colnames of Y_co_pre are used to name the returned unit_weights vector.
.fit_sdid_matrices <- function(Y_co_pre, Y_co_post,
                               Y_tr_pre_mean, Y_tr_post_mean,
                               N_tr = 1L, zeta2 = NULL, zeta_t = 1e-6) {
  T_pre  <- nrow(Y_co_pre)
  T_post <- nrow(Y_co_post)

  if (T_pre < 2L)
    stop("at least 2 pre-treatment periods required to estimate sigma_hat.",
         call. = FALSE)

  # sigma_hat^2 (Eq. 2.2)
  delta      <- diff(Y_co_pre)
  sigma2_hat <- mean((delta - mean(delta))^2)
  sigma_hat  <- sqrt(max(sigma2_hat, .Machine$double.eps))

  if (is.null(zeta2))
    zeta2 <- sqrt(N_tr * T_post) * sigma2_hat

  # Unit weights: column-demeaning concentrates out omega_0
  mu_co      <- colMeans(Y_co_pre)
  mu_tr_pre  <- mean(Y_tr_pre_mean)
  Y_co_pre_c <- sweep(Y_co_pre, 2L, mu_co, "-")
  Y_tr_c     <- Y_tr_pre_mean - mu_tr_pre

  omega  <- sdid_unit_weights_cpp(Y_co_pre_c, Y_tr_c, zeta2)
  omega0 <- mu_tr_pre - sum(omega * mu_co)

  # Time weights: row-demeaning concentrates out lambda_0
  Ybar_co_post <- colMeans(Y_co_post)
  mu_t         <- rowMeans(Y_co_pre)
  mu_post      <- mean(Ybar_co_post)
  Y_co_pre_rc  <- sweep(Y_co_pre, 1L, mu_t, "-")
  b_tilde      <- Ybar_co_post - mu_post
  zeta_t_eff   <- zeta_t * sigma_hat

  lambda  <- sdid_time_weights_cpp(Y_co_pre_rc, b_tilde, zeta_t_eff)
  lambda0 <- mu_post - sum(lambda * mu_t)

  # SDID estimate
  tau_hat <- sdid_estimate_cpp(
    Y_co_pre, Y_co_post,
    Y_tr_pre_mean, Y_tr_post_mean,
    omega, lambda
  )

  names(omega) <- colnames(Y_co_pre)

  list(
    estimate     = tau_hat,
    unit_weights = omega,
    omega0       = omega0,
    time_weights = drop(lambda),
    lambda0      = lambda0,
    zeta2        = zeta2,
    sigma2_hat   = sigma2_hat
  )
}

# ── Internal: re-estimate one SDID cohort with given control indices ──────────
# Used by sdid_inference() staggered bootstrap/jackknife.
# zeta2 = NULL: recompute from resampled data (bootstrap).
# zeta2 = cf$zeta2: fixed from original fit (jackknife).
.refit_sdid_cohort <- function(Y_mat, cf, zeta2 = NULL) {
  TT        <- nrow(Y_mat)
  pre_rows  <- seq_len(cf$T_pre)
  post_rows <- (cf$T_pre + 1L):TT
  res <- tryCatch(
    .fit_sdid_matrices(
      Y_co_pre       = Y_mat[pre_rows,  cf$idx_co, drop = FALSE],
      Y_co_post      = Y_mat[post_rows, cf$idx_co, drop = FALSE],
      Y_tr_pre_mean  = rowMeans(Y_mat[pre_rows,  cf$idx_tr, drop = FALSE]),
      Y_tr_post_mean = rowMeans(Y_mat[post_rows, cf$idx_tr, drop = FALSE]),
      N_tr  = length(cf$idx_tr),
      zeta2 = zeta2
    ),
    error = function(e) NULL
  )
  if (is.null(res)) NA_real_ else res$estimate
}

# ── Internal: cohort-by-cohort SDID for staggered adoption ───────────────────
# Arkhangelsky et al. (2021) Appendix S.8; Clarke et al. (2024).
# Y: optional pre-residualised outcome matrix. If NULL, uses pan$Y.
.fit_sdid_staggered <- function(pan, Y = NULL, zeta2 = NULL, zeta_t = 1e-6,
                                control_group = "never_treated") {
  Y       <- if (!is.null(Y)) Y else pan$Y
  T_adopt <- pan$T_adopt
  idx_tr  <- pan$idx_treat
  idx_co  <- pan$idx_control
  TT      <- nrow(Y)

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

    if (length(idx_co_g) < 2L || T_pre_g < 2L) {
      warning(sprintf(
        "SDID staggered: cohort g=%d skipped (T_pre=%d, N_co=%d).",
        g, T_pre_g, length(idx_co_g)
      ), call. = FALSE)
      return(NULL)
    }

    pre_rows  <- seq_len(T_pre_g)
    post_rows <- (T_pre_g + 1L):TT

    fit <- tryCatch(
      .fit_sdid_matrices(
        Y_co_pre       = Y[pre_rows,  idx_co_g, drop = FALSE],
        Y_co_post      = Y[post_rows, idx_co_g, drop = FALSE],
        Y_tr_pre_mean  = rowMeans(Y[pre_rows,  idx_tr_g, drop = FALSE]),
        Y_tr_post_mean = rowMeans(Y[post_rows, idx_tr_g, drop = FALSE]),
        N_tr   = length(idx_tr_g),
        zeta2  = zeta2,
        zeta_t = zeta_t
      ),
      error = function(e) {
        warning(sprintf("SDID staggered: cohort g=%d failed: %s",
                        g, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(fit)) return(NULL)

    fit$cohort    <- g
    fit$n_treated <- length(idx_tr_g)
    fit$T_pre     <- T_pre_g
    fit$T_post    <- T_post_g
    fit$weight    <- as.numeric(length(idx_tr_g)) * T_post_g
    fit$idx_co    <- idx_co_g   # column indices in Y_mat (for inference resampling)
    fit$idx_tr    <- idx_tr_g   # column indices in Y_mat
    fit
  })

  valid <- !vapply(cohort_list, is.null, logical(1L))
  if (!any(valid)) {
    stop("All cohort-level SDID fits failed; see the preceding warnings for ",
         "per-cohort reasons.", call. = FALSE)
  }
  .aggregate_cohorts(cohort_list[valid])
}

# ── Public wrapper ─────────────────────────────────────────────────────────────
fit_sdid_cpp <- function(y, d, id, time,
                         zeta2 = NULL, zeta_t = 1e-6,
                         covariates = NULL,
                         data = NULL, id_var = NULL, time_var = NULL,
                         control_group = c("never_treated", "clean"),
                         ...) {
  control_group <- match.arg(control_group)
  pan    <- panel_to_matrices(y, d, id, time)
  .check_panel_complete(pan$Y, "SDID")
  .check_control_group(control_group, pan, "SDID")
  Y      <- pan$Y
  T_pre  <- pan$T_pre
  idx_tr <- pan$idx_treat
  idx_co <- pan$idx_control

  use_cov <- !is.null(covariates) && length(covariates) > 0L && !is.null(data)

  # ── Staggered path (Arkhangelsky 2021 Appendix S.8) ──────────────────────
  if (!pan$is_sharp) {
    Y_work   <- pan$Y
    beta_hat <- numeric(0)
    if (use_cov) {
      X_arr    <- build_covariate_array(data, id_var, time_var,
                                       covariates, pan$units, pan$times)
      po       <- .sdid_partial_out_staggered(pan$Y, X_arr, pan)
      Y_work   <- po$Y_tilde
      beta_hat <- po$beta_hat
    }
    stag <- .fit_sdid_staggered(pan,
                                Y             = Y_work,
                                zeta2         = zeta2,
                                zeta_t        = zeta_t,
                                control_group = control_group)
    return(list(
      method           = "sdid",
      staggered        = TRUE,
      control_group    = control_group,
      estimate         = stag$estimate,
      cohort_estimates = stag$cohort_estimates,
      cohort_fits      = stag$cohort_fits,
      # On the partialled-out scale, which is the one every other series in
      # the fit is on -- the cohort fits, Y_mat, and the sharp path's own
      # Y_treat. Identical to the raw outcome when no covariates are given.
      Y_treat          = Y_work[, pan$idx_treat, drop = FALSE],
      Y_mat            = Y_work,
      Y_synth          = NULL,
      gap              = NULL,
      unit_weights     = NULL,
      time_weights     = NULL,
      times            = pan$times,
      T_pre            = pan$T_pre,
      beta_hat         = beta_hat,
      covariates       = covariates
    ))
  }

  # ── Sharp path ───────────────────────────────────────────────────────────
  if (length(idx_co) < 2)
    stop("SDID requires at least two control units.", call. = FALSE)

  TT     <- nrow(Y)
  T_post <- TT - T_pre
  N_tr   <- length(idx_tr)

  if (T_post < 1)
    stop("No post-treatment periods found.", call. = FALSE)

  # Partial-out covariates (Arkhangelsky 2021 S.4; Clarke et al. 2024)
  if (use_cov) {
    X_arr    <- build_covariate_array(data, id_var, time_var,
                                     covariates, pan$units, pan$times)
    po       <- .sdid_partial_out_staggered(Y, X_arr, pan)
    Y_work   <- po$Y_tilde
    beta_hat <- po$beta_hat
  } else {
    X_arr    <- NULL
    Y_work   <- Y
    beta_hat <- numeric(0)
  }

  pre_rows  <- seq_len(T_pre)
  post_rows <- (T_pre + 1L):TT

  mfit <- .fit_sdid_matrices(
    Y_co_pre       = Y_work[pre_rows,  idx_co, drop = FALSE],
    Y_co_post      = Y_work[post_rows, idx_co, drop = FALSE],
    Y_tr_pre_mean  = rowMeans(Y_work[pre_rows,  idx_tr, drop = FALSE]),
    Y_tr_post_mean = rowMeans(Y_work[post_rows, idx_tr, drop = FALSE]),
    N_tr   = N_tr,
    zeta2  = zeta2,
    zeta_t = zeta_t
  )

  omega   <- mfit$unit_weights
  Y_synth <- drop(Y_work[, idx_co, drop = FALSE] %*% omega)
  Y_treat <- rowMeans(Y_work[, idx_tr, drop = FALSE])

  list(
    method          = "sdid",
    staggered       = FALSE,
    unit_weights    = omega,
    omega0          = mfit$omega0,
    time_weights    = mfit$time_weights,
    lambda0         = mfit$lambda0,
    estimate        = mfit$estimate,
    zeta2           = mfit$zeta2,
    zeta2_user      = zeta2,
    zeta_t          = zeta_t,
    sigma2_hat      = mfit$sigma2_hat,
    Y_synth         = Y_synth,
    Y_treat         = Y_treat,
    gap             = Y_treat - Y_synth,
    times           = pan$times,
    T_pre           = T_pre,
    covariates      = covariates,
    beta_hat        = beta_hat,
    Y_co_pre        = Y_work[pre_rows,  idx_co, drop = FALSE],
    Y_co_post       = Y_work[post_rows, idx_co, drop = FALSE],
    Y_tr_pre        = Y_work[pre_rows,  idx_tr, drop = FALSE],
    Y_tr_post       = Y_work[post_rows, idx_tr, drop = FALSE],
    Y_tr_pre_mean   = rowMeans(Y_work[pre_rows,  idx_tr, drop = FALSE]),
    Y_tr_post_mean  = rowMeans(Y_work[post_rows, idx_tr, drop = FALSE]),
    Y_raw           = if (use_cov) Y,
    X_arr           = X_arr,
    units           = pan$units,
    id_var          = id_var,
    idx_tr          = idx_tr,
    idx_co          = idx_co,
    N_tr            = N_tr,
    n_treated       = N_tr
  )
}

# ── SDID Inference (Clarke et al. 2024, Algorithms 2-4) ────────────────────────

#' Inference for Synthetic Difference-in-Differences
#'
#' Computes standard errors and p-values for a SDID estimate using one of
#' four methods: permutation placebo test (Clarke et al. 2024, Algorithm 4),
#' control-unit bootstrap, leave-one-control-out jackknife (Algorithm 3), or
#' (staggered fits only) a global jackknife across the unique control units of
#' all cohorts.
#'
#' **What `method = "bootstrap"` resamples:** `sdid_inference()` resamples
#' **control (donor) units only** with replacement, holding the treated mean
#' series (`Y_tr_pre_mean`, `Y_tr_post_mean`) and the regularisation penalty
#' `zeta2` fixed at their original-fit values. To resample **both treated and
#' control units** (at the unit or cluster level, e.g. prefectures or states),
#' re-estimating \eqn{\hat\sigma^2}, \eqn{\zeta^2}, \eqn{\hat\omega}, and \eqn{\hat\lambda} on
#' every draw and obtaining period-by-period event-study confidence intervals
#' alongside the ATT standard error, use [sdid_event_study()] instead.
#'
#' For `method = "placebo"`, the p-value is the permutation p-value, while
#' the standard error is the dispersion of the placebo distribution
#' (Clarke et al. 2024, Algorithm 4) and the confidence interval is the
#' normal approximation around the estimate with that SE. The placebo SE
#' assumes the treated unit's noise is comparable to the control units';
#' interpret it with caution when the donor pool is small.
#'
#' @param fit A `coresynth` object with `method = "sdid"` (sharp or staggered).
#' @param method Inference method: `"placebo"` (permutation), `"bootstrap"`
#'   (resamples control units only; see Details and [sdid_event_study()] for
#'   unit- or cluster-level resampling of all units), `"jackknife"`, or
#'   `"jackknife_global"` (staggered only).
#' @param n_boot Number of bootstrap replications (only for `method = "bootstrap"`).
#' @param level Confidence level for the interval (all methods).
#' @param alternative Direction of the alternative hypothesis: `"two.sided"`,
#'   `"greater"`, or `"less"`.
#' @param seed Integer seed for reproducibility (only for `method = "bootstrap"`).
#' @param ... Not used. Passing `cluster` raises an error that points to
#'   [sdid_event_study()]; any other argument is reported as unused.
#'
#' @return A list with:
#' * `estimate`: The SDID point estimate.
#' * `se`: Standard error (placebo: placebo-distribution SD, Algorithm 4).
#' * `p_value`: Permutation or normal-approximation p-value.
#' * `ci_lower`, `ci_upper`: Confidence interval bounds.
#' * `method`: The inference method used.
#' * `n_controls`: Number of control units.
#' * `alternative`: The alternative hypothesis direction.
#' * `placebo_effects`: Named vector of LOO placebo effects (placebo only).
#' * `boot_ests`: Bootstrap estimate distribution (bootstrap only).
#'
#' @seealso [sdid_event_study()] for period-by-period event-study effects and
#'   cluster (or unit) bootstrap resampling over both treated and control units.
#' @export
sdid_inference <- function(
  fit,
  method      = c("placebo", "bootstrap", "jackknife", "jackknife_global"),
  n_boot      = 200L,
  level       = 0.95,
  alternative = c("two.sided", "greater", "less"),
  seed        = NULL,
  ...
) {
  dots <- list(...)
  if ("cluster" %in% names(dots)) {
    if (!identical(tryCatch(match.arg(method), error = function(e) NULL),
                   "bootstrap")) {
      stop("`cluster` applies only to bootstrap resampling, not method = \"",
           method[1L], "\". Use sdid_event_study(fit, cluster = ...) for ",
           "cluster bootstrap inference.", call. = FALSE)
    }
    stop("sdid_inference() does not accept `cluster` (its bootstrap resamples ",
         "control units only). Use sdid_event_study(fit, cluster = ...) to ",
         "resample clusters of treated and control units.", call. = FALSE)
  }
  if (length(dots) > 0L) {
    stop("unused argument(s) in sdid_inference(): ",
         paste(names(dots), collapse = ", "), ".", call. = FALSE)
  }
  method      <- match.arg(method)
  alternative <- match.arg(alternative)

  # ── Input validation ────────────────────────────────────────────────────────
  if (!inherits(fit, "coresynth_sdid"))
    stop("sdid_inference() requires a coresynth fit with method = 'sdid'.",
         call. = FALSE)
  tau_hat   <- fit$estimate
  alpha     <- 1 - level
  staggered <- inherits(fit, "coresynth_staggered")

  # ── Staggered path ────────────────────────────────────────────────────────────
  if (staggered) {
    if (is.null(fit$Y_mat))
      stop("fit$Y_mat not found. Re-estimate with the current version of coresynth.",
           call. = FALSE)

    cohort_fits <- fit$cohort_fits
    W    <- sum(vapply(cohort_fits, `[[`, numeric(1L), "weight"))
    w_g  <- vapply(cohort_fits, `[[`, numeric(1L), "weight")
    Y_mat <- fit$Y_mat

    if (!is.null(seed)) set.seed(seed)

    if (method == "placebo") {
      # never-treated controls: intersection of all cohort control pools.
      # future-treated units appear in early cohorts only; the intersection
      # yields the set that is untreated throughout all cohorts.
      never_co <- Reduce(intersect, lapply(cohort_fits, `[[`, "idx_co"))
      if (length(never_co) < 2L)
        stop("sdid_inference() staggered placebo: fewer than 2 never-treated ",
             "controls available. Use method='bootstrap' or 'jackknife'.",
             call. = FALSE)

      TT <- nrow(Y_mat)

      # Call sdid_placebo_cpp once per cohort (K calls total).
      # lambda (time_weights) and zeta2 are fixed from the original cohort fit.
      plac_by_cohort <- lapply(cohort_fits, function(cf) {
        pre_rows  <- seq_len(cf$T_pre)
        post_rows <- (cf$T_pre + 1L):TT
        sdid_placebo_cpp(
          Y_mat[pre_rows,  cf$idx_co, drop = FALSE],
          Y_mat[post_rows, cf$idx_co, drop = FALSE],
          cf$time_weights,
          cf$zeta2
        )
      })

      # For each never-treated unit j, aggregate per-cohort placebos with
      # the same cohort weights used for the staggered ATT.
      placebo_effects <- vapply(never_co, function(j) {
        att_total <- 0
        for (k in seq_along(cohort_fits)) {
          cf    <- cohort_fits[[k]]
          pos_j <- which(cf$idx_co == j)   # length 1: j always in idx_co_k
          att_total <- att_total + cf$weight * plac_by_cohort[[k]][pos_j]
        }
        att_total / W
      }, numeric(1L))

      nm <- colnames(Y_mat)
      names(placebo_effects) <-
        if (!is.null(nm)) nm[never_co] else as.character(never_co)

      p_value <- switch(alternative,
        two.sided = (1 + sum(abs(placebo_effects) >= abs(tau_hat))) / (length(never_co) + 1),
        greater   = (1 + sum(placebo_effects >= tau_hat))           / (length(never_co) + 1),
        less      = (1 + sum(placebo_effects <= tau_hat))           / (length(never_co) + 1)
      )

      # Placebo variance (Clarke et al. 2024 Alg. 4): dispersion of the
      # placebo distribution, 1/B normalization
      se       <- sqrt(mean((placebo_effects - mean(placebo_effects))^2))
      ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
      ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se

      return(structure(list(
        estimate        = tau_hat,
        se              = se,
        p_value         = p_value,
        ci_lower        = ci_lower,
        ci_upper        = ci_upper,
        method          = method,
        n_controls      = length(never_co),
        alternative     = alternative,
        placebo_effects = placebo_effects,
        boot_ests       = NULL,
        staggered       = TRUE
      ), class = c("sdid_inference", "coresynth_inference")))
    }

    if (method == "bootstrap") {
      boot_ests <- replicate(n_boot, {
        att_total <- 0
        for (cf in cohort_fits) {
          N_co_g <- length(cf$idx_co)
          idx_b  <- cf$idx_co[sample(N_co_g, replace = TRUE)]
          cf_b   <- cf; cf_b$idx_co <- idx_b
          att_g  <- tryCatch(.refit_sdid_cohort(Y_mat, cf_b),
                             error = function(e) NA_real_)
          att_total <- att_total + cf$weight * att_g
        }
        att_total / W
      })
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
        estimate        = tau_hat,
        se              = se,
        p_value         = p_value,
        ci_lower        = ci_lower,
        ci_upper        = ci_upper,
        method          = method,
        n_controls      = mean(vapply(cohort_fits, function(cf) length(cf$idx_co), numeric(1L))),
        alternative     = alternative,
        placebo_effects = NULL,
        boot_ests       = boot_ests,
        staggered       = TRUE
      ), class = c("sdid_inference", "coresynth_inference")))
    }

    if (method == "jackknife_global") {
      # Global jackknife: drop one unique control across ALL cohorts simultaneously.
      # Captures cross-cohort correlation ignored by per-cohort LOO.
      all_co    <- sort(unique(unlist(lapply(cohort_fits, `[[`, "idx_co"))))
      orig_ests <- vapply(cohort_fits, `[[`, numeric(1L), "estimate")
      jack_ests <- vapply(all_co, function(i) {
        att_total <- 0
        for (k in seq_along(cohort_fits)) {
          cf   <- cohort_fits[[k]]
          if (i %in% cf$idx_co) {
            keep <- cf$idx_co[cf$idx_co != i]
            if (length(keep) < 2L) return(NA_real_)
            cf_loo        <- cf
            cf_loo$idx_co <- keep
            att_g <- .refit_sdid_cohort(Y_mat, cf_loo, zeta2 = cf$zeta2)
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
        estimate        = tau_hat,
        se              = se,
        p_value         = p_value,
        ci_lower        = ci_lower,
        ci_upper        = ci_upper,
        method          = method,
        n_controls      = length(all_co),
        alternative     = alternative,
        placebo_effects = NULL,
        boot_ests       = NULL,
        staggered       = TRUE
      ), class = c("sdid_inference", "coresynth_inference")))
    }

    # jackknife (per-cohort LOO + delta-method; Clarke et al. 2024 Alg 3)
    var_components <- vapply(cohort_fits, function(cf) {
      N_co_g <- length(cf$idx_co)
      if (N_co_g < 2L) return(0)
      jack_g <- vapply(seq_len(N_co_g), function(i) {
        cf_loo        <- cf
        cf_loo$idx_co <- cf$idx_co[-i]
        .refit_sdid_cohort(Y_mat, cf_loo, zeta2 = cf$zeta2)
      }, numeric(1L))
      valid_g <- jack_g[!is.na(jack_g)]
      if (length(valid_g) < 2L) return(0)
      (N_co_g - 1L) / N_co_g * sum((valid_g - mean(valid_g))^2)
    }, numeric(1L))
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
    return(structure(list(
      estimate        = tau_hat,
      se              = se,
      p_value         = p_value,
      ci_lower        = ci_lower,
      ci_upper        = ci_upper,
      method          = method,
      n_controls      = mean(vapply(cohort_fits, function(cf) length(cf$idx_co), numeric(1L))),
      alternative     = alternative,
      placebo_effects = NULL,
      boot_ests       = NULL,
      staggered       = TRUE
    ), class = c("sdid_inference", "coresynth_inference")))
  }

  # ── Sharp path ───────────────────────────────────────────────────────────────
  if (method == "jackknife_global")
    stop("sdid_inference() method='jackknife_global' requires a staggered fit.",
         call. = FALSE)
  if (is.null(fit$Y_co_pre))
    stop("fit does not contain Y_co_pre. Re-estimate with the current version of ",
         "coresynth to enable inference.", call. = FALSE)

  N_co    <- ncol(fit$Y_co_pre)

  # ── Method dispatch ─────────────────────────────────────────────────────────
  if (method == "placebo") {
    effects <- drop(sdid_placebo_cpp(fit$Y_co_pre, fit$Y_co_post,
                                     fit$time_weights, fit$zeta2))
    names(effects) <- colnames(fit$Y_co_pre)

    p_value <- switch(alternative,
      two.sided = (1 + sum(abs(effects) >= abs(tau_hat))) / (N_co + 1),
      greater   = (1 + sum(effects >= tau_hat))           / (N_co + 1),
      less      = (1 + sum(effects <= tau_hat))           / (N_co + 1)
    )

    # Placebo variance (Clarke et al. 2024 Alg. 4): dispersion of the
    # placebo distribution, 1/B normalization
    se       <- sqrt(mean((effects - mean(effects))^2))
    ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
    ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se

    return(structure(list(
      estimate        = tau_hat,
      se              = se,
      p_value         = p_value,
      ci_lower        = ci_lower,
      ci_upper        = ci_upper,
      method          = method,
      n_controls      = N_co,
      alternative     = alternative,
      placebo_effects = effects,
      boot_ests       = NULL
    ), class = c("sdid_inference", "coresynth_inference")))
  }

  # bootstrap / jackknife share re-estimation via .fit_sdid_matrices()
  .refit <- function(Y_co_pre, Y_co_post) {
    .fit_sdid_matrices(
      Y_co_pre       = Y_co_pre,
      Y_co_post      = Y_co_post,
      Y_tr_pre_mean  = fit$Y_tr_pre_mean,
      Y_tr_post_mean = fit$Y_tr_post_mean,
      N_tr           = fit$N_tr,
      zeta2          = fit$zeta2
    )$estimate
  }

  if (method == "bootstrap") {
    if (!is.null(seed)) set.seed(seed)
    boot_ests <- replicate(n_boot, {
      idx <- sample(N_co, replace = TRUE)
      .refit(fit$Y_co_pre[, idx, drop = FALSE],
             fit$Y_co_post[, idx, drop = FALSE])
    })
    se       <- sd(boot_ests)
    ci_lower <- unname(quantile(boot_ests, alpha / 2))
    ci_upper <- unname(quantile(boot_ests, 1 - alpha / 2))
    z        <- tau_hat / max(se, .Machine$double.eps)
    p_value  <- switch(alternative,
      two.sided = 2 * pnorm(-abs(z)),
      greater   = pnorm(-z),
      less      = pnorm(z)
    )
    return(structure(list(
      estimate   = tau_hat,
      se         = se,
      p_value    = p_value,
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      method     = method,
      n_controls = N_co,
      alternative = alternative,
      placebo_effects = NULL,
      boot_ests  = boot_ests
    ), class = c("sdid_inference", "coresynth_inference")))
  }

  # jackknife
  jack_ests <- vapply(seq_len(N_co), function(i) {
    .refit(fit$Y_co_pre[, -i, drop = FALSE],
           fit$Y_co_post[, -i, drop = FALSE])
  }, numeric(1))
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
  structure(list(
    estimate    = tau_hat,
    se          = se,
    p_value     = p_value,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    method      = method,
    n_controls  = N_co,
    alternative = alternative,
    placebo_effects = NULL,
    boot_ests   = NULL
  ), class = c("sdid_inference", "coresynth_inference"))
}

#' @export
print.sdid_inference <- function(x, digits = 4, ...) {
  cat("SDID Inference (", x$method, ", alternative = ", x$alternative, ")\n",
      sep = "")
  cat("  Estimate  :", round(x$estimate, digits), "\n")
  if (!is.null(x$se))
    cat("  SE        :", round(x$se, digits), "\n")
  cat("  p-value   :", format.pval(x$p_value, digits = digits,
                                   eps = 10^(-digits)), "\n")
  if (!is.null(x$ci_lower))
    cat("  CI        : [", round(x$ci_lower, digits), ", ",
        round(x$ci_upper, digits), "]\n", sep = "")
  cat("  Controls  :", x$n_controls, "\n")
  invisible(x)
}

# ── SDID Event Study with Cluster / Unit Bootstrap ───────────────────────────

# Resolve cluster membership for each panel unit (in fit$units order).
# Accepts ONLY:
#   - NULL: unit-level resampling (each unit is its own cluster)
#   - A single column name in `data` (requires `data` to be supplied), or a
#     one-sided formula `~ col` naming that column
#   - A named vector (atomic) keyed by unit ID
# Positional (unnamed) vectors are rejected explicitly so unit-to-cluster
# alignment never depends on row or column order.
.resolve_sdid_clusters <- function(cluster, data, fit) {
  units  <- as.character(fit$units)
  id_var <- fit$id_var

  if (is.null(cluster)) {
    return(stats::setNames(units, units))
  }

  # Allow one-sided formula `~ pref` as syntactic sugar for column name "pref"
  if (inherits(cluster, "formula")) {
    vars <- all.vars(cluster)
    if (length(vars) != 1L || length(cluster) != 2L || !is.name(cluster[[2L]]))
      stop("`cluster` formula must be a single column name, e.g. `~ pref`.",
           call. = FALSE)
    cluster <- vars[1L]
  }

  # Case 1: Named vector keyed by unit ID (must have names; positional rejected)
  if (!is.null(names(cluster))) {
    if (any(!nzchar(names(cluster))) || anyNA(names(cluster)))
      stop("`cluster` has missing or empty names. Every entry in a named ",
           "`cluster` vector must be keyed by unit ID.", call. = FALSE)
    if (anyDuplicated(names(cluster))) {
      dup_u <- names(cluster)[duplicated(names(cluster))]
      # Allow duplicates only if they map to the same cluster; error if conflicting
      for (u in unique(dup_u)) {
        vals <- unique(cluster[names(cluster) == u])
        if (length(vals) > 1L)
          stop(sprintf("Unit '%s' is assigned to multiple clusters in `cluster` (%s).",
                       u, paste(vals, collapse = ", ")), call. = FALSE)
      }
      cluster <- cluster[!duplicated(names(cluster))]
    }
    miss <- setdiff(units, names(cluster))
    if (length(miss) > 0L)
      stop(sprintf(
        "`cluster` is missing %d unit(s) present in `fit` (first: '%s').",
        length(miss), miss[1L]), call. = FALSE)
    cl_vec <- as.character(cluster[units])
    if (anyNA(cl_vec) || any(!nzchar(cl_vec)))
      stop("`cluster` contains NA or empty cluster IDs for one or more units in `fit`.",
           call. = FALSE)
    names(cl_vec) <- units
    return(cl_vec)
  }

  # Case 2: Single character string naming a column in `data`
  if (is.character(cluster) && length(cluster) == 1L) {
    if (is.null(data)) {
      stop(sprintf(
        "`cluster = \"%s\"` names a column, so `data` must also be supplied ",
        cluster),
        "(or pass `cluster` as a named vector keyed by unit ID).",
        call. = FALSE)
    }
    if (!is.data.frame(data))
      stop("`data` must be a data.frame.", call. = FALSE)
    if (!cluster %in% names(data))
      stop(sprintf("Cluster column '%s' not found in `data`.", cluster),
           call. = FALSE)
    if (is.null(id_var) || !id_var %in% names(data))
      stop(sprintf("Unit identifier column '%s' not found in `data`.",
                   id_var %||% "<unknown>"), call. = FALSE)

    u_col  <- as.character(data[[id_var]])
    cl_col <- data[[cluster]]
    keep   <- !is.na(u_col) & u_col %in% units
    u_sub  <- u_col[keep]
    c_sub  <- cl_col[keep]

    if (anyNA(c_sub))
      stop(sprintf("Cluster column '%s' contains NA values for units in `fit`.",
                   cluster), call. = FALSE)
    c_sub <- as.character(c_sub)

    # Check that each unit belongs to a single cluster across all periods
    first_idx <- match(units, u_sub)
    if (anyNA(first_idx)) {
      miss <- units[is.na(first_idx)]
      stop(sprintf("`data` is missing %d unit(s) present in `fit` (first: '%s').",
                   length(miss), miss[1L]), call. = FALSE)
    }
    cl_map <- stats::setNames(c_sub[first_idx], units)
    if (any(c_sub != cl_map[u_sub])) {
      bad_u <- unique(u_sub[c_sub != cl_map[u_sub]])[1L]
      stop(sprintf(
        "Unit '%s' is assigned to more than one cluster in column '%s'; ",
        bad_u, cluster),
        "cluster membership must be time-invariant per unit.", call. = FALSE)
    }
    # Same level order as split(data, data[[cluster]]): factor levels for a
    # factor, numeric order for a numeric column, so a given seed draws the
    # same cluster sequence as splitting `data`.
    cl_levels <- levels(factor(cl_col[keep]))
    attr(cl_map, "cl_levels") <- cl_levels
    return(cl_map)
  }

  stop("`cluster` must be either a single column name in `data` (with `data` ",
       "supplied) or a named vector keyed by unit ID. Unnamed positional ",
       "vectors are not accepted.", call. = FALSE)
}

#' Event Study with Cluster or Unit Bootstrap for SDID
#'
#' Computes period-by-period treatment effects (event-study series) for a sharp
#' Synthetic Difference-in-Differences fit, aligned on the
#' \eqn{\hat\lambda}-weighted pre-treatment gap so that the post-treatment mean
#' of the period effects equals the SDID ATT (`fit$estimate`) exactly, together
#' with bootstrap standard errors and confidence intervals (both normal
#' approximation and percentile) for every period and for the overall ATT from
#' the same bootstrap draws.
#'
#' **Resampling scheme:** Unlike [sdid_inference()] with `method = "bootstrap"`,
#' which resamples control units only while holding the treated mean and
#' \eqn{\zeta^2} fixed, `sdid_event_study()` resamples **both treated and
#' control units** with replacement (at the cluster level when `cluster` is
#' supplied, or at the unit level when `cluster = NULL`). Duplicated units in a
#' bootstrap draw enter as distinct unit columns (equivalent to assigning fresh
#' unit IDs per draw), and \eqn{\hat\sigma^2}, \eqn{\zeta^2} (unless user-fixed in
#' `scm_fit()`), unit weights \eqn{\hat\omega}, and time weights \eqn{\hat\lambda} are
#' all re-estimated on each draw. Draws in which either the resampled treated
#' group or the resampled control group has zero units (or fewer than two
#' control units) are discarded; the counts of valid and discarded draws are
#' reported in `n_boot_valid` and `n_boot_dropped`.
#'
#' Pre-treatment period coefficients are returned and plotted for visual
#' inspection of pre-trends only; no joint Wald test is performed because the
#' pre-treatment gaps are constrained by the \eqn{\hat\omega} and
#' \eqn{\hat\lambda} weights and their covariance matrix can be singular.
#'
#' @param fit A sharp `coresynth` object fitted with `method = "sdid"`.
#'   Staggered fits are not supported and raise an error.
#' @param cluster Cluster specification for block/cluster bootstrap resampling
#'   (e.g. prefecture or state):
#'   * `NULL` (default): resample individual units (both treated and control)
#'     with replacement.
#'   * A single character string naming a column in `data` (or a one-sided
#'     formula `~ col`), in which case `data` must also be supplied.
#'   * A named vector whose names are unit IDs matching the units in `fit` and
#'     whose values give each unit's cluster ID.
#'
#'   Unnamed positional vectors are rejected to prevent silent misalignment
#'   between row/column order and unit IDs.
#' @param data Optional `data.frame` containing the unit identifier column and
#'   the `cluster` column when `cluster` is given as a column name.
#' @param n_boot Integer number of bootstrap replications (default `200L`). Set
#'   `n_boot = 0L` (or `NULL`) to skip the bootstrap and return point estimates
#'   only.
#' @param level Confidence level for the intervals (default `0.95`).
#' @param ci_type Which confidence interval type (`"normal"` or `"percentile"`)
#'   to populate `conf.low` and `conf.high` with in [tidy()], [glance()], and
#'   [plot()] by default. Both normal-approximation (`ci_norm_lower`,
#'   `ci_norm_upper`) and percentile (`ci_pct_lower`, `ci_pct_upper`) intervals
#'   are always computed and stored.
#' @param seed Optional integer seed for reproducible bootstrap resampling.
#'
#' @return An object of class `"sdid_event_study"` (a named list) containing:
#' * `terms`: A `data.frame` with one row per period (`time`, `event_time`,
#'   `period` (`"pre"` / `"post"`), `estimate`, `se`, `ci_norm_lower`,
#'   `ci_norm_upper`, `ci_pct_lower`, `ci_pct_upper`, `conf.low`, `conf.high`).
#'   `mean(terms$estimate[terms$period == "post"])` equals `estimate` (ATT).
#' * `estimate`: Overall SDID ATT (`fit$estimate`).
#' * `se`: Bootstrap standard error of the ATT (`sd(boot_att)`).
#' * `p_value`: Two-sided normal-approximation p-value for the ATT.
#' * `ci_norm_lower`, `ci_norm_upper`: Normal-approximation CI for the ATT.
#' * `ci_pct_lower`, `ci_pct_upper`: Percentile CI for the ATT.
#' * `ci_lower`, `ci_upper`: ATT CI bounds matching `ci_type`.
#' * `ci_type`: Resolved default CI type (`"normal"` or `"percentile"`).
#' * `level`: Confidence level.
#' * `n_boot`: Requested number of bootstrap draws.
#' * `n_boot_valid`: Number of valid bootstrap draws completed.
#' * `n_boot_dropped`: Number of draws discarded (empty treated/control group
#'   or solver failure).
#' * `clustered`: Logical; `TRUE` when `cluster` was supplied.
#' * `n_clusters`: Number of unique clusters (or total units when unclustered).
#' * `n_treated`: Number of treated units in `fit`.
#' * `n_controls`: Number of control units in `fit`.
#' * `times`: Time index vector (`fit$times`).
#' * `T_pre`: Number of pre-treatment periods (`fit$T_pre`).
#' * `boot_gaps`: Numeric matrix of dimension `n_boot_valid x T` holding the
#'   aligned period-by-period gap series for each valid bootstrap draw
#'   (colnames are `as.character(times)`), or `NULL` when `n_boot == 0`.
#' * `boot_att`: Numeric vector of length `n_boot_valid` holding the ATT for
#'   each valid draw (identical to `rowMeans(boot_gaps[, (T_pre + 1):T])`), or
#'   `NULL` when `n_boot == 0`.
#'
#' @seealso [scm_fit()], [sdid_inference()], [plot.sdid_event_study()]
#' @examples
#' set.seed(1)
#' panel <- expand.grid(unit = 1:12, year = 2001:2010)
#' panel$pref <- paste0("pref_", (panel$unit - 1) %% 4 + 1)
#' panel$d <- as.integer(panel$unit <= 3 & panel$year >= 2007)
#' panel$y <- panel$unit + 0.5 * (panel$year - 2000) +
#'   rnorm(nrow(panel), sd = 0.3) - 1.5 * panel$d
#'
#' fit <- scm_fit(y ~ d | unit + year, data = panel, method = "sdid")
#' es  <- sdid_event_study(fit, cluster = "pref", data = panel,
#'                         n_boot = 50, seed = 42)
#' es
#' broom::tidy(es)
#' \donttest{
#' plot(es)
#' }
#' @export
sdid_event_study <- function(
  fit,
  cluster = NULL,
  data    = NULL,
  n_boot  = 200L,
  level   = 0.95,
  ci_type = c("normal", "percentile"),
  seed    = NULL
) {
  if (!inherits(fit, "coresynth_sdid"))
    stop("sdid_event_study() requires a coresynth fit with method = 'sdid'.",
         call. = FALSE)
  if (inherits(fit, "coresynth_staggered") || isTRUE(fit$staggered))
    stop("sdid_event_study() supports sharp (single-adoption-date) SDID fits ",
         "only; staggered fits are not supported.", call. = FALSE)
  if (is.null(fit$Y_co_pre) || is.null(fit$Y_tr_pre))
    stop("`fit` is missing unit-level outcome matrices (`Y_tr_pre`). ",
         "Re-estimate the model with the current version of coresynth.",
         call. = FALSE)

  ci_type <- match.arg(ci_type)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1)
    stop("`level` must be a single number in (0, 1).", call. = FALSE)

  if (is.null(n_boot)) n_boot <- 0L
  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) ||
      n_boot < 0 || n_boot != trunc(n_boot))
    stop("`n_boot` must be a non-negative integer.", call. = FALSE)
  n_boot <- as.integer(n_boot)

  # Resolve cluster assignment per unit (in fit$units order, length N_all)
  clustered <- !is.null(cluster)
  cl_map    <- .resolve_sdid_clusters(cluster, data, fit)

  T_pre     <- fit$T_pre
  times     <- fit$times
  TT        <- length(times)
  T_post    <- TT - T_pre
  pre_rows  <- seq_len(T_pre)
  post_rows <- (T_pre + 1L):TT

  # Point estimate of period-by-period lambda-aligned gaps:
  # delta_t = (Y_tr_t - Y_synth_t) - sum_{s <= T_pre} lambda_s (Y_tr_s - Y_synth_s)
  raw_gap    <- as.numeric(fit$Y_treat) - as.numeric(fit$Y_synth)
  pre_offset <- sum(fit$time_weights * raw_gap[pre_rows])
  est_gaps   <- raw_gap - pre_offset
  tau_hat    <- fit$estimate

  is_post    <- c(rep(FALSE, T_pre), rep(TRUE, T_post))
  period_chr <- ifelse(is_post, "post", "pre")
  event_time <- seq_len(TT) - (T_pre + 1L)

  N_tr       <- ncol(fit$Y_tr_pre)
  N_co       <- ncol(fit$Y_co_pre)
  cl_levels  <- attr(cl_map, "cl_levels") %||%
                (if (clustered) sort(unique(unname(cl_map))) else unique(unname(cl_map)))
  n_clusters <- length(cl_levels)

  if (n_boot == 0L) {
    terms_df <- data.frame(
      time          = times,
      event_time    = event_time,
      period        = period_chr,
      estimate      = est_gaps,
      se            = NA_real_,
      ci_norm_lower = NA_real_,
      ci_norm_upper = NA_real_,
      ci_pct_lower  = NA_real_,
      ci_pct_upper  = NA_real_,
      conf.low      = NA_real_,
      conf.high     = NA_real_,
      stringsAsFactors = FALSE
    )
    return(structure(list(
      terms          = terms_df,
      estimate       = tau_hat,
      se             = NA_real_,
      p_value        = NA_real_,
      ci_norm_lower  = NA_real_,
      ci_norm_upper  = NA_real_,
      ci_pct_lower   = NA_real_,
      ci_pct_upper   = NA_real_,
      ci_lower       = NA_real_,
      ci_upper       = NA_real_,
      ci_type        = ci_type,
      level          = level,
      n_boot         = 0L,
      n_boot_valid   = 0L,
      n_boot_dropped = 0L,
      clustered      = clustered,
      n_clusters     = n_clusters,
      n_treated      = N_tr,
      n_controls     = N_co,
      times          = times,
      T_pre          = T_pre,
      boot_gaps      = NULL,
      boot_att       = NULL
    ), class = "sdid_event_study"))
  }

  # Pre-index each cluster's treated and control unit positions for fast draws
  idx_tr_all <- fit$idx_tr
  idx_co_all <- fit$idx_co
  units_tr   <- as.character(fit$units[idx_tr_all])
  units_co   <- as.character(fit$units[idx_co_all])
  cl_tr      <- unname(cl_map[units_tr])
  cl_co      <- unname(cl_map[units_co])

  # Map each cluster k in 1..G to:
  # - local column indices in Y_co_pre / Y_co_post (1..N_co)
  # - local column indices in Y_tr_pre / Y_tr_post (1..N_tr)
  # - global column indices in Y_raw / X_arr (for the covariate path)
  co_by_cl     <- lapply(cl_levels, function(cl) which(cl_co == cl))
  tr_by_cl     <- lapply(cl_levels, function(cl) which(cl_tr == cl))
  raw_co_by_cl <- lapply(cl_levels, function(cl) idx_co_all[cl_co == cl])
  raw_tr_by_cl <- lapply(cl_levels, function(cl) idx_tr_all[cl_tr == cl])

  has_cov    <- !is.null(fit$covariates) && length(fit$covariates) > 0L &&
                !is.null(fit$X_arr) && !is.null(fit$Y_raw)
  Y_co_pre   <- fit$Y_co_pre
  Y_co_post  <- fit$Y_co_post
  Y_tr_all   <- rbind(fit$Y_tr_pre, fit$Y_tr_post)
  zeta2_user <- fit$zeta2_user
  zeta_t_val <- fit$zeta_t %||% 1e-6

  # Pre-sum treated outcomes per cluster (T x 1 vector and count per cluster)
  # for the standard no-covariate path so each draw computes Y_tr_mean in O(G)
  tr_sum_by_cl <- lapply(tr_by_cl, function(cols) {
    if (length(cols) == 0L) NULL
    else rowSums(Y_tr_all[, cols, drop = FALSE])
  })
  tr_cnt_by_cl <- vapply(tr_by_cl, length, integer(1L))
  co_cnt_by_cl <- vapply(co_by_cl, length, integer(1L))

  if (!is.null(seed)) set.seed(seed)

  boot_One <- function() {
    draw <- sample.int(n_clusters, n_clusters, replace = TRUE)
    n_tr_b <- sum(tr_cnt_by_cl[draw])
    n_co_b <- sum(co_cnt_by_cl[draw])
    # Discard draws where treated or control group is empty (or < 2 controls)
    if (n_tr_b < 1L || n_co_b < 2L) return(NULL)

    if (!has_cov) {
      co_cols <- unlist(co_by_cl[draw], use.names = FALSE)
      # Weighted sum of cluster treated sums by draw frequency
      freq    <- tabulate(draw, nbins = n_clusters)
      active  <- which(freq > 0L & tr_cnt_by_cl > 0L)
      Y_tr_sum <- numeric(TT)
      for (k in active) {
        Y_tr_sum <- Y_tr_sum + freq[k] * tr_sum_by_cl[[k]]
      }
      Y_tr_b      <- Y_tr_sum / n_tr_b
      Y_co_pre_b  <- Y_co_pre[, co_cols, drop = FALSE]
      Y_co_post_b <- Y_co_post[, co_cols, drop = FALSE]
    } else {
      # With covariates: re-run two-way FE partial-out on the resampled panel
      # so matches re-running scm_fit() on the resampled data frame.
      co_raw <- unlist(raw_co_by_cl[draw], use.names = FALSE)
      tr_raw <- unlist(raw_tr_by_cl[draw], use.names = FALSE)
      all_b  <- c(co_raw, tr_raw)
      Y_b    <- fit$Y_raw[, all_b, drop = FALSE]
      X_b    <- fit$X_arr[, all_b, , drop = FALSE]
      pan_b  <- list(
        idx_control = seq_len(n_co_b),
        idx_treat   = n_co_b + seq_len(n_tr_b),
        T_adopt     = c(rep(NA_integer_, n_co_b), rep(T_pre + 1L, n_tr_b))
      )
      po_b <- tryCatch(
        suppressWarnings(.sdid_partial_out_staggered(Y_b, X_b, pan_b)),
        error = function(e) NULL
      )
      if (is.null(po_b)) return(NULL)
      Y_work_b    <- po_b$Y_tilde
      Y_co_pre_b  <- Y_work_b[pre_rows,  pan_b$idx_control, drop = FALSE]
      Y_co_post_b <- Y_work_b[post_rows, pan_b$idx_control, drop = FALSE]
      Y_tr_b      <- rowMeans(Y_work_b[, pan_b$idx_treat, drop = FALSE])
    }

    mfit_b <- tryCatch(
      .fit_sdid_matrices(
        Y_co_pre       = Y_co_pre_b,
        Y_co_post      = Y_co_post_b,
        Y_tr_pre_mean  = Y_tr_b[pre_rows],
        Y_tr_post_mean = Y_tr_b[post_rows],
        N_tr           = n_tr_b,
        zeta2          = zeta2_user,
        zeta_t         = zeta_t_val
      ),
      error = function(e) NULL
    )
    if (is.null(mfit_b)) return(NULL)

    omega_b  <- mfit_b$unit_weights
    lambda_b <- mfit_b$time_weights
    synth_pre_b  <- drop(Y_co_pre_b  %*% omega_b)
    synth_post_b <- drop(Y_co_post_b %*% omega_b)
    raw_gap_b    <- Y_tr_b - c(synth_pre_b, synth_post_b)
    offset_b     <- sum(lambda_b * raw_gap_b[pre_rows])
    gap_b        <- raw_gap_b - offset_b
    if (length(gap_b) != TT || any(!is.finite(gap_b))) return(NULL)
    gap_b
  }

  draws <- replicate(n_boot, boot_One(), simplify = FALSE)
  valid_mask     <- !vapply(draws, is.null, logical(1L))
  n_boot_valid   <- sum(valid_mask)
  n_boot_dropped <- n_boot - n_boot_valid

  if (n_boot_valid < 2L)
    stop(sprintf(
      "Too few valid bootstrap replications (%d of %d succeeded; %d discarded).",
      n_boot_valid, n_boot, n_boot_dropped), call. = FALSE)

  boot_gaps <- do.call(rbind, draws[valid_mask])
  colnames(boot_gaps) <- as.character(times)
  boot_att  <- rowMeans(boot_gaps[, post_rows, drop = FALSE])

  alpha <- 1 - level
  z_val <- qnorm(1 - alpha / 2)

  # Period-by-period SE and CIs (both normal and percentile)
  se_t      <- apply(boot_gaps, 2L, sd)
  ci_norm_l <- est_gaps - z_val * se_t
  ci_norm_u <- est_gaps + z_val * se_t
  ci_pct_mat <- apply(boot_gaps, 2L, quantile,
                      probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  ci_pct_l  <- as.numeric(ci_pct_mat[1L, ])
  ci_pct_u  <- as.numeric(ci_pct_mat[2L, ])

  # ATT SE and CIs from the exact same bootstrap draws
  att_se     <- sd(boot_att)
  att_norm_l <- tau_hat - z_val * att_se
  att_norm_u <- tau_hat + z_val * att_se
  att_pct    <- quantile(boot_att, probs = c(alpha / 2, 1 - alpha / 2),
                         names = FALSE)
  att_pct_l  <- att_pct[1L]
  att_pct_u  <- att_pct[2L]
  att_pval   <- 2 * pnorm(-abs(tau_hat / max(att_se, .Machine$double.eps)))

  use_norm <- identical(ci_type, "normal")
  terms_df <- data.frame(
    time          = times,
    event_time    = event_time,
    period        = period_chr,
    estimate      = est_gaps,
    se            = as.numeric(se_t),
    ci_norm_lower = as.numeric(ci_norm_l),
    ci_norm_upper = as.numeric(ci_norm_u),
    ci_pct_lower  = as.numeric(ci_pct_l),
    ci_pct_upper  = as.numeric(ci_pct_u),
    conf.low      = as.numeric(if (use_norm) ci_norm_l else ci_pct_l),
    conf.high     = as.numeric(if (use_norm) ci_norm_u else ci_pct_u),
    stringsAsFactors = FALSE
  )
  rownames(terms_df) <- NULL

  structure(list(
    terms          = terms_df,
    estimate       = tau_hat,
    se             = att_se,
    p_value        = att_pval,
    ci_norm_lower  = att_norm_l,
    ci_norm_upper  = att_norm_u,
    ci_pct_lower   = att_pct_l,
    ci_pct_upper   = att_pct_u,
    ci_lower       = if (use_norm) att_norm_l else att_pct_l,
    ci_upper       = if (use_norm) att_norm_u else att_pct_u,
    ci_type        = ci_type,
    level          = level,
    n_boot         = n_boot,
    n_boot_valid   = n_boot_valid,
    n_boot_dropped = n_boot_dropped,
    clustered      = clustered,
    n_clusters     = n_clusters,
    n_treated      = N_tr,
    n_controls     = N_co,
    times          = times,
    T_pre          = T_pre,
    boot_gaps      = boot_gaps,
    boot_att       = boot_att
  ), class = "sdid_event_study")
}

#' @export
print.sdid_event_study <- function(x, digits = 4, ...) {
  boot_lbl <- if (x$n_boot == 0L) {
    "point estimates only"
  } else if (isTRUE(x$clustered)) {
    sprintf("cluster bootstrap, %d clusters", x$n_clusters)
  } else {
    sprintf("unit bootstrap, %d units", x$n_clusters)
  }
  cat(sprintf("=== SDID Event Study (%s) ===\n", boot_lbl))
  cat("  ATT Estimate :", round(x$estimate, digits), "\n")
  if (is.finite(x$se)) {
    pct_str <- sprintf("%g%%", 100 * x$level)
    cat("  ATT SE       :", round(x$se, digits), "\n")
    cat("  p-value      :", format.pval(x$p_value, digits = digits,
                                       eps = 10^(-digits)), "\n")
    cat(sprintf("  %s CI (norm): [%.*f, %.*f]\n",
                pct_str, digits, x$ci_norm_lower, digits, x$ci_norm_upper))
    cat(sprintf("  %s CI (pct) : [%.*f, %.*f]\n",
                pct_str, digits, x$ci_pct_lower, digits, x$ci_pct_upper))
    cat(sprintf("  Bootstrap    : %d valid / %d requested (%d discarded)\n",
                x$n_boot_valid, x$n_boot, x$n_boot_dropped))
  }
  cat("\nPeriod-by-period effects (lambda-aligned; post-period mean = ATT):\n")
  df <- x$terms[, c("time", "period", "estimate", "se",
                    "ci_norm_lower", "ci_norm_upper",
                    "ci_pct_lower", "ci_pct_upper"), drop = FALSE]
  for (col in c("estimate", "se", "ci_norm_lower", "ci_norm_upper",
                "ci_pct_lower", "ci_pct_upper")) {
    df[[col]] <- round(df[[col]], digits)
  }
  print(df, row.names = FALSE)
  invisible(x)
}

