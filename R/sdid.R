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
#'   `data`. When provided, covariates are partialled out of Y via OLS on the
#'   control-unit pre-period block before computing SDID weights and the ATT
#'   (Arkhangelsky et al. 2021 S.4; Clarke et al. 2023). Requires `data`,
#'   `id_var`, and `time_var` to be supplied (done automatically when called
#'   via [scm_fit()]). For staggered adoption, a global partial-out using all
#'   Wit=0 observations is applied before cohort-level SDID (Clarke et al. 2023
#'   S.2.2 "projected" method, Kranz 2022).
#' @param data        Full long-format data frame (passed by [scm_fit()]).
#' @param id_var      Name of the unit identifier column (passed by [scm_fit()]).
#' @param time_var    Name of the time identifier column (passed by [scm_fit()]).
#' @param control_group For staggered adoption: which units serve as controls
#'   for each cohort. `"clean"` (default) uses never-treated units plus units
#'   whose adoption date is strictly later than cohort g (Clarke et al. 2023).
#'   `"never_treated"` uses only units that are never treated.
#' @return A list of class `coresynth`.
#' @noRd

# ── Internal: covariate partial-out (sharp) ──────────────────────────────────
# OLS of Y on X using control units x pre-period block.
# Returns Y_tilde = Y - X * beta_hat (all units, all periods).
.sdid_partial_out <- function(Y, X_arr, idx_co, T_pre) {
  p        <- dim(X_arr)[3]
  pre_rows <- seq_len(T_pre)
  N_co     <- length(idx_co)

  X_mat <- matrix(NA_real_, T_pre * N_co, p)
  y_vec <- as.vector(Y[pre_rows, idx_co, drop = FALSE])
  for (j in seq_len(p))
    X_mat[, j] <- as.vector(X_arr[pre_rows, idx_co, j])

  beta_hat <- drop(solve(crossprod(X_mat) + diag(1e-8, p),
                         crossprod(X_mat, y_vec)))

  Y_tilde <- Y
  for (j in seq_len(p))
    Y_tilde <- Y_tilde - X_arr[, , j] * beta_hat[j]

  list(Y_tilde = Y_tilde, beta_hat = beta_hat)
}

# ── Internal: covariate partial-out (staggered) ───────────────────────────────
# OLS of Y on X using ALL Wit=0 observations:
#   - control units x all T periods
#   - treated unit j x periods 1..(T_adopt[j]-1)
# This is Clarke et al. (2023) S.2.2 "projected" approach (Kranz 2022).
# Returns Y_tilde = Y - X * beta_hat (all units, all periods).
.sdid_partial_out_staggered <- function(Y, X_arr, pan) {
  p       <- dim(X_arr)[3]
  TT      <- nrow(Y)
  idx_co  <- pan$idx_control
  idx_tr  <- pan$idx_treat
  T_adopt <- pan$T_adopt

  # Collect (row, col) pairs for Wit = 0
  co_pairs <- do.call(rbind, lapply(idx_co, function(j)
    data.frame(row = seq_len(TT), col = j, stringsAsFactors = FALSE)))

  tr_pairs <- do.call(rbind, lapply(idx_tr, function(j) {
    t0 <- T_adopt[j]
    if (is.na(t0) || t0 <= 1L) return(NULL)
    data.frame(row = seq_len(t0 - 1L), col = j, stringsAsFactors = FALSE)
  }))

  cells <- rbind(co_pairs, tr_pairs)
  K     <- nrow(cells)

  idx_mat <- cbind(cells$row, cells$col)
  y_vec   <- Y[idx_mat]

  X_mat <- matrix(NA_real_, K, p)
  for (k in seq_len(p))
    X_mat[, k] <- X_arr[cbind(cells$row, cells$col, rep(k, K))]

  # Drop rows with any NA (missing covariate values)
  ok <- complete.cases(X_mat)
  if (sum(ok) < p + 1L)
    stop("Too few complete cases for covariate partial-out in staggered SDID.",
         call. = FALSE)
  X_mat <- X_mat[ok, , drop = FALSE]
  y_vec <- y_vec[ok]

  beta_hat <- drop(solve(crossprod(X_mat) + diag(1e-8, p),
                         crossprod(X_mat, y_vec)))

  Y_tilde <- Y
  for (k in seq_len(p))
    Y_tilde <- Y_tilde - X_arr[, , k] * beta_hat[k]

  list(Y_tilde = Y_tilde, beta_hat = beta_hat)
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
# Arkhangelsky et al. (2021) Appendix S.8; Clarke et al. (2023).
# Y: optional pre-residualised outcome matrix. If NULL, uses pan$Y.
.fit_sdid_staggered <- function(pan, Y = NULL, zeta2 = NULL, zeta_t = 1e-6,
                                control_group = "clean") {
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
  if (!any(valid)) stop("All cohort-level SDID fits failed.", call. = FALSE)
  r <- cohort_list[valid]

  w   <- vapply(r, `[[`, numeric(1L), "weight")
  tau <- vapply(r, `[[`, numeric(1L), "estimate")
  att <- sum(w * tau) / sum(w)

  cohort_df <- data.frame(
    cohort    = vapply(r, `[[`, integer(1L), "cohort"),
    n_treated = vapply(r, `[[`, integer(1L), "n_treated"),
    T_pre     = vapply(r, `[[`, integer(1L), "T_pre"),
    T_post    = vapply(r, `[[`, integer(1L), "T_post"),
    estimate  = tau,
    weight    = w / sum(w),
    stringsAsFactors = FALSE
  )

  list(estimate = att, cohort_estimates = cohort_df, cohort_fits = r)
}

# ── Public wrapper ─────────────────────────────────────────────────────────────
fit_sdid_cpp <- function(y, d, id, time,
                         zeta2 = NULL, zeta_t = 1e-6,
                         covariates = NULL,
                         data = NULL, id_var = NULL, time_var = NULL,
                         control_group = c("clean", "never_treated"),
                         ...) {
  control_group <- match.arg(control_group)
  pan    <- panel_to_matrices(y, d, id, time)
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
      Y_treat          = pan$Y[, pan$idx_treat, drop = FALSE],
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

  # Partial-out covariates (Arkhangelsky 2021 S.4; Clarke et al. 2023)
  if (use_cov) {
    X_arr    <- build_covariate_array(data, id_var, time_var,
                                     covariates, pan$units, pan$times)
    po       <- .sdid_partial_out(Y, X_arr, idx_co, T_pre)
    Y_work   <- po$Y_tilde
    beta_hat <- po$beta_hat
  } else {
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
    method         = "sdid",
    staggered      = FALSE,
    unit_weights   = omega,
    omega0         = mfit$omega0,
    time_weights   = mfit$time_weights,
    lambda0        = mfit$lambda0,
    estimate       = mfit$estimate,
    zeta2          = mfit$zeta2,
    sigma2_hat     = mfit$sigma2_hat,
    Y_synth        = Y_synth,
    Y_treat        = Y_treat,
    gap            = Y_treat - Y_synth,
    times          = pan$times,
    T_pre          = T_pre,
    covariates     = covariates,
    beta_hat       = beta_hat,
    Y_co_pre       = Y_work[pre_rows,  idx_co, drop = FALSE],
    Y_co_post      = Y_work[post_rows, idx_co, drop = FALSE],
    Y_tr_pre_mean  = rowMeans(Y_work[pre_rows,  idx_tr, drop = FALSE]),
    Y_tr_post_mean = rowMeans(Y_work[post_rows, idx_tr, drop = FALSE]),
    N_tr           = N_tr
  )
}

# ── SDID Inference (Clarke et al. 2023, Algorithms 2-4) ────────────────────────

#' Inference for Synthetic Difference-in-Differences
#'
#' Computes standard errors and p-values for a SDID estimate using one of
#' three methods: permutation placebo test (Algorithm 4), cluster bootstrap
#' (Algorithm 2), or leave-one-out jackknife (Algorithm 3), following
#' Clarke et al. (2023).
#'
#' @param fit A `coresynth` object with `method = "sdid"` (sharp adoption only).
#' @param method Inference method: `"placebo"` (permutation), `"bootstrap"`, or
#'   `"jackknife"`.
#' @param n_boot Number of bootstrap replications (only for `method = "bootstrap"`).
#' @param level Confidence level for the interval (only for `method = "bootstrap"`
#'   or `"jackknife"`).
#' @param alternative Direction of the alternative hypothesis: `"two.sided"`,
#'   `"greater"`, or `"less"`.
#' @param seed Integer seed for reproducibility (only for `method = "bootstrap"`).
#'
#' @return A list with:
#' * `estimate`: The SDID point estimate.
#' * `se`: Standard error (bootstrap / jackknife only).
#' * `p_value`: Permutation or normal-approximation p-value.
#' * `ci_lower`, `ci_upper`: Confidence interval bounds (bootstrap / jackknife).
#' * `method`: The inference method used.
#' * `n_controls`: Number of control units.
#' * `alternative`: The alternative hypothesis direction.
#' * `placebo_effects`: Named vector of LOO placebo effects (placebo only).
#' * `boot_ests`: Bootstrap estimate distribution (bootstrap only).
#'
#' @export
sdid_inference <- function(
  fit,
  method      = c("placebo", "bootstrap", "jackknife", "jackknife_global"),
  n_boot      = 200L,
  level       = 0.95,
  alternative = c("two.sided", "greater", "less"),
  seed        = NULL
) {
  method      <- match.arg(method)
  alternative <- match.arg(alternative)

  # ── Input validation ────────────────────────────────────────────────────────
  if (!inherits(fit, "coresynth") || !identical(fit$method, "sdid"))
    stop("sdid_inference() requires a coresynth fit with method = 'sdid'.",
         call. = FALSE)
  tau_hat   <- fit$estimate
  alpha     <- 1 - level
  staggered <- isTRUE(fit$staggered)

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

      return(structure(list(
        estimate        = tau_hat,
        se              = NULL,
        p_value         = p_value,
        ci_lower        = NULL,
        ci_upper        = NULL,
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

    # jackknife (per-cohort LOO + delta-method; Clarke et al. 2023 Alg 3)
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

    return(structure(list(
      estimate        = tau_hat,
      se              = NULL,
      p_value         = p_value,
      ci_lower        = NULL,
      ci_upper        = NULL,
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
  cat("  p-value   :", round(x$p_value, digits), "\n")
  if (!is.null(x$ci_lower))
    cat("  CI        : [", round(x$ci_lower, digits), ", ",
        round(x$ci_upper, digits), "]\n", sep = "")
  cat("  Controls  :", x$n_controls, "\n")
  invisible(x)
}
