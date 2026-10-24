#' Fit Time-Aware Synthetic Control (Rho et al. 2026)
#'
#' Uses a state-space model with Kalman Filter + RTS Smoother to estimate
#' counterfactual outcomes accounting for temporal autocorrelation.
#'
#' Model:
#'   z_{t+1} = A z_t + C + eta_t,  eta_t ~ N(0, Q)
#'   y_t     = mu + W z_t + eps_t,  eps_t ~ N(0, R)
#'
#' Parameters (mu, W, A, C, Q, R, and the initial state z0, P0) are estimated
#' by the exact EM algorithm of Rho et al. (2026, Algorithm 7) with the
#' missing cells of Shumway & Stoffer (1982), initialised from PCA of control
#' units.
#'
#' @param y      Outcome vector (long format)
#' @param d      Treatment indicator (long format)
#' @param id     Unit identifier (long format)
#' @param time   Time identifier (long format)
#' @param r      Number of latent state dimensions (default 2)
#' @param em_iter Maximum number of EM iterations (default 2000)
#' @param em_tol Convergence tolerance: EM stops once the fitted values
#'   change by less than this, relative to their variation, between
#'   iterations.
#' @param fix_A  Logical. If TRUE, keep A = I_r throughout EM (random-walk dynamics).
#'   Default FALSE learns A from data via OLS M-step.
#' @param intercept Logical. Give each unit its own intercept mu_i in the
#'   observation equation (default TRUE); FALSE fits the paper's model, whose
#'   observations have no intercept.
#' @return A list of class `coresynth`.
#' @noRd
fit_tasc_cpp <- function(y, d, id, time, r = 2, em_iter = 2000L, em_tol = 1e-6,
                         fix_A = FALSE, intercept = TRUE, ...) {
  pan <- panel_to_matrices(y, d, id, time)
  Y   <- pan$Y; T_pre <- pan$T_pre
  idx_tr <- pan$idx_treat; idx_co <- pan$idx_control

  TT <- nrow(Y); N <- ncol(Y)

  # Build observation matrix Y_obs: NA for each treated unit post its own adoption
  Y_obs <- Y
  for (j in idx_tr) {
    t0 <- pan$T_adopt[j]
    if (!is.na(t0) && t0 <= TT) Y_obs[t0:TT, j] <- NA_real_
  }

  # The EM (Kalman smoother + per-unit M-steps) handles missing cells in
  # Y_obs, but the SVD/OLS initialisation below cannot: impute missing cells
  # with the column mean for the starting values only.
  impute_col_mean <- function(M) {
    if (!anyNA(M)) return(M)
    cm <- colMeans(M, na.rm = TRUE)
    cm[!is.finite(cm)] <- mean(M, na.rm = TRUE)
    idx <- which(is.na(M), arr.ind = TRUE)
    M[idx] <- cm[idx[, 2L]]
    M
  }

  # Without an intercept the level of each unit has to be carried by the
  # states, and the likelihood is nearly flat along the directions that trade
  # it between z0, the drift C and a near-unit root of A: on a panel at level
  # 100 the EM drifted along them for tens of thousands of iterations, moving
  # the effect by up to 0.4. With mu_i in the observation equation a shift of
  # the states is exactly offset by the intercepts, so the fitted values -- and
  # the counterfactual -- no longer depend on where the EM stops. It is the
  # time-invariant part of the latent state of Rho et al. (2026, S.3.1.1),
  # freed from the column space of W. The intercepts start at each unit's
  # mean over the common pre-treatment window.
  Y_pre <- impute_col_mean(Y[seq_len(T_pre), , drop = FALSE])
  mu    <- if (isTRUE(intercept)) colMeans(Y_pre) else rep(0.0, N)
  Y_pre <- sweep(Y_pre, 2L, mu)

  # Initialise parameters from SVD of control-unit data (T_pre x N_co)
  # u: T_pre x r (time factors), v: N_co x r (unit loadings)
  Y_co_pre <- Y_pre[, idx_co, drop = FALSE]
  svd_res  <- svd(Y_co_pre, nu = r, nv = r)

  # W_full: N x r observation/loading matrix
  # Control rows: right singular vectors (unit loadings, N_co x r)
  W_full <- matrix(0.0, nrow = N, ncol = r)
  W_full[idx_co, ] <- svd_res$v[, seq_len(r), drop = FALSE]

  # Treated rows: least-squares loadings from pre-treatment fit
  # F_pre = U * diag(D): T_pre x r time-factor matrix
  if(T_pre >= r) {
    F_pre <- svd_res$u %*% diag(svd_res$d[seq_len(r)], r, r)   # T_pre x r
    L_tr  <- solve(t(F_pre) %*% F_pre + diag(1e-8, r),
                   t(F_pre) %*% Y_pre[, idx_tr, drop = FALSE])
    W_full[idx_tr, ] <- t(L_tr)
  }

  par <- list(
    mu = mu,
    W  = W_full,
    A  = diag(r),
    C  = rep(0.0, r),
    Q  = diag(0.1, r),
    R  = diag(var(as.vector(Y_co_pre)) * 0.5, N),
    z0 = rep(0.0, r),
    P0 = diag(1.0, r)
  )

  # EM loop: E-step = Kalman filter + RTS smoother, M-step = exact maximiser
  # of the expected complete-data log-likelihood. Convergence is judged on
  # the fitted values mu + W z_t, the quantity the counterfactual is built
  # from, relative to their variation around each unit's mean: measured
  # against the raw norm, a panel far from zero makes every change look
  # negligible.
  fitted_prev <- NULL
  converged   <- FALSE
  n_iter      <- 0L
  for (em in seq_len(em_iter)) {
    ks <- kalman_smoother_cpp(t(sweep(Y_obs, 2L, par$mu)), par$W, par$A,
                              par$C, par$Q, par$R, par$z0, par$P0)
    fitted <- par$W %*% ks$z_smooth + par$mu
    if (!is.null(fitted_prev)) {
      rel <- sqrt(sum((fitted - fitted_prev)^2)) /
        max(sqrt(sum((fitted_prev - rowMeans(fitted_prev))^2)), 1e-12)
      if (rel < em_tol) {
        converged <- TRUE
        break
      }
    }
    fitted_prev <- fitted
    par    <- .tasc_mstep(Y_obs, ks, par, fix_A, isTRUE(intercept))
    n_iter <- em
  }
  if (!converged)
    warning(sprintf(paste0(
      "TASC: EM did not converge in %d iterations (em_tol = %g); raise ",
      "em_iter or loosen em_tol."), em_iter, em_tol), call. = FALSE)

  # Final smooth
  ks_final <- kalman_smoother_cpp(t(sweep(Y_obs, 2L, par$mu)), par$W, par$A,
                                  par$C, par$Q, par$R, par$z0, par$P0)
  Z_smooth <- ks_final$z_smooth   # r x T

  Y_hat_all <- sweep(t(par$W %*% Z_smooth), 2L, par$mu, "+")   # T x N
  Y_treat   <- Y[, idx_tr, drop = FALSE]
  Y_hat_tr  <- Y_hat_all[, idx_tr, drop = FALSE]

  gap <- Y_treat - Y_hat_tr
  # ATT: average gap over each unit's own post-treatment period
  gap_post_vals <- vapply(seq_along(idx_tr), function(k) {
    j  <- idx_tr[k]
    t0 <- pan$T_adopt[j]
    if (is.na(t0) || t0 > TT) return(NA_real_)
    mean(gap[t0:TT, k, drop = FALSE], na.rm = TRUE)
  }, numeric(1L))
  att <- mean(gap_post_vals, na.rm = TRUE)

  # Staggered adoption needs no cohort-by-cohort re-fitting here: Y_obs is
  # already NA from each treated unit's own adoption date, so the single
  # state-space fit is the staggered fit. Only the aggregation changes -- the
  # cohort ATTs combine with the N_treated,g * T_post,g weights of Clarke et
  # al. (2024), as in every other estimator, instead of weighting each treated
  # unit equally.
  stag <- if (!pan$is_sharp) .impute_cohort_fits(pan, gap, Y_treat, Y_hat_tr)

  c(list(
    method       = "tasc",
    r            = r,
    A            = par$A,
    W            = par$W,
    mu           = stats::setNames(par$mu, colnames(Y)),
    C            = par$C,
    Q            = par$Q,
    R            = par$R,
    z0           = par$z0,
    P0           = par$P0,
    intercept    = isTRUE(intercept),
    n_iter       = n_iter,
    converged    = converged,
    Z_smooth     = Z_smooth,
    Y_hat        = Y_hat_all,
    Y_treat      = Y_treat,
    gap          = gap,
    times        = pan$times,
    T_pre        = T_pre,
    idx_co       = idx_co,
    idx_tr       = idx_tr,
    unit_weights = NULL,
    estimate     = stag$estimate %||% att
  ), if (!is.null(stag)) list(
    staggered        = TRUE,
    cohort_estimates = stag$cohort_estimates,
    cohort_fits      = stag$cohort_fits
  ))
}

# ── Internal: TASC M-step (Rho et al. 2026, Algorithm 7) ─────────────────────
# Every update maximises the expected complete-data log-likelihood, so each
# one carries the smoothed covariances P_t alongside the smoothed means m_t.
# Plugging in the means alone (a "hard" EM) understates Q -- smoothed states
# vary less than the states themselves -- and the understatement feeds on
# itself until Q collapses to zero and the state equation turns
# deterministic.
#
# Missing cells (treated units after adoption) follow Shumway & Stoffer
# (1982): a unit's observation row -- its intercept and loadings -- and its
# noise variance use only the periods it is observed in. Units observed
# throughout share one second-moment sum.
#
# The transition sums run over t = 0..T-1 with the smoothed initial state x_0,
# which also gives the paper's updates of z0 and P0. The drift C is estimated
# jointly with A by regressing z_{t+1} on [z_t, 1]; Q and R are diagonal, as
# in the paper.
.tasc_mstep <- function(Y_obs, ks, par, fix_A, intercept = TRUE) {
  r  <- length(par$z0)
  TT <- nrow(Y_obs)
  N  <- ncol(Y_obs)
  m  <- ks$z_smooth                                   # r x T
  P  <- array(ks$P_smooth, c(r, r, TT))
  Pc <- array(ks$P_cross,  c(r, r, TT - 1L))
  sum_slices <- function(a) matrix(rowSums(matrix(a, r * r)), r, r)

  # Smoothed initial state, one RTS step back from t = 1
  P10 <- par$A %*% par$P0 %*% t(par$A) + par$Q
  J0  <- par$P0 %*% t(par$A) %*% solve(P10)
  m0s <- drop(par$z0 + J0 %*% (m[, 1L] - (par$A %*% par$z0 + par$C)))
  P0s <- par$P0 + J0 %*% (P[, , 1L] - P10) %*% t(J0)
  Pc10 <- P[, , 1L] %*% t(J0)                         # Cov(z_1, z_0 | Y)

  # Second moments of the transitions z_t -> z_{t+1}, t = 0..T-1
  P_sum   <- sum_slices(P)
  prev    <- cbind(m0s, m[, -TT, drop = FALSE])       # z_0..z_{T-1}
  E_prev  <- P0s + P_sum - P[, , TT] + prev %*% t(prev)
  E_cross <- Pc10 + sum_slices(Pc) + m %*% t(prev)
  E_next  <- P_sum + m %*% t(m)
  s_prev  <- rowSums(prev)
  s_next  <- rowSums(m)

  if (fix_A) {
    A <- diag(r)
    C <- (s_next - s_prev) / TT
  } else {
    S00 <- rbind(cbind(E_prev, s_prev), c(s_prev, TT))
    S10 <- cbind(E_cross, s_next)
    B   <- S10 %*% solve(S00 + diag(1e-10, r + 1L))
    A   <- B[, seq_len(r), drop = FALSE]
    C   <- B[, r + 1L]
  }
  # Q = E[(z_{t+1} - A z_t - C)(.)'] averaged over the T transitions
  Q_full <- (E_next - A %*% t(E_cross) - E_cross %*% t(A) +
               A %*% E_prev %*% t(A) -
               outer(C, s_next) - outer(s_next, C) +
               A %*% outer(s_prev, C) + outer(C, s_prev) %*% t(A)) / TT +
    outer(C, C)
  Q <- diag(pmax(diag(Q_full), 1e-8), r)

  # Observation rows, per unit over its observed periods: regress y_jt on
  # [m_t, 1] (or on m_t alone without an intercept) with the smoothed second
  # moments of z_t.
  obs  <- !is.na(Y_obs)
  Yz   <- Y_obs; Yz[!obs] <- 0
  k    <- r + intercept
  S_yu <- t(Yz) %*% t(m)                              # N x r: sum y m'
  if (intercept) S_yu <- cbind(S_yu, colSums(Yz))     # ... and sum y
  E_uu <- function(t_idx, P_part) {
    m_t <- m[, t_idx, drop = FALSE]
    E   <- P_part + m_t %*% t(m_t)
    if (intercept) {
      s <- rowSums(m_t)
      E <- rbind(cbind(E, s), c(s, length(t_idx)))
    }
    E + diag(1e-8, k)
  }
  W     <- par$W
  mu    <- par$mu
  r_var <- numeric(N)
  full  <- colSums(!obs) == 0L
  if (any(full)) {
    Bf <- S_yu[full, , drop = FALSE] %*% solve(E_uu(seq_len(TT), P_sum))
    W[full, ] <- Bf[, seq_len(r), drop = FALSE]
    if (intercept) mu[full] <- Bf[, k]
    res_f <- sweep(Yz[, full, drop = FALSE] - t(W[full, , drop = FALSE] %*% m),
                   2L, mu[full])
    r_var[full] <- (colSums(res_f^2) +
      rowSums((W[full, , drop = FALSE] %*% P_sum) * W[full, , drop = FALSE])) / TT
  }
  for (j in which(!full)) {
    t_obs <- which(obs[, j])
    if (length(t_obs) == 0L) { r_var[j] <- 1e-6; next }
    P_j <- sum_slices(P[, , t_obs, drop = FALSE])
    b_j <- drop(S_yu[j, ] %*% solve(E_uu(t_obs, P_j)))
    w_j <- b_j[seq_len(r)]
    W[j, ] <- w_j
    if (intercept) mu[j] <- b_j[k]
    res_j <- Y_obs[t_obs, j] - mu[j] - drop(w_j %*% m[, t_obs, drop = FALSE])
    r_var[j] <- (sum(res_j^2) + drop(w_j %*% P_j %*% w_j)) / length(t_obs)
  }

  list(
    mu = mu,
    W  = W,
    A  = A,
    C  = C,
    Q  = Q,
    R  = diag(pmax(r_var, 1e-6), N),
    z0 = m0s,
    P0 = (P0s + t(P0s)) / 2 + diag(1e-8, r)
  )
}
