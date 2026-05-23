#' Fit Time-Aware Synthetic Control (Rho et al. 2026)
#'
#' Uses a state-space model with Kalman Filter + RTS Smoother to estimate
#' counterfactual outcomes accounting for temporal autocorrelation.
#'
#' Model:
#'   z_{t+1} = A z_t + C + eta_t,  eta_t ~ N(0, Q)
#'   y_t     = W z_t + eps_t,       eps_t ~ N(0, R)
#'
#' Parameters (W, A, C, Q, R) are estimated via the Shumway-Stoffer exact EM
#' algorithm, initialised from PCA of control units.
#'
#' @param y      Outcome vector (long format)
#' @param d      Treatment indicator (long format)
#' @param id     Unit identifier (long format)
#' @param time   Time identifier (long format)
#' @param r      Number of latent state dimensions (default 2)
#' @param em_iter Number of EM iterations (default 20)
#' @param fix_A  Logical. If TRUE, keep A = I_r throughout EM (random-walk dynamics).
#'   Default FALSE learns A from data via OLS M-step.
#' @return A list of class `coresynth`.
#' @noRd
fit_tasc_cpp <- function(y, d, id, time, r = 2, em_iter = 20, fix_A = FALSE, ...) {
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

  # Initialise parameters from SVD of control-unit data (T_pre x N_co)
  # u: T_pre x r (time factors), v: N_co x r (unit loadings)
  Y_co_pre <- Y[seq_len(T_pre), idx_co, drop = FALSE]
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
                   t(F_pre) %*% Y[seq_len(T_pre), idx_tr, drop = FALSE])
    W_full[idx_tr, ] <- t(L_tr)
  }

  C  <- rep(0.0, r)
  A  <- diag(r)
  Q  <- diag(0.1, r)
  R  <- diag(var(as.vector(Y_co_pre)) * 0.5, N)
  z0 <- rep(0.0, r)
  P0 <- diag(1.0, r)

  # EM loop (E-step = Kalman smoother, M-step = Shumway-Stoffer OLS updates)
  for(em in seq_len(em_iter)) {
    ks <- kalman_smoother_cpp(t(Y_obs), W_full, A, C, Q, R, z0, P0)
    Z    <- ks$z_smooth   # r x T
    P_sm <- ks$P_smooth   # r x r x T (3D array)

    # M-step 1: W update (OLS, per-unit to handle missing observations)
    for (j in seq_len(N)) {
      t_obs <- which(!is.na(Y_obs[, j]))
      if (length(t_obs) == 0L) next
      Z_j <- Z[, t_obs, drop = FALSE]
      y_j <- Y_obs[t_obs, j]
      W_full[j, ] <- drop(t(y_j) %*% t(Z_j) %*%
                       solve(Z_j %*% t(Z_j) + diag(1e-8, r)))
    }

    # M-step 2: A update (Shumway-Stoffer closed-form OLS)
    if (!fix_A) {
      # P_cross: r x r x (T-1) cube
      # C++ slice t (0-indexed) = P_{t+2, t+1 | T} in 1-indexed notation
      P_cross <- ks$P_cross
      num_A <- matrix(0.0, r, r)
      den_A <- matrix(0.0, r, r)
      for (tt in seq_len(TT - 1)) {
        num_A <- num_A + P_cross[,, tt] + outer(Z[, tt + 1], Z[, tt])
        den_A <- den_A + P_sm[,, tt]    + outer(Z[, tt],     Z[, tt])
      }
      A <- num_A %*% solve(den_A + diag(1e-6, r))
    }

    # M-step 3: C update (after A — C = mean(z_{t+1} - A z_t))
    diffs_AC <- Z[, -1] - A %*% Z[, -TT]   # r x (T-1)
    C <- rowMeans(diffs_AC)

    # M-step 4: Q update (process noise covariance from corrected residuals)
    resid_proc <- diffs_AC - C
    Q <- resid_proc %*% t(resid_proc) / (TT - 1) + diag(1e-8, r)

    # M-step 5: R update (diagonal only, per-unit to handle missing observations)
    r_diag <- numeric(N)
    for (j in seq_len(N)) {
      t_obs <- which(!is.na(Y_obs[, j]))
      if (length(t_obs) == 0L) { r_diag[j] <- 1e-6; next }
      resid_j <- Y_obs[t_obs, j] -
                 drop(W_full[j, , drop = FALSE] %*% Z[, t_obs])
      r_diag[j] <- max(mean(resid_j^2), 1e-6)
    }
    R <- diag(r_diag)
  }

  # Final smooth
  ks_final <- kalman_smoother_cpp(t(Y_obs), W_full, A, C, Q, R, z0, P0)
  Z_smooth <- ks_final$z_smooth   # r x T

  Y_hat_all <- t(W_full %*% Z_smooth)   # T x N
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

  list(
    method       = "tasc",
    r            = r,
    A            = A,
    W            = W_full,
    C            = C,
    Q            = Q,
    R            = R,
    Z_smooth     = Z_smooth,
    Y_hat        = Y_hat_all,
    Y_treat      = Y_treat,
    gap          = gap,
    times        = pan$times,
    T_pre        = T_pre,
    idx_co       = idx_co,
    idx_tr       = idx_tr,
    unit_weights = NULL,
    estimate     = att
  )
}
