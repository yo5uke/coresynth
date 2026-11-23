# ── Unit-Level Conformal Inference (Liu & Xu 2026, arXiv:2611.19209) ────────
#
# Finite-sample period-wise confidence sets, overall ATT intervals, and
# simultaneous bands built from across-unit exchangeability of prediction
# errors rather than temporal exchangeability of residuals.
#
# Two procedures:
#   1. Full Conformal (target-date estimators with affine dependence on the
#      treated outcome: outcomes-only SCM and SI), solved analytically in closed
#      form (SM A1.2, A2.6).
#   2. Jackknife+ (all six estimators, holding out {0, j} together per fold:
#      target-date SCM / predictor SCM / SDID / SI, and full-path GSC / MC /
#      TASC) with interval-dual p-values (Proposition 5) and simultaneous bands
#      (Theorem 4).

# ── Calibration helpers (Liu & Xu 2026, Eq. A24, A28-A29, Table 4) ───────────

# Largest integer q in 0..floor(N/2) such that (2q - 1)_+ <= gamma * N.
# The +1e-9 guard keeps IEEE-754 rounding in 1 - level (e.g. 1 - 0.90 =
# 0.099999999999999978) from pushing integer gamma * N down by one step.
.uc_q_gamma <- function(N, gamma) {
  gn <- gamma * N + 1e-9
  if (gn < 1) return(0L)
  min(as.integer(floor((gn + 1) / 2)), as.integer(floor(N / 2)))
}

.uc_guaranteed_coverage <- function(N, gamma, procedure) {
  if (identical(procedure, "full_conformal")) {
    k <- as.integer(floor(gamma * N + 1e-9))
    1 - k / N
  } else {
    q <- .uc_q_gamma(N, gamma)
    1 - max(2L * q - 1L, 0L) / N
  }
}

# Interval-dual p-value for Jackknife+ (Liu & Xu 2026, Eq. A28-A29).
# Reproduces the rejection decision of the Jackknife+ interval at level 1 - gamma:
# c is outside C_hat(gamma) iff p <= gamma (Proposition 5).
.uc_jackknife_pval <- function(y0, a, b, c = 0) {
  N <- length(a) + 1L
  L <- sum(y0 - c < a)
  U <- sum(y0 - c > b)
  m <- max(L, U)
  min(1, (2 * (N - m) - 1) / N)
}

# ── Closed-form Full Conformal quadratic solver & interval-cover (SM A1.2, A2.6) ─

# Solve { c in R : wj * |vj + kj * c| >= w0 * |tau0 - c| } as a matrix with
# columns (lo, hi) representing 1 or 2 closed intervals (may include -Inf/Inf).
.uc_quad_set <- function(wj, vj, kj, w0, tau0) {
  p1  <- wj * kj
  q1  <- wj * vj
  p0  <- -w0
  q0  <- w0 * tau0

  # Special case: leading coefficient (p1^2 - p0^2) == 0 => linear inequality
  if (abs(abs(p1) - w0) <= 1e-12 * w0) {
    # Boundary lines: (p1 - p0) c + (q1 - q0) = 0 and (p1 + p0) c + (q1 + q0) = 0
    b_lin   <- 2 * (p1 * q1 - p0 * q0)
    d_lin   <- q1^2 - q0^2
    scale_b <- 2 * w0 * max(abs(q1), abs(q0), 1e-30)
    scale_d <- max(q1^2, q0^2, 1e-30)
    if (abs(b_lin) <= 1e-12 * scale_b) {
      return(if (d_lin >= -1e-12 * scale_d) matrix(c(-Inf, Inf), nrow = 1L)
             else matrix(numeric(0L), ncol = 2L))
    }
    root <- -d_lin / b_lin
    if (b_lin > 0) return(matrix(c(root, Inf), nrow = 1L))
    return(matrix(c(-Inf, root), nrow = 1L))
  }

  # Non-degenerate quadratic: factors into two exact linear equations
  # r1 from (p1 - p0) c + (q1 - q0) = 0; r2 from (p1 + p0) c + (q1 + q0) = 0
  r1 <- -(q1 - q0) / (p1 - p0)
  r2 <- -(q1 + q0) / (p1 + p0)
  lo <- min(r1, r2)
  hi <- max(r1, r2)

  if (abs(p1) < w0) {
    # Downward parabola: >= 0 inside [lo, hi]
    matrix(c(lo, hi), nrow = 1L)
  } else {
    # Upward parabola: >= 0 outside (lo, hi), i.e. (-Inf, lo] U [hi, Inf)
    if (hi - lo <= 1e-12 * max(1, abs(lo), abs(hi))) {
      matrix(c(-Inf, Inf), nrow = 1L)
    } else {
      rbind(c(-Inf, lo), c(hi, Inf))
    }
  }
}

# Union of two 2-column closed-interval matrices (rows are [lo, hi]).
.uc_union_sets <- function(A, B) {
  if (nrow(A) == 0L) return(B)
  if (nrow(B) == 0L) return(A)
  U <- rbind(A, B)
  U <- U[order(U[, 1L], U[, 2L]), , drop = FALSE]
  out <- U[1L, , drop = FALSE]
  if (nrow(U) > 1L) {
    for (r in 2:nrow(U)) {
      cur_hi <- out[nrow(out), 2L]
      if (U[r, 1L] <= cur_hi) {
        out[nrow(out), 2L] <- max(cur_hi, U[r, 2L])
      } else {
        out <- rbind(out, U[r, , drop = FALSE])
      }
    }
  }
  out
}

# Given a list of J interval-matrices (each representing control j's accepted
# closed set E_j), return c(lower, upper) spanning {c : #{j : c in E_j} >= k}.
.uc_interval_cover <- function(sets_list, k) {
  if (k <= 0L) return(c(-Inf, Inf))
  J <- length(sets_list)
  if (k > J) return(c(NA_real_, NA_real_))

  pts <- unlist(lapply(sets_list, function(m) {
    if (nrow(m) == 0L) numeric(0L) else as.numeric(m)
  }), use.names = FALSE)
  fin <- sort(unique(pts[is.finite(pts)]))

  in_set <- function(c_val) {
    cnt <- 0L
    for (m in sets_list) {
      if (nrow(m) == 0L) next
      hit <- FALSE
      for (r in seq_len(nrow(m))) {
        lo_r  <- m[r, 1L]
        hi_r  <- m[r, 2L]
        tol_l <- if (is.finite(lo_r)) 1e-11 * max(1, abs(lo_r)) else 0
        tol_h <- if (is.finite(hi_r)) 1e-11 * max(1, abs(hi_r)) else 0
        if (c_val >= lo_r - tol_l && c_val <= hi_r + tol_h) {
          hit <- TRUE
          break
        }
      }
      if (hit) cnt <- cnt + 1L
    }
    cnt >= k
  }

  if (length(fin) == 0L) {
    return(if (in_set(0)) c(-Inf, Inf) else c(NA_real_, NA_real_))
  }

  # Left and right open rays
  left_unb  <- in_set(fin[1L] - max(1, abs(fin[1L])))
  right_unb <- in_set(fin[length(fin)] + max(1, abs(fin[length(fin)])))

  # Because every E_j is closed, the k-cover is closed, so any finite minimum
  # or maximum is attained at one of the finite endpoints `fin`.
  if (left_unb && right_unb) return(c(-Inf, Inf))

  acc <- vapply(fin, in_set, logical(1L))
  if (!any(acc) && !left_unb && !right_unb) return(c(NA_real_, NA_real_))

  lo <- if (left_unb) -Inf else min(fin[acc])
  hi <- if (right_unb) Inf else max(fin[acc])
  c(lo, hi)
}

# ── Method-specific single-unit / fold estimators ────────────────────────────

# Replicate sharp SCM weight fitting on training controls Y_B_pre (T_pre x N_B)
# and target y_i_pre (length T_pre), honouring fit's solver, v_window, oos,
# outer optimiser, predictor matrix (when present), and fixed lambda_pen.
.uc_scm_fit_weights <- function(Y_B_pre, y_i_pre, X_B = NULL, x_i = NULL, fit) {
  use_cov   <- !is.null(X_B)
  z_rows    <- fit[["z_rows"]]
  use_wolfe <- identical(fit[["qp_solver"]], "wolfe")
  eff_outer <- fit[["v_optim_effective"]] %||%
    if (use_cov) "multistart" else "coord_descent"
  lam_pen   <- fit[["lambda_pen"]] %||% NA_real_
  t_train   <- if (use_cov) as.integer(fit[["t_train"]] %||% -1L) else -1L
  is_oos    <- !use_cov && !is.null(fit[["v_rows"]]) &&
    length(fit[["v_rows"]]) < nrow(Y_B_pre)

  if (is_oos) {
    oos <- .scm_oos_outcomes(Y_B_pre, y_i_pre, eff_outer)
    if (!is.na(lam_pen)) {
      V_star   <- drop(oos$V)
      X0_u     <- oos$X0_final
      X1_u     <- oos$X1_final
      diff_mat <- X1_u - X0_u
      d_pen    <- colSums(V_star * diff_mat^2)
      V_sqrt   <- sqrt(pmax(V_star, 0))
      X0_sc    <- sweep(X0_u, 1L, V_sqrt, `*`)
      Q_mat    <- crossprod(X0_sc)
      c_base   <- as.vector(t(X0_sc) %*% (V_sqrt * X1_u))
      w        <- drop(solve_simplex_qp(Q_mat, c_base - lam_pen / 2 * d_pen))
    } else {
      w <- drop(oos$W)
    }
    return(w)
  }

  X0_use <- if (use_cov) X_B else Y_B_pre
  X1_use <- if (use_cov) x_i else y_i_pre

  # The penalised path takes only V from this step, which the main fit
  # computes without the Wolfe inner solver.
  res <- if (identical(eff_outer, "bfgs")) {
    .scm_bfgs_outer(X0_use, X1_use, Y_B_pre, y_i_pre, t_train = t_train,
                    z_rows = z_rows)
  } else {
    scm_weights_cpp(X0_use, X1_use, Y_B_pre, y_i_pre,
                    max_iter = 100L, tol = 1e-4, t_train = t_train,
                    z_rows = z_rows,
                    multistart = identical(eff_outer, "multistart"),
                    cheap_face = !use_cov,
                    wolfe = use_wolfe && is.na(lam_pen))
  }

  if (!is.na(lam_pen)) {
    V_star   <- drop(res$V)
    diff_mat <- X1_use - X0_use
    d_pen    <- colSums(V_star * diff_mat^2)
    V_sqrt   <- sqrt(pmax(V_star, 0))
    X0_sc    <- sweep(X0_use, 1L, V_sqrt, `*`)
    Q_mat    <- crossprod(X0_sc)
    c_base   <- as.vector(t(X0_sc) %*% (V_sqrt * X1_use))
    return(drop(solve_simplex_qp(Q_mat, c_base - lam_pen / 2 * d_pen)))
  }

  drop(res$W)
}

# Time-weight ridge of the main fit, held fixed across every refit like the
# other tuning values.
.uc_sdid_zeta_t_eff <- function(fit) {
  (fit[["zeta_t"]] %||% 1e-6) *
    sqrt(max(fit[["sigma2_hat"]], .Machine$double.eps))
}

# SDID time weights targeting a single control post-date vector `b_post`
# (length N_co) from Y_co_pre (T_pre x N_co).
.uc_sdid_time_weights <- function(Y_co_pre, b_post, zeta_t_eff) {
  mu_t       <- rowMeans(Y_co_pre)
  mu_post    <- mean(b_post)
  Y_co_pre_rc <- sweep(Y_co_pre, 1L, mu_t, "-")
  b_tilde    <- b_post - mu_post
  drop(sdid_time_weights_cpp(Y_co_pre_rc, b_tilde, zeta_t_eff))
}

# SDID unit weights (omega and concentrated intercept omega0) from Y_co_pre
# targeting y_tr_pre.
.uc_sdid_unit_weights <- function(Y_co_pre, y_tr_pre, zeta2) {
  mu_co      <- colMeans(Y_co_pre)
  mu_tr      <- mean(y_tr_pre)
  Y_co_pre_c <- sweep(Y_co_pre, 2L, mu_co, "-")
  y_tr_c     <- y_tr_pre - mu_tr
  omega      <- drop(sdid_unit_weights_cpp(Y_co_pre_c, y_tr_c, zeta2))
  omega0     <- mu_tr - sum(omega * mu_co)
  list(omega = omega, omega0 = omega0)
}

# Counterfactual prediction for target-date SDID given unit weights omega and
# time weights lambda.
.uc_sdid_predict <- function(Y_co_pre, y_co_post, y_tr_pre, omega, lambda) {
  y_tr_pre_l <- sum(lambda * y_tr_pre)
  y_co_pre_wl <- sum(lambda * drop(Y_co_pre %*% omega))
  y_co_post_w <- sum(omega * y_co_post)
  y_tr_pre_l + (y_co_post_w - y_co_pre_wl)
}

# ── Closed-form Full Conformal engine (target-date SCM / SI) ─────────────────

.uc_full_conformal <- function(fit, Yco, y1, T_pre, level, score, band) {
  TT     <- length(y1)
  N_co   <- ncol(Yco)
  N      <- N_co + 1L
  T_post <- TT - T_pre
  pre    <- seq_len(T_pre)
  post   <- (T_pre + 1L):TT
  gamma  <- 1 - level
  k_cov  <- as.integer(floor(gamma * N + 1e-9))

  Y_all      <- cbind(y1, Yco)
  Y_all_pre  <- Y_all[pre, , drop = FALSE]
  Y_all_post <- Y_all[post, , drop = FALSE]

  method <- fit$method
  k_si   <- if (identical(method, "si")) {
    as.integer(min(.conformal_tuning(fit[["k"]], "k", "si"), T_pre, N_co))
  } else NULL

  # Fit each unit i in {1, ..., N} (1 = treated, 2..N = controls) using the
  # remaining N - 1 units on T_pre.
  w_score <- numeric(N)
  kappa   <- numeric(N_co)          # weight placed on unit 0 by control j
  V_co    <- matrix(0.0, T_post, N_co) # v_{t, j} at c = 0
  tau_hat <- numeric(T_post)

  for (idx in seq_len(N)) {
    Y_train_pre <- Y_all_pre[, -idx, drop = FALSE]
    y_i_pre     <- Y_all_pre[, idx]
    w_i <- if (identical(method, "scm")) {
      .uc_scm_fit_weights(Y_train_pre, y_i_pre, fit = fit)
    } else {
      drop(si_pcr_cpp(Y_train_pre, Y_all_post[, -idx, drop = FALSE],
                      matrix(y_i_pre, ncol = 1L), k_si)$W)
    }

    if (identical(score, "studentized")) {
      resid_pre <- y_i_pre - drop(Y_train_pre %*% w_i)
      sig_i     <- sqrt(mean(resid_pre^2))
      if (!is.finite(sig_i) || sig_i <= 1e-12)
        stop(paste0(
          "unit_conformal(score = \"studentized\"): unit ", idx,
          " has near-zero pre-treatment residual RMSPE (", format(sig_i),
          "), so studentized weights are undefined. Use score = \"abs\"."),
          call. = FALSE)
      w_score[idx] <- 1 / sig_i
    } else {
      w_score[idx] <- 1
    }

    Y_train_post <- Y_all_post[, -idx, drop = FALSE]
    pred_post    <- drop(Y_train_post %*% w_i)
    err_post     <- Y_all_post[, idx] - pred_post

    if (idx == 1L) {
      tau_hat <- err_post
    } else {
      j <- idx - 1L
      kappa[j]  <- w_i[1L]          # unit 0 is always column 1 in -idx
      V_co[, j] <- err_post
    }
  }

  w0   <- w_score[1L]
  w_co <- w_score[-1L]

  # Period-wise sets and p-values at c = 0
  sets_by_t <- vector("list", T_post)
  ci_lo_t   <- numeric(T_post)
  ci_hi_t   <- numeric(T_post)
  p_t       <- numeric(T_post)

  for (s in seq_len(T_post)) {
    tau_s  <- tau_hat[s]
    s0     <- w0 * abs(tau_s)
    sj     <- w_co * abs(V_co[s, ])
    p_t[s] <- (1 + sum(sj >= s0 * (1 - 1e-12))) / N

    sets_s <- lapply(seq_len(N_co), function(j)
      .uc_quad_set(w_co[j], V_co[s, j], kappa[j], w0, tau_s))
    sets_by_t[[s]] <- sets_s
    bounds <- .uc_interval_cover(sets_s, k_cov)
    ci_lo_t[s] <- bounds[1L]
    ci_hi_t[s] <- bounds[2L]
  }

  # Simultaneous band (SM A2.6 exact outer set)
  band_lo_t <- NULL
  band_hi_t <- NULL
  if (isTRUE(band)) {
    if (k_cov == 0L || T_post == 1L) {
      band_lo_t <- ci_lo_t
      band_hi_t <- ci_hi_t
    } else {
      # Compute g_{js} for each control j and post date s
      g_mat <- matrix(Inf, nrow = T_post, ncol = N_co)
      for (j in seq_len(N_co)) {
        if (w_co[j] * abs(kappa[j]) < w0 - 1e-12 * w0) {
          for (s in seq_len(T_post)) {
            int_js <- sets_by_t[[s]][[j]] # single bounded interval [a, b]
            a_js   <- int_js[1L, 1L]
            b_js   <- int_js[1L, 2L]
            g_mat[s, j] <- max(
              w_co[j] * abs(V_co[s, j] + kappa[j] * a_js),
              w_co[j] * abs(V_co[s, j] + kappa[j] * b_js)
            )
          }
        }
      }

      band_lo_t <- numeric(T_post)
      band_hi_t <- numeric(T_post)
      for (t_idx in seq_len(T_post)) {
        tau_s <- tau_hat[t_idx]
        tilde_sets <- lapply(seq_len(N_co), function(j) {
          g_other <- max(g_mat[-t_idx, j])
          if (!is.finite(g_other)) return(matrix(c(-Inf, Inf), nrow = 1L))
          box <- matrix(c(tau_s - g_other / w0, tau_s + g_other / w0),
                        nrow = 1L)
          .uc_union_sets(sets_by_t[[t_idx]][[j]], box)
        })
        b_bounds <- .uc_interval_cover(tilde_sets, k_cov)
        band_lo_t[t_idx] <- b_bounds[1L]
        band_hi_t[t_idx] <- b_bounds[2L]
      }
    }
  }

  # Overall ATT interval and p-value (post-treatment averages, Section 6.2)
  tau_avg  <- mean(tau_hat)
  v_co_avg <- colMeans(V_co)
  s0_avg   <- w0 * abs(tau_avg)
  sj_avg   <- w_co * abs(v_co_avg)
  p_avg    <- (1 + sum(sj_avg >= s0_avg * (1 - 1e-12))) / N

  sets_avg <- lapply(seq_len(N_co), function(j)
    .uc_quad_set(w_co[j], v_co_avg[j], kappa[j], w0, tau_avg))
  avg_bounds <- .uc_interval_cover(sets_avg, k_cov)

  terms <- data.frame(
    time      = fit$times[post],
    estimate  = tau_hat,
    conf.low  = ci_lo_t,
    conf.high = ci_hi_t,
    p.value   = p_t,
    stringsAsFactors = FALSE
  )
  if (isTRUE(band)) {
    terms$band.low  <- band_lo_t
    terms$band.high <- band_hi_t
  }

  list(
    terms               = terms,
    p_value             = p_avg,
    ci_lower            = avg_bounds[1L],
    ci_upper            = avg_bounds[2L],
    guaranteed_coverage = .uc_guaranteed_coverage(N, gamma, "full_conformal"),
    q_gamma             = NULL,
    n_failed            = 0L
  )
}

# ── Jackknife+ engine (all 6 estimators) ─────────────────────────────────────

.uc_jackknife_fold <- function(j, fit, Yco, y1, T_pre, score) {
  TT     <- length(y1)
  N_co   <- ncol(Yco)
  pre    <- seq_len(T_pre)
  post   <- (T_pre + 1L):TT
  T_post <- length(post)
  method <- fit$method

  Y_B    <- Yco[, -j, drop = FALSE]
  Y_B_pre  <- Y_B[pre, , drop = FALSE]
  Y_B_post <- Y_B[post, , drop = FALSE]
  y0_pre   <- y1[pre]
  yj_pre   <- Yco[pre, j]
  yj_post  <- Yco[post, j]

  unconv <- 0L

  if (identical(method, "scm")) {
    X0_mat <- fit[["X0_mat"]]
    X1_vec <- fit[["X1_vec"]]
    X_B    <- if (!is.null(X0_mat)) X0_mat[, -j, drop = FALSE] else NULL
    x_j    <- if (!is.null(X0_mat)) X0_mat[, j] else NULL
    w0_f <- .uc_scm_fit_weights(Y_B_pre, y0_pre, X_B = X_B, x_i = X1_vec,
                                fit = fit)
    wj_f <- .uc_scm_fit_weights(Y_B_pre, yj_pre, X_B = X_B, x_i = x_j,
                                fit = fit)
    m0_pre  <- drop(Y_B_pre %*% w0_f)
    mj_pre  <- drop(Y_B_pre %*% wj_f)
    m0_post <- drop(Y_B_post %*% w0_f)
    mj_post <- drop(Y_B_post %*% wj_f)
    m0_avg  <- mean(m0_post)
    mj_avg  <- mean(mj_post)
  } else if (identical(method, "sdid")) {
    # When covariates were partialled out, re-run the TWFE partial-out
    # symmetrically holding both unit 0 and control j out from post dates.
    # The fold's coefficients differ from the main fit's, while the interval
    # is centred on y1, which carries the main fit's partial-out: shift unit
    # 0's predictions by the difference X_0t (beta_fold - beta_fit) so that
    # y1 - prediction equals the fold's own effect estimate.
    shift0 <- numeric(T_post)
    if (!is.null(fit[["Y_raw"]]) && !is.null(fit[["X_arr"]])) {
      B_idx  <- fit$idx_co[-j]
      tr_idx <- fit$idx_tr
      j_idx  <- fit$idx_co[j]
      cols   <- c(B_idx, tr_idx, j_idx)
      n_B    <- length(B_idx)
      pan_f  <- list(
        idx_control = seq_len(n_B),
        idx_treat   = c(n_B + 1L, n_B + 2L),
        T_adopt     = c(rep(NA_integer_, n_B), T_pre + 1L, T_pre + 1L)
      )
      po <- suppressWarnings(.sdid_partial_out_staggered(
        fit$Y_raw[, cols, drop = FALSE],
        fit$X_arr[, cols, , drop = FALSE], pan_f))
      Y_B      <- po$Y_tilde[, seq_len(n_B), drop = FALSE]
      Y_B_pre  <- Y_B[pre, , drop = FALSE]
      Y_B_post <- Y_B[post, , drop = FALSE]
      y0_pre   <- po$Y_tilde[pre, n_B + 1L]
      yj_pre   <- po$Y_tilde[pre, n_B + 2L]
      yj_post  <- po$Y_tilde[post, n_B + 2L]
      shift0   <- y1[post] - po$Y_tilde[post, n_B + 1L]
    }
    zeta2  <- fit[["zeta2"]]
    zt_eff <- .uc_sdid_zeta_t_eff(fit)

    uw0 <- .uc_sdid_unit_weights(Y_B_pre, y0_pre, zeta2)
    uwj <- .uc_sdid_unit_weights(Y_B_pre, yj_pre, zeta2)

    m0_pre <- uw0$omega0 + drop(Y_B_pre %*% uw0$omega)
    mj_pre <- uwj$omega0 + drop(Y_B_pre %*% uwj$omega)

    m0_post <- numeric(T_post)
    mj_post <- numeric(T_post)
    for (s in seq_len(T_post)) {
      lam_s <- .uc_sdid_time_weights(Y_B_pre, Y_B_post[s, ], zt_eff)
      m0_post[s] <- .uc_sdid_predict(Y_B_pre, Y_B_post[s, ], y0_pre,
                                     uw0$omega, lam_s)
      mj_post[s] <- .uc_sdid_predict(Y_B_pre, Y_B_post[s, ], yj_pre,
                                     uwj$omega, lam_s)
    }

    lam_avg <- .uc_sdid_time_weights(Y_B_pre, colMeans(Y_B_post), zt_eff)
    m0_avg  <- .uc_sdid_predict(Y_B_pre, colMeans(Y_B_post), y0_pre,
                                uw0$omega, lam_avg)
    mj_avg  <- .uc_sdid_predict(Y_B_pre, colMeans(Y_B_post), yj_pre,
                                uwj$omega, lam_avg)
    m0_post <- m0_post + shift0
    m0_avg  <- m0_avg + mean(shift0)
  } else if (identical(method, "si")) {
    k_si <- as.integer(min(.conformal_tuning(fit[["k"]], "k", "si"),
                           T_pre, N_co - 1L))
    res <- si_pcr_cpp(Y_B_pre, Y_B_post, cbind(y0_pre, yj_pre), k_si)
    m0_pre  <- drop(Y_B_pre %*% res$W[, 1L])
    mj_pre  <- drop(Y_B_pre %*% res$W[, 2L])
    m0_post <- res$Y_hat[, 1L]
    mj_post <- res$Y_hat[, 2L]
    m0_avg  <- mean(m0_post)
    mj_avg  <- mean(mj_post)
  } else if (identical(method, "gsc")) {
    r_gsc <- as.integer(.conformal_tuning(fit[["r"]], "r", "gsc"))
    force <- fit[["force"]] %||% "none"
    X_co_arr <- fit[["covariate_array_co"]]
    X_tr_arr <- fit[["covariate_array_tr"]]
    X_B  <- if (is.null(X_co_arr)) NULL else X_co_arr[, -j, , drop = FALSE]
    X_j  <- if (is.null(X_co_arr)) NULL else X_co_arr[, j, , drop = FALSE]

    ctrl <- .gsc_fit_controls(Y_B, X_B, r_gsc, force)
    if (!isTRUE(ctrl$converged)) unconv <- 1L
    cf0 <- .gsc_counterfactual(ctrl, matrix(y1, ncol = 1L), X_tr_arr, T_pre)
    cfj <- .gsc_counterfactual(ctrl, Yco[, j, drop = FALSE], X_j, T_pre)
    p0  <- drop(cf0$Y_hat)
    pj  <- drop(cfj$Y_hat)
    m0_pre  <- p0[pre]
    mj_pre  <- pj[pre]
    m0_post <- p0[post]
    mj_post <- pj[post]
    m0_avg  <- mean(m0_post)
    mj_avg  <- mean(mj_post)
  } else if (identical(method, "mc")) {
    lam      <- .conformal_tuning(fit[["lambda"]], "lambda", "mc")
    fixedeff <- isTRUE(fit[["fixedeff"]] %||% FALSE)
    mc_predict <- function(y_hold) {
      M <- cbind(Y_B, y_hold)
      O <- matrix(1.0, nrow = TT, ncol = N_co)
      O[post, N_co] <- 0.0
      na_mask <- is.na(M)
      if (any(na_mask)) {
        O[na_mask] <- 0.0
        M[na_mask] <- 0.0
      }
      soft_impute_fe_cpp(M, O, lam, estimate_fe = fixedeff)$M[, N_co]
    }
    p0 <- mc_predict(y1)
    pj <- mc_predict(Yco[, j])
    m0_pre  <- p0[pre]
    mj_pre  <- pj[pre]
    m0_post <- p0[post]
    mj_post <- pj[post]
    m0_avg  <- mean(m0_post)
    mj_avg  <- mean(mj_post)
  } else if (identical(method, "tasc")) {
    r_tasc   <- as.integer(fit[["r"]] %||% 2L)
    em_iter  <- as.integer(fit[["em_iter"]] %||% 2000L)
    em_tol   <- fit[["em_tol"]] %||% 1e-6
    fix_A    <- isTRUE(fit[["fix_A"]])
    interc   <- if (is.null(fit[["intercept"]])) TRUE else isTRUE(fit[["intercept"]])
    tasc_predict <- function(y_hold) {
      Y_obs <- cbind(Y_B, y_hold)
      Y_obs[post, N_co] <- NA_real_
      res <- .tasc_em_fit(Y_obs, T_pre, idx_tr = N_co,
                          idx_co = seq_len(N_co - 1L),
                          r = r_tasc, em_iter = em_iter, em_tol = em_tol,
                          fix_A = fix_A, intercept = interc,
                          warn_unconverged = FALSE)
      if (!isTRUE(res$converged)) unconv <<- 1L
      res$Y_hat_all[, N_co]
    }
    p0 <- tasc_predict(y1)
    pj <- tasc_predict(Yco[, j])
    m0_pre  <- p0[pre]
    mj_pre  <- pj[pre]
    m0_post <- p0[post]
    mj_post <- pj[post]
    m0_avg  <- mean(m0_post)
    mj_avg  <- mean(mj_post)
  }

  if (identical(score, "studentized")) {
    sig0 <- sqrt(mean((y0_pre - m0_pre)^2))
    sigj <- sqrt(mean((yj_pre - mj_pre)^2))
    if (!is.finite(sig0) || sig0 <= 1e-12 || !is.finite(sigj) || sigj <= 1e-12) {
      return(list(ok = FALSE, unconverged = unconv,
                  reason = paste0("near-zero pre-treatment residual RMSPE ",
                                  "(studentized score)")))
    }
    w0 <- 1 / sig0
    wj <- 1 / sigj
  } else {
    w0 <- 1
    wj <- 1
  }

  R_j0      <- wj * abs(yj_post - mj_post)
  R_j0_avg  <- wj * abs(mean(yj_post) - mj_avg)
  R_j0_path <- max(R_j0)

  list(
    ok          = TRUE,
    unconverged = unconv,
    a_t         = m0_post - R_j0 / w0,
    b_t         = m0_post + R_j0 / w0,
    a_avg       = m0_avg  - R_j0_avg / w0,
    b_avg       = m0_avg  + R_j0_avg / w0,
    a_band      = m0_post - R_j0_path / w0,
    b_band      = m0_post + R_j0_path / w0
  )
}

.uc_jackknife_plus <- function(fit, Yco, y1, T_pre, level, score, band) {
  TT     <- length(y1)
  N_co   <- ncol(Yco)
  N      <- N_co + 1L
  T_post <- TT - T_pre
  pre    <- seq_len(T_pre)
  post   <- (T_pre + 1L):TT
  gamma  <- 1 - level
  q_g    <- .uc_q_gamma(N, gamma)

  method <- fit$method
  if (identical(method, "gsc")) {
    r_gsc <- as.integer(.conformal_tuning(fit[["r"]], "r", "gsc"))
    force <- fit[["force"]] %||% "none"
    if (N_co - 1L < r_gsc)
      stop(sprintf("GSC Jackknife+ needs at least r + 1 = %d control units (found %d).",
                   r_gsc + 1L, N_co), call. = FALSE)
    if (T_pre < r_gsc + .gsc_unit_fe(force))
      stop("GSC Jackknife+: T_pre is too small for the fitted factor count.",
           call. = FALSE)
  } else if (identical(method, "tasc")) {
    r_tasc <- as.integer(fit[["r"]] %||% 2L)
    if (N_co < r_tasc)
      stop(sprintf("TASC Jackknife+ needs at least r = %d control units (found %d).",
                   r_tasc, N_co), call. = FALSE)
  } else if (identical(method, "sdid")) {
    if (T_pre < 2L)
      stop("SDID Jackknife+ requires at least 2 pre-treatment periods.",
           call. = FALSE)
  }

  a_mat    <- matrix(-Inf, nrow = T_post, ncol = N_co)
  b_mat    <- matrix( Inf, nrow = T_post, ncol = N_co)
  aband_mat <- matrix(-Inf, nrow = T_post, ncol = N_co)
  bband_mat <- matrix( Inf, nrow = T_post, ncol = N_co)
  a_avg    <- rep(-Inf, N_co)
  b_avg    <- rep( Inf, N_co)
  n_failed <- 0L
  n_unconv <- 0L
  fail_msg <- character(0L)

  for (j in seq_len(N_co)) {
    res <- tryCatch(
      .uc_jackknife_fold(j, fit, Yco, y1, T_pre, score),
      error = function(e) list(ok = FALSE, unconverged = 0L,
                               reason = conditionMessage(e))
    )
    if (!isTRUE(res$ok)) fail_msg <- c(fail_msg, res[["reason"]])
    n_unconv <- n_unconv + (res[["unconverged"]] %||% 0L)
    if (!isTRUE(res$ok)) {
      # Proposition 6 (Liu & Xu 2026): failed folds use [-Inf, +Inf] while
      # keeping N and q_gamma unchanged.
      n_failed <- n_failed + 1L
    } else {
      a_mat[, j]     <- res$a_t
      b_mat[, j]     <- res$b_t
      aband_mat[, j] <- res$a_band
      bband_mat[, j] <- res$b_band
      a_avg[j]       <- res$a_avg
      b_avg[j]       <- res$b_avg
    }
  }

  if (n_failed > 0L) {
    reasons <- table(fail_msg)
    warning(sprintf(paste0(
      "unit_conformal(): %d of %d Jackknife+ fold(s) failed or had ",
      "degenerate pre-treatment fit; treated as [-Inf, Inf] per ",
      "Proposition 6. Reasons:\n%s"),
      n_failed, N_co,
      paste(sprintf("  [%d fold(s)] %s", as.integer(reasons), names(reasons)),
            collapse = "\n")), call. = FALSE)
  }
  if (n_unconv > 0L)
    warning(sprintf(
      "unit_conformal(): %d Jackknife+ fold fit(s) stopped before iterative convergence.",
      n_unconv), call. = FALSE)

  y0_post <- y1[post]
  ci_lo_t <- numeric(T_post)
  ci_hi_t <- numeric(T_post)
  p_t     <- numeric(T_post)
  band_lo_t <- numeric(T_post)
  band_hi_t <- numeric(T_post)

  for (s in seq_len(T_post)) {
    p_t[s] <- .uc_jackknife_pval(y0_post[s], a_mat[s, ], b_mat[s, ], c = 0)
    if (q_g == 0L) {
      ci_lo_t[s]   <- -Inf
      ci_hi_t[s]   <-  Inf
      band_lo_t[s] <- -Inf
      band_hi_t[s] <-  Inf
    } else {
      u <- N - q_g
      a_ord <- sort(a_mat[s, ])
      b_ord <- sort(b_mat[s, ])
      ci_lo_t[s] <- y0_post[s] - b_ord[u]
      ci_hi_t[s] <- y0_post[s] - a_ord[q_g]

      if (isTRUE(band)) {
        ab_ord <- sort(aband_mat[s, ])
        bb_ord <- sort(bband_mat[s, ])
        band_lo_t[s] <- y0_post[s] - bb_ord[u]
        band_hi_t[s] <- y0_post[s] - ab_ord[q_g]
      }
    }
  }

  y0_avg <- mean(y0_post)
  p_avg  <- .uc_jackknife_pval(y0_avg, a_avg, b_avg, c = 0)
  if (q_g == 0L) {
    ci_avg <- c(-Inf, Inf)
  } else {
    u <- N - q_g
    ci_avg <- c(y0_avg - sort(b_avg)[u], y0_avg - sort(a_avg)[q_g])
  }

  # Period-wise point estimates from the original fit
  if (identical(method, "sdid")) {
    # Target-specific SDID point estimate per post date on the full donor pool
    Y_co_pre  <- Yco[pre, , drop = FALSE]
    Y_co_post <- Yco[post, , drop = FALSE]
    zt_eff    <- .uc_sdid_zeta_t_eff(fit)
    y0_pre    <- y1[pre]
    omega     <- if (!is.null(fit[["unit_weights"]]) &&
                     length(fit[["unit_weights"]]) == ncol(Y_co_pre)) {
      unname(fit[["unit_weights"]])
    } else {
      .uc_sdid_unit_weights(Y_co_pre, y0_pre, fit[["zeta2"]])$omega
    }
    est_t <- vapply(seq_len(T_post), function(s) {
      lam_s <- .uc_sdid_time_weights(Y_co_pre, Y_co_post[s, ], zt_eff)
      y0_post[s] - .uc_sdid_predict(Y_co_pre, Y_co_post[s, ], y0_pre,
                                    omega, lam_s)
    }, numeric(1L))
  } else {
    y_hat_full <- synthetic_outcomes(fit)
    est_t      <- y0_post - y_hat_full[post]
  }

  terms <- data.frame(
    time      = fit$times[post],
    estimate  = est_t,
    conf.low  = ci_lo_t,
    conf.high = ci_hi_t,
    p.value   = p_t,
    stringsAsFactors = FALSE
  )
  if (isTRUE(band)) {
    terms$band.low  <- band_lo_t
    terms$band.high <- band_hi_t
  }

  list(
    terms               = terms,
    p_value             = p_avg,
    ci_lower            = ci_avg[1L],
    ci_upper            = ci_avg[2L],
    guaranteed_coverage = .uc_guaranteed_coverage(N, gamma, "jackknife_plus"),
    q_gamma             = q_g,
    n_failed            = n_failed
  )
}

# ── Rolling Pre-Treatment Diagnostics (Liu & Xu 2026, S.3.3 & SM A1.7) ───────

.uc_diagnostic_slice <- function(fit, Yco, y1, T_pre, procedure, level, score,
                                 m) {
  TT     <- length(y1)
  N_co   <- ncol(Yco)
  method <- fit$method
  H_min  <- T_pre - m

  if (H_min < 1L)
    stop(sprintf(paste0(
      "pre_diagnostics = %d leaves %d period(s) before the earliest ",
      "diagnostic date; at least 1 is required (T_pre = %d)."),
      m, H_min, T_pre), call. = FALSE)
  if (identical(method, "sdid") && H_min < 2L)
    stop(sprintf(paste0(
      "pre_diagnostics = %d for SDID requires at least 2 periods before ",
      "the earliest diagnostic date (T_pre - pre_diagnostics = %d)."),
      m, H_min), call. = FALSE)
  if (identical(method, "gsc")) {
    r_gsc <- as.integer(.conformal_tuning(fit[["r"]], "r", "gsc"))
    force <- fit[["force"]] %||% "none"
    need  <- r_gsc + .gsc_unit_fe(force)
    if (H_min < need)
      stop(sprintf(paste0(
        "pre_diagnostics = %d for GSC requires at least %d periods before ",
        "the earliest diagnostic date (found %d)."),
        m, need, H_min), call. = FALSE)
  }
  if (identical(method, "scm") && !is.null(fit[["X0_mat"]]))
    stop(paste0(
      "pre_diagnostics is not supported for predictor-based SCM fits ",
      "because predictor matrices are aggregated over the full original ",
      "pre-treatment window."), call. = FALSE)
  if (identical(method, "scm") && !is.null(fit[["z_rows"]]) &&
      max(fit[["z_rows"]]) > H_min)
    stop(sprintf(paste0(
      "pre_diagnostics = %d shortens the history before diagnostic dates to ",
      "%d period(s), which excludes part of v_window."),
      m, H_min), call. = FALSE)
  if (identical(method, "scm") && !is.null(fit[["v_rows"]]) &&
      length(fit[["v_rows"]]) < T_pre && H_min < 2L)
    stop(sprintf(paste0(
      "pre_diagnostics = %d for v_selection = \"oos\" requires at least 2 ",
      "periods before the earliest diagnostic date (found %d)."),
      m, H_min), call. = FALSE)

  is_full_path <- method %in% c("gsc", "mc", "tasc")

  # Full-control GSC factors are identical across all diagnostic dates s
  gsc_ctrl_all <- if (identical(method, "gsc")) {
    r_gsc <- as.integer(.conformal_tuning(fit[["r"]], "r", "gsc"))
    force <- fit[["force"]] %||% "none"
    .gsc_fit_controls(Yco, fit[["covariate_array_co"]], r_gsc, force)
  } else NULL

  diag_dates <- seq.int(T_pre - m + 1L, T_pre)
  rows       <- vector("list", length(diag_dates))

  for (k in seq_along(diag_dates)) {
    s       <- diag_dates[k]
    T_pre_s <- s - 1L
    pre_s   <- seq_len(T_pre_s)
    fit_s   <- fit
    fit_s$T_pre <- T_pre_s

    if (!is_full_path) {
      # Target-date estimators (SM Figure 4a): truncate panel at diagnostic s
      fit_s$times <- fit$times[seq_len(s)]
      Yco_s <- Yco[seq_len(s), , drop = FALSE]
      y1_s  <- y1[seq_len(s)]

      if (identical(method, "scm") && !is.null(fit[["v_rows"]])) {
        fit_s$v_rows <- if (length(fit[["v_rows"]]) < T_pre) {
          seq_len(max(1L, T_pre_s %/% 2L))
        } else {
          seq_len(T_pre_s)
        }
      }

      if (identical(method, "sdid") &&
          !is.null(fit[["Y_raw"]]) && !is.null(fit[["X_arr"]])) {
        fit_s$Y_raw <- fit$Y_raw[seq_len(s), , drop = FALSE]
        fit_s$X_arr <- fit$X_arr[seq_len(s), , , drop = FALSE]
        cols_all <- c(fit$idx_co, fit$idx_tr)
        n_co_all <- length(fit$idx_co)
        pan_s <- list(
          idx_control = seq_len(n_co_all),
          idx_treat   = n_co_all + 1L,
          T_adopt     = c(rep(NA_integer_, n_co_all), s)
        )
        po_s <- suppressWarnings(.sdid_partial_out_staggered(
          fit_s$Y_raw[, cols_all, drop = FALSE],
          fit_s$X_arr[, cols_all, , drop = FALSE], pan_s))
        Yco_s <- po_s$Y_tilde[, seq_len(n_co_all), drop = FALSE]
        y1_s  <- po_s$Y_tilde[, n_co_all + 1L]
      }

      res_s <- if (identical(procedure, "full_conformal")) {
        .uc_full_conformal(fit_s, Yco_s, y1_s, T_pre_s, level, score,
                           band = FALSE)
      } else {
        .uc_jackknife_plus(fit_s, Yco_s, y1_s, T_pre_s, level, score,
                           band = FALSE)
      }
      row1 <- res_s$terms[1L, , drop = FALSE]

      if (identical(procedure, "jackknife_plus")) {
        # Compute genuine out-of-sample pseudo-effect on all controls (Alg. 3)
        est_s <- if (identical(method, "scm")) {
          w_s <- .uc_scm_fit_weights(Yco_s[pre_s, , drop = FALSE], y1_s[pre_s],
                                     fit = fit_s)
          y1_s[s] - sum(w_s * Yco_s[s, ])
        } else if (identical(method, "si")) {
          k_si <- as.integer(min(.conformal_tuning(fit[["k"]], "k", "si"),
                                 T_pre_s, N_co))
          y1_s[s] - drop(si_pcr_cpp(Yco_s[pre_s, , drop = FALSE],
                                    Yco_s[s, , drop = FALSE],
                                    matrix(y1_s[pre_s], ncol = 1L),
                                    k_si)$Y_hat)[1L]
        } else {
          # SDID
          uw_s  <- .uc_sdid_unit_weights(Yco_s[pre_s, , drop = FALSE],
                                         y1_s[pre_s], fit[["zeta2"]])
          lam_s <- .uc_sdid_time_weights(Yco_s[pre_s, , drop = FALSE],
                                         Yco_s[s, ], .uc_sdid_zeta_t_eff(fit))
          y1_s[s] - .uc_sdid_predict(Yco_s[pre_s, , drop = FALSE], Yco_s[s, ],
                                     y1_s[pre_s], uw_s$omega, lam_s)
        }
        row1$estimate <- est_s
      }
    } else {
      # Full-path estimators (SM Figure 4b, retrospective): keep complete 1:TT
      # for training controls and hide s:TT for held-out units.
      res_s <- .uc_jackknife_plus(fit_s, Yco, y1, T_pre_s, level, score,
                                  band = FALSE)
      row1  <- res_s$terms[1L, , drop = FALSE]

      est_s <- if (identical(method, "gsc")) {
        cf_all <- .gsc_counterfactual(gsc_ctrl_all, matrix(y1, ncol = 1L),
                                      fit[["covariate_array_tr"]], T_pre_s)
        y1[s] - drop(cf_all$Y_hat)[s]
      } else if (identical(method, "mc")) {
        lam      <- .conformal_tuning(fit[["lambda"]], "lambda", "mc")
        fixedeff <- isTRUE(fit[["fixedeff"]] %||% FALSE)
        M <- cbind(Yco, y1)
        O <- matrix(1.0, nrow = TT, ncol = N_co + 1L)
        O[s:TT, N_co + 1L] <- 0.0
        na_mask <- is.na(M)
        if (any(na_mask)) {
          O[na_mask] <- 0.0
          M[na_mask] <- 0.0
        }
        M_hat <- soft_impute_fe_cpp(M, O, lam, estimate_fe = fixedeff)$M
        y1[s] - M_hat[s, N_co + 1L]
      } else {
        # TASC
        r_tasc  <- as.integer(fit[["r"]] %||% 2L)
        em_iter <- as.integer(fit[["em_iter"]] %||% 2000L)
        em_tol  <- fit[["em_tol"]] %||% 1e-6
        fix_A   <- isTRUE(fit[["fix_A"]])
        interc  <- if (is.null(fit[["intercept"]])) TRUE else isTRUE(fit[["intercept"]])
        Y_obs <- cbind(Yco, y1)
        Y_obs[s:TT, N_co + 1L] <- NA_real_
        tasc_res <- .tasc_em_fit(Y_obs, T_pre_s, idx_tr = N_co + 1L,
                                 idx_co = seq_len(N_co),
                                 r = r_tasc, em_iter = em_iter,
                                 em_tol = em_tol, fix_A = fix_A,
                                 intercept = interc, warn_unconverged = FALSE)
        y1[s] - tasc_res$Y_hat_all[s, N_co + 1L]
      }
      row1$estimate <- est_s
    }

    rows[[k]] <- row1
  }

  pre_df <- do.call(rbind, rows)
  rownames(pre_df) <- NULL
  pre_df
}

#' Unit-Level Conformal Inference for Panel Counterfactuals
#'
#' Implements the unit-level conformal inference procedures of Liu & Xu (2026),
#' which construct finite-sample confidence sets and p-values from the
#' across-unit exchangeability of out-of-sample prediction errors rather than
#' the temporal exchangeability of residuals used by [conformal_inference()].
#' Unlike time-based conformal inference, unit-level procedures do not require
#' stationarity or weak temporal dependence and apply to all six estimators in
#' `coresynth`, including `"tasc"`.
#'
#' Two procedures are supported for sharp single-treated-unit fits:
#'
#' * **Full Conformal** (`procedure = "full_conformal"`): available in exact
#'   closed form for target-date estimators whose control predictions are
#'   affine in the treated outcome — outcomes-only `"scm"` (no `predictors`)
#'   and `"si"`. Guarantees period-wise coverage at least
#'   \eqn{1 - \lfloor (1 - \text{level}) N \rfloor / N}, where
#'   \eqn{N = N_{co} + 1}.
#' * **Jackknife+** (`procedure = "jackknife_plus"`): holds out the treated
#'   unit together with each control in turn, avoiding any candidate search and
#'   applying to all six estimators (`"scm"` with or without `predictors`,
#'   target-specific `"sdid"`, `"gsc"`, `"mc"`, `"tasc"`, and `"si"`).
#'   Guarantees period-wise coverage at least
#'   \eqn{1 - (2 q_\gamma - 1)_+ / N}, where
#'   \eqn{q_\gamma = \max\{q : (2q - 1)_+ \le (1 - \text{level}) N\}}.
#'
#' The default `procedure = "auto"` selects `"full_conformal"` for outcomes-only
#' `"scm"` and `"si"`, and `"jackknife_plus"` for `"gsc"`, `"mc"`, `"tasc"`,
#' `"sdid"`, and predictor-based `"scm"`.
#'
#' When `score = "studentized"`, each unit's prediction error is divided by its
#' in-sample pre-treatment root mean squared residual within the same fit,
#' improving coverage stability when error variances differ across units.
#' When `band = TRUE`, the returned `terms` table also includes simultaneous
#' confidence bands covering the entire post-treatment effect path jointly
#' (Theorem 4 and SM A2.6, requiring no additional fits for Jackknife+).
#'
#' Tuning choices (`r`, `lambda`, `k`, `zeta2`, `zeta_t`, `em_iter`, `em_tol`,
#' etc.) are held fixed at the values recorded by `fit` across all refits. When
#' a tuning parameter was selected data-drivenly (such as cross-validation), the
#' finite-sample coverage guarantee is conditional on that selection.
#'
#' When `pre_diagnostics > 0`, also runs the rolling pre-treatment diagnostics
#' of Liu & Xu (2026, Section 3.3 and SM Algorithm 3) over the last
#' `pre_diagnostics` pre-treatment dates, storing period-wise pseudo-effect
#' estimates, confidence intervals, and p-values for zero pseudo-effect in
#' `pre_terms`. Target-date diagnostics truncate the panel at each diagnostic
#' date; full-path diagnostics (`"gsc"`, `"mc"`, `"tasc"`) hold out each held-out
#' unit's outcomes from the diagnostic date onward while keeping complete panels
#' for training controls. These diagnostics are unadjusted for multiplicity and
#' should not be read as joint pre-trend tests.
#'
#' @param fit A sharp `coresynth` object from [scm_fit()] with a single treated
#'   unit.
#' @param procedure Inference procedure: `"auto"` (default),
#'   `"jackknife_plus"`, or `"full_conformal"`.
#' @param level Nominal confidence level in `(0, 1)` (default `0.90`). Notice
#'   that finite-sample guaranteed coverage moves in steps of `1 / N` (full
#'   conformal) or `2 / N` (Jackknife+), and returns `[-Inf, Inf]` when
#'   `(1 - level) < 1 / N`.
#' @param score Nonconformity score weighting: `"abs"` (default, absolute
#'   prediction error) or `"studentized"` (scaled by each unit's pre-treatment
#'   residual RMSPE).
#' @param band Logical; if `TRUE` (default), include simultaneous confidence
#'   bands (`band.low`, `band.high`) across all post-treatment periods in
#'   `terms`.
#' @param pre_diagnostics Non-negative integer giving how many of the last
#'   pre-treatment dates to evaluate with rolling pre-treatment diagnostics
#'   (default `0L`, disabled).
#' @param ... Unused.
#'
#' @return A list of class `c("unit_conformal", "coresynth_inference")`
#'   compatible with [tidy()], [glance()], and [plot()], containing:
#'   * `estimate`: overall ATT from `fit`.
#'   * `se`: `NA_real_` (conformal inference has no standard error).
#'   * `p_value`: p-value for \eqn{H_0: \bar{\tau} = 0} (overall post-treatment
#'     average effect).
#'   * `ci_lower`, `ci_upper`: finite-sample confidence interval for the
#'     overall post-treatment average effect.
#'   * `terms`: `data.frame` of period-wise post-treatment results with columns
#'     `time`, `estimate`, `conf.low`, `conf.high`, `p.value`, plus `band.low`
#'     and `band.high` when `band = TRUE`.
#'   * `pre_terms`: `data.frame` of rolling pre-treatment diagnostic results
#'     when `pre_diagnostics > 0L` (or `NULL`).
#'   * `procedure`, `score`, `level`, `guaranteed_coverage`, `q_gamma`,
#'     `n_failed`, `pre_diagnostics`, `base_method`, `method`, `n_controls`,
#'     `alternative`, `staggered`.
#'
#' @references Liu, Z., & Xu, Y. (2026). Unit-Level Conformal Inference for
#'   Panel Counterfactuals. *arXiv:2611.19209*.
#'
#' @seealso [overlap_check()] for the pre-treatment convex-hull overlap
#'   diagnostic, [conformal_inference()] for time-based conformal inference.
#' @export
unit_conformal <- function(
  fit,
  procedure       = c("auto", "jackknife_plus", "full_conformal"),
  level           = 0.90,
  score           = c("abs", "studentized"),
  band            = TRUE,
  pre_diagnostics = 0L,
  ...
) {
  procedure <- match.arg(procedure)
  score     <- match.arg(score)

  if (!inherits(fit, "coresynth"))
    stop("unit_conformal() requires a coresynth object from scm_fit().",
         call. = FALSE)
  if (inherits(fit, "coresynth_staggered") || inherits(fit, "coresynth_multiarm"))
    stop(paste0(
      "unit_conformal() supports sharp, single-treated-unit fits only. ",
      "For staggered or multi-arm fits, use sdid_inference(), ",
      "gsc_inference(), si_inference(), or scm_inference()."), call. = FALSE)

  method <- fit$method
  if (!method %in% c("scm", "sdid", "gsc", "mc", "tasc", "si"))
    stop(sprintf("unit_conformal() does not support method = '%s'.", method),
         call. = FALSE)

  n_tr <- if (identical(method, "scm")) {
    fit[["n_treated"]] %||% ncol(as.matrix(fit[["Y_treat_mat"]] %||% fit[["Y_treat"]]))
  } else if (identical(method, "sdid")) {
    fit[["N_tr"]] %||% fit[["n_treated"]] %||%
      ncol(as.matrix(fit[["Y_tr_pre"]] %||% fit[["Y_treat"]]))
  } else {
    yt <- fit[["Y_treat"]]
    if (is.matrix(yt)) ncol(yt) else 1L
  }
  if (!is.null(n_tr) && n_tr > 1L)
    stop(sprintf(paste0(
      "unit_conformal() requires a single treated unit (found %d). ",
      "Per Liu & Xu (2026, S.5.1), averaging multiple treated units' ",
      "prediction errors breaks single-unit exchangeability."),
      as.integer(n_tr)), call. = FALSE)

  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.logical(band) || length(band) != 1L || is.na(band))
    stop("`band` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(pre_diagnostics) || length(pre_diagnostics) != 1L ||
      !is.finite(pre_diagnostics) || pre_diagnostics < 0 ||
      pre_diagnostics != trunc(pre_diagnostics))
    stop("`pre_diagnostics` must be a single non-negative integer.",
         call. = FALSE)
  pre_diagnostics <- as.integer(pre_diagnostics)

  Yco <- donor_outcomes(fit)
  y1  <- treated_outcomes(fit)
  if (is.null(Yco) || is.null(y1))
    stop("fit does not contain the treated and donor outcome series required ",
         "by unit_conformal(). Re-run scm_fit() with the current version.",
         call. = FALSE)

  T_pre <- as.integer(fit[["T_pre"]])
  TT    <- length(y1)
  if (is.na(T_pre) || T_pre < 1L || TT - T_pre < 1L)
    stop("unit_conformal() requires at least 1 pre-treatment and 1 ",
         "post-treatment period.", call. = FALSE)
  if (any(!is.finite(y1)))
    stop("unit_conformal(): treated outcome series contains missing or ",
         "non-finite values.", call. = FALSE)

  ok_co <- apply(Yco, 2L, function(col) all(is.finite(col)))
  if (any(!ok_co)) {
    message(sprintf(
      "unit_conformal(): excluded %d control unit(s) with missing outcomes.",
      sum(!ok_co)))
    Yco <- Yco[, ok_co, drop = FALSE]
    # Keep control-indexed fit fields aligned with the retained donor pool
    # (SM A2.2 option 1: exclude incomplete controls before fitting).
    if (!is.null(fit[["covariate_array_co"]]) &&
        dim(fit[["covariate_array_co"]])[2L] == length(ok_co))
      fit$covariate_array_co <- fit$covariate_array_co[, ok_co, , drop = FALSE]
    if (!is.null(fit[["X0_mat"]]) && ncol(fit[["X0_mat"]]) == length(ok_co))
      fit$X0_mat <- fit$X0_mat[, ok_co, drop = FALSE]
    if (!is.null(fit[["idx_co"]]) && length(fit[["idx_co"]]) == length(ok_co))
      fit$idx_co <- fit$idx_co[ok_co]
    if (!is.null(fit[["unit_weights"]]) &&
        length(fit[["unit_weights"]]) == length(ok_co))
      fit$unit_weights <- NULL
  }
  if (ncol(Yco) < 2L)
    stop("unit_conformal() requires at least 2 completely observed control ",
         "units.", call. = FALSE)

  if (identical(method, "scm") &&
      length(fit[["excluded_donors"]] %||% character(0L)) > 0L)
    warning(paste0(
      "donor_mspe_threshold excluded donors using the treated unit's fit, ",
      "which breaks symmetry between treated and control units (Liu & Xu ",
      "2026, S.4.1)."), call. = FALSE)

  # Outcomes-only SCM (no X0_mat) and SI support closed-form full conformal.
  # Predictor-based SCM, SDID, and full-path estimators (GSC, MC, TASC) use
  # Jackknife+.
  has_predictors <- identical(method, "scm") && !is.null(fit[["X0_mat"]])
  fc_eligible    <- identical(method, "si") ||
    (identical(method, "scm") && !has_predictors)

  if (identical(procedure, "auto")) {
    procedure <- if (fc_eligible) "full_conformal" else "jackknife_plus"
  } else if (identical(procedure, "full_conformal") && !fc_eligible) {
    reason <- if (has_predictors) {
      "predictor-based SCM does not support closed-form full conformal here"
    } else if (identical(method, "sdid")) {
      "target-specific SDID is non-affine in the candidate effect and grid approximations lack Theorem 1's exact guarantee"
    } else {
      sprintf("'%s' is a full-path estimator requiring a joint search over all post-treatment dates (Liu & Xu 2026, Remark 1)", method)
    }
    stop(sprintf(
      "procedure = \"full_conformal\" is not supported because %s. Use procedure = \"jackknife_plus\" (or \"auto\").",
      reason), call. = FALSE)
  }

  pre_terms <- if (pre_diagnostics > 0L) {
    .uc_diagnostic_slice(fit, Yco, y1, T_pre, procedure, level, score,
                         pre_diagnostics)
  } else {
    NULL
  }

  out <- if (identical(procedure, "full_conformal")) {
    .uc_full_conformal(fit, Yco, y1, T_pre, level, score, band)
  } else {
    .uc_jackknife_plus(fit, Yco, y1, T_pre, level, score, band)
  }

  structure(list(
    estimate            = fit$estimate,
    se                  = NA_real_,
    p_value             = out$p_value,
    ci_lower            = out$ci_lower,
    ci_upper            = out$ci_upper,
    method              = paste0("unit_conformal_", procedure),
    n_controls          = ncol(Yco),
    alternative         = "two.sided",
    staggered           = FALSE,
    boot_ests           = NULL,
    terms               = out$terms,
    pre_terms           = pre_terms,
    pre_diagnostics     = pre_diagnostics,
    T_pre               = T_pre,
    times               = fit$times,
    procedure           = procedure,
    score               = score,
    level               = level,
    guaranteed_coverage = out$guaranteed_coverage,
    q_gamma             = out$q_gamma,
    n_failed            = out$n_failed,
    base_method         = method
  ), class = c("unit_conformal", "coresynth_inference"))
}

#' @export
print.unit_conformal <- function(x, digits = 4, ...) {
  cat(sprintf("Unit-Level Conformal Inference (Liu & Xu 2026, %s, base = %s, score = %s)\n",
              x$procedure, x$base_method, x$score))
  cat(sprintf("  Level         : %.2f (guaranteed coverage >= %.3f)\n",
              x$level, x$guaranteed_coverage))
  cat("  Estimate (ATT):", round(x$estimate, digits), "\n")
  cat("  p-value (ATT) :", round(x$p_value, digits), "\n")
  cat("  CI (ATT)      : [", round(x$ci_lower, digits), ", ",
      round(x$ci_upper, digits), "]\n", sep = "")
  cat("  Controls      :", x$n_controls, "\n")
  if (!is.null(x$n_failed) && x$n_failed > 0L)
    cat("  Failed folds  :", x$n_failed, "\n")
  if (!is.null(x$pre_terms) && nrow(x$pre_terms) > 0L) {
    cat(sprintf("\nPre-treatment diagnostics (m = %d):\n", nrow(x$pre_terms)))
    df_pre <- x$pre_terms
    for (col in setdiff(names(df_pre), "time")) df_pre[[col]] <- round(df_pre[[col]], digits)
    print(df_pre, row.names = FALSE)
  }
  cat("\nPeriod-wise results:\n")
  df <- x$terms
  for (col in setdiff(names(df), "time")) df[[col]] <- round(df[[col]], digits)
  print(df, row.names = FALSE)
  invisible(x)
}

#' @rdname unit_conformal
#' @param x A `unit_conformal` object returned by [unit_conformal()].
#' @param conf.int Logical; include confidence intervals and simultaneous bands
#'   in the returned table (default `TRUE`).
#' @param include_pre Logical; include rolling pre-treatment diagnostic rows
#'   (with a `period` column `"pre"` / `"post"`) when available (default
#'   `FALSE`, returning post-treatment rows only).
#' @importFrom broom tidy
#' @export
tidy.unit_conformal <- function(x, conf.int = TRUE, include_pre = FALSE, ...) {
  out <- x$terms
  if (isTRUE(include_pre) && !is.null(x$pre_terms) && nrow(x$pre_terms) > 0L) {
    pre <- x$pre_terms
    if ("band.low" %in% names(out) && !"band.low" %in% names(pre)) {
      pre$band.low  <- NA_real_
      pre$band.high <- NA_real_
    }
    pre$period <- "pre"
    out$period <- "post"
    out <- rbind(pre[, names(out), drop = FALSE], out)
    rownames(out) <- NULL
  }
  if (!isTRUE(conf.int)) {
    out$conf.low  <- NULL
    out$conf.high <- NULL
    out$band.low  <- NULL
    out$band.high <- NULL
  }
  out
}
