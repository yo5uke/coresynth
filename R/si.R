# ── Internal: multi-arm SI for K>1 treatment arms (Agarwal et al. 2025) ───────
# Latent unit factors v_il are arm-invariant, so the SVD basis from control arm
# (d=0) is reused for all K treatment arms. si_pcr_cpp is called once per arm.
.fit_si_multi <- function(Y, idx_by_arm, arm_levels, T_pre, k, times) {
  TT        <- nrow(Y)
  T_post    <- TT - T_pre
  pre_rows  <- seq_len(T_pre)
  post_rows <- (T_pre + 1L):TT
  idx_co    <- idx_by_arm[["0"]]
  treat_arms <- arm_levels[arm_levels != 0L]  # c(1L, ..., KL)

  N_co <- length(idx_co)
  if (N_co < 2L)
    stop("SI multi-arm: コントロールユニットが 2 個未満です。", call. = FALSE)
  if (T_pre < 1L)
    stop("SI multi-arm: 事前期間が 0 です。", call. = FALSE)
  if (T_post < 1L)
    stop("SI multi-arm: 事後期間が 0 です。", call. = FALSE)

  # Common donor matrices shared across all arms (unit factors are arm-invariant)
  Y_pre_co  <- Y[pre_rows,  idx_co, drop = FALSE]  # T_pre  x N_co
  Y_post_co <- Y[post_rows, idx_co, drop = FALSE]  # T_post x N_co

  if (is.null(k)) k <- max(1L, floor(sqrt(min(T_pre, N_co))))
  k <- as.integer(k)
  if (k > min(T_pre, N_co))
    stop(sprintf("k (%d) must be <= min(T_pre=%d, N_co=%d).", k, T_pre, N_co),
         call. = FALSE)

  arm_fits <- lapply(as.character(treat_arms), function(a) {
    idx_tr_d   <- idx_by_arm[[a]]
    Y_pre_tr_d <- Y[pre_rows,  idx_tr_d, drop = FALSE]
    Y_treat_d  <- Y[, idx_tr_d, drop = FALSE]

    res <- tryCatch(
      si_pcr_cpp(Y_pre_co, Y_post_co, Y_pre_tr_d, k),
      error = function(e)
        stop(sprintf("SI multi-arm arm %s: %s", a, conditionMessage(e)),
             call. = FALSE)
    )
    W_d       <- res$W       # N_co x N_tr_d
    Y_cf_d    <- res$Y_hat   # T_post x N_tr_d
    Y_synth_d <- rbind(Y_pre_co %*% W_d, Y_cf_d)  # T x N_tr_d
    gap_d     <- Y_treat_d - Y_synth_d
    att_d     <- mean(gap_d[post_rows, , drop = FALSE], na.rm = TRUE)

    list(
      arm       = as.integer(a),
      n_treated = length(idx_tr_d),
      T_pre     = T_pre,
      T_post    = T_post,
      estimate  = att_d,
      weight    = as.numeric(length(idx_tr_d)) * T_post,
      idx_tr    = idx_tr_d,
      idx_co    = idx_co,
      weights   = W_d,
      Y_cf      = Y_cf_d,
      Y_synth   = Y_synth_d,
      Y_treat   = Y_treat_d,
      gap       = gap_d
    )
  })
  names(arm_fits) <- as.character(treat_arms)

  # Weighted-average ATT: weight proportional to N_tr_d * T_post (Clarke et al. 2023)
  ws   <- vapply(arm_fits, `[[`, numeric(1L), "weight")
  taus <- vapply(arm_fits, `[[`, numeric(1L), "estimate")
  att  <- sum(ws * taus) / sum(ws)

  list(
    method        = "si",
    multi_arm     = TRUE,
    arm_levels    = treat_arms,
    k             = k,
    estimate      = att,
    arm_estimates = setNames(taus, as.character(treat_arms)),
    arm_fits      = arm_fits,
    Y_pre_co      = Y_pre_co,
    Y_post_co     = Y_post_co,
    idx_co        = idx_co,
    times         = times,
    T_pre         = T_pre,
    staggered     = FALSE,
    Y_treat       = NULL,
    Y_synth       = NULL,
    gap           = NULL,
    unit_weights  = NULL
  )
}

# ── Internal: cohort-by-cohort SI for staggered adoption ─────────────────────
# Clarke et al. (2023); Arkhangelsky et al. (2021) Appendix §8.
# SVD uses only Y_pre_co_g (T_pre_g × N_co_g): future-treated units in the
# "clean" control group are all pre-treatment during periods 1..T_pre_g, so
# there is no contamination from treatment in this window.
.fit_si_staggered <- function(pan, Y, k = NULL, control_group = "clean") {
  TT      <- nrow(Y)
  T_adopt <- pan$T_adopt
  idx_tr  <- pan$idx_treat
  idx_co  <- pan$idx_control

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
    if (N_co_g < 2L || T_pre_g < 1L) {
      warning(sprintf(
        "SI staggered: cohort g=%d skipped (T_pre=%d, N_co=%d).",
        g, T_pre_g, N_co_g), call. = FALSE)
      return(NULL)
    }

    pre_rows  <- seq_len(T_pre_g)
    post_rows <- (T_pre_g + 1L):TT

    Y_pre_co_g  <- Y[pre_rows,  idx_co_g, drop = FALSE]
    Y_post_co_g <- Y[post_rows, idx_co_g, drop = FALSE]
    Y_pre_tr_g  <- Y[pre_rows,  idx_tr_g, drop = FALSE]
    Y_treat_g   <- Y[, idx_tr_g, drop = FALSE]

    k_g <- if (!is.null(k)) {
      max(1L, min(as.integer(k), min(T_pre_g, N_co_g) - 1L))
    } else {
      max(1L, floor(sqrt(min(T_pre_g, N_co_g))))
    }

    res <- tryCatch(
      si_pcr_cpp(Y_pre_co_g, Y_post_co_g, Y_pre_tr_g, k_g),
      error = function(e) {
        warning(sprintf("SI staggered: cohort g=%d failed: %s",
                        g, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(res)) return(NULL)

    W_g       <- res$W
    Y_cf_g    <- res$Y_hat
    Y_synth_g <- rbind(Y_pre_co_g %*% W_g, Y_cf_g)
    gap_g     <- Y_treat_g - Y_synth_g
    att_g     <- mean(gap_g[post_rows, , drop = FALSE])

    list(
      cohort    = g,
      n_treated = length(idx_tr_g),
      T_pre     = T_pre_g,
      T_post    = T_post_g,
      estimate  = att_g,
      weight    = as.numeric(length(idx_tr_g)) * T_post_g,
      idx_tr    = idx_tr_g,
      idx_co    = idx_co_g,
      k         = k_g,
      weights   = W_g,
      Y_cf      = Y_cf_g,
      Y_synth   = Y_synth_g,
      Y_treat   = Y_treat_g,
      gap       = gap_g
    )
  })

  valid <- !vapply(cohort_list, is.null, logical(1L))
  if (!any(valid)) stop("All cohort-level SI fits failed.", call. = FALSE)
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

  list(estimate = att, cohort_estimates = cohort_df, cohort_fits = r_list, k = k)
}

# ── Internal: staggered + multi-arm SI (K>1 arms, staggered adoption) ────────
# Combines Phase 16b (staggered) and Phase 21 (multi-arm).
# Within each cohort g, the SVD basis (from Y_pre_co_g) is shared across all
# arms because unit factors v_il are arm-invariant (Agarwal et al. 2025 §2).
# cohort_fits is a flat list of (cohort, arm) cells — each compatible with
# .refit_si_cohort() since they carry T_pre, idx_co, idx_tr, k.
.fit_si_staggered_multi <- function(tensor, k = NULL, control_group = "clean") {
  Y          <- tensor$Y
  TT         <- nrow(Y)
  T_adopt    <- tensor$T_adopt
  idx_tr_all <- tensor$idx_treat
  idx_co_all <- tensor$idx_control
  treat_arms <- tensor$arm_levels[tensor$arm_levels != 0L]

  cohorts <- sort(unique(T_adopt[idx_tr_all]))
  cf_list <- list()

  for (g in cohorts) {
    idx_tr_g <- idx_tr_all[!is.na(T_adopt[idx_tr_all]) & T_adopt[idx_tr_all] == g]
    T_pre_g  <- g - 1L
    T_post_g <- TT - T_pre_g

    if (control_group == "never_treated") {
      idx_co_g <- idx_co_all
    } else {
      future_tr <- idx_tr_all[!is.na(T_adopt[idx_tr_all]) & T_adopt[idx_tr_all] > g]
      idx_co_g  <- c(idx_co_all, future_tr)
    }

    N_co_g <- length(idx_co_g)
    if (N_co_g < 2L || T_pre_g < 1L) {
      warning(sprintf(
        "SI staggered-multi: cohort g=%d skipped (T_pre=%d, N_co=%d).",
        g, T_pre_g, N_co_g), call. = FALSE)
      next
    }

    pre_rows  <- seq_len(T_pre_g)
    post_rows <- (T_pre_g + 1L):TT

    Y_pre_co_g  <- Y[pre_rows,  idx_co_g, drop = FALSE]
    Y_post_co_g <- Y[post_rows, idx_co_g, drop = FALSE]

    k_g <- if (!is.null(k)) {
      max(1L, min(as.integer(k), min(T_pre_g, N_co_g) - 1L))
    } else {
      max(1L, floor(sqrt(min(T_pre_g, N_co_g))))
    }

    for (a in as.character(treat_arms)) {
      idx_tr_gd <- intersect(idx_tr_g, tensor$idx_by_arm[[a]])
      if (length(idx_tr_gd) == 0L) next

      Y_pre_tr_gd <- Y[pre_rows,  idx_tr_gd, drop = FALSE]
      Y_treat_gd  <- Y[, idx_tr_gd, drop = FALSE]

      res <- tryCatch(
        si_pcr_cpp(Y_pre_co_g, Y_post_co_g, Y_pre_tr_gd, k_g),
        error = function(e) {
          warning(sprintf("SI staggered-multi: cohort g=%d arm %s failed: %s",
                          g, a, conditionMessage(e)), call. = FALSE)
          NULL
        }
      )
      if (is.null(res)) next

      W_gd       <- res$W
      Y_cf_gd    <- res$Y_hat
      Y_synth_gd <- rbind(Y_pre_co_g %*% W_gd, Y_cf_gd)
      gap_gd     <- Y_treat_gd - Y_synth_gd
      att_gd     <- mean(gap_gd[post_rows, , drop = FALSE], na.rm = TRUE)

      cf_list <- c(cf_list, list(list(
        cohort    = g,
        arm       = as.integer(a),
        n_treated = length(idx_tr_gd),
        T_pre     = T_pre_g,
        T_post    = T_post_g,
        estimate  = att_gd,
        weight    = as.numeric(length(idx_tr_gd)) * T_post_g,
        idx_tr    = idx_tr_gd,
        idx_co    = idx_co_g,
        k         = k_g,
        weights   = W_gd,
        Y_cf      = Y_cf_gd,
        Y_synth   = Y_synth_gd,
        Y_treat   = Y_treat_gd,
        gap       = gap_gd
      )))
    }
  }

  if (length(cf_list) == 0L)
    stop("All (cohort, arm) SI fits failed.", call. = FALSE)

  ws      <- vapply(cf_list, `[[`, numeric(1L), "weight")
  taus    <- vapply(cf_list, `[[`, numeric(1L), "estimate")
  W_total <- sum(ws)
  att     <- sum(ws * taus) / W_total

  arm_estimates <- setNames(
    vapply(as.character(treat_arms), function(a) {
      mask <- vapply(cf_list, function(cf) cf$arm == as.integer(a), logical(1L))
      if (!any(mask) || sum(ws[mask]) == 0) return(NA_real_)
      sum(ws[mask] * taus[mask]) / sum(ws[mask])
    }, numeric(1L)),
    as.character(treat_arms)
  )

  cohort_arm_df <- data.frame(
    cohort    = vapply(cf_list, `[[`, integer(1L),  "cohort"),
    arm       = vapply(cf_list, `[[`, integer(1L),  "arm"),
    n_treated = vapply(cf_list, `[[`, integer(1L),  "n_treated"),
    T_pre     = vapply(cf_list, `[[`, integer(1L),  "T_pre"),
    T_post    = vapply(cf_list, `[[`, integer(1L),  "T_post"),
    estimate  = taus,
    weight    = ws / W_total,
    stringsAsFactors = FALSE
  )

  list(
    method               = "si",
    staggered            = TRUE,
    multi_arm            = TRUE,
    arm_levels           = treat_arms,
    k                    = k,
    estimate             = att,
    arm_estimates        = arm_estimates,
    cohort_arm_estimates = cohort_arm_df,
    cohort_fits          = cf_list,
    Y_all                = Y,
    idx_co               = idx_co_all,
    times                = tensor$times,
    T_pre                = tensor$T_pre,
    Y_treat              = NULL,
    Y_synth              = NULL,
    gap                  = NULL,
    unit_weights         = NULL
  )
}

#' Fit Synthetic Interventions via PCR (Agarwal et al. 2025)
#'
#' Estimates counterfactual outcomes using principal component regression
#' on the pre-treatment control panel. Assumes a low-rank tensor factor model
#' for potential outcomes:
#'   Y_ti^(d) = sum_l u_tl * v_il * lambda_dl + eps_ti^(d)
#'
#' For the standard binary-treatment case (d in {0, 1}), this simplifies to:
#' find weights over control units in the pre-period using k-component PCR,
#' then apply those weights to control post-period outcomes as the counterfactual.
#'
#' Unlike SCM and SDID, SI-PCR weights are unconstrained (can be negative or
#' sum to values other than 1). The `plot(type = "weights")` view shows only
#' positive weights (> 1e-4).
#'
#' @param y    Outcome vector (long format)
#' @param d    Treatment arm indicator (integer, 0 = control, 1 = treated)
#' @param id   Unit identifier (long format)
#' @param time Time identifier (long format)
#' @param k    Number of SVD components for PCR. Default:
#'             `floor(sqrt(min(T_pre, N_co)))`, at least 1.
#' @return A list of class `coresynth` with:
#'   * `method`: `"si"`
#'   * `estimate`: ATT (average treatment effect on the treated)
#'   * `k`: number of SVD components used
#'   * `weights`: Donor weight matrix (N_co x N_tr)
#'   * `unit_weights`: Average donor weights across treated units (named N_co vector)
#'   * `Y_cf`: Counterfactual post-treatment outcomes (T_post x N_tr)
#'   * `Y_treat`: Observed treated outcomes (T x N_tr)
#'   * `Y_synth`: Full synthetic series (T x N_tr; pre = weighted control, post = Y_cf)
#'   * `gap`: Treatment effect series (T x N_tr)
#'   * `times`: Time index vector
#'   * `T_pre`: Number of pre-treatment periods
#' @noRd
fit_si_cpp <- function(y, d, id, time, k = NULL,
                       control_group = c("clean", "never_treated"), ...) {
  control_group <- match.arg(control_group)

  # ── Multi-arm detection: max(d) > 1 means K≥2 treatment arms ────────────
  d_int <- as.integer(d)
  if (max(d_int, na.rm = TRUE) > 1L) {
    tensor <- panel_to_tensor(y, d_int, id, time)
    k_use  <- if (is.null(k)) NULL else as.integer(k)
    if (!tensor$is_sharp) {
      return(.fit_si_staggered_multi(tensor, k_use, control_group))
    }
    return(.fit_si_multi(tensor$Y, tensor$idx_by_arm, tensor$arm_levels,
                         tensor$T_pre, k_use, tensor$times))
  }

  pan    <- panel_to_matrices(y, d_int, id, time)

  # ── Staggered adoption path ──────────────────────────────────────────────
  if (!pan$is_sharp) {
    res_st <- .fit_si_staggered(pan, pan$Y, k, control_group)
    res_st$method       <- "si"
    res_st$staggered    <- TRUE
    res_st$times        <- pan$times
    res_st$T_pre        <- pan$T_pre
    res_st$Y_treat      <- pan$Y[, pan$idx_treat, drop = FALSE]
    res_st$Y_synth      <- NULL
    res_st$gap          <- NULL
    res_st$unit_weights <- NULL
    res_st$Y_all        <- pan$Y             # T × N — needed by si_inference()
    res_st$idx_co       <- pan$idx_control
    res_st$idx_tr       <- pan$idx_treat
    return(res_st)
  }

  Y      <- pan$Y
  T_pre  <- pan$T_pre
  idx_tr <- pan$idx_treat
  idx_co <- pan$idx_control

  if (length(idx_co) < 2)
    stop("SI requires at least two control units.")
  if (T_pre < 1)
    stop("SI requires at least one pre-treatment period.")

  TT     <- nrow(Y)
  T_post <- TT - T_pre
  N_co   <- length(idx_co)

  if (T_post < 1)
    stop("No post-treatment periods found.")

  pre_rows  <- seq_len(T_pre)
  post_rows <- (T_pre + 1L):TT

  Y_pre_co  <- Y[pre_rows,  idx_co, drop = FALSE]  # T_pre  x N_co
  Y_post_co <- Y[post_rows, idx_co, drop = FALSE]  # T_post x N_co
  Y_pre_tr  <- Y[pre_rows,  idx_tr, drop = FALSE]  # T_pre  x N_tr

  if (is.null(k)) {
    k <- max(1L, floor(sqrt(min(T_pre, N_co))))
  }
  k <- as.integer(k)

  if (k > min(T_pre, N_co))
    stop(sprintf("k (%d) must be <= min(T_pre=%d, N_co=%d).", k, T_pre, N_co))

  res  <- si_pcr_cpp(Y_pre_co, Y_post_co, Y_pre_tr, k)
  W    <- res$W      # N_co x N_tr
  Y_cf <- res$Y_hat  # T_post x N_tr

  # Full synthetic series for plotting: pre = Y_pre_co %*% W, post = Y_cf
  Y_synth_full <- rbind(Y_pre_co %*% W, Y_cf)  # T x N_tr

  Y_treat <- Y[, idx_tr, drop = FALSE]          # T x N_tr
  gap     <- Y_treat - Y_synth_full
  att     <- mean(gap[post_rows, , drop = FALSE], na.rm = TRUE)

  unit_w        <- rowMeans(W)
  names(unit_w) <- colnames(Y_pre_co)

  list(
    method       = "si",
    k            = k,
    weights      = W,
    unit_weights = unit_w,
    Y_cf         = Y_cf,
    Y_treat      = Y_treat,
    Y_synth      = Y_synth_full,
    gap          = gap,
    times        = pan$times,
    T_pre        = T_pre,
    estimate     = att,
    Y_pre_co     = Y_pre_co,    # T_pre  × N_co — needed by si_inference()
    Y_post_co    = Y_post_co    # T_post × N_co — needed by si_inference()
  )
}

# ── Internal helper: re-estimate one SI cohort with given control indices ──────
.refit_si_cohort <- function(Y_all, cf) {
  TT        <- nrow(Y_all)
  T_pre_g   <- cf$T_pre
  pre_rows  <- seq_len(T_pre_g)
  post_rows <- (T_pre_g + 1L):TT

  Y_pre_co_g  <- Y_all[pre_rows,  cf$idx_co, drop = FALSE]
  Y_post_co_g <- Y_all[post_rows, cf$idx_co, drop = FALSE]
  Y_pre_tr_g  <- Y_all[pre_rows,  cf$idx_tr, drop = FALSE]

  res <- tryCatch(
    si_pcr_cpp(Y_pre_co_g, Y_post_co_g, Y_pre_tr_g, cf$k),
    error = function(e) NULL
  )
  if (is.null(res)) return(NA_real_)

  Y_treat_post <- Y_all[post_rows, cf$idx_tr, drop = FALSE]
  mean(Y_treat_post - res$Y_hat)
}

# ── Internal: inference for staggered + multi-arm SI (Phase 23) ──────────────
# Bootstrap: per-cohort resampling — same resampled controls applied to all arms
#   within a cohort (shared SVD basis requires shared donor bootstrap indices).
# Jackknife: per-cohort LOO refitting all arms simultaneously + delta-method.
# jackknife_global: drop one unique control from ALL cohorts and ALL arms.
.si_inference_staggered_multi <- function(fit, method, n_boot, level,
                                          alternative, seed) {
  if (!is.null(seed)) set.seed(seed)
  alpha      <- 1 - level
  tau_hat    <- fit$estimate
  Y_all      <- fit$Y_all
  cf_list    <- fit$cohort_fits
  ws         <- vapply(cf_list, `[[`, numeric(1L), "weight")
  W_total    <- sum(ws)
  cohort_ids <- vapply(cf_list, `[[`, integer(1L), "cohort")
  uc         <- sort(unique(cohort_ids))

  .p_from_z <- function(z) switch(alternative,
    two.sided = 2 * pnorm(-abs(z)),
    greater   = pnorm(-z),
    less      = pnorm(z))

  .n_co_mean <- round(mean(vapply(uc, function(g) {
    length(cf_list[cohort_ids == g][[1L]]$idx_co)
  }, numeric(1L))))

  if (method == "bootstrap") {
    .boot_one <- function() {
      att_total <- 0
      for (g in uc) {
        cells    <- cf_list[cohort_ids == g]
        idx_co_g <- cells[[1L]]$idx_co
        boot_idx <- idx_co_g[sample(length(idx_co_g), replace = TRUE)]
        for (cf in cells) {
          cf_b <- cf; cf_b$idx_co <- boot_idx
          att_gd <- .refit_si_cohort(Y_all, cf_b)
          att_total <- att_total + cf$weight * att_gd
        }
      }
      att_total / W_total
    }
    boot_ests <- replicate(n_boot, .boot_one())
    n_fail <- sum(is.na(boot_ests))
    if (n_fail > 0.1 * n_boot)
      warning(sprintf("%.0f%% of staggered-multi bootstrap replications failed.",
                      100 * n_fail / n_boot))
    valid    <- boot_ests[!is.na(boot_ests)]
    se       <- sd(valid)
    ci_lower <- unname(quantile(valid, alpha / 2))
    ci_upper <- unname(quantile(valid, 1 - alpha / 2))
    z        <- tau_hat / max(se, .Machine$double.eps)
    return(structure(list(
      estimate    = tau_hat, se = se, p_value = .p_from_z(z),
      ci_lower    = ci_lower, ci_upper = ci_upper,
      method      = method, staggered = TRUE, n_controls = .n_co_mean,
      alternative = alternative, boot_ests = boot_ests
    ), class = "coresynth_inference"))
  }

  if (method == "jackknife_global") {
    all_co    <- sort(unique(unlist(lapply(cf_list, `[[`, "idx_co"))))
    orig_ests <- vapply(cf_list, `[[`, numeric(1L), "estimate")
    jack_ests <- vapply(all_co, function(i) {
      att_total <- 0
      for (ki in seq_along(cf_list)) {
        cf <- cf_list[[ki]]
        if (i %in% cf$idx_co) {
          keep <- cf$idx_co[cf$idx_co != i]
          if (length(keep) < 2L) return(NA_real_)
          cf_loo <- cf; cf_loo$idx_co <- keep
          att_gd <- .refit_si_cohort(Y_all, cf_loo)
          if (is.na(att_gd)) return(NA_real_)
        } else {
          att_gd <- orig_ests[ki]
        }
        att_total <- att_total + cf$weight * att_gd
      }
      att_total / W_total
    }, numeric(1L))
    valid <- jack_ests[!is.na(jack_ests)]
    N_v   <- length(valid)
    se    <- sqrt((N_v - 1L) / N_v * sum((valid - mean(valid))^2))
    z     <- tau_hat / max(se, .Machine$double.eps)
    ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
    ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se
    return(structure(list(
      estimate    = tau_hat, se = se, p_value = .p_from_z(z),
      ci_lower    = ci_lower, ci_upper = ci_upper,
      method      = method, staggered = TRUE, n_controls = length(all_co),
      alternative = alternative, boot_ests = NULL
    ), class = "coresynth_inference"))
  }

  # jackknife: per-cohort LOO refitting all arms simultaneously + delta-method
  var_components <- vapply(uc, function(g) {
    cells    <- cf_list[cohort_ids == g]
    idx_co_g <- cells[[1L]]$idx_co
    N_co_g   <- length(idx_co_g)
    if (N_co_g < 2L) return(0)
    jack_g <- vapply(seq_len(N_co_g), function(i) {
      keep <- idx_co_g[-i]
      cohort_agg <- 0
      for (cf in cells) {
        cf_loo <- cf; cf_loo$idx_co <- keep
        att_gd_loo <- .refit_si_cohort(Y_all, cf_loo)
        if (is.na(att_gd_loo)) return(NA_real_)
        cohort_agg <- cohort_agg + cf$weight * att_gd_loo
      }
      cohort_agg
    }, numeric(1L))
    valid_g <- jack_g[!is.na(jack_g)]
    if (length(valid_g) < 2L) return(0)
    (N_co_g - 1L) / N_co_g * sum((valid_g - mean(valid_g))^2)
  }, numeric(1L))

  w_g_totals <- vapply(uc, function(g) {
    sum(vapply(cf_list[cohort_ids == g], `[[`, numeric(1L), "weight"))
  }, numeric(1L))
  var_st   <- sum((w_g_totals / W_total)^2 * var_components)
  se       <- sqrt(var_st)
  z        <- tau_hat / max(se, .Machine$double.eps)
  ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
  ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se
  structure(list(
    estimate    = tau_hat, se = se, p_value = .p_from_z(z),
    ci_lower    = ci_lower, ci_upper = ci_upper,
    method      = method, staggered = TRUE, n_controls = .n_co_mean,
    alternative = alternative, boot_ests = NULL
  ), class = "coresynth_inference")
}

# ── Internal: inference for multi-arm SI ──────────────────────────────────────
# Bootstrap: resample control columns (shared across all arms, since arms share
#   the same donor pool). Jackknife: LOO over control columns.
.si_inference_multi <- function(fit, method, n_boot, level, alternative, seed) {
  if (!is.null(seed)) set.seed(seed)
  alpha     <- 1 - level
  tau_hat   <- fit$estimate
  Y_pre_co  <- fit$Y_pre_co
  Y_post_co <- fit$Y_post_co
  k         <- fit$k
  arm_fits  <- fit$arm_fits
  T_pre     <- fit$T_pre
  N_co      <- ncol(Y_pre_co)

  .refit_multi <- function(co_idx) {
    Ypc  <- Y_pre_co[,  co_idx, drop = FALSE]
    Ypoc <- Y_post_co[, co_idx, drop = FALSE]
    post_rows_local <- (T_pre + 1L):nrow(arm_fits[[1L]]$Y_treat)
    ws <- taus <- numeric(length(arm_fits))
    for (i in seq_along(arm_fits)) {
      af  <- arm_fits[[i]]
      Ytr_pre <- af$Y_treat[seq_len(T_pre), , drop = FALSE]
      res <- tryCatch(si_pcr_cpp(Ypc, Ypoc, Ytr_pre, k),
                      error = function(e) NULL)
      ws[i] <- af$weight
      if (is.null(res)) {
        taus[i] <- af$estimate
      } else {
        taus[i] <- mean(
          af$Y_treat[post_rows_local, , drop = FALSE] - res$Y_hat,
          na.rm = TRUE
        )
      }
    }
    sum(ws * taus) / sum(ws)
  }

  if (method == "bootstrap") {
    boot_ests <- replicate(n_boot, .refit_multi(sample(N_co, replace = TRUE)))
    n_fail <- sum(is.na(boot_ests))
    if (n_fail > 0.1 * n_boot)
      warning(sprintf("%.0f%% of multi-arm bootstrap replications failed.",
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
    return(structure(list(
      estimate    = tau_hat, se = se, p_value = p_value,
      ci_lower    = ci_lower, ci_upper = ci_upper,
      method      = method, staggered = FALSE, n_controls = N_co,
      alternative = alternative, boot_ests = boot_ests
    ), class = "coresynth_inference"))
  }

  # jackknife
  jack_ests <- vapply(seq_len(N_co),
                      function(i) .refit_multi(seq_len(N_co)[-i]),
                      numeric(1L))
  jack_var <- (N_co - 1) / N_co * sum((jack_ests - mean(jack_ests, na.rm = TRUE))^2)
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
    estimate    = tau_hat, se = se, p_value = p_value,
    ci_lower    = ci_lower, ci_upper = ci_upper,
    method      = method, staggered = FALSE, n_controls = N_co,
    alternative = alternative, boot_ests = NULL
  ), class = "coresynth_inference")
}

#' Non-parametric Inference for SI (Agarwal et al. 2025)
#'
#' Estimates SE and confidence intervals for the ATT via non-parametric cluster
#' bootstrap or jackknife over control units. Works for both sharp and staggered
#' SI fits. For staggered fits, bootstrap resamples each cohort's control pool
#' independently, and jackknife uses a per-cohort LOO with delta-method variance
#' aggregation.
#'
#' @param fit   A `coresynth` object from [scm_fit()] with `method = "si"`.
#' @param method `"bootstrap"` (default) or `"jackknife"`.
#' @param n_boot Number of bootstrap replications (default 499L; ignored for jackknife).
#' @param level  Confidence level (default 0.95).
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param seed  RNG seed for reproducibility (default NULL).
#' @return A list of class `coresynth_inference`.
#' @export
si_inference <- function(
  fit,
  method      = c("bootstrap", "jackknife", "jackknife_global"),
  n_boot      = 499L,
  level       = 0.95,
  alternative = c("two.sided", "greater", "less"),
  seed        = NULL
) {
  method      <- match.arg(method)
  alternative <- match.arg(alternative)

  if (!inherits(fit, "coresynth") || !identical(fit$method, "si"))
    stop("si_inference() requires a coresynth fit with method = 'si'.",
         call. = FALSE)

  tau_hat   <- fit$estimate
  alpha     <- 1 - level
  staggered <- isTRUE(fit$staggered)

  # ── Multi-arm path ────────────────────────────────────────────────────────
  if (isTRUE(fit$multi_arm)) {
    if (isTRUE(fit$staggered)) {
      return(.si_inference_staggered_multi(fit, method, n_boot, level, alternative, seed))
    }
    if (method == "jackknife_global")
      stop("jackknife_global は multi-arm fit には未対応。jackknife を使用してください。",
           call. = FALSE)
    return(.si_inference_multi(fit, method, n_boot, level, alternative, seed))
  }

  if (!is.null(seed)) set.seed(seed)

  # ── Sharp path ──────────────────────────────────────────────────────────────
  if (!staggered) {
    if (method == "jackknife_global")
      stop("si_inference() method='jackknife_global' requires a staggered fit.",
           call. = FALSE)
    if (is.null(fit$Y_pre_co))
      stop("fit does not contain Y_pre_co. Re-estimate with the current version.",
           call. = FALSE)

    Y_pre_co  <- fit$Y_pre_co
    Y_post_co <- fit$Y_post_co
    Y_pre_tr  <- fit$Y_treat[seq_len(fit$T_pre), , drop = FALSE]
    k         <- fit$k
    N_co      <- ncol(Y_pre_co)
    post_rows <- (fit$T_pre + 1L):nrow(fit$Y_treat)

    .refit_sharp <- function(co_idx) {
      res <- tryCatch(
        si_pcr_cpp(Y_pre_co[, co_idx, drop = FALSE],
                   Y_post_co[, co_idx, drop = FALSE],
                   Y_pre_tr, k),
        error = function(e) NULL
      )
      if (is.null(res)) return(NA_real_)
      mean(fit$Y_treat[post_rows, , drop = FALSE] - res$Y_hat)
    }

    if (method == "bootstrap") {
      boot_ests <- replicate(n_boot, {
        .refit_sharp(sample(N_co, replace = TRUE))
      })
      n_fail <- sum(is.na(boot_ests))
      if (n_fail > 0.1 * n_boot)
        warning(sprintf("%.0f%% of bootstrap replications failed.", 100 * n_fail / n_boot))
      valid <- boot_ests[!is.na(boot_ests)]
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
  cohort_fits <- fit$cohort_fits
  w_g         <- vapply(cohort_fits, `[[`, numeric(1L), "weight")  # unnorm
  W           <- sum(w_g)

  .boot_one <- function() {
    att_total <- 0
    for (cf in cohort_fits) {
      N_co_g <- length(cf$idx_co)
      idx_b  <- cf$idx_co[sample(N_co_g, replace = TRUE)]
      cf_b   <- cf
      cf_b$idx_co <- idx_b
      att_g_b <- .refit_si_cohort(Y_all, cf_b)
      att_total <- att_total + cf$weight * att_g_b
    }
    att_total / W
  }

  if (method == "bootstrap") {
    boot_ests <- replicate(n_boot, .boot_one())
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
    n_co_mean <- round(mean(vapply(cohort_fits, function(cf) length(cf$idx_co), numeric(1L))))
    return(structure(list(
      estimate    = tau_hat, se = se, p_value = p_value,
      ci_lower    = ci_lower, ci_upper = ci_upper,
      method      = method, staggered = TRUE, n_controls = n_co_mean,
      alternative = alternative, boot_ests = boot_ests
    ), class = "coresynth_inference"))
  }

  if (method == "jackknife_global") {
    # Global jackknife: drop one unique control across ALL cohorts simultaneously
    # (Phase 18b). Captures cross-cohort correlation ignored by per-cohort LOO.
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
          att_g <- .refit_si_cohort(Y_all, cf_loo)
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
      .refit_si_cohort(Y_all, cf_loo)
    }, numeric(1))
    valid_g <- jack_g[!is.na(jack_g)]
    if (length(valid_g) < 2L) return(0)
    (N_co_g - 1) / N_co_g * sum((valid_g - mean(valid_g))^2)
  }, numeric(1))

  var_st  <- sum((w_g / W)^2 * var_components)
  se      <- sqrt(var_st)
  z       <- tau_hat / max(se, .Machine$double.eps)
  ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
  ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se
  p_value  <- switch(alternative,
    two.sided = 2 * pnorm(-abs(z)),
    greater   = pnorm(-z),
    less      = pnorm(z)
  )
  n_co_mean <- round(mean(vapply(cohort_fits, function(cf) length(cf$idx_co), numeric(1L))))
  structure(list(
    estimate    = tau_hat, se = se, p_value = p_value,
    ci_lower    = ci_lower, ci_upper = ci_upper,
    method      = method, staggered = TRUE, n_controls = n_co_mean,
    alternative = alternative, boot_ests = NULL
  ), class = "coresynth_inference")
}

#' @export
print.coresynth_inference <- function(x, digits = 4, ...) {
  cat(sprintf("Inference (method=%s, alternative=%s%s)\n",
              x$method, x$alternative,
              if (x$staggered) ", staggered" else ""))
  cat("  Estimate  :", round(x$estimate, digits), "\n")
  cat("  SE        :", round(x$se,       digits), "\n")
  cat("  p-value   :", round(x$p_value,  digits), "\n")
  cat("  CI        : [", round(x$ci_lower, digits), ", ",
      round(x$ci_upper, digits), "]\n", sep = "")
  cat("  Controls  :", x$n_controls, "\n")
  invisible(x)
}
