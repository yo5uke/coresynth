# -- Internal: multi-arm SI for K>1 treatment arms (Agarwal et al. 2025) -----
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
    stop("SI multi-arm: at least 2 control units required.", call. = FALSE)
  if (T_pre < 1L)
    stop("SI multi-arm: T_pre must be > 0.", call. = FALSE)
  if (T_post < 1L)
    stop("SI multi-arm: T_post must be > 0.", call. = FALSE)

  # Common donor matrices shared across all arms (unit factors are arm-invariant)
  Y_pre_co  <- Y[pre_rows,  idx_co, drop = FALSE]  # T_pre  x N_co
  Y_post_co <- Y[post_rows, idx_co, drop = FALSE]  # T_post x N_co

  if (is.null(k)) k <- .si_svht_k(Y_pre_co)
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

  # Weighted-average ATT: weight proportional to N_tr_d * T_post (Clarke et al. 2024)
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

# -- Internal: data-driven PCR rank via the optimal hard threshold -----------
# Agarwal et al. name cross-validation, universal singular value thresholding
# and the elbow rule as the data-driven ways to choose k. This is the second:
# Gavish & Donoho (2014) give the threshold that minimises asymptotic MSE when
# recovering a low-rank matrix from noise of unknown level, tau = omega(beta) *
# median(singular values) with beta = min(dim)/max(dim). Keeping the singular
# values above it is parameter-free and costs one SVD, where cross-validating k
# would cost a PCR fit per candidate.
#
# omega(beta) uses the cubic approximation of Gavish & Donoho (2014, Eq. 5),
# which they report as within 0.02 of the exact coefficient for beta in
# (0.001, 1]. The exact value has no closed form; at a coefficient near 2 the
# gap moves the threshold by about 1%, which can only change k when a
# singular value sits that close to it.
.si_svht_k <- function(Y_pre_co) {
  d <- svd(Y_pre_co, nu = 0, nv = 0)$d
  d <- d[is.finite(d)]
  if (length(d) == 0L || max(d) <= 0) return(1L)

  beta  <- min(dim(Y_pre_co)) / max(dim(Y_pre_co))
  omega <- 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta + 1.43
  k     <- sum(d > omega * stats::median(d))
  max(1L, min(as.integer(k), length(d)))
}

# The rank a staggered fit used, one entry per cohort. Unlike the GSC factor
# count and the MC penalty, which one global search fixes for the whole panel,
# the threshold is applied to each cohort's own donor block -- and even an
# explicitly supplied k is clipped to that cohort's dimensions -- so there is
# no single value to report. `cells` are the cohort fits, or the (cohort, arm)
# cells of a multi-arm fit, where a cohort's arms share the SVD basis and
# therefore its k.
.si_k_by_cohort <- function(cells) {
  g    <- vapply(cells, function(z) as.integer(z$cohort), integer(1L))
  kk   <- vapply(cells, function(z) as.integer(z$k), integer(1L))
  keep <- !duplicated(g)
  stats::setNames(kk[keep], as.character(g[keep]))
}

# -- Internal: cohort-by-cohort SI for staggered adoption --------------------
# Clarke et al. (2024); Arkhangelsky et al. (2021) Appendix S.8.
# The SVD sees only the pre-treatment window, but the counterfactual is
# projected over every post period, so a donor must stay untreated through
# the end of the panel; see .check_control_group() for the deprecated "clean".
.fit_si_staggered <- function(pan, Y, k = NULL, control_group = "never_treated") {
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
      max(1L, min(.si_svht_k(Y_pre_co_g), min(T_pre_g, N_co_g) - 1L))
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
  if (!any(valid)) {
    stop("All cohort-level SI fits failed; see the preceding warnings for ",
         "per-cohort reasons.", call. = FALSE)
  }
  cl <- cohort_list[valid]
  c(.aggregate_cohorts(cl), list(k = .si_k_by_cohort(cl)))
}

# -- Internal: staggered + multi-arm SI (K>1 arms, staggered adoption) -------
# Handles the combination of staggered adoption and multiple treatment arms.
# Within each cohort g, the SVD basis (from Y_pre_co_g) is shared across all
# arms because unit factors v_il are arm-invariant (Agarwal et al. 2025 S.2).
# cohort_fits is a flat list of (cohort, arm) cells, each carrying T_pre,
# idx_co, idx_tr and k, which si_inference() reads to refit them.
.fit_si_staggered_multi <- function(tensor, k = NULL, control_group = "never_treated") {
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
      max(1L, min(.si_svht_k(Y_pre_co_g), min(T_pre_g, N_co_g) - 1L))
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
    stop("All (cohort, arm) SI fits failed; see the preceding warnings for ",
         "per-cell reasons.", call. = FALSE)

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
    k                    = .si_k_by_cohort(cf_list),
    estimate             = att,
    arm_estimates        = arm_estimates,
    cohort_arm_estimates = cohort_arm_df,
    cohort_fits          = cf_list,
    Y_all                = Y,
    idx_co               = idx_co_all,
    control_group        = control_group,
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
#' @param k    Number of SVD components for PCR. `NULL` (default) applies the
#'             optimal hard threshold of Gavish & Donoho (2014) to the
#'             pre-treatment donor block. On a staggered panel the threshold
#'             is applied to each cohort's own donor block, and a supplied `k`
#'             is clipped to that cohort's dimensions, so the fit reports one
#'             value per cohort.
#' @return A list of class `coresynth` with:
#'   * `method`: `"si"`
#'   * `estimate`: ATT (average treatment effect on the treated)
#'   * `k`: number of SVD components used; on a staggered panel a named
#'     integer vector, one entry per cohort
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
                       control_group = c("never_treated", "clean"), ...) {
  control_group <- match.arg(control_group)

  # -- Multi-arm detection: max(d) > 1 means K>=2 treatment arms -------------
  d_int <- as.integer(d)
  if (any(d_int > 1L, na.rm = TRUE)) {
    tensor <- panel_to_tensor(y, d_int, id, time)
    .check_panel_complete(tensor$Y, "SI")
    .check_control_group(control_group, tensor, "SI")
    k_use  <- if (is.null(k)) NULL else as.integer(k)
    if (!tensor$is_sharp) {
      return(.fit_si_staggered_multi(tensor, k_use, control_group))
    }
    return(.fit_si_multi(tensor$Y, tensor$idx_by_arm, tensor$arm_levels,
                         tensor$T_pre, k_use, tensor$times))
  }

  pan    <- panel_to_matrices(y, d_int, id, time)
  .check_panel_complete(pan$Y, "SI")
  .check_control_group(control_group, pan, "SI")

  # -- Staggered adoption path -------------------------------------------------
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
    res_st$Y_all        <- pan$Y             # T x N -- needed by si_inference()
    res_st$idx_co       <- pan$idx_control
    res_st$idx_tr       <- pan$idx_treat
    res_st$control_group <- control_group
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

  if (is.null(k)) k <- .si_svht_k(Y_pre_co)
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
    Y_pre_co     = Y_pre_co,    # T_pre  x N_co -- needed by si_inference()
    Y_post_co    = Y_post_co    # T_post x N_co -- needed by si_inference()
  )
}

# -- Internal: SI effect for a draw of units ----------------------------------
# The estimator of fit_si_cpp() applied to never-treated columns `co` and
# treated columns `tr` of ctx$Y; see .unit_inference(). Cohorts are rebuilt
# from the treated units drawn. Within a cohort every arm shares the SVD
# basis of the donor block, so one PCR fit covers all of its treated
# columns; the rank of each cohort is the one the fit used. Since T_post is
# common within a cohort, its arm-weighted effect is the mean over its
# treated units, and the cohorts combine with the Clarke et al. (2024)
# weights.
.si_unit_estimate <- function(ctx, co, tr) {
  TT      <- nrow(ctx$Y)
  adopt   <- ctx$T_adopt[tr]
  cohorts <- sort(unique(adopt))
  est <- w <- numeric(length(cohorts))
  for (k in seq_along(cohorts)) {
    g      <- cohorts[k]
    tr_g   <- tr[adopt == g]
    pre    <- seq_len(g - 1L)
    post   <- g:TT
    donors <- if (identical(ctx$control_group, "never_treated")) co else
      c(co, tr[adopt > g])
    res <- tryCatch(
      si_pcr_cpp(ctx$Y[pre, donors, drop = FALSE],
                 ctx$Y[post, donors, drop = FALSE],
                 ctx$Y[pre, tr_g, drop = FALSE],
                 ctx$k[[as.character(g)]]),
      error = function(e) NULL)
    if (is.null(res)) return(NULL)
    est[k] <- mean(ctx$Y[post, tr_g, drop = FALSE] - res$Y_hat)
    w[k]   <- length(tr_g) * length(post)
  }
  list(estimate = sum(w * est) / sum(w), cohort = cohorts,
       cohort_estimate = est, cohort_weight = w)
}

# The panel an SI fit was estimated on, laid out for .si_unit_estimate().
# `strata` carries the arm of each treated column, so the bootstrap keeps
# every arm's size.
.si_inference_context <- function(fit) {
  if (inherits(fit, "coresynth_staggered")) {
    if (is.null(fit$Y_all))
      stop("fit does not contain Y_all. Re-estimate with the current version.",
           call. = FALSE)
    cells <- fit$cohort_fits
    tr    <- unlist(lapply(cells, `[[`, "idx_tr"), use.names = FALSE)
    arm   <- unlist(lapply(cells, function(cf)
      rep(cf[["arm"]] %||% 1L, length(cf$idx_tr))), use.names = FALSE)
    return(list(
      Y = fit$Y_all, k = as.list(.si_k_by_cohort(cells)),
      T_adopt = .adopt_from_cohorts(ncol(fit$Y_all), cells),
      co = fit$idx_co, tr = tr, strata = arm,
      control_group = .fit_control_group(fit)))
  }
  if (is.null(fit$Y_pre_co))
    stop("fit does not contain Y_pre_co. Re-estimate with the current version.",
         call. = FALSE)
  Y_co <- rbind(fit$Y_pre_co, fit$Y_post_co)
  if (inherits(fit, "coresynth_multiarm")) {
    Y_tr <- do.call(cbind, lapply(fit$arm_fits, `[[`, "Y_treat"))
    arm  <- unlist(lapply(fit$arm_fits, function(af)
      rep(af$arm, ncol(af$Y_treat))), use.names = FALSE)
  } else {
    Y_tr <- fit$Y_treat
    arm  <- rep(1L, ncol(Y_tr))
  }
  N_co <- ncol(Y_co)
  N_tr <- ncol(Y_tr)
  list(
    Y = cbind(Y_co, Y_tr),
    k = stats::setNames(list(fit$k), as.character(fit$T_pre + 1L)),
    T_adopt = c(rep(NA_integer_, N_co), rep(fit$T_pre + 1L, N_tr)),
    co = seq_len(N_co), tr = N_co + seq_len(N_tr), strata = arm,
    control_group = "never_treated")
}

#' Non-parametric Inference for SI (Agarwal et al. 2025)
#'
#' Estimates SE and confidence intervals for the ATT by re-estimating the
#' fit on resampled units: a bootstrap of the treated and the control units,
#' drawn separately, or a jackknife that drops one unit at a time. Works for
#' sharp, staggered and multi-arm SI fits. Each refit keeps the rank `k` the
#' fit used.
#'
#' Holding the treated units fixed would leave out their own noise, which is
#' most of the variance when few units are treated, so they are resampled
#' too; with a single treated unit (in any arm) there is nothing to resample
#' and the function stops. For a sharp fit with one treated unit,
#' [conformal_inference()] gives a valid interval.
#'
#' The controls are drawn once per replicate and shared by every arm and
#' every cohort, as the SVD basis is a property of the shared donor pool.
#' Treated units are drawn within their arm. On a staggered fit they are
#' pooled across cohorts, so cohort sizes vary across replicates, and the
#' cohorts are rebuilt from the draw and combined with their
#' `N_treated * T_post` weights. `"jackknife"` drops each cohort's own donors
#' and treated units and combines the per-cohort variances by the delta
#' method, which ignores the covariance that shared donors induce between
#' cohorts; `"jackknife_global"` (staggered fits only) drops every unit from
#' the whole panel in turn and keeps it. Both jackknives drop treated units,
#' so every cohort needs at least two (Clarke et al. 2024, Algorithm 3).
#'
#' @param fit   A `coresynth` object from [scm_fit()] with `method = "si"`.
#' @param method `"bootstrap"` (default), `"jackknife"`, or
#'   `"jackknife_global"` (staggered fits only).
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

  if (!inherits(fit, "coresynth_si"))
    stop("si_inference() requires a coresynth fit with method = 'si'.",
         call. = FALSE)
  if (method == "jackknife_global" && !inherits(fit, "coresynth_staggered"))
    stop("si_inference() method='jackknife_global' requires a staggered fit; ",
         "use jackknife instead.", call. = FALSE)

  ctx <- .si_inference_context(fit)
  if (!is.null(seed)) set.seed(seed)
  .unit_inference(
    fit, ctx, .si_unit_estimate, method, n_boot, level, alternative,
    label = "si_inference()",
    one_treated_hint = paste0(
      "For a sharp fit with one treated unit, conformal_inference() gives ",
      "a valid interval."))
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
