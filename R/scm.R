#' Fit Original Synthetic Control Method (Abadie et al. 2010)
#'
#' Finds donor unit weights W (on the simplex) minimising the V-weighted
#' pre-treatment covariate distance:
#'   min_W (X1 - X0 W)' V (X1 - X0 W)
#' where V is jointly optimised via coordinate descent on the pre-treatment
#' outcome prediction error.
#'
#' @param y    Outcome vector (long format)
#' @param d    Treatment indicator (long format)
#' @param id   Unit identifier (long format)
#' @param time Time identifier (long format)
#' @param v_selection How to select the V metric matrix. `"insample"` (default)
#'   minimises in-sample pre-treatment MSPE following Abadie et al. (2010).
#'   `"oos"` uses the out-of-sample validation procedure of Abadie (2021)
#'   S.3.2 / ADH (2015): in the outcomes-only case, candidate W(V) are fitted
#'   on training-half outcomes only, V* minimises validation-half MSPE, and
#'   W* is refit with V* on the last floor(T_pre/2) pre-treatment outcomes.
#' @param scale_predictors If `TRUE` (default) and `predictors` are supplied,
#'   each predictor row of X0/X1 is divided by its standard deviation across
#'   all units (treated + donors) before optimisation, matching the Synth
#'   reference implementation (Abadie, Diamond & Hainmueller 2011, JSS).
#'   V weights are then comparable across predictors with different units.
#' @param donor_mspe_threshold Numeric threshold for donor pool filtering
#'   (Abadie 2021 S.4). Each donor's individual pre-treatment MSPE is compared
#'   to the best donor's MSPE: donors with ratio > threshold are excluded before
#'   estimation. `Inf` (default) disables filtering.
#' @param lambda_pen Penalty for penalised SCM (Abadie & L'Hour 2021). `NULL`
#'   (default) runs standard SCM. `"auto"` selects lambda via out-of-sample
#'   pre-treatment MSPE. A non-negative number uses that value directly.
#' @return A list of class `coresynth` with standard fields plus:
#'   * `excluded_donors`: character vector of donors removed by filtering
#'   * `lambda_pen`: penalty value used (NA when standard SCM)
#' @noRd
fit_scm_cpp <- function(
  y,
  d,
  id,
  time,
  data = NULL,
  id_var = NULL,
  time_var = NULL,
  predictors = NULL,
  covariates = NULL,
  v_selection = c("insample", "oos"),
  donor_mspe_threshold = Inf,
  lambda_pen = NULL,
  v_optim = c("coord_descent", "auto", "bfgs"),
  control_group = c("clean", "never_treated"),
  scale_predictors = TRUE,
  ...
) {
  v_selection   <- match.arg(v_selection)
  v_optim       <- match.arg(v_optim)
  control_group <- match.arg(control_group)
  pan <- panel_to_matrices(y, d, id, time)
  .check_panel_complete(pan$Y, "SCM")

  # -- Staggered path (cohort-by-cohort SCM) -----------------------------------
  if (!pan$is_sharp) {
    if (!is.null(predictors) && length(predictors) > 0L)
      stop(
        "SCM staggered adoption does not support 'predictors'.\n",
        "Alternatives:\n",
        "  * Use 'covariates' (time-varying) instead -- these are partialled out before per-cohort SCM.\n",
        "  * Fit each cohort separately as sharp SCM (subset data by adoption time).\n",
        "  * Switch to method = 'sdid', 'gsc', 'mc', or 'tasc', which all support staggered + covariates.",
        call. = FALSE
      )

    # Covariate partial-out for staggered SCM (Clarke et al. 2023 S.2.2)
    use_cov  <- !is.null(covariates) && length(covariates) > 0L && !is.null(data)
    Y_work   <- pan$Y
    beta_hat <- numeric(0)
    if (use_cov) {
      X_arr    <- build_covariate_array(data, id_var, time_var,
                                       covariates, pan$units, pan$times)
      po       <- .sdid_partial_out_staggered(pan$Y, X_arr, pan)
      Y_work   <- po$Y_tilde
      beta_hat <- po$beta_hat
    }

    stag <- .fit_scm_staggered(
      pan,
      Y_matrix             = Y_work,
      v_selection          = v_selection,
      v_optim              = v_optim,
      donor_mspe_threshold = donor_mspe_threshold,
      lambda_pen           = lambda_pen,
      control_group        = control_group
    )
    res <- list(
      method           = "scm",
      staggered        = TRUE,
      estimate         = stag$estimate,
      cohort_estimates = stag$cohort_estimates,
      cohort_fits      = stag$cohort_fits,
      unit_weights     = NULL,
      v_weights        = NULL,
      Y_treat          = pan$Y[, pan$idx_treat, drop = FALSE],
      Y_synth          = NULL,
      gap              = NULL,
      times            = pan$times,
      T_pre            = pan$T_pre,
      excluded_donors  = character(0L),
      lambda_pen       = NA_real_,
      beta_hat         = beta_hat,
      covariates       = covariates
    )
    class(res) <- c("coresynth_scm", "coresynth")
    return(res)
  }

  Y <- pan$Y
  T_pre <- pan$T_pre
  idx_tr <- pan$idx_treat
  idx_co <- pan$idx_control

  if (length(idx_tr) != 1L) {
    stop("SCM currently supports exactly one treated unit.")
  }
  if (length(idx_co) < 2L) {
    stop("SCM requires at least two control units.")
  }
  if (T_pre < 1L) {
    stop("No pre-treatment periods: the treated unit is already treated in ",
         "the first period. SCM needs at least one pre-treatment period to ",
         "fit donor weights.", call. = FALSE)
  }

  Y_co_pre <- Y[seq_len(T_pre), idx_co, drop = FALSE]
  Y_co_all <- Y[, idx_co, drop = FALSE]
  Y_tr_pre <- Y[seq_len(T_pre), idx_tr]
  Y_tr_all <- Y[, idx_tr]

  # -- Donor pool filtering (Abadie 2021 S.4) ----------------------------------
  excluded_donors <- character(0L)
  if (!is.null(donor_mspe_threshold) && is.finite(donor_mspe_threshold)) {
    ind_mspe <- colMeans((as.vector(Y_tr_pre) - Y_co_pre)^2)
    keep     <- ind_mspe <= donor_mspe_threshold * min(ind_mspe)
    if (sum(keep) < 2L) {
      warning(
        "donor_mspe_threshold = ", donor_mspe_threshold,
        " would retain < 2 donors; ignoring filter.",
        call. = FALSE
      )
    } else {
      excluded_donors <- colnames(Y_co_pre)[!keep]
      idx_co    <- idx_co[keep]
      Y_co_pre  <- Y_co_pre[, keep, drop = FALSE]
      Y_co_all  <- Y_co_all[, keep, drop = FALSE]
    }
  }

  # -- Validate lambda_pen ------------------------------------------------------
  if (!is.null(lambda_pen)) {
    if (!(identical(lambda_pen, "auto") ||
          (is.numeric(lambda_pen) && length(lambda_pen) == 1L &&
           is.finite(lambda_pen) && lambda_pen >= 0))) {
      stop(
        "'lambda_pen' must be NULL, \"auto\", or a non-negative finite number.",
        call. = FALSE
      )
    }
  }

  if (is.character(predictors)) {
    stop(
      paste0(
        "'predictors' must be a list of pred() objects, not a character vector.\n",
        "  Use: predictors = list(pred(c('var1', 'var2'), times))"
      ),
      call. = FALSE
    )
  }

  use_cov <- !is.null(predictors) && length(predictors) > 0L

  if (use_cov) {
    if (is.null(data) || is.null(id_var) || is.null(time_var)) {
      stop(
        "'data', 'id_var', and 'time_var' must be provided when using predictors.",
        call. = FALSE
      )
    }
    pm <- build_predictor_matrices(
      data     = data,
      id_var   = id_var,
      time_var = time_var,
      units    = pan$units,
      idx_co   = idx_co,
      idx_tr   = idx_tr,
      predictors = predictors
    )
    X0         <- pm$X0
    X1         <- pm$X1
    pred_names <- pm$pred_names
  } else {
    X0         <- Y_co_pre
    X1         <- drop(Y_tr_pre)
    pred_names <- paste0("V", seq_len(T_pre))
  }

  # -- Synth-style predictor scaling (ADH 2011, JSS) ---------------------------
  # Each predictor row is divided by its sd across all units so that the
  # V-weighted loss (and the reported V weights) are not dominated by
  # predictors with large numeric scales. The outcomes-only case needs no
  # scaling (all rows share the outcome scale).
  X0_raw <- X0
  X1_raw <- X1
  if (use_cov && isTRUE(scale_predictors)) {
    pred_sds <- apply(cbind(X0, X1), 1L, stats::sd)
    pred_sds[!is.finite(pred_sds) | pred_sds < 1e-12] <- 1
    X0 <- X0 / pred_sds
    X1 <- X1 / pred_sds
  }

  Z0 <- Y_co_pre
  Z1 <- drop(Y_tr_pre)

  # OOS handling (Abadie 2021 S.3.2):
  #  * outcomes-only: proper train/validation split via .scm_oos_outcomes()
  #    (candidate W(V) must not see validation-period outcomes).
  #  * user predictors: the predictor matrix is fixed (cannot be re-measured
  #    per window), so only the MSPE evaluation window is restricted.
  oos_outcomes <- (v_selection == "oos") && !use_cov
  t_train      <- if (v_selection == "oos" && use_cov) T_pre %/% 2L else -1L

  # Determine effective outer optimiser
  k_dim <- if (oos_outcomes) max(T_pre %/% 2L, 1L) else nrow(X0)
  effective_outer <- switch(v_optim,
    "auto"          = if (k_dim <= 15L) "bfgs" else "coord_descent",
    "bfgs"          = "bfgs",
    "coord_descent" = "coord_descent"
  )

  v_rows_used <- if (use_cov) NULL else seq_len(T_pre)
  if (oos_outcomes) {
    oos <- .scm_oos_outcomes(Y_co_pre, drop(Y_tr_pre), effective_outer)
    X0          <- oos$X0_final
    X1          <- oos$X1_final
    pred_names  <- paste0("V", oos$v_rows)
    v_rows_used <- oos$v_rows
  }

  # -- Standard or penalised SCM -------------------------------------------------
  if (is.null(lambda_pen)) {
    # Standard SCM
    res <- if (oos_outcomes) {
      oos
    } else if (effective_outer == "bfgs") {
      .scm_bfgs_outer(X0, X1, Z0, Z1, t_train = t_train)
    } else {
      scm_weights_cpp(X0, X1, Z0, Z1, t_train = t_train)
    }
    unit_w <- drop(res$W)
    V_final <- drop(res$V)
    loss    <- res$loss
    lambda_pen_used <- NA_real_
  } else {
    # Step 1: get V* via selected outer optimiser (respects v_selection)
    res_v <- if (oos_outcomes) {
      oos
    } else if (effective_outer == "bfgs") {
      .scm_bfgs_outer(X0, X1, Z0, Z1, t_train = t_train)
    } else {
      scm_weights_cpp(X0, X1, Z0, Z1, t_train = t_train)
    }
    V_star <- drop(res_v$V)

    # Step 2: V*-weighted pairwise distances: d_j = (X1 - X0[:,j])' V* (X1 - X0[:,j])
    diff_mat <- X1 - X0  # k x N_co
    d_pen    <- colSums(V_star * diff_mat^2)  # N_co

    # Step 3: precompute Q = X0' V* X0 and c = X0' V* X1 for the inner QP
    V_sqrt <- sqrt(pmax(V_star, 0))
    X0_sc  <- sweep(X0, 1L, V_sqrt, `*`)    # k x N_co
    Q_mat  <- crossprod(X0_sc)              # N_co x N_co
    c_base <- as.vector(t(X0_sc) %*% (V_sqrt * X1))  # N_co

    # Step 4: select lambda
    if (identical(lambda_pen, "auto")) {
      lambda_pen_used <- .scm_tune_lambda_pen(
        Q_mat, c_base, d_pen, Z0, Z1, t_train, T_pre
      )
    } else {
      lambda_pen_used <- as.numeric(lambda_pen)
    }

    # Step 5: solve penalised QP: min W'Q W - 2(c - lambda/2 * d)' W
    c_pen  <- c_base - lambda_pen_used / 2 * d_pen
    unit_w <- drop(solve_simplex_qp(Q_mat, c_pen))
    V_final <- V_star
    resid  <- Z1 - as.vector(Z0 %*% unit_w)
    loss   <- sqrt(sum(resid^2))
  }

  names(unit_w) <- colnames(Y_co_all)
  Y_synth <- drop(Y_co_all %*% unit_w)
  Y_treat <- drop(Y_tr_all)

  if (use_cov) {
    # Balance table on the original (unscaled) predictor scale
    synth_pred <- drop(X0_raw %*% unit_w)
    predictor_table <- data.frame(
      predictor = pred_names,
      treated   = X1_raw,
      synthetic = synth_pred,
      row.names = NULL
    )
  } else {
    predictor_table <- NULL
  }

  list(
    method          = "scm",
    unit_weights    = unit_w,
    v_weights       = setNames(V_final, pred_names),
    Y_synth         = Y_synth,
    Y_treat         = Y_treat,
    gap             = Y_treat - Y_synth,
    times           = pan$times,
    T_pre           = T_pre,
    estimate        = mean((Y_treat - Y_synth)[-(seq_len(T_pre))]),
    loss            = loss,
    predictor_table = predictor_table,
    X0_mat          = if (use_cov) X0 else NULL,
    X1_vec          = if (use_cov) X1 else NULL,
    v_rows          = v_rows_used,
    Y_co_pre        = Y_co_pre,
    Y_co_post       = Y[-(seq_len(T_pre)), idx_co, drop = FALSE],
    excluded_donors = excluded_donors,
    lambda_pen      = lambda_pen_used
  )
}

# -- Abadie (2021) S.3.2 out-of-sample V selection (outcomes-only) -----------
# 1. Split the pre-period into a training half (1..t0) and a validation half
#    (t0+1..T_pre), with t0 = floor(T_pre / 2).
# 2. For each V, fit W on training-half outcomes ONLY -- validation-period
#    outcomes must not enter the inner QP, otherwise the V optimiser can
#    drive validation MSPE to zero by loading V on validation rows (the
#    leakage this procedure is designed to prevent).
# 3. Select V* minimising validation-half MSPE.
# 4. Refit W* with V* on the outcomes of the LAST t0 pre-treatment periods
#    (Abadie 2021, step 4: "data on the predictors for the last t0 periods").
.scm_oos_outcomes <- function(Y_co_pre, Y_tr_pre, optimizer = "coord_descent") {
  T_pre <- nrow(Y_co_pre)
  t0    <- T_pre %/% 2L
  if (t0 < 1L) {
    stop("v_selection = 'oos' requires at least 2 pre-treatment periods.",
         call. = FALSE)
  }
  train <- seq_len(t0)
  val   <- (t0 + 1L):T_pre

  X0_train <- Y_co_pre[train, , drop = FALSE]
  X1_train <- Y_tr_pre[train]
  Z0_val   <- Y_co_pre[val, , drop = FALSE]
  Z1_val   <- Y_tr_pre[val]

  res <- if (optimizer == "bfgs") {
    .scm_bfgs_outer(X0_train, X1_train, Z0_val, Z1_val, t_train = -1L)
  } else {
    scm_weights_cpp(X0_train, X1_train, Z0_val, Z1_val, t_train = -1L)
  }
  V_star <- drop(res$V)

  v_rows   <- (T_pre - t0 + 1L):T_pre
  X0_final <- Y_co_pre[v_rows, , drop = FALSE]
  X1_final <- Y_tr_pre[v_rows]
  W        <- drop(scm_inner_weights_cpp(X0_final, X1_final, V_star))
  loss     <- sqrt(sum((Y_tr_pre - drop(Y_co_pre %*% W))^2))

  list(W = W, V = V_star, loss = loss,
       v_rows = v_rows, X0_final = X0_final, X1_final = X1_final)
}

# Internal: L-BFGS-B outer V optimisation.
# ~O(k^2) inner QP calls vs coord_descent's k*11*iter calls.
# Mirrors scm_weights_cpp semantics (OOS window, full-data W refit).
.scm_bfgs_outer <- function(X0, X1, Z0, Z1, t_train = -1L) {
  k     <- nrow(X0)
  T_pre <- nrow(Z0)

  # OOS evaluation window -- same split as scm_weights_cpp
  do_oos <- (t_train > 0L && t_train < T_pre)
  if (do_oos) {
    idx     <- (t_train + 1L):T_pre
    Z0_eval <- Z0[idx, , drop = FALSE]
    Z1_eval <- Z1[idx]
  } else {
    Z0_eval <- Z0
    Z1_eval <- Z1
  }

  # Objective: pre-treatment SSE as a function of raw (unnormalised) V weights.
  # V is normalised to sum-to-1 inside the function, so the optimisation is
  # scale-free (L-BFGS-B with lower=0 handles non-negativity).
  fn_v <- function(v_raw) {
    v_raw <- pmax(v_raw, 0)
    v_s   <- sum(v_raw)
    if (v_s < 1e-14) return(1e10)
    W <- tryCatch(
      drop(scm_inner_weights_cpp(X0, X1, v_raw / v_s)),
      error = function(e) NULL
    )
    if (is.null(W) || !all(is.finite(W))) return(1e10)
    r   <- Z1_eval - drop(Z0_eval %*% W)
    val <- sum(r^2)
    if (!is.finite(val)) 1e10 else val
  }

  v0  <- rep(1 / k, k)

  opt <- tryCatch(
    optim(v0, fn_v, method = "L-BFGS-B",
          lower   = rep(0, k),
          upper   = rep(Inf, k),
          control = list(maxit = 300L, factr = 1e7, trace = 0L)),
    error = function(e) NULL
  )
  # Fallback: Nelder-Mead (gradient-free, robust to discontinuities)
  if (is.null(opt) || !is.finite(opt$value)) {
    opt <- tryCatch(
      optim(v0, fn_v, method = "Nelder-Mead",
            control = list(maxit = 500L, reltol = 1e-5)),
      error = function(e) NULL
    )
  }

  V_opt <- if (!is.null(opt) && is.finite(opt$value)) {
    v_raw <- pmax(opt$par, 0)
    s     <- sum(v_raw)
    if (s < 1e-14) rep(1 / k, k) else v_raw / s
  } else {
    rep(1 / k, k)
  }

  W_final    <- drop(scm_inner_weights_cpp(X0, X1, V_opt))
  final_loss <- sqrt(sum((Z1 - drop(Z0 %*% W_final))^2))
  list(W = W_final, V = V_opt, loss = final_loss)
}

# Internal helper: tune lambda for penalised SCM via OOS pre-treatment MSPE.
.scm_tune_lambda_pen <- function(Q_mat, c_base, d_pen, Z0, Z1, t_train, T_pre) {
  if (t_train < 0L || t_train >= T_pre) t_train <- T_pre %/% 2L
  Z0_val <- Z0[(t_train + 1L):T_pre, , drop = FALSE]
  Z1_val <- Z1[(t_train + 1L):T_pre]

  lambda_grid <- c(0, 10^seq(-4, 2, length.out = 19L))
  best_lam    <- 0
  best_mspe   <- Inf

  for (lam in lambda_grid) {
    c_pen <- c_base - lam / 2 * d_pen
    W_try <- tryCatch(
      drop(solve_simplex_qp(Q_mat, c_pen)),
      error = function(e) NULL
    )
    if (is.null(W_try) || !all(is.finite(W_try))) next
    mspe <- mean((Z1_val - as.vector(Z0_val %*% W_try))^2)
    if (is.finite(mspe) && mspe < best_mspe) {
      best_mspe <- mspe
      best_lam  <- lam
    }
  }
  best_lam
}

# -- Internal: cohort-by-cohort SCM for staggered adoption -------------------
# Cohort g ATT is averaged with weight N_treat_g * T_post_g, following the
# same aggregation used in Arkhangelsky et al. (2021) Appendix S.8.
# Treated units within a cohort are averaged into a single pseudo-unit.
# Only predictors = NULL is supported (all pre-treatment outcomes used as X).
.fit_scm_staggered <- function(pan,
                               Y_matrix             = NULL,
                               v_selection          = "insample",
                               v_optim              = "coord_descent",
                               donor_mspe_threshold = Inf,
                               lambda_pen           = NULL,
                               control_group        = "clean") {
  Y       <- if (!is.null(Y_matrix)) Y_matrix else pan$Y
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

    if (length(idx_co_g) < 2L || T_pre_g < 2L) {
      warning(sprintf(
        "SCM staggered: cohort g=%d skipped (T_pre=%d, N_co=%d).",
        g, T_pre_g, length(idx_co_g)
      ), call. = FALSE)
      return(NULL)
    }

    pre_rows  <- seq_len(T_pre_g)
    post_rows <- (T_pre_g + 1L):TT

    Y1_g     <- rowMeans(Y[, idx_tr_g, drop = FALSE])
    Y_co_all <- Y[, idx_co_g, drop = FALSE]
    Y_co_pre <- Y[pre_rows, idx_co_g, drop = FALSE]  # T_pre_g x N_co_g
    Y_tr_pre <- Y1_g[pre_rows]

    # Donor pool filtering per cohort (same rule as the sharp path)
    excl_g <- character(0L)
    if (is.finite(donor_mspe_threshold) && donor_mspe_threshold < Inf) {
      ind_mspe <- colMeans((Y_tr_pre - Y_co_pre)^2)
      keep     <- ind_mspe <= donor_mspe_threshold * min(ind_mspe)
      if (sum(keep) >= 2L) {
        excl_g   <- colnames(Y_co_pre)[!keep]
        idx_co_g <- idx_co_g[keep]
        Y_co_all <- Y_co_all[, keep, drop = FALSE]
        Y_co_pre <- Y_co_pre[, keep, drop = FALSE]
      }
    }

    t_train_g <- if (v_selection == "oos") T_pre_g %/% 2L else -1L
    k_dim_g   <- T_pre_g  # predictors = NULL => k = T_pre_g
    effective_outer_g <- switch(v_optim,
      "auto"          = if (k_dim_g <= 15L) "bfgs" else "coord_descent",
      "bfgs"          = "bfgs",
      "coord_descent" = "coord_descent"
    )

    fit_g <- tryCatch({
      oos_g <- (v_selection == "oos")
      if (is.null(lambda_pen)) {
        # Standard SCM (proper train/validation split when OOS)
        res <- if (oos_g) {
          .scm_oos_outcomes(Y_co_pre, Y_tr_pre, effective_outer_g)
        } else if (effective_outer_g == "bfgs") {
          .scm_bfgs_outer(Y_co_pre, Y_tr_pre, Y_co_pre, Y_tr_pre,
                          t_train = -1L)
        } else {
          scm_weights_cpp(Y_co_pre, Y_tr_pre, Y_co_pre, Y_tr_pre,
                          t_train = -1L)
        }
        unit_w          <- drop(res$W)
        lambda_pen_used <- NA_real_
      } else {
        # Penalised SCM, per cohort
        res_v <- if (oos_g) {
          .scm_oos_outcomes(Y_co_pre, Y_tr_pre, effective_outer_g)
        } else if (effective_outer_g == "bfgs") {
          .scm_bfgs_outer(Y_co_pre, Y_tr_pre, Y_co_pre, Y_tr_pre,
                          t_train = -1L)
        } else {
          scm_weights_cpp(Y_co_pre, Y_tr_pre, Y_co_pre, Y_tr_pre,
                          t_train = -1L)
        }
        V_star <- drop(res_v$V)
        # OOS: V* refers to the last-t0 window, so the penalised QP must be
        # built on the same rows (X0_final/X1_final).
        X0_pen <- if (oos_g) res_v$X0_final else Y_co_pre
        X1_pen <- if (oos_g) res_v$X1_final else Y_tr_pre
        diff_mat <- X1_pen - X0_pen           # k x N_co_g
        d_pen    <- colSums(V_star * diff_mat^2)
        V_sqrt   <- sqrt(pmax(V_star, 0))
        X0_sc    <- sweep(X0_pen, 1L, V_sqrt, `*`)
        Q_mat    <- crossprod(X0_sc)
        c_base   <- as.vector(t(X0_sc) %*% (V_sqrt * X1_pen))

        lambda_pen_used <- if (identical(lambda_pen, "auto")) {
          .scm_tune_lambda_pen(Q_mat, c_base, d_pen,
                               Y_co_pre, Y_tr_pre, t_train_g, T_pre_g)
        } else {
          as.numeric(lambda_pen)
        }
        c_pen  <- c_base - lambda_pen_used / 2 * d_pen
        unit_w <- drop(solve_simplex_qp(Q_mat, c_pen))
      }

      names(unit_w) <- colnames(Y_co_all)
      Y_synth_g <- drop(Y_co_all %*% unit_w)
      ATT_g     <- mean((Y1_g - Y_synth_g)[post_rows])

      list(
        cohort          = g,
        n_treated       = length(idx_tr_g),
        T_pre           = T_pre_g,
        T_post          = T_post_g,
        estimate        = ATT_g,
        weight          = as.numeric(length(idx_tr_g)) * T_post_g,
        unit_weights    = unit_w,
        Y_synth         = Y_synth_g,
        Y_treat         = Y1_g,
        excluded_donors = excl_g,
        lambda_pen      = lambda_pen_used
      )
    }, error = function(e) {
      warning(sprintf("SCM staggered: cohort g=%d failed: %s",
                      g, conditionMessage(e)), call. = FALSE)
      NULL
    })
    fit_g
  })

  valid <- !vapply(cohort_list, is.null, logical(1L))
  if (!any(valid)) {
    stop("All cohort-level SCM fits failed; see the preceding warnings for ",
         "per-cohort reasons.", call. = FALSE)
  }
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

#' Permutation Inference via MSPE Ratio for SCM
#'
#' Computes the Abadie et al. (2010) / Abadie (2021) permutation p-value.
#' For each control unit, a leave-one-out synthetic control is fitted.
#'
#' When `alternative = "two.sided"` (default), the test statistic is the
#' post/pre MSPE ratio, following Abadie et al. (2010). When
#' `alternative = "greater"` or `"less"`, the test statistic is the signed
#' average post-treatment gap (ATT), giving a one-sided permutation test as
#' recommended by Abadie (2021) S.3.5 for improved power when the direction
#' of the treatment effect is known.
#'
#' @param fit A `coresynth` object from [scm_fit()] with `method = "scm"`.
#' @param mspe_threshold Minimum pre-treatment MSPE for including a control
#'   unit in the two-sided test. Ignored for one-sided tests.
#'   Default: 0 (no filtering).
#' @param max_iter Passed to [scm_placebo_cpp()]. Default 100L.
#' @param tol      Passed to [scm_placebo_cpp()]. Default 1e-4.
#' @param use_covariates If `TRUE` and the fit used predictor covariates,
#'   applies the same covariate spec to each placebo unit (R-level loop).
#'   Default `FALSE` (faster C++ outcomes-only placebos).
#' @param alternative Direction of the alternative hypothesis:
#'   `"two.sided"` (default) uses the MSPE ratio statistic;
#'   `"greater"` tests whether the treatment increased the outcome;
#'   `"less"` tests whether the treatment decreased the outcome.
#'   One-sided tests use the signed ATT as the test statistic.
#' @return A list with:
#'   * `p_value`: Permutation p-value between 0 and 1
#'   * `mspe_ratio_treated`: MSPE_post / MSPE_pre for the treated unit (two.sided only)
#'   * `mspe_ratios_all`: Named numeric vector (treated first, then controls); two.sided only
#'   * `placebo_effects`: Named N_co-vector of placebo ATT estimates
#'   * `treated_effect`: ATT estimate for the treated unit
#'   * `n_placebo_used`: Number of control units used
#' @export
mspe_ratio_pval <- function(
  fit,
  mspe_threshold = 0,
  max_iter = 100L,
  tol = 1e-4,
  use_covariates = FALSE,
  alternative = c("two.sided", "greater", "less")
) {
  alternative <- match.arg(alternative)
  if (!inherits(fit, "coresynth") || fit$method != "scm") {
    stop("mspe_ratio_pval() requires a coresynth object with method = 'scm'.")
  }
  if (isTRUE(fit$staggered)) {
    stop(
      "mspe_ratio_pval() does not support staggered fits.\n",
      "Alternatives:\n",
      "  * Apply mspe_ratio_pval() to each cohort_fit individually after re-fitting it as sharp SCM.\n",
      "  * Use method = 'sdid' with sdid_inference(method = 'placebo') for staggered placebo inference.",
      call. = FALSE
    )
  }
  if (is.null(fit$Y_co_pre) || is.null(fit$Y_co_post)) {
    stop(
      "fit is missing Y_co_pre / Y_co_post. Re-run scm_fit() with the current version."
    )
  }

  T_pre <- fit$T_pre
  T_all <- length(fit$times)
  pre   <- seq_len(T_pre)
  post  <- (T_pre + 1L):T_all

  mspe_pre_tr  <- mean((fit$Y_treat[pre]  - fit$Y_synth[pre])^2)
  mspe_post_tr <- mean((fit$Y_treat[post] - fit$Y_synth[post])^2)
  ratio_tr <- if (mspe_pre_tr > mspe_threshold) {
    mspe_post_tr / mspe_pre_tr
  } else {
    NA_real_
  }

  if (use_covariates && is.null(fit$X0_mat)) {
    warning(
      "use_covariates=TRUE but fit has no covariates (X0_mat is NULL). ",
      "Falling back to outcomes-only placebo."
    )
    use_covariates <- FALSE
  }

  if (use_covariates) {
    N_co         <- ncol(fit$Y_co_pre)
    mspe_pre_vec <- numeric(N_co)
    mspe_post_vec <- numeric(N_co)
    effects_vec  <- numeric(N_co)

    for (i in seq_len(N_co)) {
      X0_loo <- fit$X0_mat[, -i, drop = FALSE]
      X1_loo <- fit$X0_mat[, i]
      Z0_loo <- fit$Y_co_pre[, -i, drop = FALSE]
      Z1_loo <- fit$Y_co_pre[, i]

      res_loo <- tryCatch(
        scm_weights_cpp(X0_loo, X1_loo, Z0_loo, Z1_loo,
                        max_iter = max_iter, tol = tol),
        error = function(e) NULL
      )
      if (is.null(res_loo)) {
        mspe_pre_vec[i]  <- NA_real_
        mspe_post_vec[i] <- NA_real_
        effects_vec[i]   <- NA_real_
        next
      }
      W_loo      <- res_loo$W
      synth_pre  <- fit$Y_co_pre[, -i, drop = FALSE]  %*% W_loo
      synth_post <- fit$Y_co_post[, -i, drop = FALSE] %*% W_loo

      mspe_pre_vec[i]  <- mean((Z1_loo - synth_pre)^2)
      mspe_post_vec[i] <- mean((fit$Y_co_post[, i] - synth_post)^2)
      effects_vec[i]   <- mean(fit$Y_co_post[, i] - synth_post)
    }

    keep      <- !is.na(mspe_pre_vec) & mspe_pre_vec > mspe_threshold
    ratios_co <- ifelse(keep, mspe_post_vec / mspe_pre_vec, NA_real_)
    co_effects <- effects_vec
  } else {
    plac      <- scm_placebo_cpp(fit$Y_co_pre, fit$Y_co_post,
                                  max_iter = max_iter, tol = tol)
    keep      <- plac$mspe_pre > mspe_threshold
    ratios_co <- ifelse(keep, plac$mspe_post / plac$mspe_pre, NA_real_)
    co_effects <- plac$effects
  }

  co_names <- colnames(fit$Y_co_pre)
  if (!is.null(co_names)) {
    names(ratios_co)  <- co_names
    names(co_effects) <- co_names
  }

  treated_effect <- fit$estimate

  if (alternative == "two.sided") {
    all_ratios <- c(treated = ratio_tr, ratios_co)
    valid      <- !is.na(all_ratios)
    p_value    <- mean(all_ratios[valid] >= all_ratios["treated"], na.rm = TRUE)
    n_used     <- sum(keep, na.rm = TRUE)
  } else {
    all_effects <- c(treated = treated_effect, co_effects)
    valid       <- !is.na(all_effects)
    if (alternative == "greater") {
      p_value <- mean(all_effects[valid] >= all_effects["treated"], na.rm = TRUE)
    } else {
      p_value <- mean(all_effects[valid] <= all_effects["treated"], na.rm = TRUE)
    }
    all_ratios <- c(treated = ratio_tr, ratios_co)
    n_used     <- sum(!is.na(co_effects))
  }

  list(
    p_value            = p_value,
    mspe_ratio_treated = ratio_tr,
    mspe_ratios_all    = all_ratios,
    placebo_effects    = co_effects,
    treated_effect     = treated_effect,
    n_placebo_used     = n_used
  )
}

# -- ADH (2015) validation exercises -------------------------------------------

#' In-Time Placebo (Backdating) Test for SCM
#'
#' Re-estimates the synthetic control after artificially backdating the
#' treatment to a pre-treatment period, following Abadie, Diamond &
#' Hainmueller (2015) and Abadie & Vives-i-Bastida (2022, principle 7:
#' "out-of-sample validation is key"). Only pre-treatment data enter the
#' exercise, so the placebo gap after `t0_placebo` is uncontaminated by the
#' actual intervention. A credible design shows no sizable divergence at the
#' backdated treatment time.
#'
#' The refit uses the outcomes of periods `1..t0_placebo` as predictors
#' (the `predictors = NULL` default), regardless of how the original fit was
#' specified, because user-supplied `pred()` windows cannot be lagged
#' automatically (ADH 2015 lag their predictors by hand).
#'
#' @param fit A sharp `coresynth` object from [scm_fit()] with
#'   `method = "scm"`.
#' @param t0_placebo Backdated treatment period as a 1-based position in
#'   `fit$times`; must satisfy `2 <= t0_placebo < T_pre`. Default
#'   `floor(T_pre / 2)`.
#' @return A list with:
#'   * `t0_placebo`: the backdated treatment period used
#'   * `times`: time values of the pre-treatment window
#'   * `unit_weights`: placebo donor weights
#'   * `Y_treat`, `Y_synth`, `gap`: series over the pre-treatment window
#'   * `placebo_att`: mean placebo gap over `(t0_placebo, T_pre]`
#'   * `fit_rmspe`: RMSPE over the placebo fitting window `1..t0_placebo`
#'   * `eval_rmspe`: RMSPE over the placebo post window `(t0_placebo, T_pre]`
#' @seealso [mspe_ratio_pval()] for in-space placebos, [loo_donors()] for
#'   donor-robustness checks.
#' @export
placebo_in_time <- function(fit, t0_placebo = NULL) {
  if (!inherits(fit, "coresynth") || fit$method != "scm") {
    stop("placebo_in_time() requires a coresynth object with method = 'scm'.",
         call. = FALSE)
  }
  if (isTRUE(fit$staggered)) {
    stop("placebo_in_time() supports sharp (single-cohort) fits only.",
         call. = FALSE)
  }
  if (is.null(fit$Y_co_pre)) {
    stop("fit is missing Y_co_pre. Re-run scm_fit() with the current version.",
         call. = FALSE)
  }

  T_pre <- fit$T_pre
  if (is.null(t0_placebo)) t0_placebo <- T_pre %/% 2L
  t0_placebo <- as.integer(t0_placebo)
  if (t0_placebo < 2L || t0_placebo >= T_pre) {
    stop("'t0_placebo' must satisfy 2 <= t0_placebo < T_pre (= ", T_pre, ").",
         call. = FALSE)
  }

  Y_co_pre <- fit$Y_co_pre
  Y_tr_pre <- fit$Y_treat[seq_len(T_pre)]
  fit_rows  <- seq_len(t0_placebo)
  eval_rows <- (t0_placebo + 1L):T_pre

  res <- scm_weights_cpp(
    Y_co_pre[fit_rows, , drop = FALSE], Y_tr_pre[fit_rows],
    Y_co_pre[fit_rows, , drop = FALSE], Y_tr_pre[fit_rows]
  )
  W <- drop(res$W)
  names(W) <- colnames(Y_co_pre)

  Y_synth <- drop(Y_co_pre %*% W)
  gap     <- Y_tr_pre - Y_synth

  list(
    t0_placebo   = t0_placebo,
    times        = fit$times[seq_len(T_pre)],
    unit_weights = W,
    Y_treat      = Y_tr_pre,
    Y_synth      = Y_synth,
    gap          = gap,
    placebo_att  = mean(gap[eval_rows]),
    fit_rmspe    = sqrt(mean(gap[fit_rows]^2)),
    eval_rmspe   = sqrt(mean(gap[eval_rows]^2))
  )
}

#' Leave-One-Out Donor Robustness for SCM
#'
#' Iteratively re-estimates the synthetic control excluding one contributing
#' donor at a time, holding the predictor weights V fixed at their baseline
#' values (Abadie, Diamond & Hainmueller 2015, footnote 20). The spread of
#' the leave-one-out ATT estimates shows how much the result hinges on any
#' single donor.
#'
#' For penalised fits (`lambda_pen` used), the same penalty is re-applied in
#' each leave-one-out QP.
#'
#' @param fit A sharp `coresynth` object from [scm_fit()] with
#'   `method = "scm"`.
#' @param weight_threshold Only donors whose baseline weight exceeds this
#'   value are dropped (removing a zero-weight donor cannot change the fit).
#'   Default `1e-6`.
#' @return A list with:
#'   * `att_original`: baseline ATT
#'   * `results`: data.frame with one row per excluded donor
#'     (`donor`, `weight`, `att_loo`)
#'   * `att_range`: range of the leave-one-out ATTs
#' @seealso [placebo_in_time()], [mspe_ratio_pval()]
#' @export
loo_donors <- function(fit, weight_threshold = 1e-6) {
  if (!inherits(fit, "coresynth") || fit$method != "scm") {
    stop("loo_donors() requires a coresynth object with method = 'scm'.",
         call. = FALSE)
  }
  if (isTRUE(fit$staggered)) {
    stop("loo_donors() supports sharp (single-cohort) fits only.",
         call. = FALSE)
  }
  if (is.null(fit$Y_co_pre) || is.null(fit$Y_co_post)) {
    stop("fit is missing Y_co_pre / Y_co_post. Re-run scm_fit() with the ",
         "current version.", call. = FALSE)
  }

  Y_co_pre  <- fit$Y_co_pre
  Y_co_post <- fit$Y_co_post
  T_pre     <- fit$T_pre
  W0        <- fit$unit_weights
  V         <- unname(fit$v_weights)
  N_co      <- ncol(Y_co_pre)
  if (N_co < 3L) {
    stop("loo_donors() requires at least 3 donors.", call. = FALSE)
  }

  # Predictor matrices consistent with the original V optimisation:
  # predictor-based fits store the (scaled) X0/X1; outcomes-only fits use
  # the pre-treatment outcome rows recorded in v_rows (handles OOS fits,
  # where V refers to the last floor(T_pre/2) periods only).
  if (!is.null(fit$X0_mat)) {
    X0 <- fit$X0_mat
    X1 <- fit$X1_vec
  } else {
    v_rows <- if (!is.null(fit$v_rows)) fit$v_rows else seq_len(T_pre)
    X0 <- Y_co_pre[v_rows, , drop = FALSE]
    X1 <- fit$Y_treat[v_rows]
  }
  use_pen <- is.finite(fit$lambda_pen) && fit$lambda_pen > 0

  refit_w <- function(keep) {
    X0_k <- X0[, keep, drop = FALSE]
    if (!use_pen) {
      return(drop(scm_inner_weights_cpp(X0_k, X1, V)))
    }
    d_pen  <- colSums(V * (X1 - X0_k)^2)
    V_sqrt <- sqrt(pmax(V, 0))
    X0_sc  <- sweep(X0_k, 1L, V_sqrt, `*`)
    Q_mat  <- crossprod(X0_sc)
    c_pen  <- as.vector(t(X0_sc) %*% (V_sqrt * X1)) -
      fit$lambda_pen / 2 * d_pen
    drop(solve_simplex_qp(Q_mat, c_pen))
  }

  drop_idx <- which(W0 > weight_threshold)
  if (length(drop_idx) == 0L) drop_idx <- seq_len(N_co)

  Y_tr_post <- fit$Y_treat[-seq_len(T_pre)]
  att_loo <- vapply(drop_idx, function(j) {
    keep  <- setdiff(seq_len(N_co), j)
    W_loo <- tryCatch(refit_w(keep), error = function(e) NULL)
    if (is.null(W_loo) || !all(is.finite(W_loo))) return(NA_real_)
    synth_post <- drop(Y_co_post[, keep, drop = FALSE] %*% W_loo)
    mean(Y_tr_post - synth_post)
  }, numeric(1L))

  results <- data.frame(
    donor   = colnames(Y_co_pre)[drop_idx],
    weight  = unname(W0[drop_idx]),
    att_loo = att_loo,
    stringsAsFactors = FALSE
  )

  list(
    att_original = fit$estimate,
    results      = results,
    att_range    = range(att_loo, na.rm = TRUE)
  )
}

# -- Augmented SCM (Ben-Michael, Feller & Rothstein 2021) ----------------------

# Internal helper: select ridge lambda via LOO-CV on control units.
.cv_ridge_scm <- function(A, b, T_pre,
                           lambda_grid = 10^seq(-4, 4, length.out = 30L)) {
  N_co <- nrow(A)
  Ipre <- diag(T_pre)

  cv_errors <- vapply(lambda_grid, function(lam) {
    loo <- vapply(seq_len(N_co), function(i) {
      A_i    <- A[-i, , drop = FALSE]
      b_i    <- b[-i]
      beta_i <- tryCatch(
        solve(crossprod(A_i) + lam * Ipre, drop(t(A_i) %*% b_i)),
        error = function(e) rep(NA_real_, T_pre)
      )
      if (anyNA(beta_i)) return(NA_real_)
      (b[i] - sum(A[i, ] * beta_i))^2
    }, numeric(1L))
    mean(loo, na.rm = TRUE)
  }, numeric(1L))

  lambda_grid[which.min(cv_errors)]
}

#' Augmented Synthetic Control Method (Ridge ASCM)
#'
#' Applies a ridge-regression-based bias correction to a fitted SCM object,
#' following Ben-Michael, Feller & Rothstein (2021, JASA). The corrected
#' estimator is:
#'
#'   tau_aug = tau_SCM + (m_tr_post - sum_j W_j * m_j_post)
#'
#' where m_i_post = Y_pre_i' beta_hat is the ridge outcome model prediction
#' for unit i's mean post-treatment outcome, and beta_hat is estimated by
#' ridge regression across control units.
#'
#' @param fit A `coresynth` object from [scm_fit()] with `method = "scm"`.
#' @param lambda_ridge Ridge penalty (non-negative). `NULL` (default) selects
#'   the penalty by leave-one-out cross-validation on the control units.
#' @return A list with:
#'   * `att_aug`: Augmented ATT estimate
#'   * `delta`: Bias correction term (m_tr_post - sum_j W_j m_j_post)
#'   * `att_scm`: Original SCM ATT for comparison
#'   * `lambda_ridge`: Ridge penalty used
#'   * `beta_hat`: Ridge regression coefficients (length T_pre)
#' @export
augment_scm <- function(fit, lambda_ridge = NULL) {
  if (!inherits(fit, "coresynth") || fit$method != "scm") {
    stop("augment_scm() requires a coresynth object with method = 'scm'.", call. = FALSE)
  }
  if (is.null(fit$Y_co_pre) || is.null(fit$Y_co_post)) {
    stop(
      "fit is missing Y_co_pre / Y_co_post. Re-run scm_fit() with the current version.",
      call. = FALSE
    )
  }
  if (!is.null(lambda_ridge) &&
      (!is.numeric(lambda_ridge) || length(lambda_ridge) != 1L ||
       !is.finite(lambda_ridge) || lambda_ridge < 0)) {
    stop("'lambda_ridge' must be NULL or a non-negative finite number.", call. = FALSE)
  }

  Y_co_pre  <- fit$Y_co_pre
  Y_co_post <- fit$Y_co_post
  W         <- fit$unit_weights
  T_pre     <- fit$T_pre
  N_co      <- ncol(Y_co_pre)

  # Design: A[i, t] = pre-treatment outcome of control unit i at time t
  A <- t(Y_co_pre)          # N_co x T_pre
  b <- colMeans(Y_co_post)  # N_co: mean post-treatment outcome per donor

  if (is.null(lambda_ridge)) {
    lambda_ridge <- .cv_ridge_scm(A, b, T_pre)
  }

  Ipre     <- diag(T_pre)
  AtA      <- crossprod(A)            # T_pre x T_pre
  Atb      <- drop(t(A) %*% b)       # T_pre
  beta_hat <- tryCatch(
    solve(AtA + lambda_ridge * Ipre, Atb),
    error = function(e) {
      warning("Ridge solve failed; increasing regularization.", call. = FALSE)
      lam2 <- max(lambda_ridge * 100, 1e-6)
      tryCatch(
        solve(AtA + lam2 * Ipre, Atb),
        error = function(e2) rep(0, T_pre)
      )
    }
  )

  Y_tr_pre <- fit$Y_treat[seq_len(T_pre)]
  m_hat_co <- as.vector(A %*% beta_hat)  # N_co: predicted post-mean per donor
  m_hat_tr <- sum(Y_tr_pre * beta_hat)   # scalar: predicted post-mean for treated

  delta   <- m_hat_tr - sum(W * m_hat_co)
  att_aug <- fit$estimate + delta

  list(
    att_aug      = att_aug,
    delta        = delta,
    att_scm      = fit$estimate,
    lambda_ridge = lambda_ridge,
    beta_hat     = beta_hat
  )
}
