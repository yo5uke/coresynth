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
  outcome_var = NULL,
  predictors = NULL,
  covariates = NULL,
  v_selection = c("insample", "oos"),
  donor_mspe_threshold = Inf,
  lambda_pen = NULL,
  v_optim = c("auto", "coord_descent", "bfgs", "multistart"),
  qp_solver = c("auto", "active_set", "wolfe"),
  v_window = NULL,
  control_group = c("clean", "never_treated"),
  scale_predictors = TRUE,
  nu = NULL,
  fixedeff = NULL,
  pool_cohorts = TRUE,
  weights = NULL,
  ...
) {
  v_selection   <- match.arg(v_selection)
  v_optim       <- match.arg(v_optim)
  qp_solver     <- match.arg(qp_solver)
  control_group <- match.arg(control_group)
  if (!is.null(nu) && !identical(nu, "auto") && !identical(nu, "separate") &&
      (!is.numeric(nu) || length(nu) != 1L || nu < 0 || nu > 1)) {
    stop("nu must be NULL, \"auto\", \"separate\", or a single number in ",
         "[0, 1].", call. = FALSE)
  }
  pan <- panel_to_matrices(y, d, id, time)
  .check_panel_complete(pan$Y, "SCM")

  if (pan$is_sharp && (!is.null(nu) || isTRUE(fixedeff) ||
                       !isTRUE(pool_cohorts))) {
    stop(
      "nu, fixedeff and pool_cohorts apply to staggered SCM fits only ",
      "(Ben-Michael, Feller & Rothstein 2022). ",
      "Under a single adoption date there is nothing to pool across cohorts, ",
      "and several treated units are already fully pooled into one cohort; ",
      "fit without these arguments.",
      call. = FALSE
    )
  }
  # Covariates enter SCM by being balanced, not by being partialled out:
  # Abadie et al. (2010) S.2.3 put them in the V-weighted predictor distance,
  # and BFR (2022) S.5.2 add a separate xi-weighted balance term, noting in
  # S.4 that "where Abadie et al. (2010, 2015) balance auxiliary covariates,
  # we focus exclusively on lagged outcomes". Partialling out is the SDID
  # device of Clarke et al. (2024) S.2.2, which the staggered path borrows
  # only because pred() is unavailable there. A sharp fit has pred(), so
  # silently ignoring `covariates` -- which is what happened before -- is the
  # one thing not to do.
  if (pan$is_sharp && !is.null(covariates) && length(covariates) > 0L) {
    stop(
      "'covariates' is the partial-out path, which applies to staggered SCM ",
      "fits only. Sharp SCM balances covariates instead, following Abadie, ",
      "Diamond & Hainmueller (2010) S.2.3: pass them as predictor rows, e.g. ",
      "predictors = list(pred(\"", outcome_var %||% "<outcome>",
      "\", <pre-treatment times>, \"mean\"), pred(\"",
      as.character(unlist(covariates))[1L],
      "\", <pre-treatment times>, \"mean\")). A time-varying covariate is ",
      "aggregated over its window, which is how BFR (2022) S.5.2 carry one ",
      "into the balance problem too.",
      call. = FALSE
    )
  }
  # BFR S.5.1 add an intercept shift to partially pooled SCM, and report in
  # S.6 that the intercept-shifted estimator "dominates either of the
  # alternatives in terms of both overall and unit-level estimates". augsynth,
  # their implementation, therefore defaults to it as well.
  if (is.null(fixedeff)) fixedeff <- !pan$is_sharp

  # Aggregation weights over the treated units (Abadie 2021 S.8: with regions
  # or countries a population- or GDP-weighted average effect is often the
  # quantity of interest). Resolved here so both paths see the same vector.
  agg_w <- .scm_treated_weights(weights, pan)

  # The Wolfe solver is wired through the sharp inner QP only. The paths
  # excluded here either run their own solver (the penalised QP and the
  # partially pooled block descent are solved in R via solve_simplex_qp()) or
  # would report weights from a solver the flag never reached, which is worse
  # than refusing. "auto" therefore resolves to Wolfe wherever it is wired and
  # falls back quietly elsewhere; an explicit "wolfe" still errors, because a
  # caller who named the solver should be told it cannot be honoured.
  bad <- c(
    if (!pan$is_sharp) "staggered adoption",
    if (v_selection == "oos") "v_selection = 'oos'",
    if (!is.null(lambda_pen)) "lambda_pen",
    if (v_optim == "bfgs") "v_optim = 'bfgs'"
  )
  if (identical(qp_solver, "wolfe") && length(bad)) {
    stop(
      "qp_solver = 'wolfe' is not available with ",
      paste(bad, collapse = ", "),
      ". These paths solve their own QP outside the sharp inner solver. ",
      "Use qp_solver = 'active_set', or the 'auto' default, which falls ",
      "back to it here.",
      call. = FALSE
    )
  }
  if (identical(qp_solver, "auto"))
    qp_solver <- if (length(bad)) "active_set" else "wolfe"
  use_wolfe <- identical(qp_solver, "wolfe")

  # -- Staggered path (cohort-by-cohort SCM) -----------------------------------
  if (!pan$is_sharp) {
    if (!is.null(v_window)) {
      stop(
        "v_window applies to sharp (single-cohort) SCM fits only. ",
        "Cohorts have different pre-treatment windows, so a single ",
        "evaluation window is not well defined for staggered fits.",
        call. = FALSE
      )
    }
    if (v_optim == "multistart") {
      stop(
        "v_optim = 'multistart' requires a predictor specification, and ",
        "staggered SCM is outcomes-only. Use v_optim = 'auto' or ",
        "'coord_descent'.",
        call. = FALSE
      )
    }
    if (!is.null(predictors) && length(predictors) > 0L)
      stop(
        "SCM staggered adoption does not support 'predictors'.\n",
        "Alternatives:\n",
        "  * Use 'covariates' (time-varying) instead -- these are partialled out before per-cohort SCM.\n",
        "  * Fit each cohort separately as sharp SCM (subset data by adoption time).\n",
        "  * Switch to method = 'sdid', 'gsc', 'mc', or 'tasc', which all support staggered + covariates.",
        call. = FALSE
      )

    # Covariate partial-out for staggered SCM (Clarke et al. 2024 S.2.2)
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

    # Partially pooled SCM (Ben-Michael, Feller & Rothstein 2022) is
    # outcomes-only with uniform lag weights, so the per-cohort path's
    # V-selection / donor-filtering / penalty knobs do not reach it.
    pooled_blocked <- !identical(donor_mspe_threshold, Inf) ||
      !is.null(lambda_pen) || v_selection == "oos"
    if (!is.null(nu) && !identical(nu, "separate") && pooled_blocked) {
      stop(
        "nu (partially pooled SCM) cannot be combined with ",
        "donor_mspe_threshold, lambda_pen, or v_selection = 'oos'. ",
        "Use nu = 'separate' for the per-cohort V-optimised path.",
        call. = FALSE
      )
    }
    # BFR choose nu by their own heuristic and warn against values close to 0
    # or 1; augsynth, their implementation, does the same when nu is left
    # unset. Separate SCM (nu = 0) is the practice they argue against, so it
    # is opt-in rather than the default -- except where the knobs above have
    # already ruled the pooled solver out.
    if (is.null(nu)) nu <- if (pooled_blocked) "separate" else "auto"

    if (!identical(nu, "separate")) {
      stag <- .fit_scm_staggered_pooled(
        pan,
        Y_matrix      = Y_work,
        nu            = nu,
        fixedeff      = fixedeff,
        control_group = control_group,
        pool_cohorts  = pool_cohorts,
        unit_weights  = agg_w
      )
    } else {
      stag <- .fit_scm_staggered(
        pan,
        unit_weights         = agg_w,
        Y_matrix             = Y_work,
        v_selection          = v_selection,
        v_optim              = v_optim,
        donor_mspe_threshold = donor_mspe_threshold,
        lambda_pen           = lambda_pen,
        control_group        = control_group,
        pool_cohorts         = pool_cohorts,
        fixedeff             = fixedeff
      )
    }
    res <- list(
      method           = "scm",
      staggered        = TRUE,
      estimate         = stag$estimate,
      cohort_estimates = stag$cohort_estimates,
      cohort_fits      = stag$cohort_fits,
      pooling          = stag$pooling,   # NULL for the per-cohort path
      fixedeff         = fixedeff,
      # Recorded on both paths: `pooling` carries it only for the pooled one,
      # and it decides what a cohort_fit block stands for.
      pool_cohorts     = pool_cohorts,
      n_treated        = length(pan$idx_treat),
      # Rescaled aggregation masses, kept so scm_balance_frontier() can refit
      # the same objective. NULL when the units count equally.
      treated_weights  = agg_w,
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
      covariates       = covariates,
      control_group    = control_group,
      Y_all            = Y_work   # T x N -- scm_rank_pval() / frontier refits
    )
    return(res)  # class tags are attached by new_coresynth() in scm_fit()
  }

  Y <- pan$Y
  T_pre <- pan$T_pre
  idx_tr <- pan$idx_treat
  idx_co <- pan$idx_control

  if (length(idx_co) < 2L) {
    stop("SCM requires at least two control units.")
  }
  if (T_pre < 1L) {
    stop("No pre-treatment periods: the treated unit is already treated in ",
         "the first period. SCM needs at least one pre-treatment period to ",
         "fit donor weights.", call. = FALSE)
  }

  # Simultaneous adoption by several units is one cohort. Ben-Michael, Feller
  # & Rothstein (2022) S.4.1 show the across-unit variation that makes
  # unit-level fits enter the ATT error bound is zero when every treated unit
  # shares an adoption date, so the bound is minimised by fully pooling them
  # -- that is, by fitting one synthetic control to their average. This is the
  # same step their Appendix A.2 hierarchy takes within each cohort, and the
  # same one the staggered path here already takes.
  n_treated  <- length(idx_tr)
  Y_treat_mat <- Y[, idx_tr, drop = FALSE]          # T x n_treated, raw series
  # `weights` makes the pooled target a weighted rather than simple average,
  # which is what carries a population-weighted ATT through the sharp path.
  w_tr <- if (is.null(agg_w)) rep(1, n_treated) else
    unname(agg_w[as.character(pan$units[idx_tr])])
  w_tr <- w_tr / sum(w_tr)
  Y_tr_all <- if (n_treated == 1L) Y[, idx_tr] else
    drop(Y_treat_mat %*% w_tr)
  Y_tr_pre <- Y_tr_all[seq_len(T_pre)]

  Y_co_pre <- Y[seq_len(T_pre), idx_co, drop = FALSE]
  Y_co_all <- Y[, idx_co, drop = FALSE]

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

  pre_times <- pan$times[seq_len(T_pre)]

  # A spec listing solely the outcome at each single pre-treatment period
  # builds X0/X1 equal to the pre-treatment outcome rows, i.e. it IS the
  # outcomes-only fit. Route it there so both ways of writing that model
  # return the identical (and much faster) fit. An explicit predictor-path
  # optimiser request (multistart/bfgs) overrides the routing.
  if (use_cov && v_optim %in% c("auto", "coord_descent")) {
    if (.is_outcome_only_spec(predictors, outcome_var, pre_times)) {
      message(
        "predictors specify the outcome at every pre-treatment period; ",
        "fitting through the outcomes-only path (identical to ",
        "predictors = NULL)."
      )
      use_cov <- FALSE
    } else {
      # Near miss. Coverage is judged against the periods `data` carries, so
      # the same `predictors` list lands on either path depending on how far
      # back `data` reaches -- a change of estimand that is invisible from
      # the call. Name the mismatch instead of letting it pass in silence.
      lag_times <- .outcome_lag_times(predictors, outcome_var)
      if (!is.null(lag_times)) {
        n_cov <- length(unique(lag_times[lag_times %in% pre_times]))
        message(
          "predictors are single-period outcome lags covering ", n_cov,
          " of the ", T_pre, " pre-treatment periods in 'data' (",
          pre_times[1L], " to ", pre_times[T_pre], "), so this is a ",
          "predictor-path fit, not the outcomes-only fit. The outer V fit ",
          "still evaluates the pre-treatment period in full",
          if (is.null(v_window)) {
            paste0("; restrict it with v_window, or subset 'data', to move ",
                   "the pre-treatment window")
          } else {
            ""
          },
          "."
        )
      }
    }
  }

  if (use_cov) {
    if (is.null(data) || is.null(id_var) || is.null(time_var)) {
      stop(
        "'data', 'id_var', and 'time_var' must be provided when using predictors.",
        call. = FALSE
      )
    }
    .check_pred_windows(predictors, pan$times, T_pre)
    pm <- build_predictor_matrices(
      data     = data,
      id_var   = id_var,
      time_var = time_var,
      units    = pan$units,
      idx_co   = idx_co,
      idx_tr   = idx_tr,
      predictors = predictors,
      # The same pooling weights the outcome rows use, so X1 and Z1 describe
      # one target series rather than two.
      tr_weights = w_tr
    )
    X0         <- pm$X0
    X1         <- pm$X1
    pred_names <- pm$pred_names
  } else {
    X0         <- Y_co_pre
    X1         <- drop(Y_tr_pre)
    pred_names <- paste0("V", seq_len(T_pre))
    # Announce the default model once per fit. Suppressed for the routed
    # per-period spec (its own message already fired) and for oos, whose
    # train/validation split is described under v_selection.
    if ((is.null(predictors) || length(predictors) == 0L) &&
        v_selection == "insample") {
      message(
        "predictors = NULL: using the outcome in each of the ", T_pre,
        " pre-treatment periods as predictors (outcomes-only SCM)."
      )
    }
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

  # -- Outer evaluation window (v_window) --------------------------------------
  # Restricts only the rows of Z on which the outer V optimisation evaluates
  # the pre-treatment fit; X and the reported full-window loss are untouched.
  z_rows <- NULL
  if (!is.null(v_window)) {
    if (v_selection == "oos") {
      stop(
        "v_window cannot be combined with v_selection = 'oos', which manages ",
        "its own train/validation split of the pre-treatment window.",
        call. = FALSE
      )
    }
    idx <- match(v_window, pre_times)
    if (anyNA(idx)) {
      stop(
        "v_window contains values outside the pre-treatment window: ",
        paste(utils::head(v_window[is.na(idx)], 5L), collapse = ", "),
        ". Pre-treatment times run from ", pre_times[1L], " to ",
        pre_times[T_pre], ".",
        call. = FALSE
      )
    }
    z_rows <- sort(unique(as.integer(idx)))
    if (length(z_rows) < 2L) {
      stop("v_window must cover at least 2 pre-treatment periods.",
           call. = FALSE)
    }
  }

  # OOS handling (Abadie 2021 S.3.2):
  #  * outcomes-only: proper train/validation split via .scm_oos_outcomes()
  #    (candidate W(V) must not see validation-period outcomes).
  #  * user predictors: the predictor matrix is fixed (cannot be re-measured
  #    per window), so only the MSPE evaluation window is restricted.
  oos_outcomes <- (v_selection == "oos") && !use_cov
  t_train      <- if (v_selection == "oos" && use_cov) T_pre %/% 2L else -1L

  # Determine effective outer optimiser. The outer problem is non-convex;
  # with user predictors a single start can settle in a poor basin, so
  # "auto" runs the deterministic multi-start there. Outcomes-only fits keep
  # the single-start coordinate descent (empirically matches the Synth
  # reference solution, and multi-start would multiply the hot path's cost).
  effective_outer <- switch(v_optim,
    "auto"          = if (use_cov) "multistart" else "coord_descent",
    "multistart"    = "multistart",
    "bfgs"          = "bfgs",
    "coord_descent" = "coord_descent"
  )
  if (effective_outer == "multistart" && !use_cov) {
    stop(
      "v_optim = 'multistart' requires a predictor specification ",
      "(predictors = list(pred(...))). The outcomes-only path uses the ",
      "single-start coordinate descent.",
      call. = FALSE
    )
  }

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
      .scm_bfgs_outer(X0, X1, Z0, Z1, t_train = t_train, z_rows = z_rows)
    } else {
      # Outcomes-only fits use the cheap bordered-KKT face solve; predictor
      # fits keep the scale-robust null-space solver (cheap_face = FALSE).
      scm_weights_cpp(X0, X1, Z0, Z1, t_train = t_train, z_rows = z_rows,
                      multistart = (effective_outer == "multistart"),
                      cheap_face = !use_cov, wolfe = use_wolfe)
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
      .scm_bfgs_outer(X0, X1, Z0, Z1, t_train = t_train, z_rows = z_rows)
    } else {
      scm_weights_cpp(X0, X1, Z0, Z1, t_train = t_train, z_rows = z_rows,
                      multistart = (effective_outer == "multistart"),
                      cheap_face = !use_cov)
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
    # Raw per-unit series behind the fitted average, under the same name the
    # staggered cohort_fits use. With one treated unit it is Y_treat as a
    # one-column matrix.
    Y_treat_mat     = Y_treat_mat,
    n_treated       = n_treated,
    treated_weights = if (is.null(agg_w)) NULL
                      else stats::setNames(w_tr,
                             as.character(pan$units[idx_tr])),
    idx_tr          = idx_tr,
    gap             = Y_treat - Y_synth,
    times           = pan$times,
    T_pre           = T_pre,
    estimate        = mean((Y_treat - Y_synth)[-(seq_len(T_pre))]),
    loss            = loss,
    predictor_table = predictor_table,
    X0_mat          = if (use_cov) X0 else NULL,
    X1_vec          = if (use_cov) X1 else NULL,
    v_rows          = v_rows_used,
    v_window        = if (!is.null(z_rows)) pan$times[z_rows] else NULL,
    z_rows          = z_rows,
    v_optim_effective = effective_outer,
    qp_solver       = qp_solver,
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
    # Outcomes-only OOS: cheap bordered-KKT face solve (no predictors).
    scm_weights_cpp(X0_train, X1_train, Z0_val, Z1_val, t_train = -1L,
                    cheap_face = TRUE)
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
# DEPRECATED: reachable only via the deprecated v_optim = "bfgs". Remove this
# helper and the "bfgs" branches in fit_scm_cpp()/.scm_oos_outcomes()/
# .fit_scm_staggered() when "bfgs" is dropped in the next major release.
.scm_bfgs_outer <- function(X0, X1, Z0, Z1, t_train = -1L, z_rows = NULL) {
  k     <- nrow(X0)
  T_pre <- nrow(Z0)

  # Evaluation window -- same precedence as scm_weights_cpp
  # (explicit v_window rows, else the OOS split, else the full window).
  do_oos <- is.null(z_rows) && (t_train > 0L && t_train < T_pre)
  if (!is.null(z_rows)) {
    Z0_eval <- Z0[z_rows, , drop = FALSE]
    Z1_eval <- Z1[z_rows]
  } else if (do_oos) {
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
# A block is a cohort by default -- treated units sharing an adoption date are
# averaged into a single pseudo-unit, which Ben-Michael, Feller & Rothstein
# (2022, S.4.1) justify by the DGP heterogeneity terms (S_rho, S_k) vanishing
# across units treated at the same time. pool_cohorts = FALSE makes a block a
# treated unit instead, as in their S.4 formulation; nothing below depends on
# which, so the two differ only in how the blocks are assembled.
# Only predictors = NULL is supported (all pre-treatment outcomes used as X).
# fixedeff = TRUE demeans every unit by its own pre-treatment mean within
# each cohort (intercept shift; Ben-Michael et al. 2022 S.5.1) and restores
# the level in the reported Y_synth via the closed-form alpha of their Eq. 8.
.fit_scm_staggered <- function(pan,
                               Y_matrix             = NULL,
                               unit_weights         = NULL,
                               v_selection          = "insample",
                               v_optim              = "coord_descent",
                               donor_mspe_threshold = Inf,
                               lambda_pen           = NULL,
                               control_group        = "clean",
                               pool_cohorts         = TRUE,
                               fixedeff             = FALSE) {
  Y       <- if (!is.null(Y_matrix)) Y_matrix else pan$Y
  TT      <- nrow(Y)
  T_adopt <- pan$T_adopt
  idx_tr  <- pan$idx_treat
  idx_co  <- pan$idx_control

  cohorts <- sort(unique(T_adopt[idx_tr]))

  # Same block granularity as .fit_scm_staggered_pooled(): a cohort, or a
  # single treated unit when the caller asked for the S.4 formulation.
  # unname() so cohort_fits is a plain list on both paths -- idx_tr carries
  # the panel's column names, and the pooled path drops them.
  block_ids <- unname(if (pool_cohorts) as.list(cohorts) else as.list(idx_tr))

  cohort_list <- lapply(block_ids, function(key) {
    if (pool_cohorts) {
      g        <- key
      idx_tr_g <- idx_tr[T_adopt[idx_tr] == g]
    } else {
      idx_tr_g <- key
      g        <- T_adopt[key]
    }
    T_pre_g  <- g - 1L
    T_post_g <- TT - T_pre_g

    if (control_group == "never_treated") {
      idx_co_g <- idx_co
    } else {
      future_tr <- idx_tr[!is.na(T_adopt[idx_tr]) & T_adopt[idx_tr] > g]
      idx_co_g  <- c(idx_co, future_tr)
    }

    blk <- .scm_block_target(Y, idx_tr_g, pan$units, unit_weights)
    if (length(idx_co_g) < 2L || T_pre_g < 2L || blk$mass <= 0) {
      warning(sprintf(
        "SCM staggered: %s skipped (T_pre=%d, N_co=%d, mass=%g).",
        if (pool_cohorts) sprintf("cohort g=%d", g) else
          sprintf("treated unit %s", as.character(pan$units[idx_tr_g])),
        T_pre_g, length(idx_co_g), blk$mass
      ), call. = FALSE)
      return(NULL)
    }

    pre_rows  <- seq_len(T_pre_g)
    post_rows <- (T_pre_g + 1L):TT

    # The block's target is the mass-weighted average of its treated units,
    # so the series fitted here is the one the aggregation weights.
    Y1_g     <- blk$y
    Y_co_all <- Y[, idx_co_g, drop = FALSE]

    # Intercept shift: fit on series demeaned by their own pre-window means
    # (the raw copies are kept so Y_synth can be reported on the raw scale)
    Y1_g_raw     <- Y1_g
    Y_co_all_raw <- Y_co_all
    if (fixedeff) {
      Y1_g     <- Y1_g - mean(Y1_g[pre_rows])
      Y_co_all <- sweep(Y_co_all, 2L,
                        colMeans(Y_co_all[pre_rows, , drop = FALSE]))
    }
    Y_co_pre <- Y_co_all[pre_rows, , drop = FALSE]  # T_pre_g x N_co_g
    Y_tr_pre <- Y1_g[pre_rows]

    # Donor pool filtering per cohort (same rule as the sharp path)
    excl_g <- character(0L)
    if (is.finite(donor_mspe_threshold) && donor_mspe_threshold < Inf) {
      ind_mspe <- colMeans((Y_tr_pre - Y_co_pre)^2)
      keep     <- ind_mspe <= donor_mspe_threshold * min(ind_mspe)
      if (sum(keep) >= 2L) {
        excl_g       <- colnames(Y_co_pre)[!keep]
        idx_co_g     <- idx_co_g[keep]
        Y_co_all     <- Y_co_all[, keep, drop = FALSE]
        Y_co_all_raw <- Y_co_all_raw[, keep, drop = FALSE]
        Y_co_pre     <- Y_co_pre[, keep, drop = FALSE]
      }
    }

    t_train_g <- if (v_selection == "oos") T_pre_g %/% 2L else -1L
    # Staggered SCM is outcomes-only, so "auto" resolves to the single-start
    # coordinate descent ("multistart" is rejected at the fit entry point).
    effective_outer_g <- switch(v_optim,
      "auto"          = "coord_descent",
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
          # Staggered SCM is outcomes-only: cheap bordered-KKT face solve.
          scm_weights_cpp(Y_co_pre, Y_tr_pre, Y_co_pre, Y_tr_pre,
                          t_train = -1L, cheap_face = TRUE)
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
                          t_train = -1L, cheap_face = TRUE)
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
      if (fixedeff) {
        alpha_g <- mean(Y1_g_raw[pre_rows]) -
          sum(unit_w * colMeans(Y_co_all_raw[pre_rows, , drop = FALSE]))
        Y_synth_g <- drop(Y_co_all_raw %*% unit_w) + alpha_g
      } else {
        alpha_g   <- 0
        Y_synth_g <- drop(Y_co_all_raw %*% unit_w)
      }
      ATT_g <- mean((Y1_g_raw - Y_synth_g)[post_rows])

      list(
        cohort          = g,
        # Named when the block is one treated unit, which is what
        # scm_rank_pval() dispatches on.
        unit            = if (length(idx_tr_g) == 1L)
          as.character(pan$units[idx_tr_g]) else NA_character_,
        n_treated       = length(idx_tr_g),
        T_pre           = T_pre_g,
        T_post          = T_post_g,
        estimate        = ATT_g,
        weight          = blk$mass * T_post_g,
        mass            = blk$mass,
        unit_masses     = blk$masses,
        unit_weights    = unit_w,
        Y_synth         = Y_synth_g,
        Y_treat         = Y1_g_raw,
        Y_treat_mat     = Y[, idx_tr_g, drop = FALSE],
        idx_tr          = idx_tr_g,
        idx_co          = idx_co_g,
        alpha           = alpha_g,
        excluded_donors = excl_g,
        lambda_pen      = lambda_pen_used
      )
    }, error = function(e) {
      warning(sprintf("SCM staggered: %s failed: %s",
                      if (pool_cohorts) sprintf("cohort g=%d", g) else
                        sprintf("treated unit %s",
                                as.character(pan$units[idx_tr_g])),
                      conditionMessage(e)), call. = FALSE)
      NULL
    })
    fit_g
  })

  valid <- !vapply(cohort_list, is.null, logical(1L))
  if (!any(valid)) {
    stop("All cohort-level SCM fits failed; see the preceding warnings for ",
         "per-cohort reasons.", call. = FALSE)
  }
  .aggregate_cohorts(cohort_list[valid], unit_weights, pan$units)
}

# -- Internal: partially pooled SCM for staggered adoption -------------------
# Ben-Michael, Feller & Rothstein (2022, JRSS-B) Eq. 6: choose all cohort
# weight vectors jointly to minimise
#   nu * (q_pool / q_pool_sep)^2 + (1 - nu) * (q_sep / q_sep_sep)^2
# where q_sep is the root mean squared per-cohort pre-treatment imbalance and
# q_pool the imbalance of the average placebo gap, aligned in event time with
# missing lags zero-padded (their tau_hat_j(-l) = 0 for l > L_j convention).
# Treated units within a cohort are averaged into one pseudo-unit (their
# Appendix A.2 modification), so "unit" in the paper maps to "cohort" here.
# nu = 0 is separate SCM with uniform lag weights (no V optimisation),
# nu = 1 fully pooled, nu = "auto" the paper's heuristic
#   nu_hat = ||pooled placebo gap||_2 / mean_g ||cohort placebo gap||_2.
# The joint QP is solved by cyclic block coordinate descent over cohorts:
# with the other cohorts fixed, the objective in gamma_g is a standard
# simplex QP handled by solve_simplex_qp() (convex objective, separable
# simplex constraints, so BCD converges to the global optimum).
.fit_scm_staggered_pooled <- function(pan,
                                      Y_matrix      = NULL,
                                      unit_weights  = NULL,
                                      nu            = "auto",
                                      fixedeff      = FALSE,
                                      control_group = "clean",
                                      pool_cohorts  = TRUE,
                                      max_sweeps    = 200L,
                                      tol           = 1e-7) {
  # tol is the relative objective-improvement threshold for the BCD loop.
  # A weight-change criterion would never trigger with warm-started FISTA
  # blocks, whose solutions jitter at the inner solver tolerance (~1e-6);
  # past that point extra sweeps change q_pool below the solver noise floor.
  Y       <- if (!is.null(Y_matrix)) Y_matrix else pan$Y
  TT      <- nrow(Y)
  T_adopt <- pan$T_adopt
  idx_tr  <- pan$idx_treat
  idx_co  <- pan$idx_control
  cohorts <- sort(unique(T_adopt[idx_tr]))

  # A block is what the pooling parameter trades off against the pooled fit.
  # BFR S.4 makes it a treated unit, with Gamma an N-by-J matrix; their
  # Appendix A.2 adds a layer that fully pools units sharing an adoption date,
  # justified in S.4.1 by there being no across-unit variation within a
  # cohort. pool_cohorts picks between the two, and the solver below is
  # indifferent: a block is a target series, its donors, and its lag window.
  block_ids <- unname(if (pool_cohorts) as.list(cohorts) else as.list(idx_tr))

  # ---- assemble blocks in event-time (lag) order ---------------------------
  blocks <- lapply(block_ids, function(key) {
    if (pool_cohorts) {
      g        <- key
      idx_tr_g <- idx_tr[T_adopt[idx_tr] == g]
    } else {
      idx_tr_g <- key
      g        <- T_adopt[key]
    }
    T_pre_g  <- g - 1L
    T_post_g <- TT - T_pre_g

    if (control_group == "never_treated") {
      idx_co_g <- idx_co
    } else {
      future_tr <- idx_tr[!is.na(T_adopt[idx_tr]) & T_adopt[idx_tr] > g]
      idx_co_g  <- c(idx_co, future_tr)
    }
    blk <- .scm_block_target(Y, idx_tr_g, pan$units, unit_weights)
    if (length(idx_co_g) < 2L || T_pre_g < 2L || blk$mass <= 0) {
      warning(sprintf(
        "SCM staggered: %s skipped (T_pre=%d, N_co=%d, mass=%g).",
        if (pool_cohorts) sprintf("cohort g=%d", g) else
          sprintf("treated unit %s", as.character(pan$units[idx_tr_g])),
        T_pre_g, length(idx_co_g), blk$mass
      ), call. = FALSE)
      return(NULL)
    }

    pre_rows  <- seq_len(T_pre_g)
    post_rows <- (T_pre_g + 1L):TT
    Y1_raw    <- blk$y
    Y_co_raw  <- Y[, idx_co_g, drop = FALSE]

    Y1_fit   <- Y1_raw
    Y_co_fit <- Y_co_raw
    if (fixedeff) {
      Y1_fit   <- Y1_raw - mean(Y1_raw[pre_rows])
      Y_co_fit <- sweep(Y_co_raw, 2L,
                        colMeans(Y_co_raw[pre_rows, , drop = FALSE]))
    }

    # lag l (= 1..T_pre_g) is calendar row T_pre_g + 1 - l, so that lag
    # vectors from different cohorts align in event time for the pooled fit
    lag_rows <- rev(pre_rows)
    list(
      cohort = g, idx_tr = idx_tr_g, idx_co = idx_co_g,
      mass = blk$mass, masses = blk$masses,
      T_pre = T_pre_g, T_post = T_post_g,
      pre_rows = pre_rows, post_rows = post_rows,
      Y1_raw = Y1_raw, Y_co_raw = Y_co_raw,
      b = Y1_fit[lag_rows],                       # length T_pre_g
      X = Y_co_fit[lag_rows, , drop = FALSE]      # T_pre_g x N_co_g
    )
  })
  blocks <- blocks[!vapply(blocks, is.null, logical(1L))]
  if (length(blocks) == 0L) {
    stop("All cohort-level SCM fits failed; see the preceding warnings for ",
         "per-cohort reasons.", call. = FALSE)
  }

  J    <- length(blocks)
  L_g  <- vapply(blocks, `[[`, integer(1L), "T_pre")
  L    <- max(L_g)
  # Block mass: the cohort size n_g of BFR Appendix A.2, or the mass the
  # caller gave its treated units. Appendix A.2 writes the cohort problem on
  # the SUM over the cohort's units -- its target is
  # sum_i 1{T_i = T(g)} Y_i and its weights satisfy sum_i gamma_ig = n_g --
  # so the block residual is n_g times the mean gap. Carrying the mean here
  # and the mass as a multiplier is the same objective, and keeps Y_synth and
  # ATT_g on the outcome's own scale for reporting.
  m_g  <- vapply(blocks, `[[`, numeric(1L), "mass")

  # ---- separate SCM baseline (nu = 0; also the normalisation constants) ----
  # The mass cancels out of each block's own least-squares problem, so the
  # separate solution does not depend on it.
  Gamma <- lapply(blocks, function(bl) {
    as.numeric(solve_simplex_qp(crossprod(bl$X), crossprod(bl$X, bl$b)))
  })

  block_gap <- function(Gamma, j)
    blocks[[j]]$b - as.numeric(blocks[[j]]$X %*% Gamma[[j]])

  # q^{pool,cohort} of Appendix A.2: the mass-weighted average gap at each
  # event-time lag, zero-padded past a block's own window, over max_g L_g.
  pooled_gap <- function(Gamma) {
    P <- numeric(L)
    for (j in seq_len(J)) {
      lg <- L_g[j]
      P[seq_len(lg)] <- P[seq_len(lg)] + m_g[j] * block_gap(Gamma, j) / J
    }
    P
  }
  # L2 norms of the (unweighted) per-block gaps; the mass enters where they
  # are used, so this stays the plain imbalance of each block's own fit.
  sep_norms <- function(Gamma)
    vapply(seq_len(J), function(j) sqrt(sum(block_gap(Gamma, j)^2)),
           numeric(1L))
  # q^{sep,cohort} squares the block SUM, hence the mass squared.
  q_sep_of  <- function(Gamma) sqrt(mean(m_g^2 * sep_norms(Gamma)^2 / L_g))
  q_pool_of <- function(Gamma) sqrt(mean(pooled_gap(Gamma)^2))

  C_sep  <- q_sep_of(Gamma)
  C_pool <- q_pool_of(Gamma)

  # Heuristic nu (BFR S.4.2, p.367): the ratio of the pooled fit of separate
  # SCM to its average block-level fit,
  #   sqrt(L) q_pool / ( (1/J) sum_j sqrt(L_j) q_j ) =
  #   || sum_j m_j gap_j || / sum_j m_j ||gap_j||,
  # bounded above by 1 via the triangle inequality.
  r_sep   <- m_g * sep_norms(Gamma)
  nu_heur <- if (mean(r_sep) < .Machine$double.eps) 0 else
    min(1, sqrt(sum(pooled_gap(Gamma)^2)) / mean(r_sep))
  nu_num <- if (identical(nu, "auto")) nu_heur else as.numeric(nu)

  # ---- block coordinate descent (skipped when the pooled term is inactive) --
  sweeps_used <- 0L
  if (nu_num > 0 && C_pool > .Machine$double.eps &&
      C_sep > .Machine$double.eps) {
    A_pool <- nu_num / (C_pool^2 * L)
    P      <- pooled_gap(Gamma)
    XtX_l  <- lapply(blocks, function(bl) crossprod(bl$X))
    Xtb_l  <- lapply(blocks, function(bl) as.numeric(crossprod(bl$X, bl$b)))
    obj_prev <- Inf
    for (s in seq_len(max_sweeps)) {
      obj_sep <- 0
      for (j in seq_len(J)) {
        bl <- blocks[[j]]
        lg <- L_g[j]
        mj <- m_g[j]
        gap_j <- bl$b - as.numeric(bl$X %*% Gamma[[j]])
        # this block's share of the pooled gap, removed so R is the rest
        R     <- P[seq_len(lg)] - mj * gap_j / J
        A_sep <- (1 - nu_num) * mj^2 / (C_sep^2 * J * lg)
        u     <- R + mj * bl$b / J
        Q     <- 2 * (A_pool * mj^2 / J^2 + A_sep) * XtX_l[[j]]
        cvec  <- 2 * (A_pool * mj / J * as.numeric(crossprod(bl$X, u)) +
                        A_sep * Xtb_l[[j]])
        # warm start from the previous block solution: successive sweeps
        # resolve near-identical QPs, so FISTA converges in a few iterations
        w_new <- as.numeric(solve_simplex_qp(Q, cvec, x0 = Gamma[[j]]))
        Gamma[[j]] <- w_new
        gap_new <- bl$b - as.numeric(bl$X %*% w_new)
        P[seq_len(lg)] <- R + mj * gap_new / J
        obj_sep <- obj_sep + A_sep * sum(gap_new^2)
      }
      obj <- A_pool * sum(P^2) + obj_sep
      sweeps_used <- s
      if (obj_prev - obj < tol * max(obj_prev, .Machine$double.eps)) break
      obj_prev <- obj
    }
  }

  # ---- per-cohort reporting on the raw outcome scale ------------------------
  r <- lapply(seq_len(J), function(j) {
    bl     <- blocks[[j]]
    unit_w <- Gamma[[j]]
    names(unit_w) <- colnames(bl$Y_co_raw)
    if (fixedeff) {
      alpha_g <- mean(bl$Y1_raw[bl$pre_rows]) -
        sum(unit_w * colMeans(bl$Y_co_raw[bl$pre_rows, , drop = FALSE]))
    } else {
      alpha_g <- 0
    }
    Y_synth_g <- drop(bl$Y_co_raw %*% unit_w) + alpha_g
    ATT_g     <- mean((bl$Y1_raw - Y_synth_g)[bl$post_rows])
    list(
      cohort          = bl$cohort,
      unit            = if (length(bl$idx_tr) == 1L)
        as.character(pan$units[bl$idx_tr]) else NA_character_,
      n_treated       = length(bl$idx_tr),
      T_pre           = bl$T_pre,
      T_post          = bl$T_post,
      estimate        = ATT_g,
      # The same mass the objective weighted this block by, so the fit and
      # the aggregation describe one estimand (Clarke et al. 2024).
      weight          = bl$mass * bl$T_post,
      mass            = bl$mass,
      unit_masses     = bl$masses,
      unit_weights    = unit_w,
      Y_synth         = Y_synth_g,
      Y_treat         = bl$Y1_raw,
      Y_treat_mat     = Y[, bl$idx_tr, drop = FALSE],
      idx_tr          = bl$idx_tr,
      idx_co          = bl$idx_co,
      alpha           = alpha_g,
      excluded_donors = character(0L),
      lambda_pen      = NA_real_
    )
  })

  c(.aggregate_cohorts(r, unit_weights, pan$units), list(
    pooling          = list(
      nu              = nu_num,
      nu_heuristic    = nu_heur,
      fixedeff        = fixedeff,
      pool_cohorts    = pool_cohorts,
      q_sep           = q_sep_of(Gamma),
      q_pool          = q_pool_of(Gamma),
      q_sep_separate  = C_sep,
      q_pool_separate = C_pool,
      sweeps          = sweeps_used
    )
  ))
}

#' Wild Bootstrap Inference for Staggered SCM
#'
#' Confidence intervals and p-values for the aggregate ATT of a staggered
#' SCM fit via the weighted multiplier (wild) bootstrap of Ben-Michael,
#' Feller & Rothstein (2022, Section 5.3), adapting Otsu & Rai (2017).
#' The aggregate ATT is written as a weighted average of per-treated-unit
#' effect contributions; each bootstrap draw perturbs those contributions
#' with independent golden-ratio two-point multipliers (mean 0, variance 1)
#' while donor weights and outcomes are kept fixed.
#'
#' Works with both staggered SCM paths (the partially pooled default and
#' `nu = "separate"`) and honours the cohort aggregation weights
#' `N_treated x T_post`. Each treated unit is resampled at its own mass --
#' its `weights` entry, or one -- so the contributions sum back to the ATT
#' the fit reports. For intercept-shifted fits (`fixedeff = TRUE`)
#' the per-unit contributions are computed in difference-in-differences
#' form, i.e. each treated unit is demeaned by its own pre-treatment mean.
#'
#' With very few treated units the multiplier distribution has few atoms,
#' so the bootstrap is unreliable; a warning is issued below 5 treated
#' units.
#'
#' @param fit A staggered SCM fit from [scm_fit()] (`method = "scm"` on a
#'   panel with multiple adoption times). Sharp fits are rejected: use
#'   [mspe_ratio_pval()] or [conformal_inference()] instead.
#' @param method Only `"wild_bootstrap"` is available.
#' @param n_boot Number of bootstrap draws. Default 1000.
#' @param level Confidence level. Default 0.95.
#' @param alternative Direction of the alternative hypothesis for the
#'   p-value: `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param seed Optional RNG seed.
#' @return A `coresynth_inference` object with the standard fields
#'   (`estimate`, `se`, `p_value`, `ci_lower`, `ci_upper`, `method`,
#'   `staggered`, `n_controls`, `alternative`, `boot_ests`), compatible
#'   with [tidy.coresynth_inference()] and [glance.coresynth_inference()].
#'   `n_treated` additionally records the number of treated units
#'   resampled by the multipliers.
#' @references Ben-Michael, E., Feller, A., & Rothstein, J. (2022).
#'   Synthetic controls with staggered adoption. *JRSS-B*, 84(2), 351-381.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- expand.grid(time = 1:20, id = paste0("u", 1:12))
#' dat$y <- rnorm(nrow(dat)) + as.numeric(factor(dat$id))
#' dat$d <- as.integer(
#'   (dat$id == "u1" & dat$time > 10) | (dat$id == "u2" & dat$time > 14)
#' )
#' fit <- scm_fit(y ~ d | id + time, data = dat, method = "scm")
#' scm_inference(fit, n_boot = 200, seed = 1)
#' }
#' @export
scm_inference <- function(fit,
                          method      = "wild_bootstrap",
                          n_boot      = 1000L,
                          level       = 0.95,
                          alternative = c("two.sided", "greater", "less"),
                          seed        = NULL) {
  method      <- match.arg(method)
  alternative <- match.arg(alternative)

  if (!inherits(fit, "coresynth_scm"))
    stop("scm_inference() requires a coresynth fit with method = 'scm'.",
         call. = FALSE)
  if (!inherits(fit, "coresynth_staggered"))
    stop(
      "scm_inference() supports staggered SCM fits only.\n",
      "For a sharp fit use mspe_ratio_pval() (in-space permutation test), ",
      "or, when several units are treated at once, scm_rank_pval(); ",
      "conformal_inference() works for either.",
      call. = FALSE
    )
  cf <- fit$cohort_fits
  if (is.null(cf[[1L]]$Y_treat_mat))
    stop("cohort_fits lack per-unit treated outcomes (Y_treat_mat). ",
         "Re-run scm_fit() with the current version.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  TT       <- length(fit$times)
  fixedeff <- isTRUE(fit$fixedeff)

  # per-treated-unit average post-treatment effects tau_i (BFR Eq. 12) and
  # aggregation weights omega_i reproducing sum_g N_g * T_post_g weighting
  tau_i   <- numeric(0)
  omega_i <- numeric(0)
  for (b in cf) {
    pre  <- seq_len(b$T_pre)
    post <- (b$T_pre + 1L):TT
    if (fixedeff) {
      # DiD form: each treated unit demeaned by its own pre-treatment mean
      tau_g <- (colMeans(b$Y_treat_mat[post, , drop = FALSE]) -
                  colMeans(b$Y_treat_mat[pre, , drop = FALSE])) -
        (mean(b$Y_synth[post]) - mean(b$Y_synth[pre]))
    } else {
      tau_g <- colMeans(b$Y_treat_mat[post, , drop = FALSE] -
                          b$Y_synth[post])
    }
    tau_i <- c(tau_i, tau_g)
    # Each unit carries its own mass (its `weights` entry, or 1) times the
    # block's post window, so the contributions sum back to the block weight
    # and reproduce the ATT the fit reports. Fits from before `unit_masses`
    # split the block weight evenly, which is the same thing under the equal
    # weighting those fits had.
    m_i <- b$unit_masses %||% rep(b$weight / b$T_post / b$n_treated,
                                  b$n_treated)
    omega_i <- c(omega_i, m_i * b$T_post)
  }
  omega_i <- omega_i / sum(omega_i)
  att_hat <- fit$estimate

  N_tr <- length(tau_i)
  if (N_tr < 5L)
    warning("Wild bootstrap with fewer than 5 treated units is unreliable; ",
            "interpret the interval with caution.", call. = FALSE)

  # golden-ratio two-point multipliers: E[W] = 0, E[W^2] = 1
  s5    <- sqrt(5)
  w_neg <- -(s5 - 1) / 2
  w_pos <- (s5 + 1) / 2
  p_neg <- (s5 + 1) / (2 * s5)

  dev <- tau_i - att_hat
  S   <- replicate(n_boot, {
    W <- sample(c(w_neg, w_pos), N_tr, replace = TRUE,
                prob = c(p_neg, 1 - p_neg))
    sum(omega_i * W * dev)
  })

  alpha_lv <- 1 - level
  se       <- stats::sd(S)
  ci_lower <- att_hat - unname(stats::quantile(S, 1 - alpha_lv / 2))
  ci_upper <- att_hat - unname(stats::quantile(S, alpha_lv / 2))
  p_value  <- switch(alternative,
    two.sided = (1 + sum(abs(S) >= abs(att_hat))) / (n_boot + 1),
    greater   = (1 + sum(S >= att_hat)) / (n_boot + 1),
    less      = (1 + sum(S <= att_hat)) / (n_boot + 1)
  )

  n_controls <- length(unique(unlist(lapply(cf, `[[`, "idx_co"))))

  structure(list(
    estimate    = att_hat,
    se          = se,
    p_value     = p_value,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    method      = method,
    staggered   = TRUE,
    n_controls  = n_controls,
    n_treated   = N_tr,
    alternative = alternative,
    boot_ests   = att_hat + S
  ), class = "coresynth_inference")
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
#' The placebo refits mirror the treated fit's outer optimiser and evaluation
#' window: a fit estimated with the multi-start outer search (`v_optim =
#' "multistart"`, or the `"auto"` default with predictors) or with a
#' `v_window` runs every placebo unit through the same configuration, keeping
#' the permutation statistic exchangeable across units.
#'
#' @param fit A `coresynth` object from [scm_fit()] with `method = "scm"`.
#' @param mspe_threshold Minimum pre-treatment MSPE for including a control
#'   unit in the two-sided test. Ignored for one-sided tests.
#'   Default: 0 (no filtering).
#' @param max_iter Passed to [scm_placebo_cpp()] (outcomes-only fits) or
#'   [scm_placebo_x_cpp()] (covariate fits). Default 100L.
#' @param tol      Passed to [scm_placebo_cpp()] or [scm_placebo_x_cpp()].
#'   Default 1e-4.
#' @param use_covariates Controls which predictor specification the placebo
#'   refits use. Default `NULL` (recommended) mirrors the treated fit: if the
#'   fit was estimated with a `predictors` specification, each placebo unit is
#'   refit with that same specification (the Abadie et al. 2010 / `Synth`
#'   convention -- treated and placebo statistics are computed under one
#'   common spec); outcomes-only fits use the fast C++ outcomes-only placebo.
#'   Set `TRUE`/`FALSE` to force either path; note that `FALSE` on a
#'   covariate fit compares a covariate-based treated statistic against
#'   outcomes-only placebo statistics, which breaks the exchangeability
#'   logic of the permutation test.
#' @param alternative Direction of the alternative hypothesis:
#'   `"two.sided"` (default) uses the MSPE ratio statistic;
#'   `"greater"` tests whether the treatment increased the outcome;
#'   `"less"` tests whether the treatment decreased the outcome.
#'   One-sided tests use the signed ATT as the test statistic.
#' @param mspe_prune Discard placebo runs whose pre-treatment MSPE exceeds
#'   this multiple of the treated unit's, following Abadie, Diamond &
#'   Hainmueller (2010, Figures 5-7), who use 20, 5 and 2. A run the synthetic
#'   control could not fit before treatment carries no information about how
#'   rare a large post-treatment gap is, so leaving it in the permutation
#'   distribution only adds noise. `Inf` (default) keeps every run. This is a
#'   different control from `mspe_threshold`, which is a floor guarding the
#'   ratio's denominator against near-zero pre-treatment MSPE.
#' @param statistic Which quantity to rank, for two-sided tests.
#'   `"ratio"` (default) is the post/pre MSPE ratio of Abadie, Diamond &
#'   Hainmueller (2010, Figure 8), whose denominator adjusts for how well each
#'   unit was fitted and so needs no cutoff. `"post_mspe"` is the alternative
#'   of Abadie (2021, S.3.5): rank the post-treatment MSPE itself, having
#'   discarded the ill-fitting runs through `mspe_prune`. Pair it with a
#'   finite `mspe_prune`; without one, placebo runs that were never fitted
#'   dominate the distribution.
#' @return An object of class `scm_placebo` (a list) with:
#'   * `p_value`: Permutation p-value between 0 and 1
#'   * `mspe_ratio_treated`: MSPE_post / MSPE_pre for the treated unit (two.sided only)
#'   * `mspe_ratios_all`: Named numeric vector (treated first, then controls); two.sided only
#'   * `placebo_effects`: Named N_co-vector of placebo ATT estimates
#'   * `treated_effect`: ATT estimate for the treated unit
#'   * `n_placebo_used`: Number of control units used
#'   * `n_placebo_pruned`: Number discarded by `mspe_prune`
#'   * `mspe_prune`, `statistic`: the settings the p-value was computed under
#'   * `gaps`: T x N_co matrix of placebo gap paths (unit minus its synthetic
#'     control over all periods), for the Abadie et al. (2010) Figure 4-7 plot
#'   * `treated_gap`: T-vector of the treated unit's gap path
#'   * `mspe_pre_treated`, `mspe_pre_placebo`: Pre-treatment MSPEs used for
#'     the relative pruning rule in [plot.scm_placebo()]
#'   * `times`, `T_pre`: Time axis metadata for plotting
#' @seealso [plot.scm_placebo()] for the placebo gap and MSPE ratio plots.
#' @export
mspe_ratio_pval <- function(
  fit,
  mspe_threshold = 0,
  max_iter = 100L,
  tol = 1e-4,
  use_covariates = NULL,
  alternative = c("two.sided", "greater", "less"),
  mspe_prune = Inf,
  statistic = c("ratio", "post_mspe")
) {
  alternative <- match.arg(alternative)
  statistic   <- match.arg(statistic)
  if (!is.numeric(mspe_prune) || length(mspe_prune) != 1L ||
      is.na(mspe_prune) || mspe_prune <= 0)
    stop("`mspe_prune` must be a single positive number (Inf keeps every ",
         "placebo run).", call. = FALSE)
  if (identical(statistic, "post_mspe") && !identical(alternative, "two.sided"))
    stop("statistic = 'post_mspe' measures the size of the post-treatment ",
         "gap, so it has no direction; use alternative = 'two.sided', or ",
         "statistic = 'ratio' for a one-sided test.", call. = FALSE)
  if (!inherits(fit, "coresynth_scm")) {
    stop("mspe_ratio_pval() requires a coresynth object with method = 'scm'.")
  }
  if (inherits(fit, "coresynth_staggered")) {
    stop(
      "mspe_ratio_pval() does not support staggered fits.\n",
      "Alternatives:\n",
      "  * Apply mspe_ratio_pval() to each cohort_fit individually after re-fitting it as sharp SCM.\n",
      "  * Use method = 'sdid' with sdid_inference(method = 'placebo') for staggered placebo inference.",
      call. = FALSE
    )
  }
  if (isTRUE((fit$n_treated %||% 1L) > 1L)) {
    stop(
      "mspe_ratio_pval() compares the treated unit against placebo runs on ",
      "single donors, so it is not valid for a fit pooling ",
      fit$n_treated, " treated units: their average is less noisy than any ",
      "one donor, which inflates its pre/post MSPE ratio. Simulated under no ",
      "treatment effect, the test rejected at 28% (3 treated units) and 53% ",
      "(6) against a nominal 10%.\n",
      "Alternatives:\n",
      "  * Use scm_rank_pval(), the rank-based permutation test of Abadie & ",
      "L'Hour (2021) for exactly this case: it ranks per-unit statistics, so ",
      "the averaged treated series never enters.\n",
      "  * Use conformal_inference(), which permutes over time rather than ",
      "across units.\n",
      "  * Fit each treated unit on its own (subset data) and test it there.",
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

  # Auto (default): mirror the treated fit's predictor specification so the
  # permutation compares statistics computed under one common spec.
  if (is.null(use_covariates)) {
    use_covariates <- !is.null(fit$X0_mat)
  } else if (use_covariates && is.null(fit$X0_mat)) {
    warning(
      "use_covariates=TRUE but fit has no covariates (X0_mat is NULL). ",
      "Falling back to outcomes-only placebo."
    )
    use_covariates <- FALSE
  }

  # Mirror the treated fit's optimiser and evaluation window in every placebo
  # refit -- the permutation test is only valid when all units run the same
  # algorithm. Old fit objects lack these fields and fall back to the
  # single-start, full-window path (their own fitting configuration).
  multistart_fit <- identical(fit$v_optim_effective, "multistart")
  z_rows_fit     <- fit$z_rows
  wolfe_fit      <- identical(fit$qp_solver, "wolfe")

  if (use_covariates) {
    plac <- scm_placebo_x_cpp(fit$X0_mat, fit$Y_co_pre, fit$Y_co_post,
                              max_iter = max_iter, tol = tol,
                              z_rows = z_rows_fit,
                              multistart = multistart_fit,
                              wolfe = wolfe_fit)
    # Armadillo vectors come back as N_co x 1 matrices; flatten for naming.
    # A placebo unit whose solver failed carries NaN through all fields.
    mspe_pre_co  <- as.numeric(plac$mspe_pre)
    mspe_post_co <- as.numeric(plac$mspe_post)
    keep      <- !is.na(mspe_pre_co) & mspe_pre_co > mspe_threshold
    ratios_co <- ifelse(keep, mspe_post_co / mspe_pre_co, NA_real_)
    co_effects <- as.numeric(plac$effects)
    gaps_mat  <- plac$gaps
  } else {
    plac      <- scm_placebo_cpp(fit$Y_co_pre, fit$Y_co_post,
                                  max_iter = max_iter, tol = tol,
                                  z_rows = z_rows_fit, wolfe = wolfe_fit)
    keep      <- plac$mspe_pre > mspe_threshold
    ratios_co <- ifelse(keep, plac$mspe_post / plac$mspe_pre, NA_real_)
    co_effects <- plac$effects
    # Armadillo vectors come back as N_co x 1 matrices; flatten for naming
    mspe_pre_co  <- as.numeric(plac$mspe_pre)
    mspe_post_co <- as.numeric(plac$mspe_post)
    gaps_mat     <- plac$gaps
  }

  co_names <- colnames(fit$Y_co_pre)
  if (!is.null(co_names)) {
    names(ratios_co)   <- co_names
    names(co_effects)  <- co_names
    names(mspe_pre_co) <- co_names
    colnames(gaps_mat) <- co_names
  }

  treated_effect <- fit$estimate

  # Abadie, Diamond & Hainmueller (2010) Figures 5-7: a placebo run that could
  # not be fitted before treatment says nothing about how rare a large
  # post-treatment gap is, so runs whose pre-treatment MSPE exceeds
  # `mspe_prune` times the treated unit's are dropped from the permutation
  # distribution. Their Figures use 20, 5 and 2.
  pruned <- is.finite(mspe_pre_co) &
    mspe_pre_co <= mspe_prune * mspe_pre_tr
  n_pruned <- sum(!pruned, na.rm = TRUE)
  if (is.finite(mspe_prune) && !any(pruned))
    stop("mspe_prune = ", mspe_prune, " discards every placebo run; the ",
         "treated unit's pre-treatment MSPE is ",
         formatC(mspe_pre_tr, format = "g", digits = 3),
         " and the best placebo's is ",
         formatC(min(mspe_pre_co, na.rm = TRUE), format = "g", digits = 3),
         ". Raise the cutoff.", call. = FALSE)

  if (alternative == "two.sided") {
    if (statistic == "post_mspe") {
      # Abadie (2021) S.3.5: with ill-fitting runs discarded, rank the
      # post-treatment fit itself rather than its ratio to the pre-treatment
      # fit, which is what the ratio statistic uses the denominator for.
      all_stats <- c(treated = mspe_post_tr, mspe_post_co)
      keep_stat <- c(TRUE, pruned)
    } else {
      all_stats <- c(treated = ratio_tr, ratios_co)
      keep_stat <- c(TRUE, keep & pruned)
    }
    valid   <- keep_stat & !is.na(all_stats)
    p_value <- mean(all_stats[valid] >= all_stats[["treated"]], na.rm = TRUE)
    n_used  <- sum(valid[-1L], na.rm = TRUE)
    all_ratios <- c(treated = ratio_tr, ratios_co)
  } else {
    all_effects <- c(treated = treated_effect, co_effects)
    valid       <- c(TRUE, pruned) & !is.na(all_effects)
    if (alternative == "greater") {
      p_value <- mean(all_effects[valid] >= all_effects[["treated"]], na.rm = TRUE)
    } else {
      p_value <- mean(all_effects[valid] <= all_effects[["treated"]], na.rm = TRUE)
    }
    all_ratios <- c(treated = ratio_tr, ratios_co)
    n_used     <- sum(valid[-1L], na.rm = TRUE)
  }

  structure(
    list(
      p_value            = p_value,
      mspe_ratio_treated = ratio_tr,
      mspe_ratios_all    = all_ratios,
      placebo_effects    = co_effects,
      treated_effect     = treated_effect,
      n_placebo_used     = n_used,
      gaps               = gaps_mat,
      treated_gap        = as.numeric(fit$gap %||% (fit$Y_treat - fit$Y_synth)),
      mspe_pre_treated   = mspe_pre_tr,
      mspe_pre_placebo   = mspe_pre_co,
      mspe_post_treated  = mspe_post_tr,
      mspe_post_placebo  = mspe_post_co,
      mspe_prune         = mspe_prune,
      n_placebo_pruned   = n_pruned,
      statistic          = statistic,
      times              = fit$times,
      T_pre              = T_pre,
      alternative        = alternative
    ),
    class = "scm_placebo"
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
  if (!inherits(fit, "coresynth_scm")) {
    stop("placebo_in_time() requires a coresynth object with method = 'scm'.",
         call. = FALSE)
  }
  if (inherits(fit, "coresynth_staggered")) {
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
    Y_co_pre[fit_rows, , drop = FALSE], Y_tr_pre[fit_rows],
    cheap_face = TRUE
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
  if (!inherits(fit, "coresynth_scm")) {
    stop("loo_donors() requires a coresynth object with method = 'scm'.",
         call. = FALSE)
  }
  if (inherits(fit, "coresynth_staggered")) {
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
  if (!inherits(fit, "coresynth_scm")) {
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
