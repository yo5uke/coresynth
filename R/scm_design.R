# Helper: build k x J predictor matrix for all J units using pred() specs
# Returns list(X_mat = k x J, pred_names = character(k))
build_design_X <- function(data, unit_var, time_var, units, predictors) {
  J    <- length(units)
  rows <- list()
  nms  <- character(0L)

  for (ps in predictors) {
    if (!inherits(ps, "pred_spec")) {
      stop("Each element of 'predictors' must be a pred() object.", call. = FALSE)
    }
    for (var in ps$vars) {
      if (!var %in% names(data)) {
        stop(paste0("Variable '", var, "' not found in data."), call. = FALSE)
      }
      row_vals <- vapply(units, function(u) {
        vals <- data[[var]][data[[unit_var]] == u & data[[time_var]] %in% ps$times]
        switch(ps$op,
          mean   = mean(vals, na.rm = TRUE),
          median = median(vals, na.rm = TRUE),
          sum    = sum(vals, na.rm = TRUE),
          mean(vals, na.rm = TRUE)
        )
      }, numeric(1L))
      rows[[length(rows) + 1L]] <- row_vals
      t_lbl <- if (length(ps$times) == 1L) as.character(ps$times[[1L]])
               else paste0(ps$times[[1L]], ":", ps$times[[length(ps$times)]])
      nms <- c(nms, paste0(var, "[", t_lbl, "]"))
    }
  }

  X_mat <- do.call(rbind, rows)
  colnames(X_mat) <- as.character(units)
  rownames(X_mat) <- nms
  list(X_mat = X_mat, pred_names = nms)
}

#' Experimental Synthetic Control Design
#'
#' Selects which units to assign to the treatment arm (and which to the control
#' arm) in a planned experiment, following Abadie and Zhao (2026).  Both sets
#' of units are chosen by minimising the distance between their weighted-average
#' predictor vectors and the population-average predictor vector \eqn{\bar{X}},
#' so the resulting estimates are less susceptible to post-randomisation bias
#' than pure random assignment.
#'
#' Three design formulations are available:
#'
#' * **`"base"`** (eq. 7): both the synthetic treated and the synthetic control
#'   independently target the population average \eqn{\bar{X}}.
#' * **`"weakly_targeted"`** (eq. 9): the synthetic treated targets
#'   \eqn{\bar{X}}; the synthetic control targets the synthetic treated predictor
#'   vector (controlled by `beta`).
#' * **`"unit_level"`** (eq. 10): each treated unit gets its own synthetic
#'   control; the aggregate control weight is a convex combination (controlled
#'   by `xi`).
#'
#' Inference uses "blank periods" — pre-experimental periods whose outcomes were
#' *not* used to estimate the weights.  Set `T_fit` strictly smaller than the
#' number of pre-experimental periods to enable the permutation test and split-
#' conformal confidence intervals from Section 3 of Abadie and Zhao (2026).
#'
#' @param data        Long-format data frame (one row per unit–time).
#' @param outcome     Name of the outcome column.
#' @param unit        Name of the unit identifier column.
#' @param time        Name of the time identifier column.
#' @param T0          Last pre-experimental period (a value present in the time
#'   column).  Periods after `T0` are the experimental periods.
#' @param T_fit       Number of fitting periods, counted from the **start** of
#'   the pre-experimental phase.  Defaults to `NULL`, which uses all
#'   pre-experimental periods for fitting (no blank periods; inference
#'   disabled).  When `T_fit` is smaller than the total number of
#'   pre-experimental periods, the remaining periods become blank periods used
#'   for inference.
#' @param m_min       Minimum number of units assigned to treatment (default 1).
#' @param m_max       Maximum number of units assigned to treatment (default 1).
#' @param f           Named numeric vector of population weights \eqn{f_j}.
#'   Defaults to uniform weights \eqn{1/J}.  Will be normalised to sum to 1.
#' @param predictors  A `list()` of [pred()] specifications that define the
#'   predictor matrix \eqn{X_j}.  Defaults to `NULL`, which uses all
#'   fitting-period outcome values as predictors.
#' @param design      Design formulation: `"base"` (default), `"weakly_targeted"`,
#'   or `"unit_level"`.
#' @param beta        Trade-off parameter \eqn{\beta > 0} for the Weakly
#'   targeted design (default 1).
#' @param xi          Trade-off parameter \eqn{\xi > 0} for the Unit-level
#'   design (default 1).
#' @param alpha       Significance level for confidence intervals (default 0.05).
#' @param normalize   If `TRUE` (default), each row of the predictor matrix is
#'   divided by its cross-unit standard deviation before optimisation, so
#'   predictors measured on different scales contribute equally.
#' @param max_subsets Maximum number of treatment-set candidates to evaluate
#'   before switching to random sampling (default 100 000).
#'
#' @return An object of class `"scm_design"` with components:
#'   * `treated_units`: unit identifiers selected for treatment
#'   * `control_units`: unit identifiers in the control pool
#'   * `w`: J-length weight vector for the synthetic treated unit (sums to 1)
#'   * `v`: J-length weight vector for the synthetic control unit (sums to 1)
#'   * `tau_hat`: estimated treatment effects for each experimental period
#'   * `p_value`: permutation p-value (NA when blank periods are unavailable)
#'   * `ci_lower`, `ci_upper`: per-period split-conformal confidence interval
#'   * `Y_synth_tr`, `Y_synth_co`: synthetic treated/control series (all periods)
#'   * `estimate`: ATT (mean of `tau_hat`)
#'
#' @references Abadie, A. and Zhao, J. (2026). "Synthetic Controls for
#'   Experimental Design." MIT Working Paper.
#'
#' @export
scm_design <- function(
  data,
  outcome,
  unit,
  time,
  T0,
  T_fit       = NULL,
  m_min       = 1L,
  m_max       = 1L,
  f           = NULL,
  predictors  = NULL,
  design      = c("base", "weakly_targeted", "unit_level"),
  beta        = 1,
  xi          = 1,
  alpha       = 0.05,
  normalize   = TRUE,
  max_subsets = 100000L
) {
  design <- match.arg(design)
  m_min  <- as.integer(m_min)
  m_max  <- as.integer(m_max)

  for (col in c(outcome, unit, time)) {
    if (!col %in% names(data)) {
      stop(paste0("Column '", col, "' not found in data."), call. = FALSE)
    }
  }

  units <- sort(unique(data[[unit]]))
  times <- sort(unique(data[[time]]))
  J     <- length(units)
  TT    <- length(times)

  if (!T0 %in% times) {
    stop("'T0' must be a value present in the time column.", call. = FALSE)
  }
  T0_idx <- which(times == T0)

  if (m_min < 1L || m_max < m_min || m_max >= J) {
    stop("Need 1 <= m_min <= m_max < J (must leave at least one control unit).",
         call. = FALSE)
  }

  # Fitting / blank / experimental row indices
  TE <- if (is.null(T_fit)) T0_idx else as.integer(T_fit)
  if (TE < 1L || TE > T0_idx) {
    stop("'T_fit' must be between 1 and the number of pre-experimental periods.",
         call. = FALSE)
  }
  E_idx    <- seq_len(TE)
  B_idx    <- setdiff(seq_len(T0_idx), E_idx)
  post_idx <- seq.int(T0_idx + 1L, TT)
  if (length(post_idx) == 0L) {
    stop("No experimental periods found after T0.", call. = FALSE)
  }

  # Outcome matrix (TT x J) — vectorised index fill
  Y <- matrix(NA_real_, nrow = TT, ncol = J,
              dimnames = list(as.character(times), as.character(units)))
  ri_all <- match(data[[time]], times)
  ci_all <- match(data[[unit]], units)
  Y[cbind(ri_all, ci_all)] <- data[[outcome]]

  # Population weights f
  f_vec <- if (is.null(f)) {
    setNames(rep(1 / J, J), as.character(units))
  } else {
    fv <- if (!is.null(names(f))) {
      f[as.character(units)]
    } else {
      if (length(f) != J) stop("Unnamed 'f' must have length J.", call. = FALSE)
      setNames(f, as.character(units))
    }
    if (any(is.na(fv))) stop("Some units in 'data' have no entry in 'f'.", call. = FALSE)
    fv / sum(fv)
  }

  # Predictor matrix X_mat (k x J)
  if (!is.null(predictors)) {
    if (!is.list(predictors)) {
      stop("'predictors' must be a list of pred() objects.", call. = FALSE)
    }
    xi_info    <- build_design_X(data, unit, time, units, predictors)
    X_mat      <- xi_info$X_mat
    pred_names <- xi_info$pred_names
  } else {
    X_mat      <- Y[E_idx, , drop = FALSE]   # TE x J
    pred_names <- paste0("Y[", times[E_idx], "]")
  }

  # Row-wise normalisation (divide by cross-unit SD)
  if (normalize) {
    row_sd <- apply(X_mat, 1L, sd, na.rm = TRUE)
    row_sd[row_sd < 1e-10] <- 1
    X_mat <- X_mat / row_sd
  }

  # Population predictor average (k-vector) and uniform V weights
  X_bar  <- drop(X_mat %*% f_vec)
  k      <- nrow(X_mat)
  V_unif <- rep(1, k)

  # ---- Enumerate treatment sets and optimise ----
  best_obj         <- Inf
  best_treated_idx <- NULL
  best_w_full      <- NULL
  best_v_full      <- NULL

  for (m in m_min:m_max) {
    n_total <- choose(J, m)
    if (n_total <= max_subsets) {
      subsets <- combn(J, m, simplify = FALSE)
    } else {
      message(sprintf(
        "C(%d,%d) = %g subsets > max_subsets = %d; using random sampling.",
        J, m, n_total, max_subsets
      ))
      subsets <- replicate(max_subsets, sort(sample.int(J, m)), simplify = FALSE)
    }

    for (S in subsets) {
      S_co <- setdiff(seq_len(J), S)
      if (length(S_co) == 0L) next

      X_S  <- X_mat[, S,    drop = FALSE]   # k x m
      X_co <- X_mat[, S_co, drop = FALSE]   # k x (J-m)

      # Treated weights w: min ||X_S w - X_bar||^2
      w_sub <- tryCatch(
        scm_inner_weights_cpp(X_S, X_bar, V_unif),
        error = function(e) NULL
      )
      if (is.null(w_sub)) next
      X_Sw  <- drop(X_S %*% w_sub)
      obj_w <- sum((X_bar - X_Sw)^2)

      # Control weights v and total objective
      v_sub <- NULL
      obj   <- NA_real_

      if (design == "base") {
        v_sub <- tryCatch(
          scm_inner_weights_cpp(X_co, X_bar, V_unif),
          error = function(e) NULL
        )
        if (!is.null(v_sub)) {
          obj <- obj_w + sum((X_bar - drop(X_co %*% v_sub))^2)
        }

      } else if (design == "weakly_targeted") {
        v_sub <- tryCatch(
          scm_inner_weights_cpp(X_co, X_Sw, V_unif),
          error = function(e) NULL
        )
        if (!is.null(v_sub)) {
          obj <- obj_w + beta * sum((X_Sw - drop(X_co %*% v_sub))^2)
        }

      } else {  # unit_level
        v_mat    <- matrix(0, length(S_co), length(S))
        obj_unit <- 0
        ok       <- TRUE
        for (jj in seq_along(S)) {
          X_j <- X_mat[, S[jj]]
          v_j <- tryCatch(
            scm_inner_weights_cpp(X_co, X_j, V_unif),
            error = function(e) NULL
          )
          if (is.null(v_j)) { ok <- FALSE; break }
          v_mat[, jj] <- v_j
          obj_unit <- obj_unit + w_sub[jj] * sum((X_j - drop(X_co %*% v_j))^2)
        }
        if (ok) {
          v_sub <- drop(v_mat %*% w_sub)   # aggregate per eq. (11)
          obj   <- obj_w + xi * obj_unit
        }
      }

      if (is.null(v_sub) || is.na(obj)) next

      if (obj < best_obj) {
        best_obj         <- obj
        best_treated_idx <- S
        w_full           <- numeric(J); w_full[S]    <- w_sub
        v_full           <- numeric(J); v_full[S_co] <- v_sub
        best_w_full      <- setNames(w_full, as.character(units))
        best_v_full      <- setNames(v_full, as.character(units))
      }
    }
  }

  if (is.null(best_treated_idx)) {
    stop("Optimisation failed: no feasible treatment set found.", call. = FALSE)
  }

  # Synthetic series for all periods
  Y_synth_tr <- drop(Y %*% best_w_full)
  Y_synth_co <- drop(Y %*% best_v_full)

  tau_hat <- (Y_synth_tr - Y_synth_co)[post_idx]
  u_blank  <- (Y_synth_tr - Y_synth_co)[B_idx]

  T_post_n  <- length(post_idx)
  T_blank_n <- length(B_idx)

  # ---- Inference (Section 3 of Abadie & Zhao 2026) ----
  S_obs <- mean(abs(tau_hat))

  if (T_blank_n == 0L) {
    # Warn only when T_fit was explicitly provided but leaves no blank periods;
    # when T_fit = NULL (default) the user is not attempting inference.
    if (!is.null(T_fit)) {
      warning(
        "No blank periods: p_value and CI are NA. ",
        "Set T_fit < number of pre-experimental periods to enable inference.",
        call. = FALSE
      )
    }
    p_value  <- NA_real_
    ci_lower <- rep(NA_real_, T_post_n)
    ci_upper <- rep(NA_real_, T_post_n)
  } else {
    # Permutation test: eq. (16)-(17)
    all_vals  <- c(u_blank, tau_hat)
    n_periods <- T_blank_n + T_post_n
    n_combs   <- choose(n_periods, T_post_n)

    if (n_combs <= 10000L) {
      combs   <- combn(n_periods, T_post_n, simplify = FALSE)
      S_perms <- vapply(
        combs, function(idx) mean(abs(all_vals[idx])), numeric(1L)
      )
    } else {
      set.seed(42L)
      S_perms <- replicate(
        10000L,
        mean(abs(all_vals[sample.int(n_periods, T_post_n)]))
      )
    }
    p_value <- mean(S_perms >= S_obs)

    # Split-conformal CI: eq. (18)-(19)
    q_alpha  <- as.numeric(quantile(abs(u_blank), 1 - alpha))
    ci_lower <- tau_hat - q_alpha
    ci_upper <- tau_hat + q_alpha
  }

  structure(
    list(
      method         = "scm_design",
      design         = design,
      treated_units  = units[best_treated_idx],
      control_units  = units[setdiff(seq_len(J), best_treated_idx)],
      w              = best_w_full,
      v              = best_v_full,
      objective      = best_obj,
      tau_hat        = tau_hat,
      u_blank        = u_blank,
      p_value        = p_value,
      ci_lower       = ci_lower,
      ci_upper       = ci_upper,
      Y_synth_tr     = Y_synth_tr,
      Y_synth_co     = Y_synth_co,
      times_post     = times[post_idx],
      times_blank    = times[B_idx],
      times_fit      = times[E_idx],
      times          = times,
      T0             = T0,
      T_fit          = TE,
      J              = J,
      m_min          = m_min,
      m_max          = m_max,
      f              = f_vec,
      units          = units,
      alpha          = alpha,
      pred_names     = pred_names,
      estimate       = mean(tau_hat),
      T_pre          = T0_idx,
      Y              = Y
    ),
    class = "scm_design"
  )
}

#' @export
print.scm_design <- function(x, ...) {
  cat("=== scm_design (Abadie & Zhao 2026) ===\n")
  cat("Design variant :", x$design, "\n")
  cat("Treated units  :", paste(x$treated_units, collapse = ", "), "\n")
  cat("ATT estimate   :", round(x$estimate, 4L), "\n")
  if (!is.na(x$p_value)) {
    cat(sprintf("p-value        : %.4f (alpha = %.2f)\n", x$p_value, x$alpha))
  } else {
    cat("p-value        : NA (no blank periods)\n")
  }
  invisible(x)
}

#' @export
summary.scm_design <- function(object, ...) {
  cat("=== scm_design summary (Abadie & Zhao 2026) ===\n")
  cat("Design variant      :", object$design, "\n")
  cat("Number of units : J =", object$J, "\n")
  cat("Fitting periods : TE =", object$T_fit, "\n")
  cat("Blank periods   : TB =", length(object$times_blank), "\n")
  cat("Exp. periods    : T_post =", length(object$times_post), "\n\n")
  cat("Treated units  :", paste(object$treated_units, collapse = ", "), "\n")
  w_nz <- object$w[object$w > 1e-4]
  if (length(w_nz) > 0L) {
    cat("Treated weights (non-zero):\n")
    print(round(w_nz, 4L))
  }
  cat("\nControl weights (non-zero):\n")
  v_nz <- object$v[object$v > 1e-4]
  if (length(v_nz) > 0L) print(round(v_nz, 4L))
  cat("\nATT estimate :", round(object$estimate, 6L), "\n")
  if (!is.na(object$p_value)) {
    cat("p-value      :", round(object$p_value, 4L), "\n\n")
    ci_df <- data.frame(
      time     = object$times_post,
      tau_hat  = round(object$tau_hat, 4L),
      ci_lower = round(object$ci_lower, 4L),
      ci_upper = round(object$ci_upper, 4L)
    )
    cat(sprintf(
      "Confidence intervals (1 - alpha = %.0f%%):\n",
      100 * (1 - object$alpha)
    ))
    print(ci_df, row.names = FALSE)
  }
  invisible(object)
}

#' Plot an scm_design object
#'
#' @param x    An `scm_design` object.
#' @param type `"outcome"` (default): synthetic treated vs synthetic control
#'   outcome series over all periods. `"gap"`: estimated treatment effect in the
#'   experimental periods, with split-conformal confidence intervals.
#' @param ...  Currently ignored.
#' @return A `ggplot` object: for `type = "outcome"`, the synthetic treated and
#'   synthetic control outcome series; for `type = "gap"`, the estimated
#'   treatment effect over the experimental periods with split-conformal
#'   confidence intervals. The object is returned for printing or further
#'   customisation.
#' @export
plot.scm_design <- function(x, type = c("outcome", "gap"), ...) {
  type <- match.arg(type)

  if (type == "outcome") {
    df <- data.frame(
      time    = x$times,
      Treated = x$Y_synth_tr,
      Control = x$Y_synth_co
    )
    df_long <- data.frame(
      time  = rep(x$times, 2L),
      value = c(x$Y_synth_tr, x$Y_synth_co),
      group = rep(
        c("Synthetic Treated", "Synthetic Control"),
        each = length(x$times)
      )
    )
    ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = value,
                                          linetype = group)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = x$T0, linetype = "dotted",
                          colour = "grey50") +
      ggplot2::scale_linetype_manual(
        values = c("Synthetic Treated" = "solid",
                   "Synthetic Control" = "dashed")
      ) +
      ggplot2::labs(
        title    = "Synthetic Control Design: Outcome Series",
        x        = "Time", y = "Outcome", linetype = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
  } else {
    df <- data.frame(
      time     = x$times_post,
      tau_hat  = x$tau_hat,
      ci_lower = x$ci_lower,
      ci_upper = x$ci_upper
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = tau_hat)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")
    if (!any(is.na(df$ci_lower))) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        alpha = 0.2, fill = "steelblue"
      )
    }
    p + ggplot2::labs(
      title = "Treatment Effect Estimate",
      x = "Time", y = "Estimated Effect"
    ) + ggplot2::theme_bw()
  }
}
