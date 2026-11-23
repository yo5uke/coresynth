# ── Pre-treatment Convex-Hull Overlap Check (Liu & Xu 2026, S.6.2 & SM A2.5) ─
#
# Model-free pre-treatment diagnostic comparing how closely each treated unit's
# demeaned pre-treatment trajectory lies inside the convex hull of the control
# units' demeaned paths, relative to the leave-one-out fit errors of the
# controls themselves.

# Extract a list of per-treated-unit pre-treatment inspection tasks from a
# sharp, staggered, or multi-arm coresynth fit.
.overlap_extract_tasks <- function(fit) {
  tasks <- list()

  if (!is.null(fit[["cohort_fits"]])) {
    Y_all_parent <- fit[["Y_all"]] %||% fit[["Y_mat"]]
    Y_co_parent  <- donor_outcomes(fit)
    is_multi     <- inherits(fit, "coresynth_multiarm")

    for (cf in fit[["cohort_fits"]]) {
      T_pre_g <- cf[["T_pre"]]
      Y_tr_g <- if (!is.null(cf[["Y_treat_mat"]])) {
        as.matrix(cf[["Y_treat_mat"]])
      } else if (!is.null(cf[["Y_treat"]])) {
        yt <- cf[["Y_treat"]]
        if (is.matrix(yt)) yt
        else if (!is.null(Y_all_parent) && !is.null(cf[["idx_tr"]]))
          Y_all_parent[, cf[["idx_tr"]], drop = FALSE]
        else matrix(yt, ncol = 1L)
      } else if (!is.null(Y_all_parent) && !is.null(cf[["idx_tr"]])) {
        Y_all_parent[, cf[["idx_tr"]], drop = FALSE]
      } else {
        NULL
      }

      Y_co_g <- if (!is.null(cf[["Y_co_pre"]])) {
        cf[["Y_co_pre"]]
      } else if (!is.null(Y_all_parent) && !is.null(cf[["idx_co"]])) {
        Y_all_parent[, cf[["idx_co"]], drop = FALSE]
      } else if (!is.null(Y_co_parent)) {
        Y_co_parent
      } else {
        NULL
      }

      if (is.null(Y_tr_g) || is.null(Y_co_g)) next
      u_names <- colnames(Y_tr_g) %||% cf[["unit"]] %||%
        paste0("cohort_", cf[["cohort"]],
               if (is_multi && !is.null(cf[["arm"]]))
                 paste0("_arm_", cf[["arm"]]) else "",
               "_unit_", seq_len(ncol(Y_tr_g)))
      if (length(u_names) != ncol(Y_tr_g))
        u_names <- paste0(as.character(u_names), "_", seq_len(ncol(Y_tr_g)))

      for (k in seq_len(ncol(Y_tr_g))) {
        tasks[[length(tasks) + 1L]] <- list(
          unit  = as.character(u_names[k]),
          y_tr  = as.numeric(Y_tr_g[, k]),
          Y_co  = Y_co_g,
          T_pre = as.integer(T_pre_g)
        )
      }
    }
  } else if (inherits(fit, "coresynth_multiarm") &&
             !is.null(fit[["arm_fits"]])) {
    Y_co <- fit[["Y_pre_co"]] %||% donor_outcomes(fit)
    for (af in fit[["arm_fits"]]) {
      Y_tr_a  <- as.matrix(af[["Y_treat"]])
      T_pre_a <- as.integer(af[["T_pre"]] %||% fit[["T_pre"]])
      u_names <- colnames(Y_tr_a) %||%
        paste0("arm_", af[["arm"]], "_unit_", seq_len(ncol(Y_tr_a)))
      for (k in seq_len(ncol(Y_tr_a))) {
        tasks[[length(tasks) + 1L]] <- list(
          unit  = as.character(u_names[k]),
          y_tr  = as.numeric(Y_tr_a[, k]),
          Y_co  = Y_co,
          T_pre = T_pre_a
        )
      }
    }
  } else {
    T_pre <- as.integer(fit[["T_pre"]])
    Y_tr  <- fit[["Y_treat_mat"]] %||% fit[["Y_tr_pre"]] %||% fit[["Y_treat"]]
    Y_co  <- donor_outcomes(fit)
    if (!is.null(Y_tr) && !is.null(Y_co)) {
      if (!is.matrix(Y_tr)) {
        nm   <- names(Y_tr)[1L] %||% "treated"
        Y_tr <- matrix(as.numeric(Y_tr), ncol = 1L,
                       dimnames = list(NULL, nm))
      }
      u_names <- colnames(Y_tr) %||% paste0("treated_", seq_len(ncol(Y_tr)))
      for (k in seq_len(ncol(Y_tr))) {
        tasks[[length(tasks) + 1L]] <- list(
          unit  = as.character(u_names[k]),
          y_tr  = as.numeric(Y_tr[, k]),
          Y_co  = Y_co,
          T_pre = T_pre
        )
      }
    }
  }

  tasks
}

#' Pre-Treatment Convex-Hull Overlap Check
#'
#' Implements the model-free pre-treatment overlap diagnostic of Liu & Xu
#' (2026, Section 6.2 and Online Supplementary Materials A2.5). For each
#' treated unit, the function takes its pre-treatment outcome path and the
#' control units observed throughout that window, subtracts each unit's
#' pre-treatment mean to remove level differences, and computes the root mean
#' squared gap when projecting the treated path onto the simplex convex hull of
#' the control paths. Each control is projected onto the remaining controls in
#' the same way, and the treated unit is flagged when its fit error exceeds the
#' `threshold` quantile (default 90th percentile) of the controls' fit errors.
#'
#' Units with fewer than 4 pre-treatment periods or fewer than 5 completely
#' observed control units, or panels where controls predict one another almost
#' perfectly (median control RMSE below \eqn{10^{-8}} times the control
#' outcome standard deviation), are reported as `"not_assessed"` rather than
#' passing.
#'
#' Like any pre-treatment check, this diagnostic cannot detect differences in
#' levels (removed by demeaning) or latent factors that are constant before
#' treatment, and carries no formal size guarantee.
#'
#' @param fit A `coresynth` object from [scm_fit()] (sharp, staggered, or
#'   multi-arm, across any of the six methods).
#' @param threshold Quantile of the controls' leave-one-out simplex fit errors
#'   above which a treated unit is flagged (default `0.9`, the 90th
#'   percentile).
#' @param max_controls Maximum number of control units (default `200L`). If
#'   more than `max_controls` completely observed controls exist,
#'   `max_controls` are sampled uniformly without replacement, and that subset
#'   serves both as the hull the treated unit is projected onto and as the
#'   control reference distribution.
#' @param seed Optional RNG seed used when subsampling controls above
#'   `max_controls` (default `NULL`; caller's `.Random.seed` is preserved).
#'
#' @return An object of class `"coresynth_overlap"` containing:
#'   * `units`: a `data.frame` with one row per treated unit and columns `unit`,
#'     `T_pre`, `n_controls`, `rmse`, `control_median`, `cutoff`,
#'     `percentile`, `flagged`, and `status` (`"ok"`, `"flagged"`, or
#'     `"not_assessed"`).
#'   * `any_flagged`: logical scalar; `TRUE` if at least one assessed unit is
#'     flagged.
#'   * `threshold`: the quantile threshold used.
#'   * `control_rmses`: named list (by treated unit) of numeric vectors of
#'     control leave-one-out RMSEs.
#'   * `method`: the estimation method from `fit$method`.
#'
#' @references Liu, Z., & Xu, Y. (2026). Unit-Level Conformal Inference for
#'   Panel Counterfactuals. *arXiv:2611.19209*.
#'
#' @seealso [unit_conformal()] for unit-level conformal confidence sets and
#'   p-values.
#' @export
overlap_check <- function(fit, threshold = 0.9, max_controls = 200L,
                          seed = NULL) {
  if (!inherits(fit, "coresynth"))
    stop("overlap_check() requires a coresynth object from scm_fit().",
         call. = FALSE)
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0 || threshold >= 1)
    stop("`threshold` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!is.numeric(max_controls) || length(max_controls) != 1L ||
      !is.finite(max_controls) || max_controls < 5 ||
      max_controls != trunc(max_controls))
    stop("`max_controls` must be a single integer >= 5.", call. = FALSE)
  max_controls <- as.integer(max_controls)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed)))
    stop("`seed` must be NULL or a single finite number.", call. = FALSE)

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed)
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  tasks <- .overlap_extract_tasks(fit)
  if (length(tasks) == 0L)
    stop("overlap_check(): fit does not contain the treated or donor outcome ",
         "matrices required. Re-run scm_fit() with the current version.",
         call. = FALSE)

  # Deduplicate treated units if any task listed the same unit twice
  seen  <- character(0L)
  keep  <- logical(length(tasks))
  for (i in seq_along(tasks)) {
    if (!tasks[[i]]$unit %in% seen) {
      seen <- c(seen, tasks[[i]]$unit)
      keep[i] <- TRUE
    }
  }
  tasks <- tasks[keep]

  co_cache      <- new.env(parent = emptyenv())
  sub_cache     <- new.env(parent = emptyenv())
  rows          <- vector("list", length(tasks))
  control_rmses <- vector("list", length(tasks))
  names(control_rmses) <- vapply(tasks, `[[`, character(1L), "unit")

  for (i in seq_along(tasks)) {
    tk      <- tasks[[i]]
    u_name  <- tk$unit
    T_pre_i <- tk$T_pre
    y_full  <- tk$y_tr
    Y_co    <- tk$Y_co

    not_assessed <- function(n_co = 0L, rmse = NA_real_, med = NA_real_) {
      data.frame(
        unit           = u_name,
        T_pre          = T_pre_i,
        n_controls     = as.integer(n_co),
        rmse           = rmse,
        control_median = med,
        cutoff         = NA_real_,
        percentile     = NA_real_,
        flagged        = NA,
        status         = "not_assessed",
        stringsAsFactors = FALSE
      )
    }

    if (is.na(T_pre_i) || T_pre_i < 4L || length(y_full) < T_pre_i ||
        is.null(Y_co) || nrow(Y_co) < T_pre_i) {
      rows[[i]] <- not_assessed(if (is.null(Y_co)) 0L else ncol(Y_co))
      next
    }

    pre_idx  <- seq_len(T_pre_i)
    y_pre    <- y_full[pre_idx]
    if (any(!is.finite(y_pre))) {
      rows[[i]] <- not_assessed(ncol(Y_co))
      next
    }

    Y_co_pre <- Y_co[pre_idx, , drop = FALSE]
    ok_co    <- apply(Y_co_pre, 2L, function(col) all(is.finite(col)))
    Y_co_pre <- Y_co_pre[, ok_co, drop = FALSE]
    n_co     <- ncol(Y_co_pre)

    if (n_co < 5L) {
      rows[[i]] <- not_assessed(n_co)
      next
    }

    pool_names <- colnames(Y_co_pre) %||% as.character(seq_len(n_co))
    pool_key   <- paste0(T_pre_i, "::", paste(pool_names, collapse = ","))

    if (n_co > max_controls) {
      if (exists(pool_key, envir = sub_cache, inherits = FALSE)) {
        pick <- get(pool_key, envir = sub_cache, inherits = FALSE)
      } else {
        pick <- sort(sample.int(n_co, max_controls))
        assign(pool_key, pick, envir = sub_cache)
      }
      Y_co_pre   <- Y_co_pre[, pick, drop = FALSE]
      n_co       <- ncol(Y_co_pre)
      pool_names <- colnames(Y_co_pre) %||% as.character(seq_len(n_co))
    }

    V_unif     <- rep(1.0 / T_pre_i, T_pre_i)
    y_tilde_i  <- y_pre - mean(y_pre)
    Y_tilde_co <- sweep(Y_co_pre, 2L, colMeans(Y_co_pre), "-")

    w_i <- drop(scm_inner_weights_cpp(Y_tilde_co, y_tilde_i, V_unif))
    e_i <- sqrt(mean((y_tilde_i - drop(Y_tilde_co %*% w_i))^2))

    if (exists(pool_key, envir = co_cache, inherits = FALSE)) {
      e_co <- get(pool_key, envir = co_cache, inherits = FALSE)
    } else {
      e_co <- vapply(seq_len(n_co), function(j) {
        w_j <- drop(scm_inner_weights_cpp(
          Y_tilde_co[, -j, drop = FALSE], Y_tilde_co[, j], V_unif))
        sqrt(mean((Y_tilde_co[, j] -
                     drop(Y_tilde_co[, -j, drop = FALSE] %*% w_j))^2))
      }, numeric(1L))
      names(e_co) <- pool_names
      assign(pool_key, e_co, envir = co_cache)
    }
    control_rmses[[i]] <- e_co

    y_sd   <- stats::sd(as.vector(Y_co_pre))
    med_co <- stats::median(e_co)
    # Liu & Xu (2026, SM A2.5 fn. 7): report "not_assessed" when controls
    # predict one another almost perfectly.
    if (!is.finite(y_sd) || y_sd <= 0 || !is.finite(med_co) ||
        med_co < 1e-8 * y_sd) {
      rows[[i]] <- not_assessed(n_co, rmse = e_i, med = med_co)
      next
    }

    cutoff <- unname(stats::quantile(e_co, probs = threshold))
    pct    <- mean(e_co <= e_i)
    flg    <- isTRUE(e_i > cutoff)

    rows[[i]] <- data.frame(
      unit           = u_name,
      T_pre          = T_pre_i,
      n_controls     = as.integer(n_co),
      rmse           = e_i,
      control_median = med_co,
      cutoff         = cutoff,
      percentile     = pct,
      flagged        = flg,
      status         = if (flg) "flagged" else "ok",
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  structure(list(
    units         = df,
    any_flagged   = any(df$flagged, na.rm = TRUE),
    threshold     = threshold,
    control_rmses = control_rmses,
    method        = fit$method
  ), class = "coresynth_overlap")
}

#' @export
print.coresynth_overlap <- function(x, digits = 4, ...) {
  cat(sprintf("Pre-treatment Overlap Check (Liu & Xu 2026, base = %s, threshold = %.0f%%)\n",
              x$method, 100 * x$threshold))
  df <- x$units
  df_show <- df
  num_cols <- c("rmse", "control_median", "cutoff", "percentile")
  for (col in num_cols) df_show[[col]] <- round(df_show[[col]], digits)
  print(df_show, row.names = FALSE)
  if (isTRUE(x$any_flagged)) {
    n_fl <- sum(df$flagged, na.rm = TRUE)
    cat(sprintf("Flagged: %d unit(s) exceed the control %.0fth percentile.\n",
                n_fl, 100 * x$threshold))
  } else if (all(df$status == "not_assessed")) {
    cat("No units could be assessed (insufficient pre-periods/controls or degenerate pool).\n")
  } else {
    cat("Flagged: none.\n")
  }
  invisible(x)
}

#' @rdname overlap_check
#' @param x A `coresynth_overlap` object returned by [overlap_check()].
#' @param ... Unused.
#' @importFrom broom tidy
#' @export
tidy.coresynth_overlap <- function(x, ...) {
  x$units
}
