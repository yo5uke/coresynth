#' @importFrom broom tidy
#' @export
tidy.coresynth <- function(x, ...) {
  if (is.null(x$unit_weights)) {
    return(data.frame())
  }
  data.frame(
    term     = names(x$unit_weights) %||% paste0("Unit_", seq_along(x$unit_weights)),
    estimate = as.numeric(x$unit_weights),
    type     = "unit_weight",
    stringsAsFactors = FALSE
  )
}

#' @export
tidy.coresynth_staggered <- function(x, ...) {
  ce <- x$cohort_estimates
  if (is.null(ce) || nrow(ce) == 0L) {
    return(data.frame())
  }
  data.frame(
    term     = paste0("cohort_", ce$cohort),
    estimate = ce$estimate,
    weight   = ce$weight,
    n_treated = ce$n_treated,
    T_pre    = ce$T_pre,
    T_post   = ce$T_post,
    type     = "cohort_estimate",
    stringsAsFactors = FALSE
  )
}

#' @export
tidy.coresynth_multiarm <- function(x, ...) {
  ca <- x$cohort_arm_estimates
  if (is.null(ca)) {
    # Sharp multi-arm fits carry no cohort-arm table; fall through to the
    # staggered/sharp implementation next in the class chain.
    return(NextMethod())
  }
  data.frame(
    term      = paste0("cohort_", ca$cohort, "_arm_", ca$arm),
    estimate  = ca$estimate,
    weight    = ca$weight,
    n_treated = ca$n_treated,
    T_pre     = ca$T_pre,
    T_post    = ca$T_post,
    type      = "cohort_arm_estimate",
    stringsAsFactors = FALSE
  )
}

#' @importFrom broom glance
#' @export
glance.coresynth <- function(x, ...) {
  T_all   <- length(x$times) %||% NA_integer_
  T_pre   <- x$T_pre %||% NA_integer_
  T_post  <- if (is.finite(T_all) && is.finite(T_pre)) T_all - T_pre else NA_integer_

  n_treated <- if (!is.null(x$Y_treat) && is.matrix(x$Y_treat)) {
    ncol(x$Y_treat)
  } else if (!is.null(x$Y_treat)) {
    1L
  } else {
    NA_integer_
  }

  n_controls <- if (!is.null(x$unit_weights)) {
    length(x$unit_weights)
  } else if (!is.null(x$idx_co)) {
    length(x$idx_co)
  } else if (isTRUE(x$staggered) && !is.null(x$cohort_fits)) {
    round(mean(vapply(x$cohort_fits,
                      function(cf) length(cf$idx_co), numeric(1L))))
  } else {
    NA_integer_
  }

  data.frame(
    method     = x$method,
    estimate   = x$estimate,
    n_controls = n_controls,
    n_treated  = n_treated,
    T_pre      = T_pre,
    T_post     = T_post,
    staggered  = isTRUE(x$staggered),
    multi_arm  = isTRUE(x$multi_arm),
    stringsAsFactors = FALSE
  )
}

#' Tidy an inference result
#'
#' Coerces a `coresynth_inference` or `sdid_inference` object to a one-row
#' tidy `data.frame` with broom-style column names so it can be combined
#' with regression output for paper tables.
#'
#' @param x A `coresynth_inference` (or `sdid_inference`) object returned by
#'   [sdid_inference()], [gsc_inference()], or [si_inference()].
#' @param conf.int Logical. Include `conf.low`/`conf.high` columns when CI is
#'   available (default `TRUE`). Methods without a CI report `NA`.
#' @param ... Unused.
#' @return A one-row `data.frame` with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, `conf.low`, `conf.high`, `method`, `alternative`,
#'   `n_controls`, `staggered`.
#' @importFrom broom tidy
#' @export
tidy.coresynth_inference <- function(x, conf.int = TRUE, ...) {
  se   <- if (is.null(x$se)) NA_real_ else x$se
  stat <- if (is.finite(se) && se > 0) x$estimate / se else NA_real_
  ci_l <- if (is.null(x$ci_lower)) NA_real_ else x$ci_lower
  ci_u <- if (is.null(x$ci_upper)) NA_real_ else x$ci_upper

  out <- data.frame(
    term        = "ATT",
    estimate    = x$estimate,
    std.error   = se,
    statistic   = stat,
    p.value     = x$p_value,
    conf.low    = ci_l,
    conf.high   = ci_u,
    method      = x$method,
    alternative = x$alternative %||% NA_character_,
    n_controls  = x$n_controls,
    staggered   = isTRUE(x$staggered),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(conf.int)) {
    out$conf.low  <- NULL
    out$conf.high <- NULL
  }
  out
}

#' Glance at an inference result
#'
#' One-row summary of a `coresynth_inference` (or `sdid_inference`) object.
#'
#' @param x  An inference object.
#' @param ... Unused.
#' @return A one-row `data.frame` with columns `method`, `n_controls`,
#'   `staggered`, `estimate`, `std.error`, `p.value`, `conf.low`,
#'   `conf.high`, `alternative`, `n_boot_valid`.
#' @importFrom broom glance
#' @export
glance.coresynth_inference <- function(x, ...) {
  n_boot_valid <- if (is.null(x$boot_ests)) NA_integer_ else
    sum(!is.na(x$boot_ests))
  data.frame(
    method       = x$method,
    n_controls   = x$n_controls,
    staggered    = isTRUE(x$staggered),
    estimate     = x$estimate,
    std.error    = x$se %||% NA_real_,
    p.value      = x$p_value,
    conf.low     = x$ci_lower %||% NA_real_,
    conf.high    = x$ci_upper %||% NA_real_,
    alternative  = x$alternative %||% NA_character_,
    n_boot_valid = n_boot_valid,
    stringsAsFactors = FALSE
  )
}

#' @importFrom broom augment
#' @export
augment.coresynth <- function(x, include_donors = FALSE, ...) {
  Y_treat <- treated_outcomes(x)
  if (is.null(Y_treat))
    stop("augment.coresynth: cannot find treated outcome in fit object.",
         call. = FALSE)

  Y_fitted <- synthetic_outcomes(x)
  if (is.null(Y_fitted))
    stop("augment.coresynth: cannot find fitted values in fit object.",
         call. = FALSE)

  T_total   <- length(Y_treat)
  T_pre     <- x$T_pre %||% 0L
  is_post   <- c(rep(FALSE, T_pre), rep(TRUE, T_total - T_pre))
  times     <- x$times %||% seq_len(T_total)
  unit_name <- if (!is.null(names(x$Y_treat))) names(x$Y_treat)[1L] else "treated"

  treated_df <- data.frame(
    .unit    = unit_name,
    .type    = "treated",
    .time    = times,
    .observed= Y_treat,
    .fitted  = Y_fitted,
    .resid   = Y_treat - Y_fitted,
    .treated = is_post,
    .period  = factor(ifelse(is_post, "post", "pre"), levels = c("pre", "post")),
    stringsAsFactors = FALSE
  )

  if (!include_donors) {
    # Drop .unit / .type for backwards compatibility
    treated_df$.unit <- NULL
    treated_df$.type <- NULL
    return(treated_df)
  }

  Y_co_mat <- donor_outcomes(x)

  if (is.null(Y_co_mat)) {
    warning(
      "include_donors = TRUE: control unit outcomes not available in this fit. ",
      "Returning treated unit rows only.",
      call. = FALSE
    )
    return(treated_df)
  }

  donor_names <- colnames(Y_co_mat) %||% paste0("Unit_", seq_len(ncol(Y_co_mat)))
  control_df <- do.call(rbind, lapply(seq_along(donor_names), function(j) {
    data.frame(
      .unit    = donor_names[j],
      .type    = "control",
      .time    = times,
      .observed= as.numeric(Y_co_mat[, j]),
      .fitted  = NA_real_,
      .resid   = NA_real_,
      .treated = FALSE,
      .period  = factor(ifelse(is_post, "post", "pre"), levels = c("pre", "post")),
      stringsAsFactors = FALSE
    )
  }))

  rbind(treated_df, control_df)
}

#' @export
augment.coresynth_staggered <- function(x, include_donors = FALSE, ...) {
  .augment_staggered(x, include_donors = include_donors)
}

# ── Internal: augment() for staggered fits ────────────────────────────────────
# Returns a long-format data.frame stacking the treated outcome and synthetic
# counterfactual for each cohort (and each treatment arm, if multi-arm).
.augment_staggered <- function(x, include_donors = FALSE) {
  cohort_fits <- x$cohort_fits
  if (is.null(cohort_fits) || length(cohort_fits) == 0L) {
    warning(
      "augment.coresynth: staggered fit has no cohort_fits to expand.",
      call. = FALSE
    )
    return(data.frame())
  }
  times_all <- x$times %||% seq_len(nrow(x$Y_treat %||% matrix(NA_real_, 0, 0)))
  is_multi  <- inherits(x, "coresynth_multiarm")

  Y_all_parent <- x$Y_all %||% x$Y_mat  # GSC/SI uses Y_all, SDID uses Y_mat

  rows <- lapply(cohort_fits, function(cf) {
    # Determine TT from cf or parent
    if (!is.null(cf$Y_treat)) {
      TT <- if (is.matrix(cf$Y_treat)) nrow(cf$Y_treat) else length(cf$Y_treat)
    } else if (!is.null(Y_all_parent)) {
      TT <- nrow(Y_all_parent)
    } else {
      return(NULL)
    }
    if (is.null(TT) || TT == 0L) return(NULL)
    T_pre_g  <- cf$T_pre
    is_post  <- c(rep(FALSE, T_pre_g), rep(TRUE, TT - T_pre_g))
    period   <- factor(ifelse(is_post, "post", "pre"), levels = c("pre", "post"))
    times_g  <- if (length(times_all) >= TT) times_all[seq_len(TT)] else seq_len(TT)

    # Resolve treated path
    Y_treat_g <- if (!is.null(cf$Y_treat)) {
      if (is.matrix(cf$Y_treat)) rowMeans(cf$Y_treat) else as.numeric(cf$Y_treat)
    } else if (!is.null(Y_all_parent) && !is.null(cf$idx_tr)) {
      rowMeans(Y_all_parent[, cf$idx_tr, drop = FALSE])
    } else {
      rep(NA_real_, TT)
    }

    # Resolve fitted in priority order (SCM/SI/GSC carry Y_synth or Y_tr_hat;
    # SDID staggered carries unit_weights + omega0 → reconstruct from Y_mat)
    Y_fit_g <- if (!is.null(cf$Y_synth)) {
      if (is.matrix(cf$Y_synth)) rowMeans(cf$Y_synth) else as.numeric(cf$Y_synth)
    } else if (!is.null(cf$Y_tr_hat)) {
      rowMeans(as.matrix(cf$Y_tr_hat))
    } else if (!is.null(cf$unit_weights) && !is.null(Y_all_parent) &&
               !is.null(cf$idx_co)) {
      omega0 <- cf$omega0 %||% 0
      drop(Y_all_parent[, cf$idx_co, drop = FALSE] %*% cf$unit_weights) + omega0
    } else {
      rep(NA_real_, TT)
    }

    n_treat_units <- if (!is.null(cf$Y_treat) && is.matrix(cf$Y_treat)) {
      ncol(cf$Y_treat)
    } else if (!is.null(cf$idx_tr)) {
      length(cf$idx_tr)
    } else {
      cf$n_treated %||% 1L
    }
    unit_label    <- sprintf("cohort_%s", as.character(cf$cohort))
    if (is_multi && !is.null(cf$arm))
      unit_label <- sprintf("%s_arm_%s", unit_label, as.character(cf$arm))

    treated_df <- data.frame(
      .cohort   = cf$cohort,
      .arm      = if (is_multi) cf$arm %||% NA_integer_ else NA_integer_,
      .unit     = unit_label,
      .type     = "treated",
      .time     = times_g,
      .observed = Y_treat_g,
      .fitted   = Y_fit_g,
      .resid    = Y_treat_g - Y_fit_g,
      .treated  = is_post,
      .period   = period,
      n_treated = n_treat_units,
      stringsAsFactors = FALSE
    )

    if (!include_donors) return(treated_df)

    Y_co_mat <- NULL
    if (!is.null(cf$Y_co_pre) && !is.null(cf$Y_co_post)) {
      Y_co_mat <- rbind(cf$Y_co_pre, cf$Y_co_post)
    } else if (!is.null(cf$idx_co) && !is.null(Y_all_parent)) {
      Y_co_mat <- Y_all_parent[, cf$idx_co, drop = FALSE]
    }
    if (is.null(Y_co_mat)) return(treated_df)

    donor_names <- colnames(Y_co_mat) %||% paste0("Unit_", seq_len(ncol(Y_co_mat)))
    control_df <- do.call(rbind, lapply(seq_along(donor_names), function(j) {
      data.frame(
        .cohort   = cf$cohort,
        .arm      = if (is_multi) cf$arm %||% NA_integer_ else NA_integer_,
        .unit     = donor_names[j],
        .type     = "control",
        .time     = times_g,
        .observed = as.numeric(Y_co_mat[, j]),
        .fitted   = NA_real_,
        .resid    = NA_real_,
        .treated  = FALSE,
        .period   = period,
        n_treated = n_treat_units,
        stringsAsFactors = FALSE
      )
    }))
    rbind(treated_df, control_df)
  })

  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  if (is.null(out)) return(data.frame())
  rownames(out) <- NULL
  if (!is_multi) out$.arm <- NULL
  out
}
