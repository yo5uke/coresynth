# -- BFR (2022) S.4.2 balance possibility frontier ---------------------------

# Rebuild the minimal panel a refit needs from a staggered SCM fit, so the
# frontier can be traced without asking for the data again.
.pan_from_staggered_fit <- function(fit) {
  Y <- fit$Y_all
  if (is.null(Y) || is.null(fit$cohort_fits))
    stop("this staggered SCM fit predates the frontier tool; re-run ",
         "scm_fit() with the current version.", call. = FALSE)
  N  <- ncol(Y)
  TT <- nrow(Y)
  idx_tr  <- sort(unique(unlist(lapply(fit$cohort_fits, `[[`, "idx_tr"))))
  T_adopt <- rep(NA_integer_, N)
  for (z in fit$cohort_fits) T_adopt[z$idx_tr] <- as.integer(z$T_pre + 1L)
  D <- matrix(0L, TT, N)
  for (j in idx_tr) D[T_adopt[j]:TT, j] <- 1L
  list(Y = Y, D = D,
       units = colnames(Y) %||% seq_len(N),
       times = fit$times,
       T_pre = min(T_adopt[idx_tr], na.rm = TRUE) - 1L,
       T_adopt = T_adopt, is_sharp = FALSE,
       idx_treat = idx_tr,
       idx_control = setdiff(seq_len(N), idx_tr))
}

#' Balance Possibility Frontier for Staggered SCM
#'
#' Traces how the two pre-treatment imbalances of partially pooled SCM trade
#' off as the pooling parameter `nu` varies: `q_pool`, the imbalance of the
#' average of the treated units, which the aggregate ATT rests on, and
#' `q_sep`, the root mean square of the per-block imbalances, which the
#' block-level effects rest on. Both are the cohort-level measures of
#' Ben-Michael, Feller & Rothstein (2022, Appendix A.2), in which a block
#' enters through the sum over its treated units, so a block of \eqn{n_g}
#' units counts \eqn{n_g} times towards `q_pool` and \eqn{n_g^2} times towards
#' `q_sep`. Supplying `weights` to [scm_fit()] replaces the unit count with
#' the mass it gives them, and the frontier is traced on that same objective.
#'
#' Ben-Michael, Feller & Rothstein (2022) call this plot "an important tool for
#' using partially pooled SCM in practice" (S.4.2): the curve is strongly
#' convex, so a little pooling usually buys a large reduction in one imbalance
#' for a small rise in the other, and where that bargain stops being worth it
#' is a judgement the data can inform but not settle. Their `nu` heuristic,
#' which [scm_fit()] uses by default, is one point on this curve.
#'
#' @param fit A staggered `coresynth` object from [scm_fit()] with
#'   `method = "scm"`.
#' @param nu_grid Values of `nu` to trace.
#' @return A `data.frame` of class `scm_frontier`, one row per `nu`, with
#'   `nu`, `q_sep`, `q_pool`, their values normalised by the separate solution
#'   (`q_sep_rel`, `q_pool_rel`), `estimate`, and `heuristic` marking the row
#'   closest to the fit's own `nu`. Plot it with [plot()].
#' @references Ben-Michael, E., Feller, A., & Rothstein, J. (2022). Synthetic
#'   controls with staggered adoption. *JRSS-B*, 84(2), 351-381.
#' @examples
#' set.seed(1)
#' panel <- expand.grid(unit = 1:15, year = 1:20)
#' adopt <- c(`1` = 12, `2` = 15, `3` = 17)
#' panel$treated <- as.integer(!is.na(adopt[as.character(panel$unit)]) &
#'   panel$year >= adopt[as.character(panel$unit)])
#' panel$treated[is.na(panel$treated)] <- 0L
#' panel$gdp <- panel$unit + 0.4 * panel$year +
#'   rnorm(nrow(panel)) + 2 * panel$treated
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "scm")
#' fr <- scm_balance_frontier(fit, nu_grid = c(0, 0.5, 1))
#' fr
#' @export
scm_balance_frontier <- function(fit, nu_grid = seq(0, 1, length.out = 11L)) {
  if (!inherits(fit, "coresynth_scm"))
    stop("scm_balance_frontier() requires a coresynth object with ",
         "method = 'scm'.", call. = FALSE)
  if (!inherits(fit, "coresynth_staggered"))
    stop("scm_balance_frontier() applies to staggered SCM fits. A single ",
         "adoption date leaves nothing to pool across.", call. = FALSE)
  if (is.null(fit$pooling))
    stop("this fit used nu = 'separate', so it carries no pooling ",
         "diagnostics. Re-fit with the default nu to trace the frontier.",
         call. = FALSE)
  if (!is.numeric(nu_grid) || length(nu_grid) < 2L || anyNA(nu_grid) ||
      any(nu_grid < 0) || any(nu_grid > 1))
    stop("`nu_grid` must be two or more values in [0, 1].", call. = FALSE)

  pan  <- .pan_from_staggered_fit(fit)
  nu_grid <- sort(unique(nu_grid))
  rows <- lapply(nu_grid, function(nu) {
    f <- .fit_scm_staggered_pooled(
      pan, nu = nu,
      # The block masses are part of the objective, so a frontier drawn
      # without them would not be the frontier the fit sits on.
      unit_weights  = fit$treated_weights,
      fixedeff      = isTRUE(fit$pooling$fixedeff),
      control_group = fit$control_group %||% "clean",
      pool_cohorts  = fit$pooling$pool_cohorts %||% TRUE)
    data.frame(nu = nu, q_sep = f$pooling$q_sep, q_pool = f$pooling$q_pool,
               estimate = f$estimate)
  })
  out <- do.call(rbind, rows)
  # Normalise by the separate solution, as the paper does, so the two
  # imbalances sit on comparable scales whatever the outcome's units
  out$q_sep_rel  <- out$q_sep  / fit$pooling$q_sep_separate
  out$q_pool_rel <- out$q_pool / fit$pooling$q_pool_separate
  out$heuristic  <- seq_len(nrow(out)) ==
    which.min(abs(out$nu - fit$pooling$nu))
  attr(out, "nu_fit") <- fit$pooling$nu
  class(out) <- c("scm_frontier", "data.frame")
  out
}

#' @export
print.scm_frontier <- function(x, ...) {
  cat("=== SCM balance possibility frontier ===\n")
  cat("nu used by the fit:", round(attr(x, "nu_fit"), 4), "\n")
  cat("q_pool: imbalance of the treated average (drives the ATT)\n")
  cat("q_sep : RMS of the per-block imbalances (drives block-level effects)\n\n")
  df <- as.data.frame(x)
  df$nu <- round(df$nu, 3)
  for (nm in c("q_sep", "q_pool", "q_sep_rel", "q_pool_rel", "estimate"))
    df[[nm]] <- round(df[[nm]], 4)
  df$heuristic <- ifelse(df$heuristic, "<- nu_hat", "")
  print(df, row.names = FALSE)
  invisible(x)
}

#' Plot a Balance Possibility Frontier
#'
#' Draws the trade-off traced by [scm_balance_frontier()]: the imbalance of the
#' treated average, which the aggregate ATT rests on, against the root mean
#' square of the per-block imbalances, which the block-level effects rest on.
#' Separate SCM sits at one end of the curve and fully pooled SCM at the other,
#' with the heuristic the fit used marked along it.
#'
#' The curve is normally strongly convex, which is the point of looking at it:
#' a small move away from either end usually buys a large reduction in the
#' other imbalance for very little, so the interesting values of `nu` are the
#' ones near the bend.
#'
#' @param x An `scm_frontier` object from [scm_balance_frontier()].
#' @param relative Plot the imbalances normalised by the separate solution
#'   (default), which puts both axes on comparable scales, rather than in the
#'   outcome's own units.
#' @param label_nu Annotate the drawn points with their `nu`.
#' @param colors Optional named overrides for the `"frontier"` line and the
#'   `"heuristic"` marker.
#' @param ... Ignored.
#' @return A `ggplot2` plot object.
#' @seealso [scm_balance_frontier()]
#' @examples
#' set.seed(1)
#' panel <- expand.grid(unit = 1:15, year = 1:20)
#' adopt <- c(`1` = 12, `2` = 15, `3` = 17)
#' panel$treated <- as.integer(!is.na(adopt[as.character(panel$unit)]) &
#'   panel$year >= adopt[as.character(panel$unit)])
#' panel$treated[is.na(panel$treated)] <- 0L
#' panel$gdp <- panel$unit + 0.4 * panel$year +
#'   rnorm(nrow(panel)) + 2 * panel$treated
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "scm")
#' plot(scm_balance_frontier(fit, nu_grid = c(0, 0.25, 0.5, 0.75, 1)))
#' @export
plot.scm_frontier <- function(x, relative = TRUE, label_nu = TRUE,
                              colors = NULL, ...) {
  if (!isTRUE(relative) && !isFALSE(relative))
    stop("`relative` must be TRUE or FALSE.", call. = FALSE)
  key_map <- c(frontier = "Frontier", heuristic = "Heuristic nu")
  cols <- stats::setNames(
    .merge_named_colors(c(frontier = "#2166ac", heuristic = "#b2182b"), colors),
    unname(key_map))

  df <- as.data.frame(x)
  df$xv <- if (relative) df$q_sep_rel  else df$q_sep
  df$yv <- if (relative) df$q_pool_rel else df$q_pool
  hit <- df[df$heuristic, , drop = FALSE]

  # ASCII only in the labels: a non-UTF-8 locale cannot measure the text
  axis_x <- if (relative) "Per-block imbalance (relative to separate SCM)"
            else "Per-block imbalance (q_sep)"
  axis_y <- if (relative) "Pooled imbalance (relative to separate SCM)"
            else "Pooled imbalance (q_pool)"

  p <- ggplot(df, aes(x = xv, y = yv)) +
    geom_path(aes(color = "Frontier"), linewidth = 0.8) +
    geom_point(aes(color = "Frontier"), size = 1.8) +
    {if (nrow(hit))
      geom_point(data = hit, aes(x = xv, y = yv, color = "Heuristic nu"),
                 size = 3.4)} +
    scale_color_manual(values = cols, breaks = names(cols)) +
    theme_minimal(base_size = 13) +
    labs(title = "Balance Possibility Frontier  [SCM]",
         subtitle = paste0("nu runs from ", round(min(df$nu), 2), " to ",
                           round(max(df$nu), 2),
                           "; lower-left is better on both"),
         x = axis_x, y = axis_y, color = NULL) +
    theme(legend.position = "bottom")

  if (isTRUE(label_nu))
    p <- p + geom_text(aes(label = formatC(nu, format = "g", digits = 2)),
                       vjust = -0.9, size = 3, show.legend = FALSE)
  p
}
