# -- Internal style-merging helpers ------------------------------------------
#
# geom_vline()/geom_hline() cannot be removed from a ggplot object once added
# (only overplotted), so suppression and restyling of reference lines must
# happen inside the plot method, before the layer is ever added.

# Merge user overrides onto a default aesthetic list for a reference line.
# override = NULL or FALSE suppresses the line entirely (returns NULL).
# override = list(...) is merged onto default via modifyList (unset keys keep
# their default value).
.line_style <- function(default, override) {
  if (is.null(override) || isFALSE(override)) return(NULL)
  if (!is.list(override))
    stop("vline/hline must be a list of aesthetic overrides ",
         "(e.g. list(color = \"red\")), or NULL/FALSE to hide the line.",
         call. = FALSE)
  utils::modifyList(default, override)
}

# Merge a user-supplied named vector (colors, labels, ...) onto the package
# default, keeping any series name the user didn't mention at its default
# value. Keys are always the canonical series names, so colors and labels can
# be overridden independently.
.merge_named_vec <- function(default, override, what = "colors",
                             example = "c(Treated = \"black\")") {
  if (is.null(override)) return(default)
  if (is.null(names(override)) || !all(nzchar(names(override))))
    stop(what, " must be a named vector, e.g. ", example, ". ",
         "Valid names: ", paste(names(default), collapse = ", "), ".",
         call. = FALSE)
  unknown <- setdiff(names(override), names(default))
  if (length(unknown) > 0)
    stop(what, " has unrecognized name(s): ", paste(unknown, collapse = ", "),
         ". Valid names: ", paste(names(default), collapse = ", "), ".",
         call. = FALSE)
  default[names(override)] <- override
  default
}

.merge_named_colors <- function(default, override) {
  .merge_named_vec(default, override, what = "colors")
}

# Legend labels default to the canonical series names themselves.
.merge_named_labels <- function(series, override) {
  .merge_named_vec(stats::setNames(series, series), override,
                   what = "labels", example = "c(Treated = \"California\")")
}

#' Plot a coresynth model
#'
#' @param x      A `coresynth` object.
#' @param type   One of `"trend"` (observed vs synthetic), `"gap"` (ATT over time),
#'               or `"weights"` (donor unit weight bar chart).
#' @param colors For `type = "trend"`: a named vector overriding series colors,
#'   e.g. `c(Treated = "black")` (valid names: `"Treated"`, `"Synthetic Control"`).
#'   For `type = "gap"`: a single color string for the gap line. Ignored for
#'   `type = "weights"` (use `fill` instead).
#' @param labels For `type = "trend"`: a named vector overriding the legend
#'   text of individual series, e.g. `c(Treated = "California")` (valid names:
#'   `"Treated"`, `"Synthetic Control"`). Series not mentioned keep their
#'   default label, and `colors` keys always refer to the original series
#'   names regardless of relabeling. Ignored for other types (no legend).
#' @param vline  Aesthetic overrides for the vertical treatment-time line, as a
#'   list passed to [ggplot2::geom_vline()] (e.g. `list(color = "red")`).
#'   `NULL` or `FALSE` hides the line entirely. Applies to `"trend"` and `"gap"`.
#' @param hline  Aesthetic overrides for the horizontal zero line in `type =
#'   "gap"`, as a list passed to [ggplot2::geom_hline()]. `NULL` or `FALSE`
#'   hides the line. Ignored for other types.
#' @param fill   For `type = "weights"`: a single color string overriding the
#'   bar fill. Ignored for other types.
#' @param ...    Ignored.
#' @return A `ggplot2` plot object.
#' @examples
#' set.seed(1)
#' panel <- expand.grid(unit = 1:10, year = 1:20)
#' panel$treated <- as.integer(panel$unit == 5 & panel$year > 15)
#' panel$gdp <- panel$unit + 0.5 * panel$year +
#'   rnorm(nrow(panel)) + 3 * panel$treated
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "scm")
#'
#' \donttest{
#' plot(fit, type = "trend")
#' plot(fit, type = "gap")
#' plot(fit, type = "weights")
#'
#' # Customize series colors, legend text, and reference lines
#' plot(fit, type = "trend",
#'      colors = c(Treated = "black"),
#'      labels = c(Treated = "Unit 5"),
#'      vline  = list(color = "red", linetype = "dashed"))
#' }
#' @import ggplot2
#' @export
plot.coresynth <- function(x, type = c("trend", "gap", "weights"),
                            colors = NULL, labels = NULL,
                            vline = list(), hline = list(),
                            fill = NULL, ...) {
  type <- match.arg(type)

  if(type %in% c("trend", "gap")) {
    if(is.null(x$times) || is.null(x$Y_treat))
      stop("fit object does not contain time series data for plotting.")

    times <- x$times
    # Coerce only types ggplot cannot place on a continuous/date axis; pass
    # Date/POSIXct through so the appropriate date scale is selected automatically.
    if (is.character(times) || is.factor(times))
      times <- as.numeric(as.character(times))
    # Multiple treated units are averaged per period by the accessors
    Y_treat <- treated_outcomes(x, na.rm = TRUE)
    Y_synth <- synthetic_outcomes(x, na.rm = TRUE)
    if(is.null(Y_synth))
      stop("fit object does not contain a synthetic/counterfactual series ",
           "to plot (staggered fits store their series per cohort).",
           call. = FALSE)
    treat_time  <- if(!is.null(x$T_pre)) times[x$T_pre + 1] else NA
    vline_style <- .line_style(list(color = "gray40", linetype = "dotted"), vline)

    if(type == "trend") {
      series_colors <- .merge_named_colors(
        c(Treated = "#2166ac", `Synthetic Control` = "#d73027"), colors
      )
      series_labels <- .merge_named_labels(c("Treated", "Synthetic Control"), labels)
      df <- data.frame(
        time     = c(times, times),
        value    = c(Y_treat, Y_synth),
        series   = rep(c("Treated", "Synthetic Control"), each = length(times))
      )
      # labels must be identical on both scales or the merged legend splits in two
      p <- ggplot(df, aes(x = time, y = value, color = series, linetype = series)) +
        geom_line(linewidth = 0.9) +
        scale_color_manual(values = series_colors, labels = series_labels) +
        scale_linetype_manual(values = c("Treated" = "solid", "Synthetic Control" = "dashed"),
                              labels = series_labels) +
        {if(!is.null(vline_style) && !is.na(treat_time)) do.call(geom_vline, c(list(xintercept = treat_time), vline_style))} +
        theme_minimal(base_size = 13) +
        labs(title    = paste0("Synthetic Control Trend  [", toupper(x$method), "]"),
             x = "Time", y = "Outcome", color = NULL, linetype = NULL)
      return(p)
    }

    if(type == "gap") {
      gap_color   <- if (is.null(colors)) "#1a9641" else unname(colors[[1]])
      hline_style <- .line_style(list(color = "gray50", linetype = "dashed"), hline)
      gap <- Y_treat - Y_synth
      df  <- data.frame(time = times, gap = gap)
      p <- ggplot(df, aes(x = time, y = gap)) +
        geom_line(color = gap_color, linewidth = 0.9) +
        {if(!is.null(hline_style)) do.call(geom_hline, c(list(yintercept = 0), hline_style))} +
        {if(!is.null(vline_style) && !is.na(treat_time)) do.call(geom_vline, c(list(xintercept = treat_time), vline_style))} +
        theme_minimal(base_size = 13) +
        labs(title = paste0("Treatment Effect Gap  [", toupper(x$method), "]"),
             x = "Time", y = "Y_treated - Y_synthetic")
      return(p)
    }
  }

  if(type == "weights") {
    w <- x$unit_weights
    if(is.null(w) || all(is.na(w)))
      stop("No unit weights available for this method (GSC/MC/TASC use factor loadings).")

    bar_fill <- fill %||% "#4575b4"
    df <- data.frame(
      unit   = names(w) %||% paste0("Unit_", seq_along(w)),
      weight = as.numeric(w)
    )
    df <- df[df$weight > 1e-4, ]
    if(nrow(df) == 0) stop("All unit weights are negligibly small.")

    p <- ggplot(df, aes(x = reorder(unit, weight), y = weight)) +
      geom_col(fill = bar_fill, alpha = 0.85) +
      coord_flip() +
      theme_minimal(base_size = 13) +
      labs(title = "Donor Unit Weights", x = NULL, y = "Weight")
    return(p)
  }
}

#' Plot SCM In-Space Placebo Results
#'
#' Visualizes the placebo study returned by [mspe_ratio_pval()], following
#' Abadie, Diamond & Hainmueller (2010, Section 3.4).
#'
#' `type = "gaps"` overlays the treated unit's gap path (treated minus
#' synthetic control) on the placebo gap paths obtained by reassigning the
#' intervention to each donor unit (ADH 2010, Figure 4). Placebo units whose
#' synthetic control fits poorly before treatment carry no information about
#' the rarity of a large post-treatment gap, so ADH exclude units whose
#' pre-treatment MSPE exceeds a multiple of the treated unit's: 20, 5, and 2
#' in their Figures 5-7 (`mspe_prune`).
#'
#' `type = "ratios"` shows the post/pre-treatment MSPE ratio of every unit
#' (ADH 2010, Figure 8), the statistic behind the two-sided permutation
#' p-value; it requires no pruning cutoff by construction.
#'
#' @param x A `scm_placebo` object from [mspe_ratio_pval()].
#' @param type One of `"gaps"` (ADH 2010, Figures 4-7) or `"ratios"`
#'   (ADH 2010, Figure 8).
#' @param mspe_prune Only for `type = "gaps"`: exclude placebo units whose
#'   pre-treatment MSPE exceeds `mspe_prune` times the treated unit's.
#'   Default `Inf` (no pruning). A rule stated on the RMSPE scale, such as
#'   tidysynth's "2 times the treated unit's pre-period RMSPE", corresponds
#'   to the squared multiple (`mspe_prune = 4`).
#' @param colors A named vector overriding series colors, e.g.
#'   `c(Treated = "black")`. Valid names: `"Treated"`, `"Placebo (donor pool)"`.
#' @param labels A named vector overriding the legend text of individual
#'   series, e.g. `c(Treated = "California")`. Valid names: `"Treated"`,
#'   `"Placebo (donor pool)"`. Series not mentioned keep their default label,
#'   and `colors` keys always refer to the original series names regardless
#'   of relabeling.
#' @param vline Only for `type = "gaps"`: aesthetic overrides for the vertical
#'   treatment-time line, as a list passed to [ggplot2::geom_vline()].
#'   `NULL` or `FALSE` hides the line entirely.
#' @param hline Only for `type = "gaps"`: aesthetic overrides for the
#'   horizontal zero line, as a list passed to [ggplot2::geom_hline()].
#'   `NULL` or `FALSE` hides the line entirely.
#' @param ... Ignored.
#' @return A `ggplot2` plot object.
#' @examples
#' set.seed(1)
#' panel <- expand.grid(unit = 1:10, year = 1:20)
#' panel$treated <- as.integer(panel$unit == 5 & panel$year > 15)
#' panel$gdp <- panel$unit + 0.5 * panel$year +
#'   rnorm(nrow(panel)) + 3 * panel$treated
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "scm")
#' placebo <- mspe_ratio_pval(fit)
#'
#' \donttest{
#' # Treated gap overlaid on the donor-pool placebo gaps (ADH 2010, Fig. 4)
#' plot(placebo, type = "gaps")
#'
#' # Prune poorly fitting placebos and relabel the legend
#' plot(placebo, type = "gaps", mspe_prune = 5,
#'      labels = c(Treated = "Unit 5"))
#'
#' # Post/pre-treatment MSPE ratios (ADH 2010, Fig. 8)
#' plot(placebo, type = "ratios")
#' }
#' @seealso [mspe_ratio_pval()]
#' @export
plot.scm_placebo <- function(x, type = c("gaps", "ratios"), mspe_prune = Inf,
                              colors = NULL, labels = NULL,
                              vline = list(), hline = list(), ...) {
  type <- match.arg(type)
  series_labels <- .merge_named_labels(c("Treated", "Placebo (donor pool)"), labels)
  if (!is.numeric(mspe_prune) || length(mspe_prune) != 1L || mspe_prune <= 0)
    stop("mspe_prune must be a single positive number (Inf = no pruning).")

  subtitle <- paste0(
    "Permutation p-value = ", formatC(x$p_value, digits = 3, format = "g"),
    " (", x$alternative, ", ", x$n_placebo_used, " placebo units)"
  )

  if (type == "gaps") {
    times <- x$times
    if (is.character(times) || is.factor(times))
      times <- as.numeric(as.character(times))
    treat_time <- times[x$T_pre + 1L]

    keep <- is.finite(x$mspe_pre_placebo) &
      x$mspe_pre_placebo <= mspe_prune * x$mspe_pre_treated
    n_pruned <- sum(!keep)
    if (!any(keep))
      warning("All placebo units were pruned; only the treated gap is shown. ",
              "Consider a larger mspe_prune.")

    gaps  <- x$gaps[, keep, drop = FALSE]
    # sprintf keeps zero-length input zero-length (paste0 would collapse it to "Donor ")
    units <- colnames(gaps) %||% sprintf("Donor %d", which(keep))
    df_pl <- data.frame(
      time = rep(times, times = ncol(gaps)),
      gap  = as.vector(gaps),
      unit = rep(units, each = length(times))
    )
    df_tr <- data.frame(time = times, gap = x$treated_gap)

    series_colors <- .merge_named_colors(
      c(Treated = "#2166ac", `Placebo (donor pool)` = "grey70"), colors
    )
    vline_style <- .line_style(list(color = "gray40", linetype = "dotted"), vline)
    hline_style <- .line_style(list(color = "gray50", linetype = "dashed"), hline)

    p <- ggplot() +
      geom_line(data = df_pl,
                aes(x = time, y = gap, group = unit, color = "Placebo (donor pool)"),
                linewidth = 0.4, alpha = 0.8) +
      geom_line(data = df_tr, aes(x = time, y = gap, color = "Treated"),
                linewidth = 1.0) +
      {if(!is.null(hline_style)) do.call(geom_hline, c(list(yintercept = 0), hline_style))} +
      {if(!is.null(vline_style)) do.call(geom_vline, c(list(xintercept = treat_time), vline_style))} +
      scale_color_manual(values = series_colors, breaks = names(series_colors),
                         labels = series_labels) +
      theme_minimal(base_size = 13) +
      labs(title = "Placebo Gaps in the Donor Pool  [SCM]",
           subtitle = subtitle,
           x = "Time", y = "Gap (unit - synthetic control)",
           color = NULL,
           caption = if (n_pruned > 0L) {
             paste0("Pruned ", n_pruned, " placebo unit(s) with pre-treatment MSPE > ",
                    mspe_prune, "x the treated unit's (ADH 2010, Figures 5-7).")
           })
    return(p)
  }

  # type == "ratios" (ADH 2010, Figure 8)
  r     <- x$mspe_ratios_all
  units <- names(r)
  if (is.null(units)) units <- c("Treated", sprintf("Donor %d", seq_len(length(r) - 1L)))
  # keep the axis tick consistent with the (possibly relabeled) legend entry
  units[1L] <- series_labels[["Treated"]]
  blank <- !nzchar(units)
  units[blank] <- sprintf("Donor %d", which(blank) - 1L)

  df <- data.frame(unit = units, ratio = as.numeric(r),
                   series = c("Treated", rep("Placebo (donor pool)", length(r) - 1L)))
  n_dropped <- sum(!is.finite(df$ratio))
  df <- df[is.finite(df$ratio), ]
  if (nrow(df) == 0L) stop("No finite MSPE ratios to plot.")

  series_colors <- .merge_named_colors(
    c(Treated = "#2166ac", `Placebo (donor pool)` = "grey60"), colors
  )

  p <- ggplot(df, aes(x = ratio, y = reorder(unit, ratio), color = series)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = series_colors, breaks = names(series_colors),
                       labels = series_labels) +
    theme_minimal(base_size = 13) +
    labs(title = "Post/Pre-Treatment MSPE Ratios  [SCM]",
         subtitle = subtitle,
         x = "MSPE ratio (post / pre)", y = NULL, color = NULL,
         caption = if (n_dropped > 0L) {
           paste0(n_dropped, " unit(s) without a finite ratio (mspe_threshold filter) omitted.")
         })
  p
}
