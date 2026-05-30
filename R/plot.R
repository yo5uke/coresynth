#' Plot a coresynth model
#'
#' @param x     A `coresynth` object.
#' @param type  One of `"trend"` (observed vs synthetic), `"gap"` (ATT over time),
#'              or `"weights"` (donor unit weight bar chart).
#' @param ...   Ignored.
#' @return A `ggplot2` plot object.
#' @import ggplot2
#' @export
plot.coresynth <- function(x, type = c("trend", "gap", "weights"), ...) {
  type <- match.arg(type)

  if(type %in% c("trend", "gap")) {
    if(is.null(x$times) || is.null(x$Y_treat))
      stop("fit object does not contain time series data for plotting.")

    times <- as.numeric(x$times)
    # Handle matrix Y_treat (multiple treated units) by averaging
    Y_treat <- if(is.matrix(x$Y_treat)) rowMeans(x$Y_treat, na.rm = TRUE) else x$Y_treat
    Y_synth <- if(!is.null(x$Y_synth)) {
      if(is.matrix(x$Y_synth)) rowMeans(x$Y_synth, na.rm = TRUE) else x$Y_synth
    } else {
      # For GSC/MC/TASC, compute from Y_hat
      Y_hat <- if(!is.null(x$Y_hat)) x$Y_hat else x$Y_tr_hat
      if(is.matrix(Y_hat)) rowMeans(Y_hat, na.rm = TRUE) else Y_hat
    }
    treat_time <- if(!is.null(x$T_pre)) times[x$T_pre + 1] else NA

    if(type == "trend") {
      df <- data.frame(
        time     = c(times, times),
        value    = c(Y_treat, Y_synth),
        series   = rep(c("Treated", "Synthetic Control"), each = length(times))
      )
      p <- ggplot(df, aes(x = time, y = value, color = series, linetype = series)) +
        geom_line(linewidth = 0.9) +
        scale_color_manual(values = c("Treated" = "#2166ac", "Synthetic Control" = "#d73027")) +
        scale_linetype_manual(values = c("Treated" = "solid", "Synthetic Control" = "dashed")) +
        {if(!is.na(treat_time)) geom_vline(xintercept = treat_time, linetype = "dotted", color = "gray40")} +
        theme_minimal(base_size = 13) +
        labs(title    = paste0("Synthetic Control Trend  [", toupper(x$method), "]"),
             x = "Time", y = "Outcome", color = NULL, linetype = NULL)
      return(p)
    }

    if(type == "gap") {
      gap <- Y_treat - Y_synth
      df  <- data.frame(time = times, gap = gap)
      p <- ggplot(df, aes(x = time, y = gap)) +
        geom_line(color = "#1a9641", linewidth = 0.9) +
        geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
        {if(!is.na(treat_time)) geom_vline(xintercept = treat_time, linetype = "dotted", color = "gray40")} +
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

    df <- data.frame(
      unit   = names(w) %||% paste0("Unit_", seq_along(w)),
      weight = as.numeric(w)
    )
    df <- df[df$weight > 1e-4, ]
    if(nrow(df) == 0) stop("All unit weights are negligibly small.")

    p <- ggplot(df, aes(x = reorder(unit, weight), y = weight)) +
      geom_col(fill = "#4575b4", alpha = 0.85) +
      coord_flip() +
      theme_minimal(base_size = 13) +
      labs(title = "Donor Unit Weights", x = NULL, y = "Weight")
    return(p)
  }
}
