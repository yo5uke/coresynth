# ── Accessor generics for coresynth fits ─────────────────────────────────────
#
# Each estimator stores its outcome series under method-specific field names
# (SCM/SDID: Y_co_pre/Y_co_post, SI: Y_pre_co/Y_post_co, GSC/MC: Y_co_all,
# fitted values: Y_synth/Y_tr_hat/Y_hat). These generics are the single place
# that knows the mapping, so downstream code (broom methods, conformal
# inference, plotting) can extract series without per-method field sniffing.

.rbind_or_null <- function(pre, post) {
  if (is.null(pre) || is.null(post)) return(NULL)
  rbind(pre, post)
}

#' Extract Outcome Series from a coresynth Fit
#'
#' Accessor generics that return the outcome series stored in a fitted
#' `coresynth` object under a uniform interface, regardless of the estimation
#' method:
#'
#' * `treated_outcomes()`: the treated unit's observed outcome series
#'   (length \eqn{T}). When several units are treated, their per-period mean.
#' * `synthetic_outcomes()`: the estimated counterfactual series
#'   (length \eqn{T}), i.e. the synthetic control or model-fitted outcome.
#' * `donor_outcomes()`: the \eqn{T \times N_{co}} matrix of observed donor
#'   (control unit) outcomes over all periods.
#'
#' Each accessor returns `NULL` when the requested series is not stored in
#' the fit. In particular, staggered-adoption fits keep their data per cohort
#' (in `fit$cohort_fits`), so the sharp-fit accessors return `NULL` for them.
#'
#' @param x A `coresynth` object from [scm_fit()].
#' @param na.rm Logical; passed to the per-period averaging over multiple
#'   treated units (default `FALSE`).
#' @param ... Passed to methods.
#' @return For `treated_outcomes()` and `synthetic_outcomes()`, a numeric
#'   vector of length \eqn{T}, or `NULL`. For `donor_outcomes()`, a
#'   \eqn{T \times N_{co}} numeric matrix (donors in columns, named when unit
#'   names are available), or `NULL`.
#' @examples
#' set.seed(1)
#' panel <- expand.grid(unit = 1:10, year = 1:20)
#' panel$treated <- as.integer(panel$unit == 1 & panel$year > 15)
#' panel$gdp <- panel$unit + 0.5 * panel$year +
#'   rnorm(nrow(panel)) + 3 * panel$treated
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "scm")
#'
#' y1   <- treated_outcomes(fit)    # observed treated series
#' y1_0 <- synthetic_outcomes(fit)  # synthetic counterfactual
#' Yco  <- donor_outcomes(fit)      # donor outcome matrix
#' all.equal(y1 - y1_0, unname(fit$gap))
#' @export
treated_outcomes <- function(x, ...) UseMethod("treated_outcomes")

#' @rdname treated_outcomes
#' @export
treated_outcomes.coresynth <- function(x, na.rm = FALSE, ...) {
  y <- x$Y_treat
  if (is.null(y)) return(NULL)
  if (is.matrix(y)) rowMeans(y, na.rm = na.rm) else as.numeric(y)
}

#' @rdname treated_outcomes
#' @export
synthetic_outcomes <- function(x, ...) UseMethod("synthetic_outcomes")

#' @rdname treated_outcomes
#' @export
synthetic_outcomes.coresynth <- function(x, na.rm = FALSE, ...) {
  s <- x$Y_synth %||% x$Y_tr_hat %||% x$Y_hat
  if (!is.null(s)) {
    return(if (is.matrix(s)) rowMeans(s, na.rm = na.rm) else as.numeric(s))
  }
  # Last resort: reconstruct from the stored gap (counterfactual = Y - gap)
  g <- x$gap
  y <- treated_outcomes(x, na.rm = na.rm)
  if (is.null(g) || is.null(y)) return(NULL)
  g <- if (is.matrix(g)) rowMeans(g, na.rm = na.rm) else as.numeric(g)
  y - g
}

#' @rdname treated_outcomes
#' @export
synthetic_outcomes.coresynth_tasc <- function(x, na.rm = FALSE, ...) {
  # TASC's Y_hat holds fitted values for all N units (T x N), unlike every
  # other method; the counterfactual is the treated columns, not the
  # all-unit average.
  if (!is.null(x$Y_hat) && !is.null(x$idx_tr))
    return(rowMeans(x$Y_hat[, x$idx_tr, drop = FALSE], na.rm = na.rm))
  # Fits from versions without idx_tr: reconstruct from the stored gap
  g <- x$gap
  y <- treated_outcomes(x, na.rm = na.rm)
  if (is.null(g) || is.null(y)) return(NULL)
  g <- if (is.matrix(g)) rowMeans(g, na.rm = na.rm) else as.numeric(g)
  y - g
}

#' @rdname treated_outcomes
#' @export
donor_outcomes <- function(x, ...) UseMethod("donor_outcomes")

#' @rdname treated_outcomes
#' @export
donor_outcomes.coresynth_scm <- function(x, ...) {
  .rbind_or_null(x$Y_co_pre, x$Y_co_post)
}

#' @rdname treated_outcomes
#' @export
donor_outcomes.coresynth_sdid <- function(x, ...) {
  .rbind_or_null(x$Y_co_pre, x$Y_co_post)
}

#' @rdname treated_outcomes
#' @export
donor_outcomes.coresynth_si <- function(x, ...) {
  .rbind_or_null(x$Y_pre_co, x$Y_post_co)
}

#' @rdname treated_outcomes
#' @export
donor_outcomes.coresynth_gsc <- function(x, ...) x$Y_co_all

#' @rdname treated_outcomes
#' @export
donor_outcomes.coresynth_mc <- function(x, ...) x$Y_co_all

#' @rdname treated_outcomes
#' @export
donor_outcomes.coresynth_tasc <- function(x, ...) NULL

#' @rdname treated_outcomes
#' @export
donor_outcomes.coresynth <- function(x, ...) {
  # Fallback for objects created by earlier package versions: try every
  # known field layout in turn.
  .rbind_or_null(x$Y_co_pre, x$Y_co_post) %||%
    .rbind_or_null(x$Y_pre_co, x$Y_post_co) %||%
    x$Y_co_all
}
