#' Fit a Synthetic Control Method Model
#'
#' Unified formula interface for Synthetic Control and related causal
#' inference methods.  The formula syntax is:
#'
#'   `outcome ~ treatment | unit_id + time_id`
#'
#' @param formula A `Formula` object, e.g. `y ~ D | unit + time`.
#' @param data    A `data.frame` in **long** format (one row per unit-time).
#' @param method  One of `"scm"`, `"sdid"`, `"gsc"`, `"mc"`, `"tasc"`, `"si"`.
#' @param predictors A `list()` of [pred()] specifications that define the
#'   predictor matrix for SCM (see Abadie et al. 2010, S.2.3). Each [pred()]
#'   entry aggregates one or more variables over a time window. Pass `NULL`
#'   (default) to use all pre-treatment outcome periods as predictors.
#'   Applies to `method = "scm"` only.
#' @param covariates An optional named `list` of additional time-varying
#'   covariates to partial out before estimation. Each element is a character
#'   string naming a column in `data`. Supported for `method = "sdid"`,
#'   `"scm"`, and `"gsc"`.
#' @param v_selection V matrix selection method for `method = "scm"`.
#'   `"insample"` (default) follows Abadie et al. (2010): V is chosen by
#'   minimising in-sample pre-treatment MSPE. `"oos"` follows Abadie (2021)
#'   S.3.2: the pre-treatment window is split in half; V is selected to minimise
#'   MSPE on the validation half, then W is refit on the full window.
#' @param donor_mspe_threshold Donor pool filtering threshold (Abadie 2021 S.4).
#'   For `method = "scm"` only. Each donor's individual pre-treatment MSPE
#'   (using that donor alone as the counterfactual) is divided by the minimum
#'   such MSPE across all donors. Donors whose ratio exceeds this threshold are
#'   excluded from estimation. `Inf` (default) disables filtering.
#' @param lambda_pen Penalised SCM parameter (Abadie & L'Hour 2021, JASA).
#'   For `method = "scm"` only. `NULL` (default) runs standard unpenalised SCM.
#'   `"auto"` selects the penalty via out-of-sample pre-treatment MSPE on the
#'   same validation window as `v_selection = "oos"`. A non-negative number
#'   uses that value directly.
#' @param v_optim Outer V-optimisation method for `method = "scm"`.
#'   `"coord_descent"` (default) uses the existing C++ coordinate descent with
#'   11-point grid search -- fastest when `k = T_pre` is large (outcomes-only).
#'   `"bfgs"` uses R's L-BFGS-B, which requires only O(k^2) inner QP calls and
#'   is faster when `k` is small (e.g. external predictors with k <= 15).
#'   `"auto"` selects `"bfgs"` when `k <= 15`, otherwise `"coord_descent"`.
#' @param ...     Additional arguments forwarded to the specific method
#'                (e.g. `r`, `lambda`, `zeta2`).
#'
#' @return An object of classes `c("coresynth_<method>", "coresynth")`.
#'   All methods return at minimum:
#'   * `method`: estimator name
#'   * `estimate`: average treatment effect (ATT)
#'   * `times`: time index vector
#'   * `T_pre`: number of pre-treatment periods
#'   * `Y_treat`: treated unit outcome series
#'   * `gap`: treatment effect series (Y_treat - counterfactual)
#'
#' @export
#'
#' @examples
#' # Synthetic balanced panel: 10 units over 20 periods, unit 1 treated
#' # after period 15.
#' set.seed(1)
#' panel <- expand.grid(unit = 1:10, year = 1:20)
#' panel$treated <- as.integer(panel$unit == 1 & panel$year > 15)
#' panel$gdp <- panel$unit + 0.5 * panel$year +
#'   rnorm(nrow(panel)) + 3 * panel$treated
#'
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "sdid")
#' summary(fit)
#'
#' \donttest{
#' # Visualise the estimated gap (requires ggplot2)
#' plot(fit, type = "gap")
#' }
scm_fit <- function(
  formula,
  data,
  method = c("scm", "sdid", "gsc", "mc", "tasc", "si"),
  predictors = NULL,
  covariates = NULL,
  v_selection = c("insample", "oos"),
  donor_mspe_threshold = Inf,
  lambda_pen = NULL,
  v_optim = c("coord_descent", "auto", "bfgs"),
  ...
) {
  v_selection <- match.arg(v_selection)
  v_optim     <- match.arg(v_optim)
  method <- match.arg(method)

  # Parse Formula
  f_parts <- Formula::Formula(formula)

  if (length(f_parts)[2] < 2) {
    stop(
      "Formula must specify unit and time after '|', e.g. y ~ D | unit + time"
    )
  }

  y_var <- all.vars(formula(f_parts, lhs = 1, rhs = 0))
  d_var <- all.vars(formula(f_parts, lhs = 0, rhs = 1))
  idx_vars <- all.vars(formula(f_parts, lhs = 0, rhs = 2))

  if (length(y_var) != 1) {
    stop("Exactly one outcome variable required.")
  }
  if (length(d_var) != 1) {
    stop("Exactly one treatment variable required.")
  }
  if (length(idx_vars) != 2) {
    stop("Must specify exactly two index variables: unit_id and time_id.")
  }

  id_var <- idx_vars[1]
  time_var <- idx_vars[2]

  # Validate columns exist
  for (v in c(y_var, d_var, id_var, time_var)) {
    if (!v %in% names(data)) {
      stop(paste0("Variable '", v, "' not found in data."))
    }
  }

  y <- as.numeric(data[[y_var]])
  d <- as.integer(data[[d_var]])
  id <- data[[id_var]]
  time <- data[[time_var]]

  # Dispatch
  res <- switch(
    method,
    "scm" = fit_scm_cpp(
      y,
      d,
      id,
      time,
      data = data,
      id_var = id_var,
      time_var = time_var,
      predictors = predictors,
      covariates = covariates,
      v_selection = v_selection,
      donor_mspe_threshold = donor_mspe_threshold,
      lambda_pen = lambda_pen,
      v_optim = v_optim,
      ...
    ),
    "sdid" = fit_sdid_cpp(
      y, d, id, time,
      covariates = covariates,
      data       = data,
      id_var     = id_var,
      time_var   = time_var,
      ...
    ),
    "gsc" = fit_gsc_cpp(
      y,
      d,
      id,
      time,
      data = data,
      id_var = id_var,
      time_var = time_var,
      covariates = covariates,
      ...
    ),
    "mc" = fit_mc_cpp(y, d, id, time, ...),
    "tasc" = fit_tasc_cpp(y, d, id, time, ...),
    "si" = fit_si_cpp(y, d, id, time, ...),
    stop(paste0("Unknown method: '", method, "'"))
  )

  class(res) <- c(paste0("coresynth_", method), "coresynth")
  res
}

#' @export
print.coresynth <- function(x, ...) {
  cat("=== coresynth fit ===\n")
  cat("Method :", toupper(x$method), "\n")
  if (isTRUE(x$multi_arm)) {
    stag_label <- if (isTRUE(x$staggered)) " (staggered)" else ""
    cat(sprintf("Multi-arm SI%s (K = %d arms)\n", stag_label, length(x$arm_levels)))
    cat("Per-arm ATT:",
        paste(names(x$arm_estimates), round(x$arm_estimates, 4),
              sep = "=", collapse = "  "), "\n")
  }
  cat("Estimate (ATT):", round(x$estimate, 4), "\n")
  cat("Pre-treatment periods:", x$T_pre, "\n")
  invisible(x)
}

#' @export
summary.coresynth <- function(object, ...) {
  cat("=== coresynth summary ===\n")
  cat("Method :", toupper(object$method), "\n")
  if (isTRUE(object$multi_arm)) {
    stag_label <- if (isTRUE(object$staggered)) " (staggered)" else ""
    cat(sprintf("Multi-arm SI%s (K = %d treatment arms)\n",
                stag_label, length(object$arm_levels)))
    cat("Per-arm ATT:\n")
    print(round(object$arm_estimates, 6))
  }
  cat(
    "Periods : T_pre =",
    object$T_pre,
    "| T_post =",
    length(object$times) - object$T_pre,
    "\n"
  )
  cat("ATT estimate:", round(object$estimate, 6), "\n")
  if (!is.null(object$unit_weights)) {
    cat("Unit weights (non-zero donors):\n")
    w <- object$unit_weights
    print(round(w[w > 1e-4], 4))
  }
  if (!is.null(object$predictor_table)) {
    cat("\nPredictor balance:\n")
    pt <- object$predictor_table
    pt$treated <- round(pt$treated, 4)
    pt$synthetic <- round(pt$synthetic, 4)
    print(pt, row.names = FALSE)
  }
  invisible(object)
}
