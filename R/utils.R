#' Panel Data Helper: Reshape long-format data to matrices
#'
#' Internal utility used by all fit_* functions.
#'
#' @param y      Numeric outcome vector (long format)
#' @param d      Binary treatment indicator (0/1, long format)
#' @param id     Unit identifier (long format)
#' @param time   Time identifier (long format)
#' @return A named list with:
#'   * `Y`: Outcome matrix (T x N), units as columns
#'   * `D`: Treatment matrix (T x N)
#'   * `units`: Unique unit identifiers
#'   * `times`: Unique time identifiers
#'   * `T_pre`: Number of pre-treatment periods (global minimum across treated units)
#'   * `T_adopt`: Integer vector (length N). Per-unit first treated row; NA for controls.
#'   * `is_sharp`: Logical. TRUE iff all treated units share one adoption date.
#'   * `idx_treat`: Column indices of treated units
#'   * `idx_control`: Column indices of control units
#' @noRd
panel_to_matrices <- function(y, d, id, time) {
  # Sort by id then time
  ord <- order(id, time)
  y <- y[ord]
  d <- d[ord]
  id <- id[ord]
  time <- time[ord]

  units <- unique(id)
  times <- unique(sort(time))
  N <- length(units)
  TT <- length(times)

  Y <- matrix(
    NA_real_,
    nrow = TT,
    ncol = N,
    dimnames = list(as.character(times), as.character(units))
  )
  D <- matrix(
    0L,
    nrow = TT,
    ncol = N,
    dimnames = list(as.character(times), as.character(units))
  )

  for (i in seq_along(id)) {
    ri <- which(times == time[i])
    ci <- which(units == id[i])
    Y[ri, ci] <- y[i]
    D[ri, ci] <- d[i]
  }

  # Determine treated units: any unit that is ever treated
  ever_treated <- which(colSums(D) > 0)
  ever_control <- which(colSums(D) == 0)

  # Per-unit adoption row index (first row where D[,j] > 0, i.e. any treatment arm)
  T_adopt <- rep(NA_integer_, N)
  for (j in ever_treated) {
    first_t <- which(D[, j] > 0L)
    if (length(first_t) > 0L) T_adopt[j] <- first_t[1L]
  }

  # Global T_pre = min(T_adopt) - 1 (backward-compatible semantics)
  T_pre <- min(T_adopt[ever_treated], na.rm = TRUE) - 1L
  is_sharp <- length(unique(T_adopt[ever_treated])) == 1L

  list(
    Y = Y,
    D = D,
    units = units,
    times = times,
    T_pre = T_pre,
    T_adopt = T_adopt,
    is_sharp = is_sharp,
    idx_treat = ever_treated,
    idx_control = ever_control
  )
}

check_sharp_adoption <- function(pan, method) {
  if (!isTRUE(pan$is_sharp)) {
    stop(
      sprintf(
        paste0(
          "'%s' does not support staggered adoption ",
          "(units with different treatment timing detected). ",
          "Use method = 'mc', 'tasc', or 'sdid' instead."
        ),
        method
      ),
      call. = FALSE
    )
  }
}

#' Predictor Specification for SCM
#'
#' Creates a single predictor specification for use in [scm_fit()] with
#' `method = "scm"`. Pass a `list()` of `pred()` calls as the `predictors`
#' argument to define the full covariate matrix.
#'
#' @param vars  Character vector of variable names. All variables share the
#'   same `times` window and `op` operator. Use separate `pred()` calls for
#'   variables with different time windows.
#' @param times Numeric/integer vector of time values to aggregate over.
#' @param op    Aggregation operator applied to each variable over `times`.
#'   One of `"mean"` (default), `"median"`, or `"sum"`.
#'
#' @return A `pred_spec` object (a named list with class `"pred_spec"`).
#'
#' @seealso [scm_fit()] for the `predictors` argument that consumes a `list()`
#'   of `pred_spec` objects.
#'
#' @export
#'
#' @examples
#' # Three variables averaged over the same window
#' pred(c("lnincome", "retprice", "age15to24"), 1980:1988)
#'
#' # Single variable at a specific year
#' pred("cigsale", 1975)
#'
#' # Single variable averaged over a range
#' pred("beer", 1984:1988)
#'
#' # Abadie, Diamond & Hainmueller (2010) California Prop 99 style:
#' # combine several covariates aggregated over different windows plus
#' # three outcome lags at specific years. Pass the list to scm_fit(..., predictors = ...).
#' \dontrun{
#' predictors <- list(
#'   pred(c("lnincome", "retprice", "age15to24"), 1980:1988),
#'   pred("beer",    1984:1988),
#'   pred("cigsale", 1988),
#'   pred("cigsale", 1980),
#'   pred("cigsale", 1975)
#' )
#' fit <- scm_fit(cigsale ~ treated | state + year,
#'                data = prop99, method = "scm",
#'                predictors = predictors)
#' }
pred <- function(vars, times, op = "mean") {
  if (!is.character(vars) || length(vars) == 0L) {
    stop("'vars' must be a non-empty character vector.", call. = FALSE)
  }
  if (length(times) == 0L) {
    stop("'times' must be a non-empty vector of time values.", call. = FALSE)
  }
  if (!op %in% c("mean", "median", "sum")) {
    stop("'op' must be one of \"mean\", \"median\", or \"sum\".", call. = FALSE)
  }
  structure(list(vars = vars, times = times, op = op), class = "pred_spec")
}

#' @export
print.pred_spec <- function(x, ...) {
  times_str <- if (length(x$times) == 1L) {
    as.character(x$times)
  } else {
    sprintf("%s:%s", min(x$times), max(x$times))
  }
  cat(sprintf(
    "pred(%s, %s, op = \"%s\")\n",
    paste(x$vars, collapse = ", "),
    times_str,
    x$op
  ))
  invisible(x)
}

#' Build predictor matrices X0 and X1 for SCM
#'
#' Constructs the (k x N_co) predictor matrix X0 and (k x 1) vector X1
#' from a list of [pred()] specifications, following Abadie et al. (2010)
#' Section 2.3. Each `pred()` entry expands to one row per variable.
#'
#' @param data       Full long-format data frame.
#' @param id_var     Name of the unit identifier column.
#' @param time_var   Name of the time identifier column.
#' @param units      All unit identifiers (length N), same order as Y columns.
#' @param idx_co     Integer indices of control units in `units`.
#' @param idx_tr     Integer index of the treated unit in `units`.
#' @param predictors List of `pred_spec` objects created by [pred()].
#' @return A list with:
#'   * `X0`: k x N_co numeric matrix (predictors x control units)
#'   * `X1`: numeric vector of length k (predictors for treated unit)
#'   * `pred_names`: character vector of length k with predictor labels
#' @noRd
build_predictor_matrices <- function(
  data,
  id_var,
  time_var,
  units,
  idx_co,
  idx_tr,
  predictors
) {
  co_units <- units[idx_co]

  agg_unit <- function(var, times, op) {
    fn <- match.fun(op)
    vapply(
      units,
      function(u) {
        sub <- data[[var]][data[[id_var]] == u & data[[time_var]] %in% times]
        if (length(sub) == 0L) NA_real_ else fn(sub, na.rm = TRUE)
      },
      numeric(1L)
    )
  }

  pred_label <- function(var, times) {
    if (length(times) == 1L) {
      sprintf("%s[%s]", var, times)
    } else {
      sprintf("%s[%s:%s]", var, min(times), max(times))
    }
  }

  rows_X0 <- list()
  rows_X1 <- list()
  pred_names <- character(0L)

  for (p in predictors) {
    if (!inherits(p, "pred_spec")) {
      stop(
        paste0(
          "Each element of 'predictors' must be a pred_spec object ",
          "created by pred()."
        ),
        call. = FALSE
      )
    }
    for (var in p$vars) {
      if (!var %in% names(data)) {
        stop(sprintf("Variable '%s' not found in data.", var), call. = FALSE)
      }
      vals <- agg_unit(var, p$times, p$op)
      rows_X0 <- c(rows_X0, list(vals[idx_co]))
      rows_X1 <- c(rows_X1, list(vals[idx_tr]))
      pred_names <- c(pred_names, pred_label(var, p$times))
    }
  }

  X0 <- do.call(rbind, rows_X0) # k x N_co
  X1 <- unlist(rows_X1) # length k

  colnames(X0) <- as.character(co_units)
  rownames(X0) <- pred_names

  list(X0 = X0, X1 = X1, pred_names = pred_names)
}

#' Build a time-varying covariate array for GSC
#'
#' Constructs a T × N × p 3D R array suitable for passing to gsc_ife_cpp as an
#' arma::cube. arr[t, i, j] is the value of the j-th covariate for unit i at
#' time t.
#'
#' @param data            Long-format data frame.
#' @param id_var          Name of the unit identifier column.
#' @param time_var        Name of the time identifier column.
#' @param covariate_names Character vector of covariate column names (length p).
#' @param units           Character vector of unit IDs (length N).
#' @param times           Vector of time values (length T).
#' @return A T × N × p numeric array.
#' @noRd
build_covariate_array <- function(
  data,
  id_var,
  time_var,
  covariate_names,
  units,
  times
) {
  T_all <- length(times)
  N <- length(units)
  p <- length(covariate_names)
  arr <- array(NA_real_, dim = c(T_all, N, p))
  for (j in seq_len(p)) {
    var <- covariate_names[j]
    if (!var %in% names(data)) {
      stop(sprintf("covariate '%s' not found in data.", var), call. = FALSE)
    }
    for (i in seq_len(N)) {
      idx <- data[[id_var]] == units[i]
      sub_t <- data[[time_var]][idx]
      sub_v <- data[[var]][idx]
      arr[match(sub_t, times), i, j] <- sub_v
    }
  }
  arr
}

#' Panel Data Helper: Reshape long-format multi-arm data to tensor structure
#'
#' Extends [panel_to_matrices()] for the multi-arm Synthetic Interventions
#' setting (Agarwal et al. 2025). Each unit belongs to exactly one treatment
#' arm (d = 0 for control, d = 1,...,K for treatment arms). Before the
#' treatment date, all treatment-arm units have d = 0; at the treatment date
#' they switch to their assigned arm value.
#'
#' @param y    Numeric outcome vector (long format)
#' @param d    Treatment arm indicator (integer, 0 = control, 1,...,K = arms)
#' @param id   Unit identifier (long format)
#' @param time Time identifier (long format)
#' @return All fields returned by [panel_to_matrices()] plus:
#'   * `arm_levels`: sorted integer vector of unique arm values (`c(0L, 1L, ..., KL)`)
#'   * `idx_by_arm`: named list of column indices, one entry per arm level
#' @noRd
panel_to_tensor <- function(y, d, id, time) {
  d_int <- as.integer(d)
  pan   <- panel_to_matrices(y, d_int, id, time)

  # Per-unit arm = max(D[, j]): control units always 0, treated units = their arm
  arm_of_unit <- as.integer(apply(pan$D, 2, max))
  arm_levels  <- sort(unique(arm_of_unit))

  if (!0L %in% arm_levels)
    stop("panel_to_tensor: arm 0 (control) が存在しません。", call. = FALSE)
  if (any(arm_levels < 0L))
    stop("panel_to_tensor: arm は非負整数である必要があります。", call. = FALSE)

  idx_by_arm <- setNames(
    lapply(arm_levels, function(a) which(arm_of_unit == a)),
    as.character(arm_levels)
  )

  c(pan, list(arm_levels = arm_levels, idx_by_arm = idx_by_arm))
}
