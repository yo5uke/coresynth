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
  if (length(y) == 0L) {
    stop("Panel data have 0 rows. Check that upstream filtering or cleaning ",
         "did not remove all observations before calling scm_fit().",
         call. = FALSE)
  }
  if (anyNA(id) || anyNA(time)) {
    stop("Unit or time identifiers contain NA values. Drop or recode these ",
         "rows before fitting.", call. = FALSE)
  }
  if (anyNA(d)) {
    stop(sprintf(paste0(
      "Treatment indicator contains %d NA value(s). Units with NA treatment ",
      "status cannot be classified as treated or control; recode or drop ",
      "these rows before fitting."), sum(is.na(d))), call. = FALSE)
  }
  if (any(d < 0)) {
    stop("Treatment indicator must be non-negative (0 = control, ",
         "1 = treated, 2+ = additional arms for method = 'si').",
         call. = FALSE)
  }

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

  # Vectorised fill (avoids O(n * (T + N)) which() lookups inside a loop)
  idx_mat <- cbind(match(time, times), match(id, units))

  # A balanced panel requires each (id, time) cell to be unique. Duplicates would
  # otherwise be silently overwritten by the last assignment, dropping rows.
  dup <- duplicated(idx_mat)
  if (any(dup)) {
    i <- which(dup)[1L]
    stop(sprintf(
      paste0("Duplicate (id, time) entries detected (%d rows): unit '%s' at ",
             "time '%s' appears more than once. Each unit-time cell must be ",
             "unique (balanced panel)."),
      sum(dup), as.character(id[i]), as.character(time[i])),
      call. = FALSE)
  }

  Y[idx_mat] <- y
  D[idx_mat] <- d

  # Determine treated units: any unit that is ever treated
  ever_treated <- which(colSums(D) > 0)
  ever_control <- which(colSums(D) == 0)

  if (length(ever_treated) == 0L) {
    stop("No treated units found: the treatment indicator is 0 for all ",
         "observations. Check that the treated unit(s) were not dropped ",
         "during data preparation and that the treatment column is coded ",
         "correctly.", call. = FALSE)
  }

  # Per-unit adoption row index (first row where D[,j] > 0, i.e. any treatment arm)
  T_adopt <- rep(NA_integer_, N)
  for (j in ever_treated) {
    first_t <- which(D[, j] > 0L)
    if (length(first_t) > 0L) T_adopt[j] <- first_t[1L]
  }

  # Global T_pre = min(T_adopt) - 1 (backward-compatible semantics)
  T_pre <- min(T_adopt[ever_treated], na.rm = TRUE) - 1L
  is_sharp <- length(unique(T_adopt[ever_treated])) == 1L

  # Under staggered adoption a panel with no never-treated units can still be
  # estimable (future adopters serve as clean controls), so only reject the
  # sharp case where the donor pool is structurally empty.
  if (length(ever_control) == 0L && is_sharp) {
    stop("No control units found: every unit is treated at the same time, ",
         "leaving no donor pool. Check the treatment coding and upstream ",
         "filtering.", call. = FALSE)
  }

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

#' Validate that the outcome matrix is fully observed
#'
#' SCM/SDID/GSC/SI operate on dense matrices and would otherwise fail deep in
#' the C++ solvers with cryptic messages (e.g. "eig_sym(): decomposition
#' failed") when the panel has missing cells. MC and TASC handle missing
#' outcomes by design and must not call this.
#'
#' @param Y            Outcome matrix (T x N) from `panel_to_matrices()`.
#' @param method_label Estimator name used in the error message.
#' @noRd
.check_panel_complete <- function(Y, method_label) {
  bad <- !is.finite(Y)
  n_bad <- sum(bad)
  if (n_bad == 0L) return(invisible(NULL))
  ij <- which(bad, arr.ind = TRUE)[1L, ]
  stop(sprintf(paste0(
    "%s requires a balanced panel with fully observed outcomes, but %d ",
    "unit-time cell(s) are missing or non-finite (first: unit '%s' at time ",
    "'%s'). This happens when (id, time) rows are absent from the data or ",
    "the outcome contains NA/Inf values. Complete the panel, drop the ",
    "affected units, or use method = 'mc' or 'tasc', which handle missing ",
    "outcomes."),
    method_label, n_bad, colnames(Y)[ij[2L]], rownames(Y)[ij[1L]]),
    call. = FALSE)
}

#' Aggregate per-cohort fits into a staggered ATT
#'
#' Shared tail of every staggered engine: combines the per-cohort estimates
#' into the weighted-average ATT of Clarke et al. (2024) and builds the
#' `cohort_estimates` table. Each element of `r` must carry `cohort`,
#' `n_treated`, `T_pre`, `T_post`, `estimate`, and `weight`
#' (`N_treated,g * T_post,g`); `weight` is renormalised to sum to one in the
#' returned table.
#'
#' Supplying `unit_weights` replaces the treated-unit count in each block's
#' weight, so a block contributes `sum(weights of its units) * T_post,g`. This
#' is how Abadie (2021, S.8) describes reporting a population- or GDP-weighted
#' average effect, which is often the quantity of interest when the units are
#' regions or countries.
#'
#' @param r List of per-cohort fit lists (already filtered to the valid ones).
#' @param unit_weights Optional named numeric vector over treated unit labels.
#' @param unit_labels Optional character vector naming every panel column, so
#'   each block's `idx_tr` can be matched against `unit_weights`.
#' @return A list with `estimate`, `cohort_estimates`, and `cohort_fits`.
#' @noRd
.aggregate_cohorts <- function(r, unit_weights = NULL, unit_labels = NULL) {
  if (!is.null(unit_weights)) {
    if (is.null(unit_labels))
      stop("internal: unit_weights needs unit_labels.", call. = FALSE)
    r <- lapply(r, function(z) {
      lab <- as.character(unit_labels[z$idx_tr])
      wt  <- sum(unit_weights[lab])
      z$weight    <- wt * z$T_post
      z$unit_mass <- wt
      z
    })
  }
  w   <- vapply(r, `[[`, numeric(1L), "weight")
  tau <- vapply(r, `[[`, numeric(1L), "estimate")

  cohort_df <- data.frame(
    cohort    = vapply(r, `[[`, integer(1L), "cohort"),
    n_treated = vapply(r, `[[`, integer(1L), "n_treated"),
    T_pre     = vapply(r, `[[`, integer(1L), "T_pre"),
    T_post    = vapply(r, `[[`, integer(1L), "T_post"),
    estimate  = tau,
    weight    = w / sum(w),
    stringsAsFactors = FALSE
  )

  list(
    estimate         = sum(w * tau) / sum(w),
    cohort_estimates = cohort_df,
    cohort_fits      = r
  )
}

#' Cohort table for the imputation estimators (MC, TASC)
#'
#' MC and TASC recover the counterfactual from a single global fit: treated
#' post-adoption cells are masked per unit, so staggered adoption needs no
#' cohort-by-cohort re-fitting. The cohort table is therefore a partition of
#' that one fit's gaps rather than a set of separate fits, and every cohort
#' shares the same donor set -- the never-treated units, plus the observed
#' (pre-adoption) cells of the other treated units, which the factorisation
#' uses too.
#'
#' @param pan     Result of `panel_to_matrices()`.
#' @param gap     T x N_treated matrix of observed minus imputed outcomes.
#' @param Y_treat T x N_treated matrix of observed treated outcomes.
#' @param Y_cf    T x N_treated matrix of imputed counterfactual outcomes.
#' @return A list with `estimate`, `cohort_estimates`, and `cohort_fits`.
#' @noRd
.impute_cohort_fits <- function(pan, gap, Y_treat, Y_cf) {
  TT      <- nrow(gap)
  idx_tr  <- pan$idx_treat
  T_adopt <- pan$T_adopt
  cohorts <- sort(unique(T_adopt[idx_tr]))

  r <- lapply(cohorts, function(g) {
    k_g       <- which(T_adopt[idx_tr] == g)   # gap columns for this cohort
    T_pre_g   <- g - 1L
    T_post_g  <- TT - T_pre_g
    post_rows <- (T_pre_g + 1L):TT

    list(
      cohort    = g,
      n_treated = length(k_g),
      T_pre     = T_pre_g,
      T_post    = T_post_g,
      estimate  = mean(gap[post_rows, k_g, drop = FALSE], na.rm = TRUE),
      weight    = as.numeric(length(k_g)) * T_post_g,
      idx_tr    = idx_tr[k_g],
      idx_co    = pan$idx_control,
      Y_treat   = Y_treat[, k_g, drop = FALSE],
      Y_tr_hat  = Y_cf[, k_g, drop = FALSE]
    )
  })

  .aggregate_cohorts(r)
}

#' NULL-coalescing helper
#'
#' Returns `l` unless it is `NULL`, in which case `r` is returned. Shared
#' across the package (broom, plot, fit helpers).
#' @noRd
`%||%` <- function(l, r) if (is.null(l)) r else l

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
#'   Values are matched against the time index of the panel passed to
#'   [scm_fit()], and are expected to be pre-treatment: a window naming times
#'   the panel does not carry silently aggregates the overlap only, and one
#'   reaching past the treatment date aggregates values the treatment has
#'   already moved, so [scm_fit()] warns about both.
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
#' # Abadie, Diamond & Hainmueller (2010) California Prop 99 style: combine
#' # several covariates aggregated over different windows plus three outcome
#' # lags at specific years. The resulting list is passed to
#' # scm_fit(..., predictors = predictors).
#' predictors <- list(
#'   pred(c("lnincome", "retprice", "age15to24"), 1980:1988),
#'   pred("beer",    1984:1988),
#'   pred("cigsale", 1988),
#'   pred("cigsale", 1980),
#'   pred("cigsale", 1975)
#' )
#' predictors
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

#' Single-period outcome lag times in a predictor specification
#'
#' Returns the time points when every pred() entry names the outcome
#' variable at one single period, and NULL otherwise. That shape builds
#' predictor rows which are plain pre-treatment outcome rows (any `op`
#' collapses to the value itself on a single period), so it is the only
#' shape that can coincide with the outcomes-only fit.
#'
#' @param predictors  List of pred_spec objects (non-empty).
#' @param outcome_var Name of the outcome variable, or NULL when unknown
#'   (direct internal calls) -- returns NULL then.
#' @noRd
.outcome_lag_times <- function(predictors, outcome_var) {
  if (is.null(outcome_var)) return(NULL)
  ts <- vector("list", length(predictors))
  for (i in seq_along(predictors)) {
    p <- predictors[[i]]
    if (!inherits(p, "pred_spec")) return(NULL)
    if (length(p$vars) != 1L || p$vars != outcome_var) return(NULL)
    if (length(p$times) != 1L) return(NULL)
    ts[[i]] <- p$times
  }
  unlist(ts, use.names = FALSE)
}

#' Detect the canonical outcomes-only predictor specification
#'
#' TRUE when the entries are single-period outcome lags that jointly cover
#' each pre-treatment period exactly once. X0/X1 then equal the
#' pre-treatment outcome rows, so the fit is the outcomes-only fit and can
#' be routed to that path. Coverage is judged against `pre_times`, i.e.
#' against the periods `data` actually carries.
#'
#' @param predictors  List of pred_spec objects (non-empty).
#' @param outcome_var Name of the outcome variable, or NULL when unknown
#'   (direct internal calls) -- returns FALSE then.
#' @param pre_times   The pre-treatment time values.
#' @noRd
.is_outcome_only_spec <- function(predictors, outcome_var, pre_times) {
  ts <- .outcome_lag_times(predictors, outcome_var)
  !is.null(ts) && length(ts) == length(pre_times) &&
    anyDuplicated(ts) == 0L && all(ts %in% pre_times)
}

#' Flag pred() windows that do not aggregate what they read as
#'
#' Two silent traps that leave the fit well defined but not the model that
#' was written down: a window only partly present in the panel is averaged
#' over the intersection while keeping its full label, and a window reaching
#' past the treatment date averages values the treatment has already moved.
#' A window with no overlap at all is left alone -- it produces an all-NA
#' predictor row, which build_predictor_matrices() rejects outright.
#'
#' @param predictors List of pred_spec objects.
#' @param times      Panel time index, in order.
#' @param T_pre      Number of pre-treatment periods.
#' @noRd
.check_pred_windows <- function(predictors, times, T_pre) {
  post_times <- times[-seq_len(T_pre)]
  fmt <- function(x) {
    if (length(x) > 5L) {
      paste0(paste(x[seq_len(5L)], collapse = ", "), ", ...")
    } else {
      paste(x, collapse = ", ")
    }
  }
  absent <- character(0L)
  post   <- character(0L)
  for (p in predictors) {
    if (!inherits(p, "pred_spec")) next
    lab  <- paste(p$vars, collapse = ", ")
    miss <- setdiff(p$times, times)
    if (length(miss) > 0L && length(miss) < length(p$times)) {
      absent <- c(absent, sprintf("pred(%s, ...): %s", lab, fmt(miss)))
    }
    hit <- intersect(p$times, post_times)
    if (length(hit) > 0L) {
      post <- c(post, sprintf("pred(%s, ...): %s", lab, fmt(hit)))
    }
  }
  if (length(absent) > 0L) {
    warning(
      "pred() windows name times the panel does not contain, so those ",
      "predictor rows aggregate fewer periods than the window reads as -- ",
      paste(absent, collapse = "; "), ".",
      call. = FALSE
    )
  }
  if (length(post) > 0L) {
    warning(
      "pred() windows reach into the post-treatment periods, whose values ",
      "the treatment has already moved -- ", paste(post, collapse = "; "),
      ". SCM predictors are measured before treatment.",
      call. = FALSE
    )
  }
  invisible(NULL)
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
#' @param idx_tr     Integer indices of the treated units in `units`. Several
#'   units sharing one adoption date are one cohort, fitted through their
#'   average (Ben-Michael, Feller & Rothstein 2022, S.4.1), so each predictor
#'   row of `X1` is their average too.
#' @param tr_weights Weights over `idx_tr` summing to one, giving the average
#'   that `X1` targets. `NULL` weights the treated units equally. Must be the
#'   same weights the outcome path pools with, or the predictor rows and the
#'   outcome rows would describe different targets.
#' @param predictors List of `pred_spec` objects created by [pred()].
#' @return A list with:
#'   * `X0`: k x N_co numeric matrix (predictors x control units)
#'   * `X1`: numeric vector of length k (predictors for the treated average)
#'   * `pred_names`: character vector of length k with predictor labels
#' @noRd
build_predictor_matrices <- function(
  data,
  id_var,
  time_var,
  units,
  idx_co,
  idx_tr,
  predictors,
  tr_weights = NULL
) {
  co_units <- units[idx_co]
  w_tr <- if (is.null(tr_weights)) rep(1 / length(idx_tr), length(idx_tr))
          else tr_weights

  # One grouped aggregation per predictor row instead of one full-data scan
  # per unit: id_pos maps each data row to its unit's position in `units`.
  id_pos <- match(data[[id_var]], units)

  agg_unit <- function(var, rows_in_window, op) {
    fn  <- match.fun(op)
    grp <- split(data[[var]][rows_in_window], id_pos[rows_in_window])
    out <- rep(NA_real_, length(units))
    out[as.integer(names(grp))] <- vapply(grp, fn, numeric(1L), na.rm = TRUE)
    if (is.character(units)) names(out) <- units
    out
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
    rows_in_window <- data[[time_var]] %in% p$times
    for (var in p$vars) {
      if (!var %in% names(data)) {
        stop(sprintf("Variable '%s' not found in data.", var), call. = FALSE)
      }
      vals <- agg_unit(var, rows_in_window, p$op)
      rows_X0 <- c(rows_X0, list(vals[idx_co]))
      rows_X1 <- c(rows_X1, list(sum(vals[idx_tr] * w_tr)))
      pred_names <- c(pred_names, pred_label(var, p$times))
    }
  }

  X0 <- do.call(rbind, rows_X0) # k x N_co
  X1 <- unlist(rows_X1) # length k

  colnames(X0) <- as.character(co_units)
  rownames(X0) <- pred_names

  # Non-finite predictor values (e.g. a pred() window with no data for some
  # unit) would silently corrupt the downstream QP -- fail loudly instead.
  bad <- (rowSums(!is.finite(X0)) > 0L) | !is.finite(X1)
  if (any(bad)) {
    stop(
      "Predictor rows contain missing or non-finite values: ",
      paste(pred_names[bad], collapse = ", "),
      ". Check the pred() time windows against the available data.",
      call. = FALSE
    )
  }

  list(X0 = X0, X1 = X1, pred_names = pred_names)
}

#' Build the covariate balance matrix for partially pooled SCM
#'
#' Ben-Michael, Feller & Rothstein (2022, S.5.2) add auxiliary covariates to
#' the partially pooled problem by *balancing* them alongside the lagged
#' outcomes, not by partialling them out. Each unit contributes one vector
#' X_i: the paper takes covariates measured before treatment, and says a
#' time-varying one enters through "the values at time periods before the
#' first treatment time T1", so each column here is that covariate's mean
#' over the periods preceding the earliest adoption.
#'
#' The covariates are standardised across units and then multiplied by
#' `xi_sd`. Squaring inside the imbalance norm turns that factor into the
#' hyperparameter xi of their Equation (11), which is why the paper can
#' describe xi as a variance while it is applied here as a standard
#' deviation.
#'
#' @param data        Long-format data frame.
#' @param id_var,time_var Names of the unit and time columns.
#' @param covariates  Character vector of covariate column names.
#' @param units,times Panel unit and time indices, in matrix order.
#' @param T_pre       Number of periods before the earliest adoption.
#' @param idx_ref     Column indices of the never-treated units, whose mean and
#'   standard deviation standardise each covariate.
#' @param xi_sd       Scale applied after standardising (see above).
#' @return An N x d numeric matrix, units in rows.
#' @noRd
.scm_cov_balance_matrix <- function(data, id_var, time_var, covariates,
                                    units, times, T_pre, idx_ref, xi_sd) {
  pre_times <- times[seq_len(T_pre)]
  rows      <- data[[time_var]] %in% pre_times
  id_pos    <- match(data[[id_var]], units)

  Z <- vapply(covariates, function(v) {
    if (!v %in% names(data))
      stop(sprintf("Covariate '%s' not found in data.", v), call. = FALSE)
    grp <- split(data[[v]][rows], id_pos[rows])
    out <- rep(NA_real_, length(units))
    out[as.integer(names(grp))] <- vapply(grp, mean, numeric(1L), na.rm = TRUE)
    out
  }, numeric(length(units)))
  Z <- matrix(Z, nrow = length(units),
              dimnames = list(as.character(units), covariates))

  bad <- covariates[colSums(!is.finite(Z)) > 0L]
  if (length(bad))
    stop("Covariate(s) with missing pre-treatment values for some units: ",
         paste(bad, collapse = ", "),
         ". Balancing needs one value per unit.", call. = FALSE)

  # Standardised on the never-treated units, as `augsynth` does: they are the
  # pool the weights are drawn from, and their moments do not move with which
  # units of the panel happen to be treated.
  ctr <- colMeans(Z[idx_ref, , drop = FALSE])
  sds <- apply(Z[idx_ref, , drop = FALSE], 2L, stats::sd)
  flat <- covariates[!is.finite(sds) | sds < 1e-12]
  if (length(flat))
    stop("Covariate(s) with no variation across the donor units: ",
         paste(flat, collapse = ", "),
         ". A constant cannot be balanced.", call. = FALSE)

  sweep(sweep(Z, 2L, ctr, "-"), 2L, sds / xi_sd, "/")
}

#' Build a time-varying covariate array for GSC
#'
#' Constructs a T x N x p 3D R array suitable for passing to gsc_ife_cpp as an
#' arma::cube. `arr[t, i, j]` is the value of the j-th covariate for unit i at
#' time t.
#'
#' @param data            Long-format data frame.
#' @param id_var          Name of the unit identifier column.
#' @param time_var        Name of the time identifier column.
#' @param covariate_names Character vector of covariate column names (length p).
#' @param units           Character vector of unit IDs (length N).
#' @param times           Vector of time values (length T).
#' @return A T x N x p numeric array.
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
  # Vectorised fill: one matrix-index assignment per covariate instead of a
  # full-data scan per unit. Rows whose id is not in `units` are skipped,
  # matching the per-unit subsetting behaviour.
  id_pos <- match(data[[id_var]], units)
  t_pos  <- match(data[[time_var]], times)
  keep   <- !is.na(id_pos)
  for (j in seq_len(p)) {
    var <- covariate_names[j]
    if (!var %in% names(data)) {
      stop(sprintf("covariate '%s' not found in data.", var), call. = FALSE)
    }
    arr[cbind(t_pos[keep], id_pos[keep], j)] <- data[[var]][keep]
  }
  arr
}

#' Panel Data Helper: Reshape long-format multi-arm data to tensor structure
#'
#' Extends `panel_to_matrices()` for the multi-arm Synthetic Interventions
#' setting (Agarwal et al. 2025). Each unit belongs to exactly one treatment
#' arm (d = 0 for control, d = 1,...,K for treatment arms). Before the
#' treatment date, all treatment-arm units have d = 0; at the treatment date
#' they switch to their assigned arm value.
#'
#' @param y    Numeric outcome vector (long format)
#' @param d    Treatment arm indicator (integer, 0 = control, 1,...,K = arms)
#' @param id   Unit identifier (long format)
#' @param time Time identifier (long format)
#' @return All fields returned by `panel_to_matrices()` plus:
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
    stop("panel_to_tensor: arm 0 (control) not found.", call. = FALSE)
  if (any(arm_levels < 0L))
    stop("panel_to_tensor: arm values must be non-negative integers.", call. = FALSE)

  idx_by_arm <- setNames(
    lapply(arm_levels, function(a) which(arm_of_unit == a)),
    as.character(arm_levels)
  )

  c(pan, list(arm_levels = arm_levels, idx_by_arm = idx_by_arm))
}

# Validate and align aggregation weights over the treated units. Abadie (2021)
# S.8 notes that with aggregate units a population- or GDP-weighted average
# effect is often the quantity of interest; `weights` supplies that mass.
# Returned named over the treated units' labels, or NULL for equal weights.
.scm_treated_weights <- function(weights, pan) {
  if (is.null(weights)) return(NULL)
  labs <- as.character(pan$units[pan$idx_treat])
  if (!is.numeric(weights) || anyNA(weights))
    stop("`weights` must be a numeric vector without NA values, named over ",
         "the treated units.", call. = FALSE)
  if (is.null(names(weights)))
    stop("`weights` must be named over the treated units (",
         paste(labs, collapse = ", "), ").", call. = FALSE)
  missing_lab <- setdiff(labs, names(weights))
  if (length(missing_lab))
    stop("`weights` is missing the treated unit(s): ",
         paste(missing_lab, collapse = ", "), ".", call. = FALSE)
  w <- weights[labs]
  if (any(w < 0))
    stop("`weights` must be non-negative.", call. = FALSE)
  if (sum(w) <= 0)
    stop("`weights` must have a positive sum.", call. = FALSE)
  # Rescaled to average one, so a block's mass is the sum of its units' and
  # reduces to the cohort size n_g of Ben-Michael, Feller & Rothstein (2022,
  # Appendix A.2) when the weights are equal. Only the ratios matter --
  # everything downstream normalises -- so the rescaling changes no estimate.
  w / mean(w)
}

#' Mass and target series of one staggered block
#'
#' A block is a cohort (or a single treated unit under `pool_cohorts = FALSE`).
#' Its mass is what it counts for, both in the pooled/separate imbalance
#' measures of Ben-Michael, Feller & Rothstein (2022, Appendix A.2), where it
#' is the cohort size \eqn{n_g}, and in the Clarke et al. (2024) aggregation
#' weight \eqn{n_g T_{post,g}}. Its target is the correspondingly weighted
#' average of the treated series, so the weights are fitted towards the same
#' quantity the aggregation reports.
#'
#' `unit_weights` (Abadie 2021, S.8) replaces the unit count with the mass the
#' caller gives each treated unit. Left `NULL` every unit counts once.
#'
#' @param Y            Outcome matrix (T x N).
#' @param idx_tr_g     Column indices of the block's treated units.
#' @param unit_labels  Unit labels for every panel column (`pan$units`).
#' @param unit_weights Optional named numeric vector over treated unit labels.
#' @return A list with `mass` (scalar), `masses` (one per treated unit, summing
#'   to `mass`) and `y` (length-T target series).
#' @noRd
.scm_block_target <- function(Y, idx_tr_g, unit_labels, unit_weights = NULL) {
  if (is.null(unit_weights))
    return(list(mass   = as.numeric(length(idx_tr_g)),
                masses = rep(1, length(idx_tr_g)),
                y      = rowMeans(Y[, idx_tr_g, drop = FALSE])))
  w <- unname(unit_weights[as.character(unit_labels[idx_tr_g])])
  m <- sum(w)
  # An all-zero block carries no mass in the ATT; the caller's skip path drops
  # it rather than leaving the solver an objective that does not identify it.
  if (m <= 0)
    return(list(mass = 0, masses = w,
                y = rowMeans(Y[, idx_tr_g, drop = FALSE])))
  list(mass = m, masses = w,
       y = drop(Y[, idx_tr_g, drop = FALSE] %*% (w / m)))
}
