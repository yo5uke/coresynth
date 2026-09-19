#' Fit Matrix Completion (Athey et al. 2021)
#'
#' Treats treated post-intervention outcomes as missing and recovers them
#' via nuclear-norm-regularised matrix completion (Soft-Impute algorithm).
#'
#' @param y      Outcome vector (long format)
#' @param d      Treatment indicator (long format)
#' @param id     Unit identifier (long format)
#' @param time   Time identifier (long format)
#' @param lambda Nuclear norm penalty. `NULL` (default) selects it by the
#'   K-fold cross-validation of Athey et al. (2021, S.4).
#' @param cv_folds Number of cross-validation folds used when `lambda = NULL`.
#' @param cv_nlambda Length of the candidate penalty path.
#' @param cv_seed Seed for the fold draw, so the selected penalty is
#'   reproducible without consuming the caller's RNG state.
#' @return A list of class `coresynth`.
#' @noRd
# ── Internal: Athey et al. (2021) S.4 -- cross-validate the nuclear norm penalty
# K random subsets O_k of the observed cells are drawn with |O_k| chosen so the
# observed fraction inside a fold matches that of the full panel,
# |O_k|/|O| = |O|/(NT). For each fold the whole penalty path is solved and
# scored on the cells held out (O \ O_k), and the lambda with the smallest
# average squared error across folds wins. The paper's warm start along the
# path is a speed-up rather than part of the estimator, and soft_impute_cpp
# always starts from zero, so it is not used here.
#
# The path runs from the spectral norm of the (masked) panel, which shrinks
# everything to zero, down towards zero, which fits the observed cells exactly.
.mc_cv_lambda <- function(Y_obs, O, n_folds = 5L, n_lambda = 20L, seed = 1L) {
  obs_idx <- which(O > 0)
  n_obs   <- length(obs_idx)
  n_cell  <- length(O)
  if (n_obs < 2L) return(NULL)

  lam_max <- svd(Y_obs * O, nu = 0, nv = 0)$d[1L]
  if (!is.finite(lam_max) || lam_max <= 0) return(NULL)
  # Geometric path down to a small fraction of lambda_max; the paper's
  # lambda_L = 0 is the interpolating fit, which cross-validation never wants
  # and which costs the most to compute.
  lambdas <- lam_max * 10^seq(0, -3, length.out = as.integer(n_lambda))

  fold_size <- max(1L, floor(n_obs^2 / n_cell))
  if (fold_size >= n_obs) return(NULL)

  # Own RNG stream: the penalty must not depend on, or disturb, the caller's
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  on.exit({
    if (is.null(old_seed))
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    else assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)

  err <- matrix(NA_real_, nrow = as.integer(n_folds), ncol = length(lambdas))
  for (k in seq_len(as.integer(n_folds))) {
    keep_idx <- sample(obs_idx, fold_size)
    O_k <- matrix(0.0, nrow(O), ncol(O))
    O_k[keep_idx] <- 1.0
    held <- setdiff(obs_idx, keep_idx)
    if (length(held) == 0L) next

    for (i in seq_along(lambdas)) {
      L_i <- tryCatch(soft_impute_cpp(Y_obs, O_k, lambdas[i]),
                      error = function(e) NULL)
      if (is.null(L_i)) next
      err[k, i] <- mean((Y_obs[held] - L_i[held])^2)
    }
  }

  score <- colMeans(err, na.rm = TRUE)
  if (all(is.na(score))) return(NULL)
  lambdas[which.min(score)]
}

fit_mc_cpp <- function(y, d, id, time, lambda = NULL,
                       cv_folds = 5L, cv_nlambda = 20L, cv_seed = 1L, ...) {
  pan <- panel_to_matrices(y, d, id, time)
  Y   <- pan$Y; T_pre <- pan$T_pre
  idx_tr <- pan$idx_treat; idx_co <- pan$idx_control

  TT <- nrow(Y); N <- ncol(Y)

  # Per-unit observation mask (0 = treated post-adoption, 1 = observed)
  O <- matrix(1.0, nrow = TT, ncol = N)
  for (j in idx_tr) {
    t0 <- pan$T_adopt[j]
    if (!is.na(t0) && t0 <= TT) O[t0:TT, j] <- 0.0
  }

  # Missing cells (unbalanced panel or NA outcomes) are unobserved, like
  # treated post-adoption cells: mask them out, then zero-fill so the masked
  # positions carry no weight in the C++ solver.
  O[is.na(Y)] <- 0.0
  Y_obs <- Y; Y_obs[is.na(Y_obs)] <- 0.0

  # Athey et al. (2021) select lambda by cross-validation. Fall back to the
  # old fixed fraction of the spectral norm only when the panel is too small
  # for a fold to hold anything out.
  lambda_cv <- FALSE
  if (is.null(lambda)) {
    lambda <- .mc_cv_lambda(Y_obs, O, n_folds = cv_folds,
                            n_lambda = cv_nlambda, seed = cv_seed)
    lambda_cv <- !is.null(lambda)
    if (!lambda_cv) lambda <- 0.01 * svd(Y_obs, nu = 0, nv = 0)$d[1L]
  }

  L_hat <- soft_impute_cpp(Y_obs, O, lambda)

  Y_treat <- Y[, idx_tr, drop = FALSE]
  Y_hat   <- L_hat[, idx_tr, drop = FALSE]

  gap <- Y_treat - Y_hat
  # ATT: average gap over each unit's own post-treatment period
  gap_post_vals <- vapply(seq_along(idx_tr), function(k) {
    j  <- idx_tr[k]
    t0 <- pan$T_adopt[j]
    if (is.na(t0) || t0 > TT) return(NA_real_)
    mean(gap[t0:TT, k, drop = FALSE], na.rm = TRUE)
  }, numeric(1L))
  att <- mean(gap_post_vals, na.rm = TRUE)

  # Staggered adoption needs no cohort-by-cohort re-fitting here: the mask O
  # already blanks each treated unit from its own adoption date, so the single
  # factorisation is the staggered fit. Only the aggregation changes -- the
  # cohort ATTs combine with the N_treated,g * T_post,g weights of Clarke et
  # al. (2024), as in every other estimator, instead of weighting each treated
  # unit equally.
  stag <- if (!pan$is_sharp) .impute_cohort_fits(pan, gap, Y_treat, Y_hat)

  c(list(
    method       = "mc",
    lambda       = lambda,
    lambda_cv    = lambda_cv,
    L_hat        = L_hat,
    Y_hat        = Y_hat,
    Y_treat      = Y_treat,
    gap          = gap,
    times        = pan$times,
    T_pre        = T_pre,
    idx_co       = idx_co,
    idx_tr       = idx_tr,
    Y_co_all     = Y[, idx_co, drop = FALSE],  # T x N_co -- needed by conformal_inference()
    unit_weights = NULL,
    estimate     = stag$estimate %||% att
  ), if (!is.null(stag)) list(
    staggered        = TRUE,
    cohort_estimates = stag$cohort_estimates,
    cohort_fits      = stag$cohort_fits
  ))
}
