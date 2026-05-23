#' Fit Matrix Completion (Athey et al. 2021)
#'
#' Treats treated post-intervention outcomes as missing and recovers them
#' via nuclear-norm-regularised matrix completion (Soft-Impute algorithm).
#'
#' @param y      Outcome vector (long format)
#' @param d      Treatment indicator (long format)
#' @param id     Unit identifier (long format)
#' @param time   Time identifier (long format)
#' @param lambda Nuclear norm penalty. If NULL, defaults to 1% of the
#'               spectral norm of the observed data matrix.
#' @return A list of class `coresynth`.
#' @noRd
fit_mc_cpp <- function(y, d, id, time, lambda = NULL, ...) {
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

  # Replace NAs in Y with 0 (unobserved positions will be ignored via mask)
  Y_obs <- Y; Y_obs[is.na(Y_obs)] <- 0.0

  # Default lambda: small fraction of spectral norm
  if(is.null(lambda)) {
    sv <- svd(Y_obs, nu = 0, nv = 0)$d
    lambda <- 0.01 * sv[1]
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

  list(
    method       = "mc",
    lambda       = lambda,
    L_hat        = L_hat,
    Y_hat        = Y_hat,
    Y_treat      = Y_treat,
    gap          = gap,
    times        = pan$times,
    T_pre        = T_pre,
    idx_co       = idx_co,
    idx_tr       = idx_tr,
    unit_weights = NULL,
    estimate     = att
  )
}
