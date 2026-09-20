# -- Abadie & L'Hour (2021) S.3.2 rank-based permutation inference -----------

#' Rank-Based Permutation Test for Several Treated Units
#'
#' Permutation inference on the aggregate effect when more than one unit is
#' treated, following Abadie & L'Hour (2021, S.3.2). Each treated unit
#' contributes its own test statistic; the identities of the treated units are
#' then reassigned at random among all units, and the sum of the ranks of the
#' real treated units' statistics is compared with its permutation
#' distribution.
#'
#' Ranking is what makes the test usable here. [mspe_ratio_pval()] compares one
#' treated unit against placebo runs on single donors, which breaks down as
#' soon as the treated series is an average: the average is less noisy than any
#' single donor, so its pre/post MSPE ratio is inflated and the test rejects
#' far too often. Ranks sidestep that, because every statistic entering the
#' comparison, real or permuted, is computed on a single unit.
#'
#' The fit must expose per-unit effects, so a staggered fit has to be made with
#' `pool_cohorts = FALSE`. A sharp fit with several treated units pools them by
#' design, and its per-unit statistics are computed here from the stored
#' per-unit series.
#'
#' @param fit A `coresynth` object from [scm_fit()] with `method = "scm"` and
#'   more than one treated unit.
#' @param B Number of treatment reassignments.
#' @param statistic Per-unit statistic to rank. `"mspe_ratio"` (default) is the
#'   post/pre MSPE ratio of Abadie, Diamond & Hainmueller (2010), which adjusts
#'   for how well each unit was fitted; `"effect"` is the mean post-treatment
#'   gap and `"abs_effect"` its absolute value.
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param seed Seed for the reassignment draw, so the p-value is reproducible
#'   without consuming the caller's RNG state.
#' @return An object of class `c("scm_rank_test", "coresynth_inference")` with
#'   `p_value`, `statistic` (the observed rank sum), `rank_sums` (its
#'   permutation distribution), `unit_statistics`, `n_treated`, `B` and
#'   `alternative`.
#' @references Abadie, A., & L'Hour, J. (2021). A penalized synthetic control
#'   estimator for disaggregated data. *JASA*, 116(536), 1817-1834.
#' @seealso [mspe_ratio_pval()] for the single-treated-unit test.
#' @examples
#' set.seed(1)
#' panel <- expand.grid(unit = 1:12, year = 1:20)
#' panel$treated <- as.integer(panel$unit <= 3 & panel$year > 12)
#' panel$gdp <- panel$unit + 0.4 * panel$year +
#'   rnorm(nrow(panel)) + 2 * panel$treated
#' fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = "scm")
#' scm_rank_pval(fit, B = 50)
#' @export
scm_rank_pval <- function(fit, B = 500L,
                          statistic = c("mspe_ratio", "effect", "abs_effect"),
                          alternative = c("two.sided", "greater", "less"),
                          seed = 1L) {
  statistic   <- match.arg(statistic)
  alternative <- match.arg(alternative)
  if (!inherits(fit, "coresynth_scm"))
    stop("scm_rank_pval() requires a coresynth object with method = 'scm'.",
         call. = FALSE)

  spec <- .scm_rank_spec(fit)
  n1   <- length(spec$idx_tr)
  if (n1 < 2L)
    stop("scm_rank_pval() is for several treated units; with one treated ",
         "unit use mspe_ratio_pval().", call. = FALSE)

  Y       <- spec$Y
  T_adopt <- spec$T_adopt
  pool    <- c(spec$idx_tr, spec$idx_co)

  # Statistic for one unit treated from row t0, against `donors`. Uniform lag
  # weights, as in Abadie & L'Hour.
  unit_stat <- function(j, t0, donors) {
    if (is.na(t0) || t0 < 3L || t0 > nrow(Y)) return(NA_real_)
    pre  <- seq_len(t0 - 1L)
    post <- t0:nrow(Y)
    X0   <- Y[pre, donors, drop = FALSE]
    if (ncol(X0) < 2L) return(NA_real_)
    w <- tryCatch(
      as.numeric(solve_simplex_qp(crossprod(X0), crossprod(X0, Y[pre, j]))),
      error = function(e) NULL)
    if (is.null(w)) return(NA_real_)
    gap <- Y[, j] - drop(Y[, donors, drop = FALSE] %*% w)
    switch(statistic,
      effect     = mean(gap[post]),
      abs_effect = abs(mean(gap[post])),
      mspe_ratio = {
        pre_mspe <- mean(gap[pre]^2)
        if (pre_mspe <= .Machine$double.eps) NA_real_ else
          mean(gap[post]^2) / pre_mspe
      })
  }

  # Statistics for one assignment. The units not coded treated form that
  # assignment's donor pool, so every draw sees a pool of the same size.
  assignment_stats <- function(treated, adopt) {
    donors <- setdiff(pool, treated)
    vapply(seq_along(treated),
           function(k) unit_stat(treated[k], adopt[k], donors), numeric(1L))
  }

  adopt_obs <- T_adopt[spec$idx_tr]
  obs       <- assignment_stats(spec$idx_tr, adopt_obs)

  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  on.exit({
    if (is.null(old_seed))
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    else assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)

  # Adoption dates stay attached to the slots rather than to the units, so a
  # permuted assignment has the same event-time structure as the real one.
  perm <- lapply(seq_len(B), function(b)
    assignment_stats(sample(pool, n1), adopt_obs))

  all_stats <- c(list(obs), perm)
  flat      <- unlist(all_stats)
  if (all(is.na(flat)))
    stop("scm_rank_pval(): every permutation statistic was undefined; the ",
         "donor pool is probably too small.", call. = FALSE)

  # Rank across all n1 * (B + 1) unit-level statistics at once, then sum the
  # ranks within each assignment (Abadie & L'Hour S.3.2, step 3).
  ranks   <- rank(flat, na.last = "keep", ties.method = "average")
  by_draw <- split(ranks, rep(seq_along(all_stats), lengths(all_stats)))
  SR      <- vapply(by_draw, function(r) sum(r, na.rm = TRUE), numeric(1L))
  SR_obs  <- unname(SR[1L])
  SR_perm <- unname(SR[-1L])

  p_value <- switch(alternative,
    greater   = (1 + sum(SR_perm >= SR_obs)) / (B + 1),
    less      = (1 + sum(SR_perm <= SR_obs)) / (B + 1),
    two.sided = {
      centre <- mean(SR_perm)
      (1 + sum(abs(SR_perm - centre) >= abs(SR_obs - centre))) / (B + 1)
    })

  structure(list(
    p_value         = p_value,
    statistic       = SR_obs,
    rank_sums       = SR_perm,
    unit_statistics = stats::setNames(obs,
                        as.character(spec$units[spec$idx_tr])),
    stat_type       = statistic,
    n_treated       = n1,
    n_controls      = length(spec$idx_co),
    B               = as.integer(B),
    estimate        = fit$estimate,
    staggered       = inherits(fit, "coresynth_staggered"),
    alternative     = alternative
  ), class = c("scm_rank_test", "coresynth_inference"))
}

# Recover the panel and treatment pattern the rank test needs, from either a
# sharp fit pooling several treated units or a unit-level staggered fit.
.scm_rank_spec <- function(fit) {
  if (inherits(fit, "coresynth_staggered")) {
    cf <- fit$cohort_fits
    if (is.null(cf) || any(vapply(cf, function(z) length(z$idx_tr) != 1L,
                                  logical(1L))))
      stop("scm_rank_pval() needs one synthetic control per treated unit. ",
           "Re-fit with pool_cohorts = FALSE.", call. = FALSE)
    if (is.null(fit$Y_all))
      stop("staggered SCM fit is missing Y_all; re-run scm_fit() with the ",
           "current version.", call. = FALSE)
    idx_tr  <- vapply(cf, function(z) as.integer(z$idx_tr), integer(1L))
    idx_co  <- setdiff(sort(unique(unlist(lapply(cf, `[[`, "idx_co")))), idx_tr)
    T_adopt <- rep(NA_integer_, ncol(fit$Y_all))
    for (z in cf) T_adopt[z$idx_tr] <- as.integer(z$T_pre + 1L)
    return(list(Y = fit$Y_all, idx_tr = idx_tr, idx_co = idx_co,
                T_adopt = T_adopt,
                units = colnames(fit$Y_all) %||% seq_len(ncol(fit$Y_all))))
  }
  if (is.null(fit$Y_treat_mat) || is.null(fit$Y_co_pre))
    stop("fit is missing the per-unit treated series; re-run scm_fit() with ",
         "the current version.", call. = FALSE)
  Y_co <- rbind(fit$Y_co_pre, fit$Y_co_post)
  Y    <- cbind(fit$Y_treat_mat, Y_co)
  n1   <- ncol(fit$Y_treat_mat)
  list(Y = Y, idx_tr = seq_len(n1), idx_co = n1 + seq_len(ncol(Y_co)),
       T_adopt = rep(as.integer(fit$T_pre + 1L), ncol(Y)),
       units = colnames(Y) %||% seq_len(ncol(Y)))
}

#' @export
print.scm_rank_test <- function(x, ...) {
  cat("=== SCM rank-based permutation test ===\n")
  cat("Treated units :", x$n_treated, "| donors:", x$n_controls, "\n")
  cat("Statistic     :", x$stat_type, "\n")
  cat("Reassignments :", x$B, "\n")
  cat("Rank sum      :", round(x$statistic, 2), "\n")
  cat("p-value       :", format.pval(x$p_value, digits = 3),
      paste0("(", x$alternative, ")"), "\n")
  invisible(x)
}
