# ── Unit-resampling inference shared by the GSC, SI and SDID routines ────────
#
# Each estimator is re-run on draws of units. `estimate_fn(ctx, co, tr)`
# returns the effect for the never-treated columns `co` and the treated
# columns `tr` of `ctx$Y` -- indices may repeat, as in a bootstrap draw -- as
# a list with `estimate` and the per-cohort `cohort`, `cohort_estimate` and
# `cohort_weight` (N_treated,g * T_post,g), or NULL when a refit fails.
# `jack_fn` has the same contract and serves the jackknife, for estimators
# whose jackknife keeps the fitted weights instead of re-estimating them.
#
# The treated units are resampled along with the controls. Holding them fixed
# leaves out their own noise -- their post-treatment shocks and the error in
# what their pre-treatment periods identify -- which is most of the variance
# when few units are treated: with a single treated unit a 95% interval built
# from the controls alone covered the truth 30% of the time. By default the
# bootstrap draws the treated units within `ctx$strata` (the arms of a
# multi-arm fit) and the controls separately. With `ctx$pooled` it draws N
# units from the whole panel and redraws a sample without a treated unit or
# with fewer than two controls, which is Algorithm 2 of Arkhangelsky et al.
# (2021). Either way the cohorts are pooled, so cohort sizes vary across
# draws as in Clarke et al. (2024, Algorithm 2), and the controls drawn serve
# every cohort, which keeps the correlation that shared donors induce between
# the cohort effects. The jackknife drops one unit at a time, treated or
# control.
#
# `label` names the caller in messages and `one_treated_hint` says where to go
# when a stratum has a single treated unit, whose noise no resampling of units
# can reach.
.unit_inference <- function(fit, ctx, estimate_fn, method, n_boot, level,
                            alternative, label, one_treated_hint,
                            jack_fn = estimate_fn) {
  alpha     <- 1 - level
  tau_hat   <- fit$estimate
  staggered <- inherits(fit, "coresynth_staggered")
  co        <- ctx$co
  tr        <- ctx$tr
  strata    <- ctx$strata %||% rep(1L, length(tr))

  if (any(table(strata) < 2L))
    stop(sprintf(paste0(
      "%s resamples the treated units along with the controls, so it needs ",
      "at least two treated units%s. With one, that unit's own noise never ",
      "varies across replicates and the interval comes out far too narrow. ",
      "%s"),
      label, if (length(unique(strata)) > 1L) " in every arm" else "",
      one_treated_hint), call. = FALSE)

  # Clarke et al. (2024, Algorithm 3): dropping a cohort's only treated unit
  # leaves its effect undefined, so the jackknife needs two in every cohort.
  if (staggered && method != "bootstrap") {
    n_g <- table(ctx$T_adopt[tr])
    if (any(n_g < 2L))
      stop(sprintf(paste0(
        "%s: the jackknife drops each treated unit in turn, so every adoption ",
        "cohort needs at least two treated units (Clarke et al. 2024, ",
        "Algorithm 3); %d cohort(s) have one. Use method = \"bootstrap\", ",
        "which draws the treated units from all cohorts together."),
        label, sum(n_g < 2L)), call. = FALSE)
  }

  n_unconverged <- 0L
  counted <- function(fn) function(co_b, tr_b) {
    res <- fn(ctx, co_b, tr_b)
    if (is.null(res)) return(NULL)
    n_unconverged <<- n_unconverged + (res[["unconverged"]] %||% 0L)
    res
  }
  refit     <- counted(estimate_fn)
  jack_once <- counted(jack_fn)
  value <- function(res) if (is.null(res)) NA_real_ else res$estimate

  boot_ests <- NULL
  if (method == "bootstrap") {
    if (isTRUE(ctx$pooled)) {
      pool  <- c(co, tr)
      is_tr <- rep(c(FALSE, TRUE), c(length(co), length(tr)))
      draw  <- function() {
        for (attempt in seq_len(1000L)) {
          pick <- sample.int(length(pool), length(pool), replace = TRUE)
          n_tr <- sum(is_tr[pick])
          if (n_tr >= 1L && length(pick) - n_tr >= 2L)
            return(list(co = pool[pick[!is_tr[pick]]],
                        tr = pool[pick[is_tr[pick]]]))
        }
        NULL
      }
    } else {
      by_stratum <- split(tr, strata)
      draw <- function() list(
        co = co[sample.int(length(co), length(co), replace = TRUE)],
        tr = unlist(lapply(by_stratum, function(s)
          s[sample.int(length(s), length(s), replace = TRUE)]),
          use.names = FALSE))
    }
    boot_ests <- replicate(n_boot, {
      d <- draw()
      if (is.null(d)) NA_real_ else value(refit(d$co, d$tr))
    })
    n_fail <- sum(is.na(boot_ests))
    if (n_fail > 0.1 * n_boot)
      warning(sprintf("%.0f%% of bootstrap replications failed.",
                      100 * n_fail / n_boot), call. = FALSE)
    valid    <- boot_ests[!is.na(boot_ests)]
    se       <- stats::sd(valid)
    ci_lower <- unname(quantile(valid, alpha / 2))
    ci_upper <- unname(quantile(valid, 1 - alpha / 2))
  } else {
    # One refit per dropped unit serves both jackknives: the global one reads
    # the overall effect, the per-cohort one each cohort's.
    units   <- unique(c(co, tr))
    dropped <- lapply(units, function(u) jack_once(co[co != u], tr[tr != u]))
    jack_var <- function(th) {
      th <- th[!is.na(th)]
      n  <- length(th)
      if (n < 2L) return(NA_real_)
      (n - 1) / n * sum((th - mean(th))^2)
    }

    if (!staggered || method == "jackknife_global") {
      v <- jack_var(vapply(dropped, value, numeric(1L)))
    } else {
      # Per-cohort jackknife with delta-method aggregation: each cohort's
      # variance from dropping its own donors and treated units, combined
      # with the squared cohort weights. It ignores the covariance that
      # shared donors induce between cohorts; jackknife_global keeps it.
      full <- jack_fn(ctx, co, tr)
      if (is.null(full)) stop(label, ": the fit could not be reproduced.",
                              call. = FALSE)
      adopt <- ctx$T_adopt[tr]
      var_g <- vapply(seq_along(full$cohort), function(k) {
        g      <- full$cohort[k]
        donors <- if (identical(ctx$control_group, "never_treated")) co else
          c(co, tr[adopt > g])
        own <- units %in% unique(c(donors, tr[adopt == g]))
        th  <- vapply(dropped[own], function(res) {
          if (is.null(res)) return(NA_real_)
          v <- res$cohort_estimate[res$cohort == g]
          if (length(v) == 1L) v else NA_real_
        }, numeric(1L))
        vg <- jack_var(th)
        if (is.na(vg)) 0 else vg
      }, numeric(1L))
      w <- full$cohort_weight / sum(full$cohort_weight)
      v <- sum(w^2 * var_g)
    }
    se       <- sqrt(v)
    ci_lower <- tau_hat - qnorm(1 - alpha / 2) * se
    ci_upper <- tau_hat + qnorm(1 - alpha / 2) * se
  }

  if (n_unconverged > 0L)
    warning(sprintf(
      "%s: %d refit(s) stopped before their iterative fit converged.",
      label, n_unconverged), call. = FALSE)

  z <- tau_hat / max(se, .Machine$double.eps)
  structure(list(
    estimate    = tau_hat,
    se          = se,
    p_value     = switch(alternative,
                         two.sided = 2 * pnorm(-abs(z)),
                         greater   = pnorm(-z),
                         less      = pnorm(z)),
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    method      = method,
    staggered   = staggered,
    n_controls  = length(unique(co)),
    alternative = alternative,
    boot_ests   = boot_ests
  ), class = "coresynth_inference")
}

# Adoption row of every column of a staggered fit's panel, NA for the
# never-treated columns, rebuilt from the cohort fits.
.adopt_from_cohorts <- function(N, cohort_fits) {
  T_adopt <- rep(NA_integer_, N)
  for (cf in cohort_fits) T_adopt[cf$idx_tr] <- as.integer(cf$cohort)
  T_adopt
}

# The control group a staggered fit used. Fits that do not record it are
# "never_treated" exactly when every cohort's donors are the never-treated
# units.
.fit_control_group <- function(fit) {
  fit[["control_group"]] %||% {
    same <- vapply(fit$cohort_fits, function(cf)
      setequal(cf$idx_co, fit$idx_co), logical(1L))
    if (all(same)) "never_treated" else "clean"
  }
}
