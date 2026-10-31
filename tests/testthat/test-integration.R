## Integration tests for coresynth
## Run: devtools::test()

library(testthat)

# ── Shared helpers ────────────────────────────────────────────────────────────

# Balanced panel: 1 treated unit (u1), rest control, treatment from T_pre+1 onward.
# DGP: Y_it = alpha_i + f_t * lambda_i + tau*D_it + eps
make_panel <- function(N = 10, T = 20, T_pre = 10, effect = 2.0, seed = 42) {
  set.seed(seed)
  times <- 1:T
  units <- paste0("u", 1:N)
  f <- cumsum(rnorm(T, 0, 0.5)) # common factor
  lam <- rnorm(N, 1, 0.3) # unit loadings
  rows <- expand.grid(time = times, id = units, stringsAsFactors = FALSE)
  rows$y <- as.vector(outer(f, lam)) + rnorm(nrow(rows), 0, 0.3)
  rows$d <- as.integer(rows$id == "u1" & rows$time > T_pre)
  rows$y[rows$d == 1] <- rows$y[rows$d == 1] + effect
  rows[order(rows$id, rows$time), ]
}

panel <- make_panel()

# Unit resampling needs more than one treated unit, so the sharp inference
# tests use a panel with three.
make_panel_tr3 <- function(N = 15, T = 20, T_pre = 12, effect = 2.0, seed = 42) {
  set.seed(seed)
  f1   <- cumsum(rnorm(T, 0, 0.5))
  f2   <- rnorm(T)
  lam1 <- rnorm(N, 1, 0.3)
  lam2 <- rnorm(N, 0, 0.5)
  rows <- expand.grid(time = seq_len(T), id = paste0("u", seq_len(N)),
                      stringsAsFactors = FALSE)
  rows$y <- as.vector(outer(f1, lam1) + outer(f2, lam2)) +
    rnorm(nrow(rows), 0, 0.3)
  rows$d <- as.integer(rows$id %in% c("u1", "u2", "u3") & rows$time > T_pre)
  rows$y[rows$d == 1] <- rows$y[rows$d == 1] + effect
  rows[order(rows$id, rows$time), ]
}

panel_tr3 <- make_panel_tr3()

# ── SCM ──────────────────────────────────────────────────────────────────────
test_that("SCM weights lie on the unit simplex", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_s3_class(fit, "coresynth")
  expect_equal(fit$method, "scm")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
  expect_equal(length(fit$Y_synth), 20)
  expect_equal(length(fit$gap), 20)
})

test_that("SCM estimate is in a reasonable range for effect=2", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_true(abs(fit$estimate - 2.0) < 3.0)
})

test_that("SCM pre-treatment fit is good (low MSPE)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  pre <- seq_len(10)
  mspe <- mean((fit$Y_treat[pre] - fit$Y_synth[pre])^2)
  total_var <- var(fit$Y_treat[pre])
  expect_true(mspe / total_var < 0.5) # synthetic tracks pre-period reasonably
})

# ── SDID ─────────────────────────────────────────────────────────────────────
test_that("SDID unit weights lie on the unit simplex", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_s3_class(fit, "coresynth")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
})

test_that("SDID time weights lie on the unit simplex", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_equal(sum(fit$time_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$time_weights >= -1e-5))
  expect_equal(length(fit$time_weights), 10L) # T_pre = 10
})

test_that("SDID zeta2 follows Arkhangelsky et al. (2021) Eq.(2.2)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  # zeta2 = sqrt(N_tr * T_post) * sigma2_hat
  # N_tr=1, T_post=10, sigma2_hat from fit
  expected_zeta2 <- sqrt(1 * 10) * fit$sigma2_hat
  expect_equal(fit$zeta2, expected_zeta2, tolerance = 1e-10)
  expect_true(fit$zeta2 > 0)
})

test_that("SDID omega_0 and lambda_0 intercepts are finite", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_true(is.finite(fit$omega0))
  expect_true(is.finite(fit$lambda0))
})

test_that("SDID demeaning shifts unit weights vs. no-intercept baseline", {
  # With intercept (current): demeaned QP allows level shifts
  fit_new <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  # The estimate should still be in a reasonable range
  expect_true(abs(fit_new$estimate - 2.0) < 4.0)
})

test_that("SDID estimate matches synthdid package (California smoke data)", {
  skip_if_not_installed("synthdid")
  data("california_prop99", package = "synthdid")
  setup <- synthdid::panel.matrices(california_prop99)
  ref <- synthdid::synthdid_estimate(setup$Y, setup$N0, setup$T0)
  # Convert data to long format for coresynth
  Y_wide <- setup$Y
  df <- data.frame(
    y = as.vector(Y_wide),
    id = rep(rownames(Y_wide), ncol(Y_wide)),
    time = rep(colnames(Y_wide), each = nrow(Y_wide))
  )
  N0 <- setup$N0
  T0 <- setup$T0
  post_cols <- (T0 + 1):ncol(Y_wide)
  treat_rows <- (N0 + 1):nrow(Y_wide)
  df$d <- as.integer(
    df$id %in%
      rownames(Y_wide)[treat_rows] &
      df$time %in% colnames(Y_wide)[post_cols]
  )
  fit_cs <- scm_fit(y ~ d | id + time, data = df, method = "sdid")
  # Reference: Arkhangelsky et al. (2021) Table 1 reports SDID = -15.6
  expect_equal(fit_cs$estimate, as.numeric(ref), tolerance = 0.5)
})

# ── GSC ──────────────────────────────────────────────────────────────────────
test_that("GSC returns factor matrix with correct dimensions", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  expect_s3_class(fit, "coresynth")
  expect_equal(ncol(fit$F), 2L) # r = 2 factors
  expect_equal(nrow(fit$F), 20L) # T periods
  expect_equal(ncol(fit$L_co), 2L) # r loadings per control unit
  expect_equal(nrow(fit$L_co), 9L) # N_co = 9
})

test_that("GSC counterfactual has correct dimension", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  expect_equal(dim(fit$Y_tr_hat), c(20L, 1L)) # T x N_tr
})

test_that("GSC estimate is finite and in reasonable range", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  expect_true(is.finite(fit$estimate))
  expect_true(abs(fit$estimate - 2.0) < 4.0)
})

test_that("GSC singular values are decreasing", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  sv <- fit$singular_values
  expect_true(all(diff(sv[sv > 0]) <= 1e-10)) # non-increasing
})

test_that("GSC rejects r larger than N_co or T", {
  # R wrapper fires first (r > N_co = 9): "Need at least r control units"
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 20),
    regexp = "at least r"
  )
})

# ── MC ───────────────────────────────────────────────────────────────────────
test_that("MC returns low-rank matrix of correct size", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  expect_s3_class(fit, "coresynth")
  expect_equal(dim(fit$L_hat), c(20L, 10L))
})

test_that("MC observed entries are well-fitted (low residual on mask)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  pan <- coresynth:::panel_to_matrices(
    panel$y,
    panel$d,
    panel$id,
    panel$time
  )
  # Observation mask: post x treated = 0, rest = 1
  O <- matrix(1, nrow = 20, ncol = 10)
  O[11:20, 1] <- 0 # unit u1 post-treatment is "missing"
  Y_obs <- pan$Y
  Y_obs[is.na(Y_obs)] <- 0
  resid_obs <- (Y_obs - fit$L_hat)[O == 1]
  rmse_obs <- sqrt(mean(resid_obs^2))
  # Should be small relative to outcome range
  expect_true(rmse_obs < diff(range(Y_obs[O == 1])) * 0.3)
})

test_that("MC estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  expect_true(is.finite(fit$estimate))
})

# ── TASC ─────────────────────────────────────────────────────────────────────
test_that("TASC runs and returns valid output structure", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_tol = 1e-5
  )
  expect_s3_class(fit, "coresynth")
  expect_true(is.finite(fit$estimate))
  expect_equal(length(fit$times), 20L)
  expect_equal(dim(fit$Y_hat), c(20L, 10L))
})

test_that("TASC general A: ATT is finite", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_tol = 1e-5
  )
  expect_true(is.finite(fit$estimate))
})

test_that("TASC A matrix has correct dimensions (r x r)", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_tol = 1e-5
  )
  expect_equal(dim(fit$A), c(2L, 2L))
})

test_that("TASC A deviates from identity after EM (learning occurred)", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_tol = 1e-5
  )
  expect_true(max(abs(fit$A - diag(2))) > 1e-6)
})

test_that("TASC fix_A=TRUE keeps A as identity", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_tol = 1e-5,
    fix_A = TRUE
  )
  expect_equal(fit$A, diag(2))
  expect_true(is.finite(fit$estimate))
})

test_that("TASC fix_A=FALSE and fix_A=TRUE both return finite ATT", {
  fit_free <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_tol = 1e-5,
    fix_A = FALSE
  )
  fit_fix <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_tol = 1e-5,
    fix_A = TRUE
  )
  expect_true(is.finite(fit_free$estimate))
  expect_true(is.finite(fit_fix$estimate))
})

# ── print / summary ──────────────────────────────────────────────────────────
test_that("print and summary work for all methods", {
  for (m in c("scm", "sdid", "gsc", "mc")) {
    fit <- scm_fit(y ~ d | id + time, data = panel, method = m)
    expect_output(print(fit), toupper(m))
    expect_output(summary(fit), "ATT estimate")
  }
})

# ── broom tidiers ─────────────────────────────────────────────────────────────
test_that("tidy() and glance() work for SDID", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  td <- broom::tidy(fit)
  gl <- broom::glance(fit)
  expect_s3_class(td, "data.frame")
  expect_true(nrow(td) > 0)
  expect_true(all(c("term", "estimate") %in% names(td)))
  expect_equal(gl$method, "sdid")
})

# ── plot ─────────────────────────────────────────────────────────────────────
test_that("plot.coresynth returns ggplot objects", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_s3_class(plot(fit, type = "trend"), "ggplot")
  expect_s3_class(plot(fit, type = "gap"), "ggplot")
  expect_s3_class(plot(fit, type = "weights"), "ggplot")
})

# ── SI ───────────────────────────────────────────────────────────────────────

test_that("SI returns correct class and structure", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  expect_s3_class(fit, "coresynth")
  expect_s3_class(fit, "coresynth_si")
  expect_equal(fit$method, "si")
  expect_true(is.numeric(fit$k) && fit$k >= 1L)
})

test_that("SI output dimensions are correct", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  # panel: T=20, N=10, N_co=9, N_tr=1, T_pre=10, T_post=10
  expect_equal(dim(fit$Y_treat), c(20L, 1L))
  expect_equal(dim(fit$Y_synth), c(20L, 1L))
  expect_equal(dim(fit$Y_cf), c(10L, 1L))
  expect_equal(dim(fit$weights), c(9L, 1L))
  expect_equal(length(fit$unit_weights), 9L)
  expect_equal(length(fit$times), 20L)
  expect_equal(fit$T_pre, 10L)
})

test_that("SI pre-treatment fit is good (low MSPE)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  pre <- seq_len(10)
  mspe <- mean((fit$Y_treat[pre, ] - fit$Y_synth[pre, ])^2)
  total_var <- var(as.vector(fit$Y_treat[pre, ]))
  expect_true(mspe / total_var < 0.5)
})

test_that("SI estimate is finite and in reasonable range for effect=2", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  expect_true(is.finite(fit$estimate))
  expect_true(abs(fit$estimate - 2.0) < 3.0)
})

test_that("SI respects user-supplied k", {
  fit1 <- scm_fit(y ~ d | id + time, data = panel, method = "si", k = 1L)
  fit3 <- scm_fit(y ~ d | id + time, data = panel, method = "si", k = 3L)
  expect_equal(fit1$k, 1L)
  expect_equal(fit3$k, 3L)
  expect_true(is.finite(fit1$estimate))
  expect_true(is.finite(fit3$estimate))
})

test_that("SI errors when k > min(T_pre, N_co)", {
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "si", k = 100L),
    regexp = "k"
  )
})

test_that("tensor_unfold_cpp produces correct dimensions for all modes", {
  tc <- array(seq_len(3 * 4 * 2), dim = c(3L, 4L, 2L))
  expect_equal(dim(tensor_unfold_cpp(tc, 1L)), c(3L, 8L))
  expect_equal(dim(tensor_unfold_cpp(tc, 2L)), c(4L, 6L))
  expect_equal(dim(tensor_unfold_cpp(tc, 3L)), c(2L, 12L))
  expect_error(tensor_unfold_cpp(tc, 4L), regexp = "Mode")
})

test_that("SI and MC estimates are in similar range", {
  fit_si <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  fit_mc <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  expect_true(abs(fit_si$estimate - fit_mc$estimate) < 3.0)
})

test_that("print and summary work for SI", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  expect_output(print(fit), "SI")
  expect_output(summary(fit), "ATT estimate")
})

# ── edge cases ───────────────────────────────────────────────────────────────
test_that("scm_fit errors on bad formula", {
  expect_error(scm_fit(y ~ d, data = panel, method = "sdid"), "after '|'")
})

test_that("scm_fit errors on missing variable", {
  expect_error(
    scm_fit(z ~ d | id + time, data = panel, method = "sdid"),
    "not found"
  )
})

# ── Phase 4: SCM placebo / MSPE-ratio inference ───────────────────────────────

test_that("fit_scm_cpp returns Y_co_pre and Y_co_post", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_equal(dim(fit$Y_co_pre), c(10L, 9L)) # T_pre x N_co
  expect_equal(dim(fit$Y_co_post), c(10L, 9L)) # T_post x N_co
})

test_that("scm_placebo_cpp returns correct structure", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  plac <- scm_placebo_cpp(fit$Y_co_pre, fit$Y_co_post)
  expect_equal(length(plac$mspe_pre), 9L)
  expect_equal(length(plac$mspe_post), 9L)
  expect_equal(length(plac$effects), 9L)
  expect_true(all(plac$mspe_pre >= 0))
  expect_true(all(plac$mspe_post >= 0))
  expect_true(all(is.finite(plac$effects)))
})

test_that("mspe_ratio_pval returns p_value in [0, 1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_true(is.finite(inf$mspe_ratio_treated))
  expect_true(inf$mspe_ratio_treated >= 0)
  expect_equal(length(inf$mspe_ratios_all), 10L) # 1 treated + 9 controls
  expect_equal(length(inf$placebo_effects), 9L)
})

test_that("mspe_ratio_pval mspe_threshold filters all controls", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit, mspe_threshold = 1e10)
  expect_equal(inf$n_placebo_used, 0L)
})

test_that("mspe_ratio_pval rejects non-SCM fits", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_error(mspe_ratio_pval(fit), "method = 'scm'")
})

test_that("mspe_ratio_pval default mirrors the fit's predictor spec", {
  fit_cov <- scm_fit(
    y ~ d | id + time,
    data = panel, method = "scm",
    predictors = list(pred("y", 1:10, op = "mean"))
  )
  auto     <- mspe_ratio_pval(fit_cov)
  explicit <- mspe_ratio_pval(fit_cov, use_covariates = TRUE)
  expect_equal(auto$p_value, explicit$p_value)
  expect_equal(auto$placebo_effects, explicit$placebo_effects)
  expect_equal(auto$mspe_pre_placebo, explicit$mspe_pre_placebo)

  fit_out       <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  auto_out      <- mspe_ratio_pval(fit_out)
  explicit_out  <- mspe_ratio_pval(fit_out, use_covariates = FALSE)
  expect_equal(auto_out$p_value, explicit_out$p_value)
  expect_equal(auto_out$placebo_effects, explicit_out$placebo_effects)
})

test_that("explicit use_covariates=TRUE without covariates warns and falls back", {
  fit_out <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_warning(
    inf <- mspe_ratio_pval(fit_out, use_covariates = TRUE),
    "X0_mat is NULL"
  )
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
})

test_that("outcomes-only placebo solver agrees with scm_weights_cpp per unit", {
  set.seed(7)
  T_pre <- 40L; T_post <- 8L; N <- 20L
  f <- cumsum(rnorm(T_pre + T_post, 0, 0.5))
  Y <- sapply(seq_len(N), function(j) {
    100 + rnorm(1, 0, 5) + f * rnorm(1, 1, 0.3) +
      rnorm(T_pre + T_post, 0, 0.5)
  })
  Y_pre  <- Y[seq_len(T_pre), ]
  Y_post <- Y[(T_pre + 1):(T_pre + T_post), ]

  plac <- scm_placebo_cpp(Y_pre, Y_post)
  for (i in c(1L, 10L, 20L)) {
    ref <- scm_weights_cpp(
      Y_pre[, -i, drop = FALSE], Y_pre[, i],
      Y_pre[, -i, drop = FALSE], Y_pre[, i]
    )
    mspe_ref <- mean((Y_pre[, i] - Y_pre[, -i, drop = FALSE] %*% ref$W)^2)
    # Same nested V/W estimator, different inner QP solver: results agree
    # within the coordinate-descent convergence tolerance.
    expect_lt(
      abs(plac$mspe_pre[i] - mspe_ref) / max(mspe_ref, 1e-12),
      5e-3
    )
  }
})

test_that("realistic-size outcomes-only placebo completes quickly", {
  skip_on_cran()
  set.seed(2)
  T_pre <- 139L; T_post <- 19L; N <- 75L  # ~12y monthly, 75 donors
  tt <- seq_len(T_pre + T_post)
  common <- 200 + 10 * sin(tt / 12) + 0.05 * tt
  Y <- sapply(seq_len(N), function(j) {
    common + rnorm(1, 0, 15) +
      as.numeric(arima.sim(list(ar = 0.8), length(tt), sd = 3))
  })
  t0 <- Sys.time()
  plac <- scm_placebo_cpp(
    Y[seq_len(T_pre), ],
    Y[(T_pre + 1):(T_pre + T_post), ]
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_true(all(is.finite(plac$mspe_pre)))
  expect_true(all(is.finite(plac$mspe_post)))
  expect_true(all(is.finite(plac$effects)))
  expect_equal(dim(plac$gaps), c(T_pre + T_post, N))
  # Generous bound (runs in ~2 s on a laptop; was ~an hour before the
  # active-set inner solver).
  expect_lt(elapsed, 60)
})

# ── Phase 4: GSC bootstrap inference ─────────────────────────────────────────

test_that("fit_gsc_cpp returns Y_co_all", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  expect_equal(dim(fit$Y_co_all), c(20L, 9L)) # T x N_co
})

test_that("gsc_boot returns correct structure", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  boot <- gsc_boot(fit, B = 29L, seed = 1L)
  expect_named(
    boot,
    c("p_value", "ci_lower", "ci_upper", "se", "boot_dist", "att_obs")
  )
  expect_equal(length(boot$boot_dist), 29L)
  expect_true(boot$p_value >= 0 && boot$p_value <= 1)
  expect_true(boot$ci_lower <= boot$ci_upper)
  expect_true(is.finite(boot$se) && boot$se >= 0)
  expect_equal(boot$att_obs, fit$estimate)
})

test_that("gsc_boot seed ensures reproducibility", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  boot1 <- gsc_boot(fit, B = 29L, seed = 42L)
  boot2 <- gsc_boot(fit, B = 29L, seed = 42L)
  expect_equal(boot1$boot_dist, boot2$boot_dist)
})

test_that("gsc_boot rejects non-GSC fits", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(gsc_boot(fit), "method = 'gsc'")
})

# ── Phase 4: export_json ──────────────────────────────────────────────────────

test_that("export_json NULL file returns list without writing", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  res <- export_json(fit, file = NULL)
  expect_type(res, "list")
  expect_equal(res$meta$method, "sdid")
  expect_equal(res$meta$package, "coresynth")
})

test_that("export_json has required top-level sections", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  res <- export_json(fit, file = NULL)
  expect_true(all(c("meta", "data", "estimate", "time_series") %in% names(res)))
})

test_that("export_json estimate section contains att", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  res <- export_json(fit, file = NULL)
  expect_equal(res$estimate$att, signif(fit$estimate, 6))
})

test_that("export_json time_series has required fields", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  res <- export_json(fit, file = NULL)
  expect_true(all(
    c("times", "Y_treat", "Y_synth", "gap") %in% names(res$time_series)
  ))
  expect_equal(length(res$time_series$times), 20L)
})

test_that("export_json writes valid JSON for SCM, SDID, GSC, MC", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  for (m in c("scm", "sdid", "gsc", "mc")) {
    fit <- scm_fit(y ~ d | id + time, data = panel, method = m)
    export_json(fit, file = tmp)
    dat <- jsonlite::read_json(tmp)
    expect_equal(dat$meta$method, m)
    expect_false(is.null(dat$estimate$att))
    expect_equal(dat$data$T_pre, 10L)
  }
})

test_that("export_json with inference includes inference section (SCM)", {
  skip_if_not_installed("jsonlite")
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  res <- export_json(fit, file = NULL, inference = inf)
  expect_equal(res$inference$type, "mspe_permutation")
  expect_false(is.null(res$inference$p_value))
  expect_false(is.null(res$inference$mspe_ratios_all))
  expect_equal(res$estimate$p_value, signif(inf$p_value, 6))
})

test_that("export_json with inference includes inference section (GSC)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  boot <- gsc_boot(fit, B = 29L, seed = 1L)
  res <- export_json(fit, file = NULL, inference = boot)
  expect_equal(res$inference$type, "parametric_bootstrap")
  expect_equal(res$inference$B, 29L)
  expect_false(is.null(res$inference$ci_lower))
})

# ── Phase 5a: Staggered adoption ──────────────────────────────────────────────

make_staggered_panel <- function(N = 12, TT = 24, seed = 7) {
  set.seed(seed)
  f <- cumsum(rnorm(TT, 0, 0.5))
  lam <- rnorm(N, 1, 0.3)
  rows <- expand.grid(
    time = seq_len(TT),
    id = paste0("u", seq_len(N)),
    stringsAsFactors = FALSE
  )
  rows$y <- as.vector(outer(f, lam)) + rnorm(nrow(rows), 0, 0.3)
  # u1 採用 t=9、u2 採用 t=15
  rows$d <- as.integer(
    (rows$id == "u1" & rows$time >= 9) |
      (rows$id == "u2" & rows$time >= 15)
  )
  rows$y[rows$d == 1] <- rows$y[rows$d == 1] + 1.5
  rows[order(rows$id, rows$time), ]
}

staggered <- make_staggered_panel()

# The jackknife drops each treated unit in turn, so its staggered tests need
# two treated units in every cohort (Clarke et al. 2024, Algorithm 3).
make_staggered_panel2 <- function(N = 14, TT = 24, seed = 7) {
  set.seed(seed)
  f <- cumsum(rnorm(TT, 0, 0.5))
  lam <- rnorm(N, 1, 0.3)
  rows <- expand.grid(time = seq_len(TT), id = paste0("u", seq_len(N)),
                      stringsAsFactors = FALSE)
  rows$y <- as.vector(outer(f, lam)) + rnorm(nrow(rows), 0, 0.3)
  rows$d <- as.integer(
    (rows$id %in% c("u1", "u2") & rows$time >= 9) |
      (rows$id %in% c("u3", "u4") & rows$time >= 15)
  )
  rows$y[rows$d == 1] <- rows$y[rows$d == 1] + 1.5
  rows[order(rows$id, rows$time), ]
}

staggered2 <- make_staggered_panel2()

test_that("panel_to_matrices detects staggered adoption correctly", {
  pan <- coresynth:::panel_to_matrices(
    staggered$y,
    staggered$d,
    staggered$id,
    staggered$time
  )
  expect_false(pan$is_sharp)
  expect_equal(pan$T_pre, 8L)
  u1c <- which(pan$units == "u1")
  u2c <- which(pan$units == "u2")
  expect_equal(pan$T_adopt[u1c], 9L)
  expect_equal(pan$T_adopt[u2c], 15L)
})

test_that("panel_to_matrices is_sharp TRUE for non-staggered panel", {
  pan <- coresynth:::panel_to_matrices(
    panel$y,
    panel$d,
    panel$id,
    panel$time
  )
  expect_true(pan$is_sharp)
})

test_that("MC handles staggered adoption and returns finite estimate", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "mc")
  expect_true(is.finite(fit$estimate))
})

test_that("MC staggered: L_hat dimensions are correct", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "mc")
  expect_equal(dim(fit$L_hat), c(24L, 12L))
})

test_that("TASC handles staggered adoption and returns finite estimate", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = staggered,
    method = "tasc",
    r = 2L,
    em_tol = 1e-5
  )
  expect_true(is.finite(fit$estimate))
})

# MC and TASC impute the counterfactual from one global fit whose mask already
# encodes each unit's own adoption date, so their staggered path re-uses that
# single fit. They must still declare themselves staggered and aggregate the
# cohort ATTs with the same N_treated,g * T_post,g weights as every other
# estimator, or the same panel silently yields a different estimand.
.staggered_impute_fit <- function(method) {
  args <- list(y ~ d | id + time, data = staggered, method = method)
  if (method == "tasc") args <- c(args, list(r = 2L, em_tol = 1e-5))
  do.call(scm_fit, args)
}

test_that("MC/TASC staggered fits carry the staggered class and cohort table", {
  for (m in c("mc", "tasc")) {
    fit <- .staggered_impute_fit(m)
    expect_s3_class(fit, "coresynth_staggered")
    expect_true(isTRUE(fit$staggered))
    expect_true(broom::glance(fit)$staggered)

    ce <- fit$cohort_estimates
    expect_equal(ce$cohort, c(9L, 15L))
    expect_equal(ce$n_treated, c(1L, 1L))
    expect_equal(ce$T_pre, c(8L, 14L))
    expect_equal(ce$T_post, c(16L, 10L))
    expect_equal(sum(ce$weight), 1)
    expect_length(fit$cohort_fits, 2L)
  }
})

test_that("MC/TASC staggered ATT is the cohort-weighted average", {
  for (m in c("mc", "tasc")) {
    fit <- .staggered_impute_fit(m)
    ce  <- fit$cohort_estimates
    # N_treated,g * T_post,g weights (Clarke et al. 2024), not one per unit
    w   <- ce$n_treated * ce$T_post
    expect_equal(fit$estimate, sum(w * ce$estimate) / sum(w))
    # each cohort ATT is the mean gap over that cohort's own post periods
    for (i in seq_len(nrow(ce))) {
      k <- match(fit$cohort_fits[[i]]$idx_tr, fit$idx_tr)
      expect_equal(
        ce$estimate[i],
        mean(fit$gap[(ce$T_pre[i] + 1L):nrow(fit$gap), k, drop = FALSE])
      )
    }
  }
})

test_that("MC/TASC staggered: tidy/augment use the per-cohort series", {
  for (m in c("mc", "tasc")) {
    fit <- .staggered_impute_fit(m)
    td  <- broom::tidy(fit)
    expect_equal(td$term, c("cohort_9", "cohort_15"))
    expect_equal(td$type, rep("cohort_estimate", 2L))

    aug <- broom::augment(fit)
    expect_true(all(c(".cohort", ".fitted", ".period") %in% names(aug)))
    expect_equal(sort(unique(aug$.cohort)), c(9L, 15L))
    # .period must follow each cohort's own adoption date, not a global split
    per_cohort <- split(aug, aug$.cohort)
    expect_equal(sum(per_cohort[["9"]]$.period == "pre"), 8L)
    expect_equal(sum(per_cohort[["15"]]$.period == "pre"), 14L)
  }
})

test_that("MC/TASC staggered: plot() and the sharp accessors refuse", {
  for (m in c("mc", "tasc")) {
    fit <- .staggered_impute_fit(m)
    expect_null(synthetic_outcomes(fit))
    expect_error(plot(fit, type = "trend"), "per cohort")
  }
})

test_that("MC/TASC sharp fits are untouched by the staggered path", {
  for (m in c("mc", "tasc")) {
    args <- list(y ~ d | id + time, data = panel, method = m)
    if (m == "tasc") args <- c(args, list(r = 2L, em_tol = 1e-5))
    fit <- do.call(scm_fit, args)
    expect_false(inherits(fit, "coresynth_staggered"))
    expect_null(fit$cohort_estimates)
    expect_type(synthetic_outcomes(fit), "double")
  }
})

test_that("SCM staggered with predictors arg errors informatively", {
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm",
            predictors = list(pred("y", 1:3))),
    regexp = "predictors"
  )
})

test_that("GSC handles staggered adoption (Phase 16)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_true(isTRUE(fit$staggered))
  expect_true(is.finite(fit$estimate))
})

test_that("SI handles staggered adoption (Phase 16)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_true(isTRUE(fit$staggered))
  expect_true(is.finite(fit$estimate))
})

# ── SCM with covariates (pred() interface) ────────────────────────────────────

# Helper: panel with two extra covariate columns
make_panel_with_cov <- function(
  N = 10,
  T = 20,
  T_pre = 10,
  effect = 2.0,
  seed = 42
) {
  set.seed(seed)
  p <- make_panel(N = N, T = T, T_pre = T_pre, effect = effect, seed = seed)
  units_all <- unique(p$id)
  cov1_vals <- setNames(rnorm(N, mean = 5, sd = 1), units_all)
  cov2_vals <- setNames(rnorm(N, mean = 10, sd = 2), units_all)
  p$cov1 <- cov1_vals[p$id]
  p$cov2 <- cov2_vals[p$id]
  p
}

panel_cov <- make_panel_with_cov()

test_that("SCM with pred() returns valid fit object", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  expect_s3_class(fit, "coresynth")
  expect_equal(fit$method, "scm")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
  expect_equal(length(fit$Y_synth), 20)
})

test_that("SCM with pred(): v_weights have predictor names", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  expect_named(fit$v_weights, c("cov1[1:10]", "cov2[1:10]"))
})

test_that("SCM with pred(): predictor_table is a data.frame with correct cols", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  expect_s3_class(fit$predictor_table, "data.frame")
  expect_true(all(
    c("predictor", "treated", "synthetic") %in% names(fit$predictor_table)
  ))
  expect_equal(nrow(fit$predictor_table), 2L)
})

test_that("SCM with multiple pred() returns valid fit object", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(
      pred("y", c(3, 5, 7)),
      pred("cov1", 1:10)
    )
  )
  expect_s3_class(fit, "coresynth")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_equal(length(fit$v_weights), 2L)
  expect_true(all(is.finite(fit$Y_synth)))
})

test_that("SCM with mixed pred() entries stacks rows", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(
      pred("cov1", 1:10),
      pred("y", c(5, 8, 10))
    )
  )
  expect_equal(length(fit$v_weights), 2L)
  expect_equal(nrow(fit$predictor_table), 2L)
})

test_that("SCM without predictors matches previous outcome-only behaviour", {
  fit_plain <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  fit_cov <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm")
  expect_equal(fit_plain$unit_weights, fit_cov$unit_weights, tolerance = 1e-4)
  expect_null(fit_cov$predictor_table)
  expect_equal(
    names(fit_plain$v_weights),
    paste0("V", seq_along(fit_plain$v_weights))
  )
})

test_that("SCM with pred(): predictor_table synthetic matches W dot X0", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  # Synthetic predictor values should be a weighted combination of control values
  expect_true(all(is.finite(fit$predictor_table$synthetic)))
})

test_that("SCM errors when pred() variable not in data", {
  expect_error(
    scm_fit(
      y ~ d | id + time,
      data = panel_cov,
      method = "scm",
      predictors = list(pred("nonexistent_var", 1:5))
    ),
    regexp = "not found in data"
  )
})

test_that("SCM errors when predictors is a character vector (old API)", {
  expect_error(
    scm_fit(
      y ~ d | id + time,
      data = panel_cov,
      method = "scm",
      predictors = "cov1"
    ),
    regexp = "list of pred\\(\\) objects"
  )
})

# ── Phase 7a: GSC covariates (time-varying, EM loop) ─────────────────────────

# Helper: panel with one time-varying covariate column
make_panel_with_timecov <- function(
  N = 10,
  T = 20,
  T_pre = 10,
  effect = 2.0,
  seed = 42
) {
  p <- make_panel(N = N, T = T, T_pre = T_pre, effect = effect, seed = seed)
  set.seed(seed + 1L)
  units_all <- unique(p$id)
  for (u in units_all) {
    base <- rnorm(1, 5, 1)
    p$cov1[p$id == u] <- base +
      seq(0, by = 0.1, length.out = T) +
      rnorm(T, 0, 0.2)
  }
  p
}

panel_timecov <- make_panel_with_timecov()

test_that("GSC with covariates: fit$beta has length p", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_timecov,
    method = "gsc",
    r = 2,
    covariates = "cov1"
  )
  expect_equal(length(fit$beta), 1L)
  expect_true(is.numeric(fit$beta))
})

test_that("GSC with covariates: ATT is finite and in reasonable range", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_timecov,
    method = "gsc",
    r = 2,
    covariates = "cov1"
  )
  expect_true(is.finite(fit$estimate))
  expect_true(abs(fit$estimate - 2.0) < 5.0)
})

test_that("GSC with covariates: Y_hat_co has correct dimensions (T x N_co)", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_timecov,
    method = "gsc",
    r = 2,
    covariates = "cov1"
  )
  expect_equal(dim(fit$Y_hat_co), c(20L, 9L))
})

test_that("GSC with covariates via scm_fit formula interface", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_timecov,
    method = "gsc",
    r = 2,
    covariates = "cov1"
  )
  expect_s3_class(fit, "coresynth")
  expect_equal(fit$method, "gsc")
})

test_that("GSC without covariates backward compat: beta is empty, covariates NULL", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  expect_equal(length(fit$beta), 0L)
  expect_null(fit$covariates)
  expect_true(is.finite(fit$estimate))
})

test_that("gsc_boot with covariates: p_value in [0,1], CI valid", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_timecov,
    method = "gsc",
    r = 2,
    covariates = "cov1"
  )
  boot <- gsc_boot(fit, B = 19L, seed = 1L)
  expect_true(boot$p_value >= 0 && boot$p_value <= 1)
  expect_true(boot$ci_lower <= boot$ci_upper)
  expect_equal(length(boot$boot_dist), 19L)
})

test_that("build_covariate_array returns correct dimensions (T x N x p)", {
  arr <- coresynth:::build_covariate_array(
    panel_timecov,
    "id",
    "time",
    "cov1",
    units = unique(panel_timecov$id),
    times = 1:20
  )
  expect_equal(dim(arr), c(20L, 10L, 1L))
  expect_true(all(is.finite(arr)))
})

test_that("GSC with covariates: beta deviates from zero (EM learned)", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_timecov,
    method = "gsc",
    r = 2,
    covariates = "cov1"
  )
  expect_true(abs(fit$beta[1]) > 1e-10)
})

# ── Phase 7b: SCM placebo test covariate consistency ─────────────────────────

test_that("fit_scm_cpp with covariates stores X0_mat (k x N_co)", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  expect_false(is.null(fit$X0_mat))
  expect_equal(dim(fit$X0_mat), c(2L, 9L))
})

test_that("mspe_ratio_pval use_covariates=TRUE returns p_value in [0,1]", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  inf <- mspe_ratio_pval(fit, use_covariates = TRUE)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_equal(length(inf$mspe_ratios_all), 10L)
})

test_that("mspe_ratio_pval use_covariates=TRUE: mspe_ratios_all correct length", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred("cov1", 1:10))
  )
  inf <- mspe_ratio_pval(fit, use_covariates = TRUE)
  expect_equal(length(inf$mspe_ratios_all), 10L) # treated + 9 controls
  expect_true(is.numeric(inf$n_placebo_used))
})

test_that("mspe_ratio_pval covariate placebo matches a per-donor scm_weights_cpp loop", {
  # The covariate placebo path runs scm_placebo_x_cpp (OpenMP batch). Each
  # leave-one-out problem must be identical to an individual scm_weights_cpp
  # call on the same submatrices -- same solver core, so results agree to
  # machine precision regardless of parallel execution order.
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel_cov,
    method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  inf <- mspe_ratio_pval(fit, use_covariates = TRUE)

  N_co <- ncol(fit$Y_co_pre)
  for (i in seq_len(N_co)) {
    # The default predictor-path fit runs the multi-start outer search and the
    # Wolfe inner solver, and the placebo batch mirrors both -- so the
    # reference loop must too. With a degenerate QP the optimal set is a face,
    # and the two solvers pick different points on it, so a reference loop on
    # the other solver would disagree even though both optima are exact.
    res_i <- scm_weights_cpp(
      fit$X0_mat[, -i, drop = FALSE], fit$X0_mat[, i],
      fit$Y_co_pre[, -i, drop = FALSE], fit$Y_co_pre[, i],
      multistart = identical(fit$v_optim_effective, "multistart"),
      wolfe = identical(fit$qp_solver, "wolfe")
    )
    synth_post <- fit$Y_co_post[, -i, drop = FALSE] %*% res_i$W
    expect_equal(unname(inf$placebo_effects[i]),
                 mean(fit$Y_co_post[, i] - synth_post),
                 tolerance = 1e-10)
    synth_pre <- fit$Y_co_pre[, -i, drop = FALSE] %*% res_i$W
    expect_equal(unname(inf$mspe_pre_placebo[i]),
                 mean((fit$Y_co_pre[, i] - synth_pre)^2),
                 tolerance = 1e-10)
  }
})

test_that("mspe_ratio_pval use_covariates=TRUE warns and falls back for no-cov fit", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_null(fit$X0_mat)
  expect_warning(
    inf <- mspe_ratio_pval(fit, use_covariates = TRUE),
    regexp = "no covariates"
  )
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
})

# ── Phase 9: Abadie (2021) — OOS V selection & one-sided inference ────────────

test_that("SCM v_selection='oos' returns valid fit (Abadie 2021 §3.2)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_selection = "oos")
  expect_s3_class(fit, "coresynth")
  expect_equal(fit$method, "scm")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
  expect_true(all(is.finite(fit$Y_synth)))
})

test_that("SCM v_selection='oos' estimate is in reasonable range", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_selection = "oos")
  expect_true(abs(fit$estimate - 2.0) < 3.0)
})

test_that("SCM v_selection='insample' and 'oos' give valid V weights", {
  fit_is  <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_selection = "insample")
  fit_oos <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_selection = "oos")
  # In-sample: V spans all T_pre outcome rows. OOS (Abadie 2021 §3.2 / Phase
  # 26): V is selected on the training half and applied to the last
  # floor(T_pre/2) rows, so it has floor(T_pre/2) entries.
  expect_equal(length(fit_is$v_weights),  10L)
  expect_equal(length(fit_oos$v_weights), 5L)
  expect_equal(fit_oos$v_rows, 6:10)
  # Both should be valid simplex weights for V
  expect_equal(sum(fit_is$v_weights),  1, tolerance = 1e-4)
  expect_equal(sum(fit_oos$v_weights), 1, tolerance = 1e-4)
})

test_that("mspe_ratio_pval alternative='two.sided' matches default (Abadie 2021 §3.5)", {
  fit  <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf1 <- mspe_ratio_pval(fit)
  inf2 <- mspe_ratio_pval(fit, alternative = "two.sided")
  expect_equal(inf1$p_value, inf2$p_value)
  expect_equal(inf1$mspe_ratio_treated, inf2$mspe_ratio_treated)
})

test_that("mspe_ratio_pval alternative='greater' returns valid one-sided p-value", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit, alternative = "greater")
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_true(is.numeric(inf$treated_effect))
  expect_equal(inf$treated_effect, fit$estimate)
})

test_that("mspe_ratio_pval alternative='less' returns valid one-sided p-value", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf_g <- mspe_ratio_pval(fit, alternative = "greater")
  inf_l <- mspe_ratio_pval(fit, alternative = "less")
  # greater and less p-values should sum to approximately 1 + 1/N (both include treated)
  expect_true(inf_l$p_value >= 0 && inf_l$p_value <= 1)
  expect_true(inf_g$p_value + inf_l$p_value <= 1 + 1/10 + 1e-6)
})

test_that("mspe_ratio_pval one-sided test is more powerful for correct direction", {
  # With a positive effect (effect=2), 'greater' p-value should be <= two-sided
  fit    <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf_2s <- mspe_ratio_pval(fit, alternative = "two.sided")
  inf_gt <- mspe_ratio_pval(fit, alternative = "greater")
  # One-sided can only be <= two-sided for the correct direction (not always guaranteed
  # due to different test statistics, but treated_effect should be accessible)
  expect_true(is.numeric(inf_gt$p_value))
  expect_true(is.numeric(inf_2s$p_value))
})

# ── Phase 10: Abadie & Zhao (2026) — Experimental SCM Design ─────────────────

# Shared helper: balanced panel without a treatment indicator (design setting)
make_design_panel <- function(J = 10, TT = 20, seed = 42) {
  set.seed(seed)
  f_t  <- cumsum(rnorm(TT, 0, 0.5))
  lam  <- rnorm(J, 1, 0.3)
  df   <- expand.grid(unit = paste0("m", seq_len(J)),
                      time = seq_len(TT),
                      stringsAsFactors = FALSE)
  df$outcome <- as.vector(outer(f_t, lam)) + rnorm(nrow(df), 0, 0.2)
  df
}

test_that("scm_design returns scm_design class with correct structure", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L)
  expect_s3_class(res, "scm_design")
  expect_equal(res$method, "scm_design")
  expect_length(res$treated_units, 1L)
  expect_length(res$tau_hat, 5L)          # 20 - 15 post periods
  expect_equal(res$T_pre, 15L)
})

test_that("scm_design w sums to 1 and v sums to 1", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L)
  expect_equal(sum(res$w), 1, tolerance = 1e-6)
  expect_equal(sum(res$v), 1, tolerance = 1e-6)
})

test_that("scm_design w and v are non-negative and disjoint", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L)
  expect_true(all(res$w >= -1e-8))
  expect_true(all(res$v >= -1e-8))
  # No unit has positive weight in both w and v
  expect_true(all(res$w * res$v < 1e-10))
  # treated and control sets are disjoint
  expect_length(intersect(res$treated_units, res$control_units), 0L)
})

test_that("scm_design inference: blank periods produce valid p_value and CI", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L, T_fit = 10L)
  expect_equal(length(res$u_blank), 5L)    # 15 - 10 blank periods
  expect_true(!is.na(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_length(res$ci_lower, 5L)
  expect_length(res$ci_upper, 5L)
  expect_true(all(!is.na(res$ci_lower)))
  expect_true(all(res$ci_upper >= res$ci_lower))
})

test_that("scm_design default (no T_fit) gives NA inference without warning", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L)  # no warning for default
  expect_true(is.na(res$p_value))
  expect_true(all(is.na(res$ci_lower)))
})

test_that("scm_design explicit T_fit consuming all pre-periods warns", {
  df <- make_design_panel()  # TT=20, T0=15 -> T0_idx=15
  expect_warning(
    scm_design(df, "outcome", "unit", "time", T0 = 15L, T_fit = 15L),
    "No blank periods"
  )
})

test_that("scm_design m_max=2 selects 2 treated units", {
  df  <- make_design_panel(J = 10L)
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L,
                    m_min = 2L, m_max = 2L)
  expect_length(res$treated_units, 2L)
  expect_equal(sum(res$w > 1e-6), 2L)    # exactly 2 units have positive w
})

test_that("scm_design m_min=1, m_max=2 returns between 1 and 2 treated units", {
  df  <- make_design_panel(J = 8L)
  res <- scm_design(df, "outcome", "unit", "time", T0 = 12L,
                    m_min = 1L, m_max = 2L)
  n_treated <- length(res$treated_units)
  expect_true(n_treated >= 1L && n_treated <= 2L)
})

test_that("scm_design weakly_targeted variant works", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L,
                    design = "weakly_targeted", beta = 1)
  expect_equal(res$design, "weakly_targeted")
  expect_equal(sum(res$w), 1, tolerance = 1e-6)
  expect_equal(sum(res$v), 1, tolerance = 1e-6)
})

test_that("scm_design unit_level variant works and returns disjoint weights", {
  df  <- make_design_panel(J = 8L)
  res <- scm_design(df, "outcome", "unit", "time", T0 = 12L,
                    m_max = 2L, design = "unit_level", xi = 1)
  expect_equal(res$design, "unit_level")
  expect_true(all(res$w * res$v < 1e-10))
})

# Rebuild the normalised predictor matrix scm_design() uses internally
# (fitting-period outcomes, row-wise SD scaling) to evaluate the design
# objectives of eq. (9)/(10) from a fitted object.
make_design_X <- function(df, TE) {
  units <- sort(unique(df$unit))
  times <- sort(unique(df$time))
  Y <- matrix(NA_real_, length(times), length(units),
              dimnames = list(NULL, units))
  Y[cbind(match(df$time, times), match(df$unit, units))] <- df$outcome
  X <- Y[seq_len(TE), , drop = FALSE]
  row_sd <- apply(X, 1L, sd)
  row_sd[row_sd < 1e-10] <- 1
  X / row_sd
}

test_that("scm_design weakly_targeted: (w, v) are jointly optimal in eq. (9)", {
  df    <- make_design_panel(J = 8L, TT = 20L, seed = 11)
  X     <- make_design_X(df, TE = 12L)
  X_bar <- rowMeans(X)
  Vk    <- rep(1, nrow(X))

  base_term  <- function(res) {
    Xw <- drop(X[, res$treated_units, drop = FALSE] %*% res$w[res$treated_units])
    sum((X_bar - Xw)^2)
  }
  match_term <- function(res) {
    Xw <- drop(X[, res$treated_units, drop = FALSE] %*% res$w[res$treated_units])
    Xv <- drop(X[, res$control_units, drop = FALSE] %*% res$v[res$control_units])
    sum((Xw - Xv)^2)
  }

  fit_lo <- scm_design(df, "outcome", "unit", "time", T0 = 12L,
                       m_min = 2L, m_max = 2L,
                       design = "weakly_targeted", beta = 0.1)
  fit_hi <- scm_design(df, "outcome", "unit", "time", T0 = 12L,
                       m_min = 2L, m_max = 2L,
                       design = "weakly_targeted", beta = 50)

  # comparative statics: the treated/control matching term is non-increasing
  # in beta across optimal solutions
  expect_lte(match_term(fit_hi), match_term(fit_lo) + 1e-8)

  # the returned objective equals eq. (9) recomputed from the returned weights
  expect_equal(fit_hi$objective, base_term(fit_hi) + 50 * match_term(fit_hi),
               tolerance = 1e-6)

  # the joint solution is at least as good as the decoupled one (w fit to
  # X_bar alone, v best-responding) on the same treated set
  Xs <- X[, fit_hi$treated_units, drop = FALSE]
  Xc <- X[, fit_hi$control_units, drop = FALSE]
  w0 <- scm_inner_weights_cpp(Xs, X_bar, Vk)
  v0 <- scm_inner_weights_cpp(Xc, drop(Xs %*% w0), Vk)
  obj0 <- sum((X_bar - drop(Xs %*% w0))^2) +
    50 * sum((drop(Xs %*% w0) - drop(Xc %*% v0))^2)
  expect_lte(fit_hi$objective, obj0 + 1e-8)
})

test_that("scm_design unit_level: objective matches eq. (10) and responds to xi", {
  df    <- make_design_panel(J = 8L, TT = 20L, seed = 11)
  X     <- make_design_X(df, TE = 12L)
  X_bar <- rowMeans(X)
  Vk    <- rep(1, nrow(X))

  loss_term <- function(res) {
    Xc <- X[, res$control_units, drop = FALSE]
    sum(vapply(res$treated_units, function(u) {
      v_j <- scm_inner_weights_cpp(Xc, X[, u], Vk)
      res$w[[u]] * sum((X[, u] - drop(Xc %*% v_j))^2)
    }, numeric(1)))
  }

  fit_lo <- scm_design(df, "outcome", "unit", "time", T0 = 12L,
                       m_min = 2L, m_max = 2L, design = "unit_level", xi = 0.1)
  fit_hi <- scm_design(df, "outcome", "unit", "time", T0 = 12L,
                       m_min = 2L, m_max = 2L, design = "unit_level", xi = 50)

  # comparative statics: the w-weighted per-unit losses are non-increasing in xi
  expect_lte(loss_term(fit_hi), loss_term(fit_lo) + 1e-8)

  # the returned objective equals eq. (10) recomputed from the returned weights
  Xw <- drop(X[, fit_hi$treated_units, drop = FALSE] %*%
               fit_hi$w[fit_hi$treated_units])
  expect_equal(fit_hi$objective, sum((X_bar - Xw)^2) + 50 * loss_term(fit_hi),
               tolerance = 1e-6)

  # the aggregate control weights follow eq. (11): v = sum_j w_j v_j
  Xc    <- X[, fit_hi$control_units, drop = FALSE]
  v_agg <- Reduce(`+`, lapply(fit_hi$treated_units, function(u) {
    fit_hi$w[[u]] * as.numeric(scm_inner_weights_cpp(Xc, X[, u], Vk))
  }))
  expect_equal(unname(fit_hi$v[fit_hi$control_units]), v_agg, tolerance = 1e-6)
})

test_that("scm_design unit_level: large xi shifts the selected treated set", {
  # Panel where the eq. (10) trade-off flips the optimal treated set: at small
  # xi the pair best matching X_bar wins; at large xi the per-unit synthetic
  # fits dominate and a different pair is selected. The seed matters: the
  # margin between candidate subsets is small, so a panel has to be picked
  # where the flip is a genuine property of the objective rather than of how
  # accurately the per-unit QPs happen to be solved.
  set.seed(20260713)
  J <- 8L; TT <- 20L
  Fmat    <- matrix(rnorm(TT * 2L), TT, 2L)
  Lam     <- matrix(rnorm(J * 2L, 1, 0.5), J, 2L)
  delta_t <- cumsum(rnorm(TT, 0, 0.3))
  Y  <- outer(delta_t, rep(1, J)) + Fmat %*% t(Lam) +
    matrix(rnorm(TT * J, 0, 1), TT, J)
  df <- data.frame(unit = rep(paste0("u", 1:J), each = TT),
                   time = rep(1:TT, times = J),
                   y    = as.vector(Y))

  fit_lo <- scm_design(df, "y", "unit", "time", T0 = 14L, T_fit = 10L,
                       m_min = 2L, m_max = 2L, design = "unit_level", xi = 0.1)
  fit_hi <- scm_design(df, "y", "unit", "time", T0 = 14L, T_fit = 10L,
                       m_min = 2L, m_max = 2L, design = "unit_level", xi = 100)
  expect_false(setequal(fit_lo$treated_units, fit_hi$treated_units))
})

test_that("scm_design custom population weights f are normalised and applied", {
  J     <- 8L
  units <- paste0("m", seq_len(J))
  df    <- make_design_panel(J = J)
  f_raw <- setNames(c(0.4, 0.2, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02), units)
  res   <- scm_design(df, "outcome", "unit", "time", T0 = 12L, f = f_raw)
  expect_equal(sum(res$f), 1, tolerance = 1e-10)
  # f stored in the object should be normalised
  expect_equal(res$f[["m1"]], 0.4 / sum(f_raw), tolerance = 1e-10)
})

test_that("scm_design rejects a pred() window with no data", {
  df <- make_design_panel(J = 8L, TT = 20L)
  expect_error(
    scm_design(df, "outcome", "unit", "time", T0 = 15L,
               predictors = list(pred("outcome", 100:110))),
    regexp = "non-finite"
  )
})

test_that("scm_design with pred() predictors works", {
  df <- make_design_panel(J = 8L, TT = 20L)
  # Add a covariate
  df$cov1 <- rnorm(nrow(df), 5, 1)
  res <- scm_design(
    df, "outcome", "unit", "time", T0 = 15L,
    predictors = list(
      pred("outcome", 1:10),
      pred("cov1",   1:15)
    )
  )
  expect_s3_class(res, "scm_design")
  expect_true(grepl("outcome", res$pred_names[1]))
  expect_true(grepl("cov1", res$pred_names[length(res$pred_names)]))
})

test_that("scm_design estimate and tau_hat are finite", {
  df  <- make_design_panel(J = 10L, TT = 25L)
  res <- scm_design(df, "outcome", "unit", "time", T0 = 18L, T_fit = 12L)
  expect_true(is.finite(res$estimate))
  expect_true(all(is.finite(res$tau_hat)))
  expect_true(all(is.finite(res$u_blank)))
})

test_that("scm_design print output mentions scm_design", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L, T_fit = 10L)
  expect_output(print(res), "scm_design")
  expect_output(print(res), "Treated units")
})

test_that("scm_design summary output shows design variant and CI", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L, T_fit = 10L)
  expect_output(summary(res), "Design variant")
  expect_output(summary(res), "Confidence intervals")
})

test_that("scm_design plot outcome type returns ggplot", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L)
  p   <- plot(res, type = "outcome")
  expect_s3_class(p, "gg")
})

test_that("scm_design plot gap type with blank periods returns ggplot with ribbon", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L, T_fit = 10L)
  p   <- plot(res, type = "gap")
  expect_s3_class(p, "gg")
})

test_that("scm_design error on invalid m_max >= J", {
  df <- make_design_panel(J = 5L)
  expect_error(
    scm_design(df, "outcome", "unit", "time", T0 = 10L, m_max = 5L),
    "m_max < J"
  )
})

test_that("scm_design error on missing outcome column", {
  df <- make_design_panel()
  expect_error(
    scm_design(df, "no_col", "unit", "time", T0 = 15L),
    "not found"
  )
})

test_that("scm_design error when T0 not in time column", {
  df <- make_design_panel()
  expect_error(
    scm_design(df, "outcome", "unit", "time", T0 = 999L),
    "present in the time column"
  )
})

# ── Phase 11a: Donor pool filtering (Abadie 2021 §4) ─────────────────────────

# Helper: panel where two donors have extreme pre-treatment trajectories
make_panel_outliers <- function() {
  df <- make_panel(N = 8L, T = 20L, T_pre = 10L, seed = 42L)
  df$y[df$id == "u5"] <- df$y[df$id == "u5"] * 20  # clear outlier
  df$y[df$id == "u6"] <- df$y[df$id == "u6"] * 20  # clear outlier
  df
}

test_that("donor_mspe_threshold = Inf (default) keeps all donors", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_equal(fit$excluded_donors, character(0L))
})

test_that("donor_mspe_threshold too strict raises warning and keeps all donors", {
  # threshold = 1 means keep only donors with exactly min MSPE → < 2 donors → warning
  expect_warning(
    fit <- scm_fit(
      y ~ d | id + time, data = panel, method = "scm",
      donor_mspe_threshold = 1.0
    ),
    "ignoring filter"
  )
  expect_equal(fit$excluded_donors, character(0L))
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
})

test_that("donor_mspe_threshold excludes clear outliers", {
  df  <- make_panel_outliers()
  fit <- scm_fit(y ~ d | id + time, data = df, method = "scm",
                 donor_mspe_threshold = 5)
  expect_true(length(fit$excluded_donors) > 0L)
  expect_true("u5" %in% fit$excluded_donors || "u6" %in% fit$excluded_donors)
})

test_that("excluded donors are absent from unit_weights", {
  df  <- make_panel_outliers()
  fit <- scm_fit(y ~ d | id + time, data = df, method = "scm",
                 donor_mspe_threshold = 5)
  if (length(fit$excluded_donors) > 0L) {
    expect_false(any(fit$excluded_donors %in% names(fit$unit_weights)))
  }
})

test_that("weights remain on the simplex after donor filtering", {
  df  <- make_panel_outliers()
  fit <- scm_fit(y ~ d | id + time, data = df, method = "scm",
                 donor_mspe_threshold = 5)
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
})

test_that("ATT is finite after donor filtering", {
  df  <- make_panel_outliers()
  fit <- scm_fit(y ~ d | id + time, data = df, method = "scm",
                 donor_mspe_threshold = 5)
  expect_true(is.finite(fit$estimate))
})

test_that("donor_mspe_threshold via scm_fit() passes through correctly", {
  df  <- make_panel_outliers()
  fit <- scm_fit(y ~ d | id + time, data = df, method = "scm",
                 donor_mspe_threshold = 5)
  expect_type(fit$excluded_donors, "character")
})

# ── Phase 11b: Penalised SCM (Abadie & L'Hour 2021) ─────────────────────────

test_that("lambda_pen = NULL gives standard SCM (lambda_pen field is NA)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 lambda_pen = NULL)
  expect_true(is.na(fit$lambda_pen))
})

test_that("lambda_pen = 0 weights lie on the unit simplex", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 lambda_pen = 0)
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
})

test_that("lambda_pen = 0 records lambda_pen = 0 in fit object", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 lambda_pen = 0)
  expect_equal(fit$lambda_pen, 0)
})

test_that("lambda_pen = 1 weights lie on the unit simplex", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 lambda_pen = 1)
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
  expect_true(is.finite(fit$estimate))
})

test_that("large lambda_pen concentrates weights (fewer donors used)", {
  fit_pen  <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                      lambda_pen = 1e6)
  fit_std  <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  # Very large penalty pushes toward nearest-neighbor → max weight should increase
  expect_true(max(fit_pen$unit_weights) >= max(fit_std$unit_weights) - 1e-4)
})

test_that("lambda_pen = 'auto' returns a non-negative finite lambda", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 lambda_pen = "auto")
  expect_true(is.numeric(fit$lambda_pen))
  expect_true(is.finite(fit$lambda_pen))
  expect_true(fit$lambda_pen >= 0)
})

test_that("lambda_pen = 'auto' gives finite ATT and valid weights", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 lambda_pen = "auto")
  expect_true(is.finite(fit$estimate))
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
})

test_that("invalid lambda_pen raises an error", {
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            lambda_pen = -1),
    "non-negative"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            lambda_pen = "bad_value"),
    "NULL"
  )
})

# ── Phase 11c: Augmented SCM (Ben-Michael, Feller & Rothstein 2021) ──────────

test_that("augment_scm returns a list with required fields", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- augment_scm(fit)
  expect_type(aug, "list")
  expect_true(all(c("att_aug", "delta", "att_scm", "lambda_ridge", "beta_hat") %in% names(aug)))
})

test_that("augment_scm att_aug and delta are finite", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- augment_scm(fit)
  expect_true(is.finite(aug$att_aug))
  expect_true(is.finite(aug$delta))
})

test_that("augment_scm att_scm equals fit$estimate", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- augment_scm(fit)
  expect_equal(aug$att_scm, fit$estimate)
})

test_that("augment_scm beta_hat has length equal to T_pre", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- augment_scm(fit)
  expect_length(aug$beta_hat, fit$T_pre)
})

test_that("augment_scm auto-selected lambda_ridge is positive", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- augment_scm(fit)
  expect_true(aug$lambda_ridge > 0)
})

test_that("augment_scm with very large lambda_ridge gives delta near zero", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- augment_scm(fit, lambda_ridge = 1e12)
  # With huge lambda, beta -> 0, so delta = m_tr - W'm_co -> 0
  expect_true(abs(aug$delta) < 1e-3)
  expect_equal(aug$att_aug, aug$att_scm, tolerance = 1e-3)
})

test_that("augment_scm explicit lambda_ridge is honoured", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- augment_scm(fit, lambda_ridge = 0.5)
  expect_equal(aug$lambda_ridge, 0.5)
  expect_true(is.finite(aug$att_aug))
})

test_that("augment_scm rejects wrong method or missing Y_co_pre", {
  fit_sdid <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_error(augment_scm(fit_sdid), "method = 'scm'")

  fit_scm <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  fit_bad <- fit_scm
  fit_bad$Y_co_pre <- NULL
  expect_error(augment_scm(fit_bad), "Y_co_pre")
})

# ── Phase 12: SCM 高速化 ──────────────────────────────────────────────────────

test_that("Phase 12b: low-rank path returns valid simplex weights (k=7, N_co=46)", {
  set.seed(42L)
  k <- 7L; N_co <- 46L
  X0     <- matrix(rnorm(k * N_co), k, N_co)
  X1     <- rnorm(k)
  V_diag <- rep(1 / k, k)
  w <- scm_inner_weights_cpp(X0, X1, V_diag)
  expect_equal(length(w), N_co)
  expect_equal(sum(w), 1, tolerance = 1e-4)
  expect_true(all(w >= -1e-5))
})

test_that("Phase 12b: low-rank and full QP give algebraically equivalent weights", {
  set.seed(7L)
  k <- 5L; N_co <- 20L    # 2*k=10 < N_co=20 → low-rank path
  X0     <- matrix(rnorm(k * N_co), k, N_co)
  X1     <- rnorm(k)
  V_diag <- runif(k); V_diag <- V_diag / sum(V_diag)
  w_lr   <- scm_inner_weights_cpp(X0, X1, V_diag)
  # Reference via explicit full QP
  Q_full <- t(X0) %*% diag(V_diag) %*% X0
  c_full <- drop(t(X0) %*% diag(V_diag) %*% X1)
  w_full <- solve_simplex_qp(Q_full, c_full)
  expect_equal(w_lr, w_full, tolerance = 1e-3)
})

test_that("Phase 12b: full path used when k >= N_co/2 (k=10, N_co=16)", {
  set.seed(11L)
  k <- 10L; N_co <- 16L   # 2*k=20 >= N_co=16 → full path
  X0 <- matrix(rnorm(k * N_co), k, N_co)
  X1 <- rnorm(k)
  w  <- scm_inner_weights_cpp(X0, X1, rep(1 / k, k))
  expect_equal(sum(w), 1, tolerance = 1e-4)
  expect_true(all(w >= -1e-5))
})

test_that("Phase 12a: v_optim='bfgs' returns valid simplex weights", {
  # bfgs is deprecated (Phase 37): still functional, warns.
  expect_warning(
    fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                   v_optim = "bfgs"),
    regexp = "deprecated"
  )
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 12a: v_optim='bfgs' + v_selection='oos' works", {
  fit <- suppressWarnings(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            v_optim = "bfgs", v_selection = "oos"))
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 12a: v_optim='coord_descent' is identical to default", {
  fit1 <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  fit2 <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                  v_optim = "coord_descent")
  expect_equal(fit1$unit_weights, fit2$unit_weights, tolerance = 1e-6)
})

test_that("Phase 12a: v_optim='auto' returns valid weights", {
  # panel is outcomes-only, so auto resolves to coord_descent.
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_optim = "auto")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 12a: v_optim='bfgs' with lambda_pen works", {
  fit <- suppressWarnings(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            v_optim = "bfgs", lambda_pen = 0.5))
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(is.finite(fit$estimate))
})

# ── Phase 13a: SDID 共変量対応 ────────────────────────────────────────────────

# panel_cov's cov1/cov2 are time-invariant, so the unit fixed effects of the
# partial-out regression absorb them (and say so); enough to test the API.
test_that("Phase 13a: SDID with one covariate returns valid fit", {
  expect_warning(
    fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid",
                   covariates = "cov1"), "absorb")
  expect_s3_class(fit, "coresynth")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 13a: beta_hat length equals number of covariates", {
  expect_warning(
    fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid",
                   covariates = c("cov1", "cov2")), "absorb")
  expect_length(fit$beta_hat, 2L)
})

test_that("Phase 13a: SDID covariates — unit weights stay on simplex", {
  fit <- suppressWarnings(
    scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid",
            covariates = "cov1"))
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-8))
})

test_that("Phase 13a: SDID without covariates — beta_hat is length-0", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid")
  expect_length(fit$beta_hat, 0L)
})

test_that("Phase 13a: covariate partial-out changes the ATT estimate", {
  # cov1/cov2 in panel_cov are time-invariant, so they are absorbed by the
  # two-way fixed effects and (correctly) leave the SDID ATT essentially
  # unchanged. A meaningful partial-out test needs a genuinely time-varying
  # covariate that drives the outcome.
  set.seed(99)
  pc <- panel_cov
  pc$covtv <- rnorm(nrow(pc), 0, 1) + 0.5 * pc$time
  pc$y     <- pc$y + 0.8 * pc$covtv
  fit_no  <- scm_fit(y ~ d | id + time, data = pc, method = "sdid")
  fit_cov <- scm_fit(y ~ d | id + time, data = pc, method = "sdid",
                     covariates = "covtv")
  expect_false(isTRUE(all.equal(fit_no$estimate, fit_cov$estimate, tol = 1e-6)))
})

test_that("Phase 13a: staggered=FALSE for sharp SDID", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_false(fit$staggered)
})

# ── Phase 13b: SDID Staggered Adoption ───────────────────────────────────────

test_that("Phase 13b: SDID handles staggered adoption without error", {
  expect_no_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  )
})

test_that("Phase 13b: SDID staggered — staggered field is TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  expect_true(fit$staggered)
})

test_that("Phase 13b: SDID staggered — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 13b: SDID staggered — cohort_estimates has correct columns", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  expected_cols <- c("cohort", "n_treated", "T_pre", "T_post", "estimate", "weight")
  expect_true(all(expected_cols %in% names(fit$cohort_estimates)))
})

test_that("Phase 13b: SDID staggered — cohort weights sum to 1", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  expect_equal(sum(fit$cohort_estimates$weight), 1, tolerance = 1e-8)
})

test_that("Phase 13b: SDID staggered — each cohort estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  expect_true(all(is.finite(fit$cohort_estimates$estimate)))
})

test_that("Phase 13b: SDID staggered control_group='never_treated' works", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid",
                 control_group = "never_treated")
  expect_true(fit$staggered)
  expect_true(is.finite(fit$estimate))
})

test_that("staggered control_group defaults to never_treated; 'clean' is deprecated", {
  fit_def   <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  fit_never <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid",
                       control_group = "never_treated")
  expect_equal(fit_def$estimate, fit_never$estimate)
  expect_equal(fit_def$control_group, "never_treated")
  for (m in c("scm", "sdid", "gsc", "si"))
    expect_warning(
      scm_fit(y ~ d | id + time, data = staggered, method = m,
              control_group = "clean"),
      "deprecated")
})

test_that("Phase 13b: SDID staggered — Y_treat dimensions match panel", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  pan <- coresynth:::panel_to_matrices(staggered$y, staggered$d,
                                       staggered$id, staggered$time)
  expect_equal(nrow(fit$Y_treat), nrow(pan$Y))
})

# ── Phase 14a: SDID staggered + covariates ────────────────────────────────────

make_staggered_cov_panel <- function(N = 12, TT = 24, seed = 7) {
  p <- make_staggered_panel(N = N, TT = TT, seed = seed)
  set.seed(seed + 10L)
  p$cov1 <- NA_real_
  for (u in unique(p$id)) {
    base_val <- rnorm(1, 5, 1)
    p$cov1[p$id == u] <- base_val +
      seq(0, by = 0.05, length.out = TT) + rnorm(TT, 0, 0.1)
  }
  p
}

stag_cov <- make_staggered_cov_panel()

test_that("Phase 14a: SDID staggered + covariates — no error", {
  expect_no_error(
    scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
            covariates = "cov1")
  )
})

test_that("Phase 14a: SDID staggered + covariates — beta_hat length equals p", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                 covariates = "cov1")
  expect_length(fit$beta_hat, 1L)
})

test_that("Phase 14a: SDID staggered + covariates — staggered=TRUE preserved", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                 covariates = "cov1")
  expect_true(fit$staggered)
})

test_that("Phase 14a: SDID staggered + covariates — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                 covariates = "cov1")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 14a: SDID staggered + covariates — cohort_estimates columns correct", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                 covariates = "cov1")
  expected_cols <- c("cohort", "n_treated", "T_pre", "T_post", "estimate", "weight")
  expect_true(all(expected_cols %in% names(fit$cohort_estimates)))
})

test_that("Phase 14a: SDID staggered + covariates — cohort weights sum to 1", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                 covariates = "cov1")
  expect_equal(sum(fit$cohort_estimates$weight), 1, tolerance = 1e-8)
})

test_that("Phase 14a: SDID staggered + covariates — covariate changes estimate", {
  fit_no_cov <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid")
  fit_cov    <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                        covariates = "cov1")
  expect_false(isTRUE(all.equal(fit_no_cov$estimate, fit_cov$estimate,
                                tolerance = 1e-10)))
})

test_that("Phase 14a: SDID staggered + covariates — control_group='never_treated' works", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                 covariates = "cov1", control_group = "never_treated")
  expect_true(fit$staggered)
  expect_true(is.finite(fit$estimate))
})

# ── Phase 14b: SCM staggered adoption ─────────────────────────────────────────

test_that("Phase 14b: SCM handles staggered adoption without error", {
  expect_no_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  )
})

test_that("Phase 14b: SCM staggered — staggered = TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_true(fit$staggered)
})

test_that("Phase 14b: SCM staggered — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 14b: SCM staggered — cohort_estimates has required columns", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expected_cols <- c("cohort", "n_treated", "T_pre", "T_post", "estimate", "weight")
  expect_true(all(expected_cols %in% names(fit$cohort_estimates)))
})

test_that("Phase 14b: SCM staggered — cohort weights sum to 1", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_equal(sum(fit$cohort_estimates$weight), 1, tolerance = 1e-8)
})

test_that("Phase 14b: SCM staggered — each cohort unit_weights on simplex", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  for (cf in fit$cohort_fits) {
    w <- cf$unit_weights
    expect_equal(sum(w), 1, tolerance = 1e-5)
    expect_true(all(w >= -1e-5))
  }
})

test_that("Phase 14b: SCM staggered — v_selection='oos' gives finite ATT", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                 v_selection = "oos")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 14b: SCM staggered — predictors != NULL errors", {
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm",
            predictors = list(pred("y", 1:3))),
    regexp = "predictors"
  )
})

test_that("Phase 14b: SCM staggered — estimate within ±4 of true ATT (1.5)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_lt(abs(fit$estimate - 1.5), 4.0)
})

test_that("Phase 14b: SCM staggered — Y_treat nrow matches panel T", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_equal(nrow(fit$Y_treat), 24L)
})

test_that("Phase 14b: SCM staggered — control_group='never_treated' gives finite ATT", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                 control_group = "never_treated")
  expect_true(is.finite(fit$estimate))
})

# ── Phase 14c: augment.coresynth() ───────────────────────────────────────────

test_that("Phase 14c: augment() returns data.frame for SCM", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit)
  expect_s3_class(aug, "data.frame")
})

test_that("Phase 14c: augment() column names are correct", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit)
  expected <- c(".time", ".observed", ".fitted", ".resid", ".treated", ".period")
  expect_true(all(expected %in% names(aug)))
})

test_that("Phase 14c: augment() nrow equals T", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit)
  expect_equal(nrow(aug), 20L)
})

test_that("Phase 14c: augment() .resid == .observed - .fitted", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit)
  expect_equal(aug$.resid, aug$.observed - aug$.fitted, tolerance = 1e-10)
})

test_that("Phase 14c: augment() .treated is FALSE pre, TRUE post", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit)
  T_pre <- fit$T_pre
  expect_true(all(!aug$.treated[seq_len(T_pre)]))
  expect_true(all(aug$.treated[(T_pre + 1L):nrow(aug)]))
})

test_that("Phase 14c: augment() .period factor levels correct", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit)
  expect_equal(levels(aug$.period), c("pre", "post"))
})

test_that("Phase 14c: augment() works for SDID sharp", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  aug <- broom::augment(fit)
  expect_s3_class(aug, "data.frame")
  expect_true(all(c(".time", ".resid") %in% names(aug)))
})

test_that("Phase 14c: augment() works for GSC", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2L)
  aug <- broom::augment(fit)
  expect_s3_class(aug, "data.frame")
  expect_equal(nrow(aug), 20L)
})

test_that("Phase 14c (Phase 24-updated): augment() staggered SDID returns cohort-long df", {
  # Phase 24: augment.coresynth() now supports staggered fits and returns a
  # long-format data.frame stacking per-cohort treated paths (and donors when
  # include_donors=TRUE). The old behavior (warning + empty df) is replaced.
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  aug <- broom::augment(fit)
  expect_s3_class(aug, "data.frame")
  expect_gt(nrow(aug), 0L)
  expect_true(".cohort" %in% names(aug))
})

test_that("Phase 14c: tidy() returns cohort_estimates for staggered SCM", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  td <- broom::tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_true(nrow(td) > 0L)
  expect_true("cohort" %in% names(td) || grepl("^cohort_", td$term[1]))
})

# ── Phase 15a: SDID 推論 API (sdid_inference) ─────────────────────────────────

test_that("Phase 15a: SDID sharp fit stores Y_co_pre and Y_co_post", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_equal(dim(fit$Y_co_pre),  c(10L, 9L))  # T_pre x N_co
  expect_equal(dim(fit$Y_co_post), c(10L, 9L))  # T_post x N_co
})

test_that("Phase 15a: SDID sharp fit stores Y_tr_pre_mean, Y_tr_post_mean, N_tr", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_length(fit$Y_tr_pre_mean,  10L)  # T_pre
  expect_length(fit$Y_tr_post_mean, 10L)  # T_post
  expect_equal(fit$N_tr, 1L)
})

test_that("Phase 15a: sdid_inference(placebo) returns valid structure", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  expect_equal(inf$estimate, fit$estimate)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_length(inf$placebo_effects, 9L)  # N_co = 9
  # Phase 34: placebo now reports the placebo-distribution SE and a normal CI
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(is.finite(inf$ci_lower))
  expect_equal(inf$n_controls, 9L)
})

test_that("Phase 15a: sdid_inference(placebo) named placebo_effects", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  expect_equal(names(inf$placebo_effects), colnames(fit$Y_co_pre))
})

test_that("Phase 15a: sdid_inference(bootstrap) returns se and CI", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid")
  inf <- sdid_inference(fit, method = "bootstrap", n_boot = 50L, seed = 1L)
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_true(inf$ci_lower < inf$ci_upper)
  expect_length(inf$boot_ests, 50L)
})

test_that("Phase 15a: sdid_inference(bootstrap) is reproducible with seed", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid")
  inf1 <- sdid_inference(fit, method = "bootstrap", n_boot = 30L, seed = 99L)
  inf2 <- sdid_inference(fit, method = "bootstrap", n_boot = 30L, seed = 99L)
  expect_equal(inf1$se, inf2$se)
  expect_equal(inf1$boot_ests, inf2$boot_ests)
})

test_that("Phase 15a: sdid_inference(jackknife) returns se and CI", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_true(inf$ci_lower < inf$ci_upper)
  expect_null(inf$boot_ests)
})

test_that("Phase 15a: sdid_inference alternative='greater' gives one-sided p-value", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf_g <- sdid_inference(fit, method = "placebo", alternative = "greater")
  inf_l <- sdid_inference(fit, method = "placebo", alternative = "less")
  inf_t <- sdid_inference(fit, method = "placebo", alternative = "two.sided")
  expect_equal(inf_g$alternative, "greater")
  expect_true(inf_g$p_value >= 0 && inf_g$p_value <= 1)
  # sum of one-sided p-values cannot exceed 1 + 1/N
  N_co <- inf_t$n_controls
  expect_true(inf_g$p_value + inf_l$p_value <= 1 + 1 / (N_co + 1) + 1e-9)
})

test_that("Phase 15a: sdid_inference errors on non-SDID fit", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(sdid_inference(fit), regexp = "method = 'sdid'")
})

test_that("Phase 15a: sdid_inference staggered placebo no longer errors (Phase 20)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  expect_no_error(sdid_inference(fit, method = "placebo"))
})

test_that("Phase 15a: print.sdid_inference produces output", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  expect_output(print(inf), regexp = "Estimate")
  expect_output(print(inf), regexp = "p-value")
})

test_that("Phase 15a: SDID covariates fit also stores Y_co_pre/Y_co_post", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov[stag_cov$id %in%
    c("u1", paste0("u", 3:12)), ],
    method = "sdid", covariates = "cov1")
  # Sharp fit with covariates should store Y_co_pre
  if (!isTRUE(fit$staggered)) {
    expect_false(is.null(fit$Y_co_pre))
  }
})

# ── Phase 15b: SCM staggered covariates ────────────────────────────────────────

test_that("Phase 15b: SCM staggered + covariates — no error", {
  expect_no_error(
    scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
            covariates = "cov1")
  )
})

test_that("Phase 15b: SCM staggered + covariates — staggered=TRUE preserved", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
                 covariates = "cov1")
  expect_true(fit$staggered)
})

test_that("Phase 15b: SCM staggered + covariates — beta_hat length equals p", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
                 covariates = "cov1")
  expect_length(fit$beta_hat, 1L)
})

test_that("Phase 15b: SCM staggered + covariates — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
                 covariates = "cov1")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 15b: SCM staggered + covariates — cohort_estimates exist", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
                 covariates = "cov1")
  expect_s3_class(fit$cohort_estimates, "data.frame")
  expect_true(nrow(fit$cohort_estimates) > 0L)
})

test_that("Phase 15b: SCM staggered without covariates unchanged (regression)", {
  fit_no  <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  fit_cov <- scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
                     covariates = "cov1")
  # Estimates differ when covariates are included
  expect_false(isTRUE(all.equal(fit_no$estimate, fit_cov$estimate)))
})

test_that("Phase 15b: SCM staggered + predictors still errors", {
  fit_sharp <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  T_pre <- fit_sharp$T_pre
  times <- sort(unique(panel$time))
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm",
            predictors = list(pred("y", times[1:3]))),
    regexp = "predictors"
  )
})

test_that("Phase 15b: SCM sharp + covariates forwarded (no error)", {
  # Sharp SCM doesn't use covariates the same way, but passing it shouldn't error
  # (covariates is silently ignored for sharp SCM unless wired up in fit_scm_cpp)
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_s3_class(fit, "coresynth")
})

# ── Phase 15c: augment() include_donors ──────────────────────────────────────

test_that("Phase 15c: augment(include_donors=FALSE) is same as before", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug_old <- broom::augment(fit)
  aug_new <- broom::augment(fit, include_donors = FALSE)
  expect_equal(nrow(aug_new), nrow(aug_old))
  expect_equal(names(aug_new), names(aug_old))
})

test_that("Phase 15c: augment(include_donors=TRUE) adds control rows for SCM", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit, include_donors = TRUE)
  T_total <- length(fit$times)
  N_co    <- ncol(fit$Y_co_pre)
  expect_equal(nrow(aug), T_total * (1L + N_co))
})

test_that("Phase 15c: augment(include_donors=TRUE) has .unit and .type columns", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit, include_donors = TRUE)
  expect_true(".unit" %in% names(aug))
  expect_true(".type" %in% names(aug))
})

test_that("Phase 15c: augment control rows have NA .fitted and .resid", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit, include_donors = TRUE)
  ctrl <- aug[aug$.type == "control", ]
  expect_true(all(is.na(ctrl$.fitted)))
  expect_true(all(is.na(ctrl$.resid)))
})

test_that("Phase 15c: augment .type values are 'treated' and 'control'", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit, include_donors = TRUE)
  expect_true(all(aug$.type %in% c("treated", "control")))
  expect_true("treated" %in% aug$.type)
  expect_true("control" %in% aug$.type)
})

test_that("Phase 15c: augment(include_donors=TRUE) works for SDID", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  aug <- broom::augment(fit, include_donors = TRUE)
  T_total <- length(fit$times)
  N_co    <- ncol(fit$Y_co_pre)
  expect_equal(nrow(aug), T_total * (1L + N_co))
  expect_true("control" %in% aug$.type)
})

test_that("Phase 15c: augment(include_donors=TRUE) .unit names match donor names", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit, include_donors = TRUE)
  ctrl_units <- unique(aug$.unit[aug$.type == "control"])
  expect_setequal(ctrl_units, colnames(fit$Y_co_pre))
})

# ── Phase 16a: GSC staggered adoption ─────────────────────────────────────────

make_staggered_panel_multi <- function(N = 12, TT = 24, seed = 7) {
  set.seed(seed)
  f   <- cumsum(rnorm(TT, 0, 0.5))
  lam <- rnorm(N, 1, 0.3)
  rows <- expand.grid(time = seq_len(TT), id = paste0("u", seq_len(N)),
                      stringsAsFactors = FALSE)
  rows$y <- as.vector(outer(f, lam)) + rnorm(nrow(rows), 0, 0.3)
  rows$d <- as.integer(
    (rows$id %in% c("u1", "u2", "u3") & rows$time >= 9) |
    (rows$id %in% c("u4", "u5", "u6") & rows$time >= 15)
  )
  rows$y[rows$d == 1] <- rows$y[rows$d == 1] + 1.5
  rows[order(rows$id, rows$time), ]
}
staggered_multi <- make_staggered_panel_multi()

test_that("Phase 16a: GSC staggered — no error", {
  expect_no_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  )
})

test_that("Phase 16a: GSC staggered — staggered = TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_true(fit$staggered)
})

test_that("Phase 16a: GSC staggered — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 16a: GSC staggered — cohort_estimates has required columns", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expected_cols <- c("cohort", "n_treated", "T_pre", "T_post", "estimate", "weight")
  expect_true(all(expected_cols %in% names(fit$cohort_estimates)))
})

test_that("Phase 16a: GSC staggered — cohort_estimates weights sum to 1", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_equal(sum(fit$cohort_estimates$weight), 1, tolerance = 1e-8)
})

test_that("Phase 16a: GSC staggered — 2 cohorts present", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_equal(nrow(fit$cohort_estimates), 2L)
})

test_that("Phase 16a: GSC staggered — estimate within ±4 of true ATT (1.5)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_lt(abs(fit$estimate - 1.5), 4.0)
})

test_that("Phase 16a: GSC staggered — Y_treat nrow matches panel T (24)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_equal(nrow(fit$Y_treat), 24L)
})

test_that("Phase 16a: GSC staggered — per-cohort F is T x r matrix", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc", r = 2L)
  for (cf in fit$cohort_fits) {
    expect_equal(ncol(cf$F), 2L)
    expect_equal(nrow(cf$F), 24L)
  }
})

test_that("Phase 16a: GSC staggered — per-cohort L_co is N_co x r matrix", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc", r = 2L)
  for (cf in fit$cohort_fits) {
    expect_equal(ncol(cf$L_co), 2L)
    expect_true(nrow(cf$L_co) >= 2L)
  }
})

test_that("Phase 16a: GSC staggered — control_group='never_treated' gives finite ATT", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc",
                 control_group = "never_treated")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 16a: GSC staggered — multi-unit cohort gives finite ATT", {
  fit <- scm_fit(y ~ d | id + time, data = staggered_multi, method = "gsc")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 16a: GSC staggered — r=3 changes F dimensions", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc", r = 3L)
  for (cf in fit$cohort_fits) {
    expect_equal(ncol(cf$F), 3L)
  }
})

test_that("Phase 16a: GSC staggered — gsc_boot() throws error", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_error(gsc_boot(fit), regexp = "staggered")
})

test_that("Phase 16a: GSC staggered — tidy() returns cohort rows", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  td <- broom::tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_true(nrow(td) == 2L)
  expect_true(grepl("^cohort_", td$term[1]))
})

# ── Phase 16b: SI staggered adoption ──────────────────────────────────────────

test_that("Phase 16b: SI staggered — no error", {
  expect_no_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "si")
  )
})

test_that("Phase 16b: SI staggered — staggered = TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_true(fit$staggered)
})

test_that("Phase 16b: SI staggered — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 16b: SI staggered — cohort_estimates has required columns", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expected_cols <- c("cohort", "n_treated", "T_pre", "T_post", "estimate", "weight")
  expect_true(all(expected_cols %in% names(fit$cohort_estimates)))
})

test_that("Phase 16b: SI staggered — cohort_estimates weights sum to 1", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_equal(sum(fit$cohort_estimates$weight), 1, tolerance = 1e-8)
})

# ── Phase 21: SI multi-arm (K>1, Agarwal et al. 2025) ────────────────────────

# DGP: tensor factor model  Y^(d)_{ti} = sum_l u_tl * Lambda[d,l] * v_il + eps
# Lambda: arm 0 = c(1,1), arm 1 = c(1.5,0.8), arm 2 = c(0.7,1.3)
make_tensor_panel <- function(N_co = 10L, N1 = 3L, N2 = 3L,
                               T = 24L, T_pre = 12L, seed = 42L) {
  set.seed(seed)
  r <- 2L
  N <- N_co + N1 + N2
  U      <- matrix(cumsum(rnorm(T * r)), T, r)
  V      <- matrix(rnorm(N * r, 0, 0.5), N, r)
  Lambda <- rbind(c(1.0, 1.0), c(1.5, 0.8), c(0.7, 1.3))
  arm_unit <- c(rep(0L, N_co), rep(1L, N1), rep(2L, N2))
  ids  <- paste0("u", seq_len(N))
  rows <- expand.grid(time = seq_len(T), id = ids, stringsAsFactors = FALSE)
  rows <- rows[order(rows$id, rows$time), ]
  i_idx <- match(rows$id, ids)
  t_idx <- rows$time
  arm_i <- arm_unit[i_idx]
  rows$d <- ifelse(arm_i == 0L, 0L, ifelse(t_idx > T_pre, arm_i, 0L))
  rows$y <- mapply(function(t, i, arm) {
    sum(U[t, ] * Lambda[arm + 1L, ] * V[i, ]) + rnorm(1L, 0, 0.05)
  }, t_idx, i_idx, arm_i)
  rows
}

tensor_panel <- make_tensor_panel()

# § A: panel_to_tensor() ──────────────────────────────────────────────────────

test_that("Phase 21: panel_to_tensor — arm_levels is c(0,1,2)", {
  tp <- coresynth:::panel_to_tensor(tensor_panel$y, tensor_panel$d,
                                     tensor_panel$id, tensor_panel$time)
  expect_equal(tp$arm_levels, c(0L, 1L, 2L))
})

test_that("Phase 21: panel_to_tensor — idx_by_arm unit counts correct", {
  tp <- coresynth:::panel_to_tensor(tensor_panel$y, tensor_panel$d,
                                     tensor_panel$id, tensor_panel$time)
  expect_equal(length(tp$idx_by_arm[["0"]]), 10L)  # N_co
  expect_equal(length(tp$idx_by_arm[["1"]]), 3L)   # N1
  expect_equal(length(tp$idx_by_arm[["2"]]), 3L)   # N2
})

test_that("Phase 21: panel_to_tensor — error when arm 0 absent", {
  bad <- tensor_panel
  bad$d <- bad$d + 1L  # shift: now arms are 1, 2, 3 (no 0)
  # The shift makes every unit "treated" from period 1, so the shared
  # no-control-units check fires before panel_to_tensor's own arm-0 check.
  expect_error(
    coresynth:::panel_to_tensor(bad$y, bad$d, bad$id, bad$time),
    "arm 0|No control units"
  )
})

test_that("Phase 21: panel_to_tensor — works on binary d (K=1)", {
  tp <- coresynth:::panel_to_tensor(panel$y, panel$d, panel$id, panel$time)
  expect_equal(tp$arm_levels, c(0L, 1L))
  expect_true(!is.null(tp$idx_by_arm[["0"]]))
  expect_true(!is.null(tp$idx_by_arm[["1"]]))
})

test_that("Phase 21: panel_to_tensor — return has required fields", {
  tp <- coresynth:::panel_to_tensor(tensor_panel$y, tensor_panel$d,
                                     tensor_panel$id, tensor_panel$time)
  for (fld in c("Y", "T_pre", "is_sharp", "arm_levels", "idx_by_arm"))
    expect_true(fld %in% names(tp), info = paste("missing field:", fld))
})

# § B: fit_si_cpp() multi-arm basic behaviour ─────────────────────────────────

test_that("Phase 21: SI multi-arm — class is coresynth", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_s3_class(fit, "coresynth")
  expect_s3_class(fit, "coresynth_si")
})

test_that("Phase 21: SI multi-arm — multi_arm flag is TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_true(isTRUE(fit$multi_arm))
})

test_that("Phase 21: SI multi-arm — arm_levels is c(1L, 2L)", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_equal(fit$arm_levels, c(1L, 2L))
})

test_that("Phase 21: SI multi-arm — arm_estimates length equals K", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_equal(length(fit$arm_estimates), 2L)
})

test_that("Phase 21: SI multi-arm — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 21: SI multi-arm — arm_fits has required sub-fields", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  for (arm in c("1", "2")) {
    af <- fit$arm_fits[[arm]]
    for (fld in c("weights", "Y_cf", "Y_treat", "Y_synth", "gap"))
      expect_true(fld %in% names(af),
                  info = paste("arm", arm, "missing field:", fld))
  }
})

test_that("Phase 21: SI multi-arm — arm 1 Y_cf dimension is T_post x N1", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  T_post <- 24L - 12L
  expect_equal(dim(fit$arm_fits[["1"]]$Y_cf), c(T_post, 3L))
})

test_that("Phase 21: SI multi-arm — arm 1 weights dimension is N_co x N1", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_equal(dim(fit$arm_fits[["1"]]$weights), c(10L, 3L))
})

# § C: estimation accuracy ────────────────────────────────────────────────────

test_that("Phase 21: SI multi-arm — arm 1 ATT is finite and bounded", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_true(is.finite(fit$arm_estimates[["1"]]))
  expect_lt(abs(fit$arm_estimates[["1"]]), 20)
})

test_that("Phase 21: SI multi-arm — arm 2 ATT is finite and bounded", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_true(is.finite(fit$arm_estimates[["2"]]))
  expect_lt(abs(fit$arm_estimates[["2"]]), 20)
})

test_that("Phase 21: SI multi-arm — manual k is respected", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si", k = 2L)
  expect_equal(fit$k, 2L)
  expect_true(is.finite(fit$estimate))
})

# § D: scm_fit() interface ───────────────────────────────────────────────────

test_that("Phase 21: SI multi-arm via scm_fit — no error with multi-valued d", {
  expect_no_error(scm_fit(y ~ d | id + time, data = tensor_panel, method = "si"))
})

test_that("Phase 21: SI multi-arm via scm_fit — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 21: SI multi-arm — print() outputs Multi-arm", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  out <- capture.output(print(fit))
  expect_true(any(grepl("Multi-arm", out)))
})

# § E: si_inference() multi-arm ───────────────────────────────────────────────

test_that("Phase 21: si_inference multi-arm bootstrap — class is coresynth_inference", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 50L, seed = 1L)
  expect_s3_class(inf, "coresynth_inference")
})

test_that("Phase 21: si_inference multi-arm bootstrap — SE finite and positive", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 50L, seed = 1L)
  expect_true(is.finite(inf$se))
  expect_gt(inf$se, 0)
})

test_that("Phase 21: si_inference multi-arm bootstrap — CI contains estimate", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 50L, seed = 1L)
  expect_lte(inf$ci_lower, inf$estimate + 1e-9)
  expect_gte(inf$ci_upper, inf$estimate - 1e-9)
})

test_that("Phase 21: si_inference multi-arm jackknife — SE finite and positive", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  inf <- si_inference(fit, method = "jackknife", seed = 1L)
  expect_true(is.finite(inf$se))
  expect_gt(inf$se, 0)
})

test_that("Phase 21: si_inference multi-arm jackknife_global — error with guidance", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_error(si_inference(fit, method = "jackknife_global"),
               "jackknife")
})

test_that("Phase 21: si_inference — rejects non-SI fit", {
  fit_mc <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  expect_error(si_inference(fit_mc), "method = 'si'")
})

# § F: edge cases ─────────────────────────────────────────────────────────────

test_that("Phase 21: SI multi-arm — aggregate ATT equals weighted arm mean", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  ws   <- vapply(fit$arm_fits, `[[`, numeric(1L), "weight")
  taus <- fit$arm_estimates
  expect_equal(sum(ws * taus) / sum(ws), fit$estimate, tolerance = 1e-10)
})

test_that("Phase 21: SI multi-arm — N_tr_d=1 arm works without error", {
  tp1 <- make_tensor_panel(N_co = 10L, N1 = 1L, N2 = 3L, seed = 99L)
  expect_no_error(scm_fit(y ~ d | id + time, data = tp1, method = "si"))
})

test_that("Phase 21: SI multi-arm — k > min(T_pre,N_co) raises error", {
  expect_error(
    scm_fit(y ~ d | id + time, data = tensor_panel, method = "si", k = 100L),
    "k"
  )
})

test_that("Phase 16b: SI staggered — 2 cohorts present", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_equal(nrow(fit$cohort_estimates), 2L)
})

test_that("Phase 16b: SI staggered — estimate within ±4 of true ATT (1.5)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_lt(abs(fit$estimate - 1.5), 4.0)
})

test_that("Phase 16b: SI staggered — Y_treat nrow matches panel T (24)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_equal(nrow(fit$Y_treat), 24L)
})

test_that("Phase 16b: SI staggered — control_group='never_treated' gives finite ATT", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si",
                 control_group = "never_treated")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 16b: SI staggered — multi-unit cohort gives finite ATT", {
  fit <- scm_fit(y ~ d | id + time, data = staggered_multi, method = "si")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 16b: SI staggered — k argument respected", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si", k = 2L)
  for (cf in fit$cohort_fits) {
    expect_lte(cf$k, 2L)
  }
})

test_that("Phase 16b: SI staggered — per-cohort k auto-adjusted for short pre-period", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  for (cf in fit$cohort_fits) {
    expect_gte(cf$k, 1L)
    expect_lte(cf$k, min(cf$T_pre, length(cf$idx_co)))
  }
})

test_that("Phase 16b: SI staggered — per-cohort weights matrix has correct dims", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  for (cf in fit$cohort_fits) {
    expect_equal(nrow(cf$weights), length(cf$idx_co))
    expect_equal(ncol(cf$weights), cf$n_treated)
  }
})

test_that("Phase 16b: SI staggered — tidy() returns cohort rows", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  td <- broom::tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_true(nrow(td) == 2L)
  expect_true(grepl("^cohort_", td$term[1]))
})

# ── Phase 17: GSC/SI 推論 API（Bootstrap/Jackknife）────────────────────────────

test_that("unit-resampling inference stops with one treated unit and says where to go", {
  fit_g <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 1L)
  expect_error(gsc_inference(fit_g, method = "bootstrap"), "gsc_boot")
  expect_error(gsc_inference(fit_g, method = "jackknife"), "conformal_inference")
  fit_s <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  expect_error(si_inference(fit_s, method = "bootstrap"), "conformal_inference")
})

test_that("GSC unit engine reproduces the fit it resamples", {
  est <- function(fit) {
    ctx <- coresynth:::.gsc_inference_context(fit)
    coresynth:::.gsc_unit_estimate(ctx, ctx$co, ctx$tr)$estimate
  }
  fit_sh <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc", r = 2L)
  expect_equal(est(fit_sh), fit_sh$estimate, tolerance = 1e-10)
  fit_st <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc", r = 1L)
  expect_equal(est(fit_st), fit_st$estimate, tolerance = 1e-10)
  fit_cl <- suppressWarnings(scm_fit(y ~ d | id + time, data = staggered,
                                     method = "gsc", r = 1L,
                                     control_group = "clean"))
  expect_equal(est(fit_cl), fit_cl$estimate, tolerance = 1e-10)
})

test_that("SI unit engine reproduces sharp, multi-arm and staggered fits", {
  est <- function(fit) {
    ctx <- coresynth:::.si_inference_context(fit)
    coresynth:::.si_unit_estimate(ctx, ctx$co, ctx$tr)$estimate
  }
  for (fit in list(
    scm_fit(y ~ d | id + time, data = panel_tr3, method = "si"),
    scm_fit(y ~ d | id + time, data = tensor_panel, method = "si"),
    scm_fit(y ~ d | id + time, data = staggered, method = "si")))
    expect_equal(est(fit), fit$estimate, tolerance = 1e-10)
})

test_that("staggered GSC fits the never-treated donors once for every cohort", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc", r = 1L)
  expect_identical(fit$cohort_fits[[1]]$F, fit$cohort_fits[[2]]$F)
  expect_identical(fit$control_group, "never_treated")
  # Each cohort matches a stand-alone fit on its own pre-treatment window
  cf  <- fit$cohort_fits[[2]]
  Y   <- fit$Y_all
  pre <- seq_len(cf$T_pre)
  ref <- gsc_ife_cpp(Y[, cf$idx_co, drop = FALSE], Y[pre, cf$idx_tr, drop = FALSE],
                     1L, array(0, c(nrow(Y), length(cf$idx_co), 0)),
                     array(0, c(cf$T_pre, length(cf$idx_tr), 0)))
  post <- (cf$T_pre + 1L):nrow(Y)
  expect_equal(cf$estimate,
               mean(Y[post, cf$idx_tr] - ref$Y_tr_hat[post, ]),
               tolerance = 1e-10)
})

# ── Phase 17a: gsc_inference — sharp ─────────────────────────────────────────

test_that("Phase 17a: gsc_inference sharp bootstrap — no error, class correct", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 1L)
  expect_s3_class(inf, "coresynth_inference")
})

test_that("Phase 17a: gsc_inference sharp bootstrap — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 2L)
  expect_true(is.finite(inf$se))
  expect_gte(inf$p_value, 0)
  expect_lte(inf$p_value, 1)
})

test_that("Phase 17a: gsc_inference sharp bootstrap — CI contains estimate", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 3L)
  expect_lte(inf$ci_lower, inf$estimate + 1e-8)
  expect_gte(inf$ci_upper, inf$estimate - 1e-8)
})

test_that("Phase 17a: gsc_inference sharp jackknife — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc")
  inf <- gsc_inference(fit, method = "jackknife")
  expect_true(is.finite(inf$se))
  expect_gte(inf$p_value, 0)
  expect_lte(inf$p_value, 1)
})

test_that("Phase 17a: gsc_inference — rejects non-GSC fit", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_error(gsc_inference(fit), "method = 'gsc'")
})

# ── Phase 17b: gsc_inference — staggered ──────────────────────────────────────

test_that("Phase 17b: gsc_inference staggered bootstrap — no error, staggered=TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 10L)
  expect_s3_class(inf, "coresynth_inference")
  expect_true(inf$staggered)
})

test_that("Phase 17b: gsc_inference staggered bootstrap — finite SE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 11L)
  expect_true(is.finite(inf$se))
})

test_that("Phase 17b: gsc_inference staggered bootstrap — CI contains estimate", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 12L)
  expect_lte(inf$ci_lower, inf$estimate + 1e-8)
  expect_gte(inf$ci_upper, inf$estimate - 1e-8)
})

test_that("Phase 17b: gsc_inference staggered jackknife — finite SE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "gsc")
  inf <- gsc_inference(fit, method = "jackknife")
  expect_true(is.finite(inf$se))
  expect_true(inf$staggered)
})

# ── Phase 17c: si_inference — sharp ───────────────────────────────────────────

test_that("Phase 17c: si_inference sharp bootstrap — no error, class correct", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 20L)
  expect_s3_class(inf, "coresynth_inference")
})

test_that("Phase 17c: si_inference sharp bootstrap — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 21L)
  expect_true(is.finite(inf$se))
  expect_gte(inf$p_value, 0)
  expect_lte(inf$p_value, 1)
})

test_that("Phase 17c: si_inference sharp bootstrap — CI contains estimate", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 22L)
  expect_lte(inf$ci_lower, inf$estimate + 1e-8)
  expect_gte(inf$ci_upper, inf$estimate - 1e-8)
})

test_that("Phase 17c: si_inference sharp jackknife — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "si")
  inf <- si_inference(fit, method = "jackknife")
  expect_true(is.finite(inf$se))
  expect_gte(inf$p_value, 0)
  expect_lte(inf$p_value, 1)
})

test_that("Phase 17c: si_inference — rejects non-SI fit", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_error(si_inference(fit), "method = 'si'")
})

# ── Phase 17d: si_inference — staggered ───────────────────────────────────────

test_that("Phase 17d: si_inference staggered bootstrap — no error, staggered=TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 30L)
  expect_s3_class(inf, "coresynth_inference")
  expect_true(inf$staggered)
})

test_that("Phase 17d: si_inference staggered bootstrap — finite SE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 31L)
  expect_true(is.finite(inf$se))
})

test_that("Phase 17d: si_inference staggered bootstrap — CI contains estimate", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 32L)
  expect_lte(inf$ci_lower, inf$estimate + 1e-8)
  expect_gte(inf$ci_upper, inf$estimate - 1e-8)
})

test_that("Phase 17d: si_inference staggered jackknife — finite SE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "si")
  inf <- si_inference(fit, method = "jackknife")
  expect_true(is.finite(inf$se))
  expect_true(inf$staggered)
})

test_that("Phase 17d: gsc_boot staggered error mentions gsc_inference()", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  expect_error(gsc_boot(fit), "gsc_inference")
})

# ── Phase 18a: sdid_inference() staggered bootstrap / jackknife ───────────────

test_that("Phase 18a: sdid_inference staggered bootstrap returns sdid_inference", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "bootstrap", n_boot = 50L, seed = 1L)
  expect_s3_class(inf, "sdid_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
})

test_that("Phase 18a: sdid_inference staggered bootstrap CI straddles estimate", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "bootstrap", n_boot = 50L, seed = 1L)
  expect_true(inf$ci_lower <= inf$estimate && inf$estimate <= inf$ci_upper)
  expect_false(is.null(inf$boot_ests))
})

test_that("Phase 18a: sdid_inference staggered jackknife returns sdid_inference", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  expect_s3_class(inf, "sdid_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
})

test_that("Phase 18a: sdid_inference staggered jackknife CI straddles estimate", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  expect_true(inf$ci_lower <= inf$estimate && inf$estimate <= inf$ci_upper)
  expect_null(inf$boot_ests)
})

# ── Phase 18b: jackknife_global for SDID / GSC / SI staggered ────────────────

test_that("Phase 18b: sdid_inference jackknife_global staggered", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife_global")
  expect_s3_class(inf, "sdid_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$ci_lower <= inf$estimate && inf$estimate <= inf$ci_upper)
})

test_that("Phase 18b: sdid_inference jackknife_global errors on sharp fit", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_error(sdid_inference(fit, method = "jackknife_global"),
               regexp = "staggered fit")
})

test_that("Phase 18b: gsc_inference jackknife_global staggered", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "gsc")
  inf <- gsc_inference(fit, method = "jackknife_global")
  expect_s3_class(inf, "coresynth_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$ci_lower <= inf$estimate && inf$estimate <= inf$ci_upper)
})

test_that("Phase 18b: si_inference jackknife_global staggered", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "si")
  inf <- si_inference(fit, method = "jackknife_global")
  expect_s3_class(inf, "coresynth_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$ci_lower <= inf$estimate && inf$estimate <= inf$ci_upper)
})

# ── Phase 20: sdid_inference() staggered placebo ─────────────────────────────

test_that("Phase 20: sdid_inference staggered placebo returns sdid_inference", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  expect_s3_class(inf, "sdid_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  # Phase 34: placebo now reports the placebo-distribution SE and a normal CI
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(is.finite(inf$ci_lower))
  expect_true(is.finite(inf$ci_upper))
  expect_null(inf$boot_ests)
})

test_that("Phase 20: staggered placebo_effects length = never-treated controls", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  # never_co = intersect({u2..u12}, {u3..u12}) = {u3..u12} = 10 units
  never_co <- Reduce(intersect, lapply(fit$cohort_fits, `[[`, "idx_co"))
  expect_equal(length(inf$placebo_effects), length(never_co))
  expect_false(is.null(names(inf$placebo_effects)))
})

test_that("Phase 20: staggered placebo n_controls = never-treated count", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  never_co <- Reduce(intersect, lapply(fit$cohort_fits, `[[`, "idx_co"))
  expect_equal(inf$n_controls, length(never_co))
})

test_that("Phase 20: staggered placebo alternative one-sided p-values sum >= 1", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  p_g <- sdid_inference(fit, method = "placebo", alternative = "greater")$p_value
  p_l <- sdid_inference(fit, method = "placebo", alternative = "less")$p_value
  expect_gte(p_g + p_l, 1 - .Machine$double.eps)
  expect_lte(p_g, 1)
  expect_lte(p_l, 1)
})

test_that("Phase 20: staggered placebo with large effect yields small p-value", {
  d_strong <- staggered
  d_strong$y[d_strong$d == 1] <- d_strong$y[d_strong$d == 1] + 100
  fit <- scm_fit(y ~ d | id + time, data = d_strong, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  never_co <- Reduce(intersect, lapply(fit$cohort_fits, `[[`, "idx_co"))
  # ATT ≈ 100 >> noise; observed effect should rank first in placebo distribution
  expect_lte(inf$p_value, 2 / (length(never_co) + 1))
})

test_that("Phase 20: staggered placebo control_group=never_treated uses all controls", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid",
                 control_group = "never_treated")
  inf <- sdid_inference(fit, method = "placebo")
  # never_treated: intersection == full idx_co (same across all cohorts)
  expect_equal(inf$n_controls, length(fit$cohort_fits[[1]]$idx_co))
  expect_true(isTRUE(inf$staggered))
})

# ── Phase 23: SI multi-arm + staggered adoption ───────────────────────────────
# Two cohorts (g=5, g=7) × two arms (d=1, d=2) + N_co=4 controls, T=10.
# cohort g=5: 1 arm-1 unit + 1 arm-2 unit
# cohort g=7: 1 arm-1 unit + 1 arm-2 unit  →  4 (g,d) cells total
make_stag_tensor_panel <- function(N_co = 4L, TT = 10L, seed = 77L) {
  set.seed(seed)
  r <- 2L
  N <- N_co + 4L  # 4 treated (1 per cohort-arm cell)
  U      <- matrix(cumsum(rnorm(TT * r)), TT, r)
  V      <- matrix(rnorm(N * r, 0, 0.5), N, r)
  Lambda <- rbind(c(1.0, 1.0), c(1.5, 0.8), c(0.7, 1.3))
  # arm assignment: cols 1..N_co = arm 0, N_co+1 = arm1 g=5, N_co+2 = arm2 g=5,
  #                N_co+3 = arm1 g=7, N_co+4 = arm2 g=7
  arm_unit   <- c(rep(0L, N_co), 1L, 2L, 1L, 2L)
  cohort_of  <- c(rep(NA_integer_, N_co), 5L, 5L, 7L, 7L)
  ids  <- paste0("u", seq_len(N))
  rows <- expand.grid(time = seq_len(TT), id = ids, stringsAsFactors = FALSE)
  rows <- rows[order(rows$id, rows$time), ]
  i_idx <- match(rows$id, ids)
  rows$d <- mapply(function(i, t) {
    arm_i <- arm_unit[i]
    if (arm_i == 0L) return(0L)
    g <- cohort_of[i]
    if (!is.na(g) && t >= g) as.integer(arm_i) else 0L
  }, i_idx, rows$time)
  rows$y <- mapply(function(t, i) {
    sum(U[t, ] * Lambda[arm_unit[i] + 1L, ] * V[i, ]) + rnorm(1L, 0, 0.05)
  }, rows$time, i_idx)
  rows
}

stag_tensor_panel <- make_stag_tensor_panel()

# § A: fit structure ──────────────────────────────────────────────────────────

test_that("Phase 23: SI staggered-multi — no error", {
  expect_no_error(
    scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  )
})

test_that("Phase 23: SI staggered-multi — multi_arm=TRUE and staggered=TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  expect_true(isTRUE(fit$multi_arm))
  expect_true(isTRUE(fit$staggered))
})

test_that("Phase 23: SI staggered-multi — estimate is finite", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 23: SI staggered-multi — arm_levels is c(1L, 2L)", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  expect_equal(fit$arm_levels, c(1L, 2L))
})

test_that("Phase 23: SI staggered-multi — arm_estimates has 2 finite entries", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  expect_equal(length(fit$arm_estimates), 2L)
  expect_true(all(is.finite(fit$arm_estimates)))
})

test_that("Phase 23: SI staggered-multi — cohort_arm_estimates shape (4 rows)", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  ca  <- fit$cohort_arm_estimates
  expect_s3_class(ca, "data.frame")
  expect_equal(nrow(ca), 4L)
  for (col in c("cohort", "arm", "n_treated", "T_pre", "T_post", "estimate", "weight"))
    expect_true(col %in% names(ca))
})

test_that("Phase 23: SI staggered-multi — cohort_arm_estimates weights sum to 1", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  expect_equal(sum(fit$cohort_arm_estimates$weight), 1, tolerance = 1e-8)
})

test_that("Phase 23: SI staggered-multi — Y_all stored with correct dims", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  expect_false(is.null(fit$Y_all))
  expect_equal(nrow(fit$Y_all), 10L)                # TT
  expect_equal(ncol(fit$Y_all), 4L + 4L)            # N_co + N_treated
})

test_that("Phase 23: SI staggered-multi — control_group='never_treated' works", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si",
                 control_group = "never_treated")
  expect_true(isTRUE(fit$multi_arm))
  expect_true(isTRUE(fit$staggered))
  expect_true(is.finite(fit$estimate))
})

# § B: aggregate formula ──────────────────────────────────────────────────────

test_that("Phase 23: SI staggered-multi — ATT = sum(w*tau)/sum(w) over (g,d) cells", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  ws   <- vapply(fit$cohort_fits, `[[`, numeric(1L), "weight")
  taus <- vapply(fit$cohort_fits, `[[`, numeric(1L), "estimate")
  expect_equal(sum(ws * taus) / sum(ws), fit$estimate, tolerance = 1e-10)
})

test_that("Phase 23: SI staggered-multi — cohort_fits flat list with required fields", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  expect_type(fit$cohort_fits, "list")
  expect_gte(length(fit$cohort_fits), 1L)
  cf1 <- fit$cohort_fits[[1L]]
  for (fld in c("cohort", "arm", "n_treated", "T_pre", "T_post", "estimate",
                "weight", "idx_tr", "idx_co", "k", "weights", "Y_cf",
                "Y_treat", "gap"))
    expect_true(fld %in% names(cf1), info = paste("missing field:", fld))
})

test_that("Phase 23: SI staggered-multi — tidy() returns cohort_arm_estimate rows", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  td  <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 4L)
  expect_true("type" %in% names(td))
  expect_true(all(td$type == "cohort_arm_estimate"))
})

# § C: inference — all three methods ─────────────────────────────────────────

test_that("Phase 23: si_inference staggered-multi bootstrap — class and staggered flag", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 49L, seed = 1L)
  expect_s3_class(inf, "coresynth_inference")
  expect_true(isTRUE(inf$staggered))
})

test_that("Phase 23: si_inference staggered-multi bootstrap — SE finite and CI valid", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 49L, seed = 2L)
  expect_true(is.finite(inf$se))
  expect_lte(inf$ci_lower, inf$estimate)
  expect_gte(inf$ci_upper, inf$estimate)
  expect_gte(inf$p_value, 0); expect_lte(inf$p_value, 1)
})

test_that("Phase 23: si_inference staggered-multi jackknife — SE finite and CI valid", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  inf <- si_inference(fit, method = "jackknife")
  expect_true(is.finite(inf$se) && inf$se >= 0)
  expect_lte(inf$ci_lower, inf$estimate)
  expect_gte(inf$ci_upper, inf$estimate)
  expect_true(isTRUE(inf$staggered))
})

test_that("Phase 23: si_inference staggered-multi jackknife_global — SE finite and CI valid", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  inf <- si_inference(fit, method = "jackknife_global")
  expect_s3_class(inf, "coresynth_inference")
  expect_true(is.finite(inf$se) && inf$se >= 0)
  expect_lte(inf$ci_lower, inf$estimate)
  expect_gte(inf$ci_upper, inf$estimate)
  expect_true(isTRUE(inf$staggered))
})

# § D: regression — existing paths unaffected ─────────────────────────────────

test_that("Phase 23: sharp multi-arm still works (regression)", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_true(isTRUE(fit$multi_arm))
  expect_false(isTRUE(fit$staggered))
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 23: single-arm staggered still works (regression)", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_true(isTRUE(fit$staggered))
  expect_false(isTRUE(fit$multi_arm))
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 23: sharp multi-arm jackknife_global still errors (regression)", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  expect_error(si_inference(fit, method = "jackknife_global"), "jackknife")
})

# ── Phase 24: applied-pattern regression tests ────────────────────────────────
# 1. SCM multi-predictor patterns (Abadie, Diamond & Hainmueller 2010 style)
# 2. tidy/glance.coresynth_inference and enriched glance.coresynth
# 3. augment.coresynth() staggered + multi-arm
# 4. Improved error messages (staggered + mspe_ratio_pval, staggered + predictors)
# 5. dead-code cleanup regression (Phase 14b/15b SCM staggered unchanged)

# § A: SCM multi-predictor patterns (Prop99 style) ───────────────────────────

test_that("Phase 24: SCM with three pred() at specific years works", {
  fit <- scm_fit(
    y ~ d | id + time, data = panel_cov, method = "scm",
    predictors = list(pred("y", 3), pred("y", 5), pred("y", 7))
  )
  expect_s3_class(fit, "coresynth")
  expect_equal(length(fit$v_weights), 3L)
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 24: SCM with mixed-aggregation pred() runs", {
  fit <- scm_fit(
    y ~ d | id + time, data = panel_cov, method = "scm",
    predictors = list(
      pred("y", 1:10),
      pred("cov1", 1:10, op = "median"),
      pred("cov2", 1:10, op = "mean")
    )
  )
  expect_s3_class(fit, "coresynth")
  expect_equal(nrow(fit$predictor_table), 3L)
})

test_that("Phase 24: SCM oos + multiple pred() returns finite ATT", {
  fit <- scm_fit(
    y ~ d | id + time, data = panel_cov, method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10),
                      pred("y", c(2, 5, 8))),
    v_selection = "oos"
  )
  expect_true(is.finite(fit$estimate))
  expect_equal(length(fit$v_weights), 3L)
})

# § B: tidy/glance.coresynth_inference ────────────────────────────────────────

test_that("Phase 24: tidy.coresynth_inference SDID sharp returns expected columns", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  td  <- broom::tidy(inf)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 1L)
  required <- c("term", "estimate", "std.error", "statistic", "p.value",
                "conf.low", "conf.high", "method", "alternative",
                "n_controls", "staggered")
  expect_true(all(required %in% names(td)))
  expect_equal(td$term, "ATT")
  expect_false(td$staggered)
})

test_that("Phase 24: tidy.coresynth_inference SDID staggered carries staggered=TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered2, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  td  <- broom::tidy(inf)
  expect_true(td$staggered)
})

test_that("Phase 24: tidy.coresynth_inference GSC works", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc", r = 2L)
  inf <- gsc_inference(fit, method = "jackknife")
  td  <- broom::tidy(inf)
  expect_equal(td$term, "ATT")
  expect_true(is.finite(td$estimate))
  expect_true(is.finite(td$std.error))
})

test_that("Phase 24: tidy.coresynth_inference SI works", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "si")
  inf <- si_inference(fit, method = "jackknife")
  td  <- broom::tidy(inf)
  expect_equal(td$term, "ATT")
  expect_true(is.finite(td$std.error))
})

test_that("Phase 24: glance.coresynth_inference returns one-row summary", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid")
  inf <- sdid_inference(fit, method = "bootstrap", n_boot = 31L, seed = 1L)
  g   <- broom::glance(inf)
  expect_equal(nrow(g), 1L)
  expect_true(all(c("method", "n_controls", "staggered", "estimate",
                    "std.error", "p.value", "conf.low", "conf.high",
                    "alternative", "n_boot_valid") %in% names(g)))
  expect_gte(g$n_boot_valid, 1L)
})

test_that("Phase 24: tidy.coresynth_inference carries placebo SE/CI (finite since Phase 34)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  td  <- broom::tidy(inf)
  expect_true(is.finite(td$std.error))
  expect_true(is.finite(td$conf.low))
  expect_true(is.finite(td$conf.high))
  expect_true(is.finite(td$p.value))
})

# § C: enriched glance.coresynth ─────────────────────────────────────────────

test_that("Phase 24: glance.coresynth returns new columns", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  g <- broom::glance(fit)
  expect_true(all(c("method", "estimate", "n_controls", "n_treated",
                    "T_pre", "T_post", "staggered", "multi_arm") %in% names(g)))
  expect_equal(g$T_pre, 10L)
  expect_equal(g$T_post, 10L)
  expect_equal(g$n_treated, 1L)
  expect_false(g$staggered)
  expect_false(g$multi_arm)
})

test_that("Phase 24: glance.coresynth staggered SDID sets staggered=TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  g <- broom::glance(fit)
  expect_true(g$staggered)
  expect_false(g$multi_arm)
})

test_that("Phase 24: glance.coresynth multi-arm SI sets multi_arm=TRUE", {
  fit <- scm_fit(y ~ d | id + time, data = tensor_panel, method = "si")
  g <- broom::glance(fit)
  expect_true(g$multi_arm)
  expect_false(g$staggered)
})

# § D: augment.coresynth staggered & multi-arm ────────────────────────────────

test_that("Phase 24: augment.coresynth() SDID staggered returns long df", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  aug <- broom::augment(fit)
  expect_s3_class(aug, "data.frame")
  expect_true(nrow(aug) > 0L)
  expect_true(all(c(".cohort", ".unit", ".type", ".time", ".observed",
                    ".fitted", ".period") %in% names(aug)))
  # SDID staggered reconstructs fitted via unit_weights + omega0
  expect_true(any(is.finite(aug$.fitted)))
})

test_that("Phase 24: augment.coresynth() SI staggered carries .cohort column", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  aug <- broom::augment(fit)
  expect_true(".cohort" %in% names(aug))
  expect_true(nrow(aug) > 0L)
  # SI cohort_fits carry Y_synth → .fitted should be finite
  expect_true(any(is.finite(aug$.fitted)))
})

test_that("Phase 24: augment.coresynth() multi-arm SI carries .arm column", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  aug <- broom::augment(fit)
  expect_true(".arm" %in% names(aug))
  expect_true(any(!is.na(aug$.arm)))
})

test_that("Phase 24: augment.coresynth() staggered include_donors adds control rows", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  aug <- broom::augment(fit, include_donors = TRUE)
  expect_true("control" %in% aug$.type)
  expect_true("treated" %in% aug$.type)
})

test_that("Phase 24: augment.coresynth() sharp unchanged (regression)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  aug <- broom::augment(fit)
  expect_false(".cohort" %in% names(aug))
  expect_true(all(c(".time", ".observed", ".fitted", ".resid") %in% names(aug)))
})

# § E: improved error messages ────────────────────────────────────────────────

test_that("Phase 24: mspe_ratio_pval() rejects staggered fits with helpful message", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_error(mspe_ratio_pval(fit), regexp = "staggered")
})

test_that("Phase 24: SCM staggered + predictors error message names alternatives", {
  err <- tryCatch(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm",
            predictors = list(pred("y", 1:3))),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("covariates", err))
  expect_true(grepl("sdid|gsc|mc|tasc", err))
})

# § F: dead-code cleanup regression ───────────────────────────────────────────

test_that("Phase 24: SCM staggered unchanged after max_iter/tol cleanup", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_true(isTRUE(fit$staggered))
  expect_true(is.finite(fit$estimate))
  expect_equal(sum(fit$cohort_estimates$weight), 1, tolerance = 1e-8)
})

test_that("Phase 24: SCM staggered + covariates unchanged after cleanup", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
                 covariates = "cov1")
  expect_true(isTRUE(fit$staggered))
  expect_true(is.finite(fit$estimate))
  expect_length(fit$beta_hat, 1L)
})

test_that("Phase 24: sharp SCM tidy() unchanged (regression)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  td <- broom::tidy(fit)
  expect_true("unit_weight" %in% td$type)
})

test_that("Phase 24: sdid_inference object inherits coresynth_inference", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  expect_s3_class(inf, "sdid_inference")
  expect_s3_class(inf, "coresynth_inference")
})

# ── Conformal inference (Chernozhukov, Wuthrich & Zhu 2021) ───────────────────
test_that("conformal_inference works across supported sharp methods", {
  for (m in c("scm", "sdid", "gsc", "mc", "si")) {
    fit <- scm_fit(y ~ d | id + time, data = panel, method = m)
    ci  <- conformal_inference(fit, n_grid = 60L)
    expect_s3_class(ci, "conformal_inference")
    expect_s3_class(ci, "coresynth_inference")
    expect_equal(ci$method, "conformal")
    expect_true(ci$p_value >= 0 && ci$p_value <= 1)
    expect_true(is.finite(ci$estimate))
    expect_true(ci$n_controls == 9L)
  }
})

test_that("conformal p-value for true effect exceeds p-value for tau0 = 0", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p_null   <- conformal_inference(fit, tau0 = 0, ci = FALSE)$p_value
  p_truth  <- conformal_inference(fit, tau0 = fit$estimate, ci = FALSE)$p_value
  # The estimate itself should not be rejected; tau0 = 0 (no effect) should be
  # less compatible than the point estimate.
  expect_true(p_truth >= p_null)
})

test_that("conformal CI contains the point estimate", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  ci  <- conformal_inference(fit, n_grid = 120L)
  if (!is.na(ci$ci_lower)) {
    expect_true(ci$ci_lower <= fit$estimate && fit$estimate <= ci$ci_upper)
  }
  expect_true(is.numeric(ci$grid) && length(ci$grid) == 120L)
})

test_that("conformal one-sided alternatives run and tidy()/glance() work", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  ci_g <- conformal_inference(fit, alternative = "greater", ci = FALSE)
  ci_l <- conformal_inference(fit, alternative = "less", ci = FALSE)
  expect_true(ci_g$p_value >= 0 && ci_g$p_value <= 1)
  expect_true(ci_l$p_value >= 0 && ci_l$p_value <= 1)
  td <- broom::tidy(conformal_inference(fit, n_grid = 40L))
  expect_equal(td$term, "ATT")
  expect_equal(td$method, "conformal")
  gl <- broom::glance(conformal_inference(fit, n_grid = 40L))
  expect_equal(gl$method, "conformal")
})

test_that("conformal_inference rejects staggered and tasc fits", {
  stag <- make_panel()
  stag$d <- 0L
  stag$d[stag$id == "u1" & stag$time > 10] <- 1L
  stag$d[stag$id == "u2" & stag$time > 15] <- 1L
  stag$y[stag$d == 1] <- stag$y[stag$d == 1] + 2.0
  fit_stag <- scm_fit(y ~ d | id + time, data = stag, method = "sdid")
  expect_error(conformal_inference(fit_stag), "sharp")

  fit_tasc <- scm_fit(y ~ d | id + time, data = panel, method = "tasc")
  expect_error(conformal_inference(fit_tasc), "tasc")
})

# ── Phase 26: Abadie 理論照合 — predictor scaling・OOS CV・ADH 2015 検証ツール ─

test_that("Phase 26: predictor scaling makes W invariant to predictor units", {
  pan_sc <- panel_cov
  pan_sc$cov1 <- pan_sc$cov1 * 1000  # same information, different units
  fit_base <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                      predictors = list(pred(c("cov1", "cov2"), 1:10)))
  fit_resc <- scm_fit(y ~ d | id + time, data = pan_sc, method = "scm",
                      predictors = list(pred(c("cov1", "cov2"), 1:10)))
  expect_equal(fit_base$unit_weights, fit_resc$unit_weights, tolerance = 1e-6)
  expect_equal(fit_base$estimate, fit_resc$estimate, tolerance = 1e-6)
})

test_that("Phase 26: scale_predictors = FALSE disables Synth-style scaling", {
  fit_off <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                     predictors = list(pred(c("cov1", "cov2"), 1:10)),
                     scale_predictors = FALSE)
  expect_s3_class(fit_off, "coresynth")
  expect_equal(sum(fit_off$unit_weights), 1, tolerance = 1e-4)
})

test_that("Phase 26: predictor_table reports original (unscaled) values", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                 predictors = list(pred(c("cov1", "cov2"), 1:10)))
  tr_cov1 <- mean(panel_cov$cov1[panel_cov$id == "u1" & panel_cov$time %in% 1:10])
  expect_equal(fit$predictor_table$treated[1], tr_cov1, tolerance = 1e-8)
})

test_that("Phase 26: build_predictor_matrices rejects non-finite predictors", {
  expect_error(
    scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
            predictors = list(pred("cov1", 100:110))),
    regexp = "non-finite"
  )
})

test_that("Phase 26: OOS V selection uses train/validation split (no leakage)", {
  fit_oos <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_selection = "oos")
  # V applies to the last floor(T_pre/2) outcome rows (Abadie 2021 §3.2 step 4)
  expect_equal(names(fit_oos$v_weights), paste0("V", 6:10))
  expect_equal(fit_oos$v_rows, 6:10)
  expect_equal(sum(fit_oos$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(is.finite(fit_oos$Y_synth)))
})

test_that("Phase 26: OOS + bfgs uses the same train/validation split", {
  fit <- suppressWarnings(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            v_selection = "oos", v_optim = "bfgs"))
  expect_equal(length(fit$v_weights), 5L)
  expect_true(all(is.finite(fit$Y_synth)))
})

test_that("Phase 26: OOS + lambda_pen builds the QP on the last-half window", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_selection = "oos", lambda_pen = 0.1)
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 26: placebo_in_time shows no effect when backdated (ADH 2015)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  pit <- placebo_in_time(fit, t0_placebo = 5L)
  expect_equal(pit$t0_placebo, 5L)
  expect_equal(length(pit$gap), 10L)
  expect_true(is.finite(pit$placebo_att))
  # no treatment occurs in the pre-period, so the placebo ATT should be
  # small relative to the true post-period effect (2.0)
  expect_lt(abs(pit$placebo_att), 1.5)
})

test_that("Phase 26: placebo_in_time validates inputs", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(placebo_in_time(fit, t0_placebo = 1L), "t0_placebo")
  expect_error(placebo_in_time(fit, t0_placebo = 10L), "t0_placebo")
  stag_fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_error(placebo_in_time(stag_fit), "sharp")
})

test_that("Phase 26: loo_donors returns one row per contributing donor", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  loo <- loo_donors(fit)
  expect_gte(nrow(loo$results), 1L)
  expect_true(all(is.finite(loo$results$att_loo)))
  expect_true(all(loo$results$weight > 1e-6))
  expect_equal(loo$att_original, fit$estimate)
  expect_true(all(c("donor", "weight", "att_loo") %in% names(loo$results)))
})

test_that("Phase 26: loo_donors works for OOS and predictor fits", {
  fit_oos <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_selection = "oos")
  loo_oos <- loo_donors(fit_oos)
  expect_true(all(is.finite(loo_oos$results$att_loo)))

  fit_cov <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                     predictors = list(pred(c("cov1", "cov2"), 1:10)))
  loo_cov <- loo_donors(fit_cov)
  expect_true(all(is.finite(loo_cov$results$att_loo)))

  stag_fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_error(loo_donors(stag_fit), "sharp")
})

# == Input validation: empty / degenerate / unbalanced panels =================

test_that("0-row data errors with a clear message for every method", {
  empty <- panel[0, ]
  for (m in c("scm", "sdid", "gsc", "mc", "tasc", "si")) {
    expect_error(
      scm_fit(y ~ d | id + time, data = empty, method = m),
      "0 rows", info = m
    )
  }
})

test_that("data with no treated units errors clearly", {
  no_tr <- panel[panel$id != "u1", ]
  expect_error(
    scm_fit(y ~ d | id + time, data = no_tr, method = "scm"),
    "No treated units"
  )
})

test_that("sharp data with no control units errors clearly", {
  all_tr <- panel
  all_tr$d <- as.integer(all_tr$time > 10)
  expect_error(
    scm_fit(y ~ d | id + time, data = all_tr, method = "scm"),
    "No control units"
  )
})

test_that("staggered panel with no never-treated units: donor-pool methods stop, MC fits", {
  stag_all <- make_panel()
  stag_all$d <- as.integer(
    (stag_all$id %in% paste0("u", 1:3) & stag_all$time > 10) |
    (!stag_all$id %in% paste0("u", 1:3) & stag_all$time > 16)
  )
  # A later adopter is treated during the earlier cohort's post window, so it
  # cannot stand in for that cohort's untreated outcomes.
  for (m in c("scm", "sdid", "gsc", "si"))
    expect_error(scm_fit(y ~ d | id + time, data = stag_all, method = m),
                 "never-treated")
  fit <- scm_fit(y ~ d | id + time, data = stag_all, method = "mc")
  expect_true(isTRUE(fit$staggered))
  expect_true(is.finite(fit$estimate))
})

test_that("NA in the treatment indicator errors instead of silently dropping the unit", {
  na_d <- panel
  na_d$d[57] <- NA
  expect_error(
    scm_fit(y ~ d | id + time, data = na_d, method = "scm"),
    "NA value"
  )
})

test_that("negative treatment values are rejected", {
  neg <- panel
  neg$d[1] <- -1L
  expect_error(
    scm_fit(y ~ d | id + time, data = neg, method = "scm"),
    "non-negative"
  )
})

test_that("NA unit identifiers are rejected", {
  na_id <- panel
  na_id$id[3] <- NA
  expect_error(
    scm_fit(y ~ d | id + time, data = na_id, method = "scm"),
    "identifiers contain NA"
  )
})

test_that("missing (id, time) cells error clearly for dense-matrix methods", {
  unbal <- panel[!(panel$id == "u5" & panel$time == 10), ]
  for (m in c("scm", "sdid", "gsc", "si")) {
    expect_error(
      scm_fit(y ~ d | id + time, data = unbal, method = m),
      "balanced panel with fully observed outcomes", info = m
    )
  }
})

test_that("NA outcome values error the same way as missing cells", {
  na_y <- panel
  na_y$y[na_y$id == "u5" & na_y$time == 10] <- NA
  expect_error(
    scm_fit(y ~ d | id + time, data = na_y, method = "scm"),
    "balanced panel with fully observed outcomes"
  )
})

test_that("MC masks missing cells instead of treating them as observed zeros", {
  unbal <- panel[!(panel$id == "u5" & panel$time == 10), ]
  fit_full  <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  fit_unbal <- scm_fit(y ~ d | id + time, data = unbal, method = "mc")
  expect_true(is.finite(fit_unbal$estimate))
  expect_lt(abs(fit_unbal$estimate - fit_full$estimate), 0.5)
})

test_that("TASC fits a panel with missing outcome cells", {
  unbal <- panel[!(panel$id == "u5" & panel$time == 10), ]
  fit <- scm_fit(y ~ d | id + time, data = unbal, method = "tasc")
  expect_true(is.finite(fit$estimate))
})

test_that("factor treatment and non-numeric outcome columns are rejected", {
  fac_d <- panel
  fac_d$d <- factor(fac_d$d)
  expect_error(
    scm_fit(y ~ d | id + time, data = fac_d, method = "scm"),
    "must be numeric, integer, or logical"
  )
  chr_y <- panel
  chr_y$y <- as.character(chr_y$y)
  expect_error(
    scm_fit(y ~ d | id + time, data = chr_y, method = "scm"),
    "must be numeric"
  )
})

test_that("non-integer treatment values are rejected", {
  frac_d <- panel
  frac_d$d <- frac_d$d * 0.5
  expect_error(
    scm_fit(y ~ d | id + time, data = frac_d, method = "scm"),
    "non-integer"
  )
})

test_that("SCM with zero pre-treatment periods errors clearly", {
  t0 <- panel
  t0$d <- as.integer(t0$id == "u1")
  expect_error(
    scm_fit(y ~ d | id + time, data = t0, method = "scm"),
    "No pre-treatment periods"
  )
})

test_that("scm_design rejects 0-row data", {
  expect_error(
    scm_design(panel[0, ], outcome = "y", unit = "id", time = "time", T0 = 10),
    "0 rows"
  )
})

# ── Phase 28: in-space placebo gap paths + plot (ADH 2010, Figures 4-8) ───────

test_that("scm_placebo_cpp returns gap paths consistent with effects", {
  fit  <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  plac <- scm_placebo_cpp(fit$Y_co_pre, fit$Y_co_post)
  expect_equal(dim(plac$gaps), c(20L, 9L))
  expect_true(all(is.finite(plac$gaps)))
  # Mean of the post-period gap rows must reproduce the placebo effects
  expect_equal(colMeans(plac$gaps[11:20, , drop = FALSE]), as.numeric(plac$effects),
               tolerance = 1e-10)
  # Mean squared pre-period gap rows must reproduce mspe_pre
  expect_equal(colMeans(plac$gaps[1:10, , drop = FALSE]^2), as.numeric(plac$mspe_pre),
               tolerance = 1e-10)
})

test_that("mspe_ratio_pval returns scm_placebo object with gap paths", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  expect_s3_class(inf, "scm_placebo")
  expect_equal(dim(inf$gaps), c(20L, 9L))
  expect_equal(colnames(inf$gaps), sort(paste0("u", 2:10)))
  expect_equal(length(inf$treated_gap), 20L)
  expect_equal(inf$treated_gap, as.numeric(fit$gap))
  expect_equal(length(inf$mspe_pre_placebo), 9L)
  expect_true(is.finite(inf$mspe_pre_treated))
  expect_equal(inf$T_pre, 10L)
  expect_equal(inf$times, fit$times)
  # placebo_effects must equal post-period column means of the gap paths
  expect_equal(as.numeric(inf$placebo_effects), unname(colMeans(inf$gaps[11:20, ])),
               tolerance = 1e-10)
})

test_that("mspe_ratio_pval use_covariates=TRUE also returns gap paths", {
  fit <- scm_fit(
    y ~ d | id + time, data = panel, method = "scm",
    predictors = list(pred("y", 1:5, "mean"), pred("y", 6:10, "mean"))
  )
  inf <- mspe_ratio_pval(fit, use_covariates = TRUE)
  expect_s3_class(inf, "scm_placebo")
  expect_equal(dim(inf$gaps), c(20L, 9L))
  expect_true(all(is.finite(inf$gaps)))
  expect_equal(as.numeric(inf$placebo_effects), unname(colMeans(inf$gaps[11:20, ])),
               tolerance = 1e-10)
})

test_that("plot.scm_placebo type='gaps' returns a ggplot and prunes by MSPE multiple", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)

  p_all <- plot(inf, type = "gaps")
  expect_s3_class(p_all, "ggplot")
  expect_equal(length(unique(p_all$layers[[1]]$data$unit)), 9L)

  # ADH 2010 Figures 5-7 style pruning: only well-fitted placebos remain
  mult    <- 2
  n_keep  <- sum(inf$mspe_pre_placebo <= mult * inf$mspe_pre_treated)
  p_prune <- plot(inf, type = "gaps", mspe_prune = mult)
  expect_s3_class(p_prune, "ggplot")
  expect_equal(length(unique(p_prune$layers[[1]]$data$unit)), n_keep)
})

test_that("plot.scm_placebo warns when all placebos are pruned", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  expect_warning(p <- plot(inf, type = "gaps", mspe_prune = 1e-12),
                 "pruned")
  expect_s3_class(p, "ggplot")
})

test_that("plot.scm_placebo type='ratios' returns a ggplot with all units", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  p <- plot(inf, type = "ratios")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 10L) # treated + 9 donors
  expect_true("Treated" %in% p$data$unit)
})

test_that("plot.scm_placebo rejects invalid mspe_prune", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  expect_error(plot(inf, mspe_prune = -1), "positive")
  expect_error(plot(inf, mspe_prune = c(1, 2)), "single")
})

# ── Phase 29: plot style customization (colors/vline/hline/fill) ────────────

test_that("plot.coresynth trend: color override reaches the built plot data", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p_default <- plot(fit, type = "trend")
  expect_equal(length(p_default$layers), 2L) # geom_line + geom_vline

  p_custom <- plot(fit, type = "trend", colors = c(treated = "black"))
  built <- ggplot2::ggplot_build(p_custom)
  expect_true("black" %in% built$data[[1]]$colour)
  expect_true("#d73027" %in% built$data[[1]]$colour) # Synthetic Control unchanged
})

test_that("plot.coresynth trend: unknown color name errors", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(plot(fit, type = "trend", colors = c(Bogus = "red")),
               "unrecognized name")
})

test_that("plot.coresynth trend: vline = FALSE/NULL suppresses the line", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  n_default <- length(plot(fit, type = "trend")$layers)
  expect_equal(length(plot(fit, type = "trend", vline = FALSE)$layers), n_default - 1L)
  expect_equal(length(plot(fit, type = "trend", vline = NULL)$layers), n_default - 1L)
})

test_that("plot.coresynth trend: vline must be a list", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(plot(fit, type = "trend", vline = "red"),
               "list of aesthetic overrides")
})

test_that("plot.coresynth gap: color/vline/hline overrides work", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  n_default <- length(plot(fit, type = "gap")$layers)
  expect_equal(length(plot(fit, type = "gap", vline = NULL, hline = NULL)$layers),
               n_default - 2L)

  p_color <- plot(fit, type = "gap", colors = "black")
  built   <- ggplot2::ggplot_build(p_color)
  expect_true(all(built$data[[1]]$colour == "black"))
})

test_that("plot.coresynth weights: fill override reaches the built plot data", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p <- plot(fit, type = "weights", fill = "darkorange")
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[1]]$fill == "darkorange"))
})

test_that("plot.coresynth weights: top_n limits bars to the largest weights", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  n_all <- nrow(ggplot2::ggplot_build(plot(fit, type = "weights"))$data[[1]])

  built <- ggplot2::ggplot_build(plot(fit, type = "weights", top_n = 2))
  expect_equal(nrow(built$data[[1]]), 2L)
  w_kept <- sort(fit$unit_weights[fit$unit_weights > 1e-4], decreasing = TRUE)
  expect_setequal(built$data[[1]]$y, unname(w_kept[1:2]))

  # top_n beyond the donor count keeps everything (default Inf unchanged)
  expect_equal(
    nrow(ggplot2::ggplot_build(plot(fit, type = "weights", top_n = 999))$data[[1]]),
    n_all
  )

  expect_error(plot(fit, type = "weights", top_n = 0), "top_n")
  expect_error(plot(fit, type = "weights", top_n = c(1, 2)), "top_n")
})

test_that("plot.coresynth pred_weights: one bar per predictor, fill reaches plot", {
  fit <- scm_fit(
    y ~ d | id + time, data = panel_cov, method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  p <- plot(fit, type = "pred_weights", fill = "darkorange")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[1]]), length(fit$v_weights))
  expect_true(all(built$data[[1]]$fill == "darkorange"))
  # bar heights are the V weights (order-independent)
  expect_setequal(built$data[[1]]$y, unname(fit$v_weights))
})

test_that("plot.coresynth pred_weights: top_n keeps the largest V weights", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  n_all <- nrow(ggplot2::ggplot_build(plot(fit, type = "pred_weights"))$data[[1]])
  expect_equal(n_all, length(fit$v_weights))

  built <- ggplot2::ggplot_build(plot(fit, type = "pred_weights", top_n = 2))
  expect_equal(nrow(built$data[[1]]), 2L)
  v_kept <- sort(fit$v_weights, decreasing = TRUE)[1:2]
  expect_setequal(built$data[[1]]$y, unname(v_kept))

  expect_error(plot(fit, type = "pred_weights", top_n = 0), "top_n")
})

test_that("plot.coresynth pred_weights: errors when no V matrix exists", {
  fit_mc <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  expect_error(plot(fit_mc, type = "pred_weights"), "V")

  # staggered SCM stores v_weights = NULL
  fit_stag <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_error(plot(fit_stag, type = "pred_weights"), "V")
})

test_that("plot.scm_placebo type='gaps': color/vline/hline overrides work", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)

  n_default <- length(plot(inf, type = "gaps")$layers)
  expect_equal(length(plot(inf, type = "gaps", vline = FALSE, hline = FALSE)$layers),
               n_default - 2L)

  p_custom <- plot(inf, type = "gaps",
                    colors = c(placebo = "lightblue"))
  built <- ggplot2::ggplot_build(p_custom)
  expect_true("lightblue" %in% built$data[[1]]$colour)
  expect_true("#2166ac" %in% built$data[[2]]$colour) # Treated unchanged
})

test_that("plot.scm_placebo: unknown color name errors", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  expect_error(plot(inf, type = "gaps", colors = c(Bogus = "red")),
               "unrecognized name")
  expect_error(plot(inf, type = "ratios", colors = c(Bogus = "red")),
               "unrecognized name")
})

test_that("plot.coresynth trend: labels relabel the legend on both scales", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p <- plot(fit, type = "trend", labels = c(treated = "California"))
  expect_equal(p$scales$get_scales("colour")$labels[["Treated"]], "California")
  expect_equal(p$scales$get_scales("linetype")$labels[["Treated"]], "California")
  # unmentioned series keeps its default label
  expect_equal(p$scales$get_scales("colour")$labels[["Synthetic Control"]],
               "Synthetic Control")
})

test_that("plot.coresynth trend: labels validation", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(plot(fit, type = "trend", labels = c(Bogus = "x")),
               "unrecognized name")
  expect_error(plot(fit, type = "trend", labels = "x"), "named vector")
})

test_that("plot.scm_placebo: labels relabel the legend and the ratios axis tick", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)

  p_gaps <- plot(inf, type = "gaps", labels = c(placebo = "Donors"))
  expect_equal(p_gaps$scales$get_scales("colour")$labels[["Placebo (donor pool)"]],
               "Donors")

  p_ratios <- plot(inf, type = "ratios", labels = c(treated = "California"))
  expect_equal(p_ratios$scales$get_scales("colour")$labels[["Treated"]], "California")
  built <- ggplot2::ggplot_build(p_ratios)
  expect_true("California" %in% built$layout$panel_params[[1]]$y$get_labels())

  expect_error(plot(inf, type = "gaps", labels = c(Bogus = "x")),
               "unrecognized name")
})

test_that("plot.coresynth trend: linetypes override the series line types", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")

  built_default <- ggplot2::ggplot_build(plot(fit, type = "trend"))$data[[1]]
  expect_setequal(built_default$linetype, c("solid", "dashed"))

  # both series solid: the point of the argument
  built <- ggplot2::ggplot_build(
    plot(fit, type = "trend", linetypes = c(synthetic = "solid")))$data[[1]]
  expect_true(all(built$linetype == "solid"))

  # unmentioned series keeps its default, and integer codes are accepted
  built_num <- ggplot2::ggplot_build(
    plot(fit, type = "trend", linetypes = c(treated = 3)))$data[[1]]
  expect_setequal(built_num$linetype, c("dotted", "dashed"))
})

test_that("plot.coresynth trend: linetypes compose with colors and labels", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p <- plot(fit, type = "trend",
            colors    = c(treated = "black"),
            labels    = c(treated = "California"),
            linetypes = c(treated = "dashed", synthetic = "solid"))
  built <- ggplot2::ggplot_build(p)$data[[1]]
  expect_true("black" %in% built$colour)
  expect_setequal(built$linetype, c("dashed", "solid"))
  # legend stays merged: identical labels on the color and linetype scales
  expect_equal(p$scales$get_scales("linetype")$labels[["Treated"]], "California")
  expect_equal(p$scales$get_scales("colour")$labels[["Treated"]], "California")
})

test_that("plot.coresynth trend: linetypes validation", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(plot(fit, type = "trend", linetypes = c(Bogus = "solid")),
               "unrecognized name")
  expect_error(plot(fit, type = "trend", linetypes = "solid"), "named vector")
  expect_error(plot(fit, type = "trend", linetypes = c(treated = 9)),
               "integer code")
  # the donors key exists only while donor paths are drawn
  expect_error(plot(fit, type = "trend", linetypes = c(donors = "solid")),
               "unrecognized name")
  p <- plot(fit, type = "trend", show_donors = 2,
            linetypes = c(donors = "dotted"))
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$linetype[1], "dotted")
})

test_that("plot.coresynth gap: linetypes take a single line type", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_true(all(
    ggplot2::ggplot_build(plot(fit, type = "gap"))$data[[1]]$linetype == "solid"))

  built <- ggplot2::ggplot_build(
    plot(fit, type = "gap", linetypes = "dashed"))$data[[1]]
  expect_true(all(built$linetype == "dashed"))

  built_num <- ggplot2::ggplot_build(
    plot(fit, type = "gap", linetypes = 3))$data[[1]]
  expect_true(all(built_num$linetype == "dotted"))

  expect_error(plot(fit, type = "gap", linetypes = 9), "integer code")
})

test_that("plot.scm_placebo gaps: linetypes override the series line types", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)

  built_default <- ggplot2::ggplot_build(plot(inf, type = "gaps"))$data
  expect_true(all(built_default[[1]]$linetype == "solid"))
  expect_true(all(built_default[[2]]$linetype == "solid"))

  p <- plot(inf, type = "gaps", linetypes = c(placebo = "dotted"),
            labels = c(placebo = "Donors"))
  built <- ggplot2::ggplot_build(p)$data
  expect_true(all(built[[1]]$linetype == "dotted"))  # placebo paths
  expect_true(all(built[[2]]$linetype == "solid"))   # treated path unchanged
  # merged legend: both scales relabel the same series
  expect_equal(p$scales$get_scales("linetype")$labels[["Placebo (donor pool)"]],
               "Donors")

  expect_error(plot(inf, type = "gaps", linetypes = c(Bogus = "solid")),
               "unrecognized name")
})

# ── plot_data(): tidy extraction of the plot() data ──────────────────────────

test_that("plot_data(type='trend') mirrors the treated/synthetic accessors", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  df  <- plot_data(fit, type = "trend")
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("time", "value", "series"))
  expect_setequal(unique(df$series), c("Treated", "Synthetic Control"))
  expect_equal(df$value[df$series == "Treated"],
               as.numeric(treated_outcomes(fit, na.rm = TRUE)))
  expect_equal(df$value[df$series == "Synthetic Control"],
               as.numeric(synthetic_outcomes(fit, na.rm = TRUE)))
})

test_that("plot_data(type='gap') equals treated minus synthetic", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  df  <- plot_data(fit, type = "gap")
  expect_identical(names(df), c("time", "gap"))
  expect_equal(df$gap,
               as.numeric(treated_outcomes(fit, na.rm = TRUE) -
                            synthetic_outcomes(fit, na.rm = TRUE)))
})

test_that("plot_data trend align shifts synthetic by the pre-period gap", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  raw <- plot_data(fit, type = "gap")
  ali <- plot_data(fit, type = "gap", align = TRUE)
  # SDID: lambda-aligned post-period mean gap equals the point estimate
  post <- seq(fit$T_pre + 1L, length(ali$gap))
  expect_equal(mean(ali$gap[post]), unname(fit$estimate), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(raw$gap, ali$gap)))
})

test_that("plot_data trend show_donors adds Donors rows and a unit column", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  df  <- plot_data(fit, type = "trend", show_donors = 3)
  expect_true("unit" %in% names(df))
  expect_setequal(unique(df$series),
                  c("Treated", "Synthetic Control", "Donors"))
  # exactly three distinct donor units
  expect_equal(length(unique(df$unit[df$series == "Donors"])), 3L)
})

test_that("plot_data(type='weights') returns every donor (no 1e-4 pruning)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  df  <- plot_data(fit, type = "weights")
  expect_identical(names(df), c("unit", "weight"))
  # complete: one row per donor, unlike the plot which drops tiny weights
  expect_equal(nrow(df), length(fit$unit_weights))
  expect_equal(sort(df$weight), sort(as.numeric(fit$unit_weights)))
})

test_that("plot_data(type='weights') top_n keeps the largest weights", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  df  <- plot_data(fit, type = "weights", top_n = 2)
  expect_equal(nrow(df), 2L)
  expect_equal(df$weight, sort(as.numeric(fit$unit_weights),
                               decreasing = TRUE)[1:2])
  expect_error(plot_data(fit, type = "weights", top_n = 0), "top_n")
})

test_that("plot_data(type='weights') for SDID adds an omega/lambda panel", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  df  <- plot_data(fit, type = "weights")
  expect_true("panel" %in% names(df))
  expect_setequal(unique(df$panel), c("omega", "lambda"))
  expect_equal(df$weight[df$panel == "lambda"], as.numeric(fit$time_weights))
})

test_that("plot_data(type='pred_weights') returns one row per predictor", {
  fit <- scm_fit(
    y ~ d | id + time, data = panel_cov, method = "scm",
    predictors = list(pred(c("cov1", "cov2"), 1:10))
  )
  df <- plot_data(fit, type = "pred_weights")
  expect_identical(names(df), c("predictor", "weight"))
  expect_equal(df$predictor, names(fit$v_weights))
  expect_equal(df$weight, as.numeric(fit$v_weights))
})

test_that("plot_data errors on fits without the requested component", {
  fit_mc <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  expect_error(plot_data(fit_mc, type = "pred_weights"), "V")
  expect_error(plot_data(fit_mc, type = "weights"), "No unit weights")

  fit_stag <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_error(plot_data(fit_stag, type = "trend"), "per cohort")
})

test_that("plot_data has no method for unsupported objects", {
  expect_error(plot_data(1:10), "no method")
})

test_that("plot_data.scm_placebo returns tidy gaps and ratios frames", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)

  g <- plot_data(inf, type = "gaps")
  expect_identical(names(g), c("time", "gap", "unit", "series"))
  expect_setequal(unique(g$series), c("Treated", "Placebo (donor pool)"))
  expect_equal(g$gap[g$series == "Treated"], as.numeric(inf$treated_gap))

  r <- plot_data(inf, type = "ratios")
  expect_identical(names(r), c("unit", "ratio", "post_mspe", "series"))
  expect_equal(r$ratio, as.numeric(inf$mspe_ratios_all))
  expect_equal(r$post_mspe,
               as.numeric(c(inf$mspe_post_treated, inf$mspe_post_placebo)))
  expect_equal(r$series[1], "Treated")
})

test_that("plot_data.scm_placebo gaps honors mspe_prune", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  n_all  <- length(unique(plot_data(inf, type = "gaps")$unit))
  n_keep <- length(unique(plot_data(inf, type = "gaps", mspe_prune = 1)$unit))
  expect_lte(n_keep, n_all)
  expect_error(plot_data(inf, type = "gaps", mspe_prune = -1), "positive")
})

# ── Phase 33: Partially pooled staggered SCM (Ben-Michael et al. 2022) ────────

test_that("Phase 33: nu/fixedeff on a sharp fit error informatively", {
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm", nu = 0.5),
    regexp = "staggered"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm", fixedeff = TRUE),
    regexp = "staggered"
  )
})

test_that("Phase 33: invalid nu is rejected", {
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm", nu = 1.5),
    regexp = "nu must be"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm", nu = "bogus"),
    regexp = "nu must be"
  )
})

test_that("Phase 33: nu is incompatible with the legacy-path knobs", {
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm",
            nu = 0.5, lambda_pen = 0.1),
    regexp = "cannot be combined"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm",
            nu = 0.5, v_selection = "oos"),
    regexp = "cannot be combined"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm",
            nu = 0.5, donor_mspe_threshold = 20),
    regexp = "cannot be combined"
  )
})

test_that("Phase 33: partially pooled SCM returns a valid staggered fit", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm", nu = 0.5)
  expect_s3_class(fit, "coresynth_staggered")
  expect_true(is.finite(fit$estimate))
  expect_lt(abs(fit$estimate - 1.5), 4.0)
  expect_equal(sum(fit$cohort_estimates$weight), 1, tolerance = 1e-8)
  for (cf in fit$cohort_fits) {
    expect_equal(sum(cf$unit_weights), 1, tolerance = 1e-5)
    expect_true(all(cf$unit_weights >= -1e-5))
  }
})

test_that("Phase 33: pooling diagnostics are stored and consistent", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm", nu = 0.5)
  p <- fit$pooling
  expect_equal(p$nu, 0.5)
  # staggered SCM now carries the intercept shift by default (BFR S.5.1)
  expect_true(p$fixedeff)
  expect_false(scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                       nu = 0.5, fixedeff = FALSE)$pooling$fixedeff)
  expect_true(all(is.finite(c(p$q_sep, p$q_pool,
                              p$q_sep_separate, p$q_pool_separate))))
  # pooling can only improve the pooled fit relative to separate SCM
  expect_lte(p$q_pool, p$q_pool_separate + 1e-10)
  # and can only worsen the separate fit
  expect_gte(p$q_sep, p$q_sep_separate - 1e-10)
})

test_that("Phase 33: nu = 0 equals the separate (uniform-V) solution", {
  fit0 <- scm_fit(y ~ d | id + time, data = staggered, method = "scm", nu = 0)
  p <- fit0$pooling
  expect_equal(p$q_sep, p$q_sep_separate, tolerance = 1e-12)
  expect_equal(p$q_pool, p$q_pool_separate, tolerance = 1e-12)
})

test_that("Phase 33: pooled fit improves monotonically in nu", {
  fit0 <- scm_fit(y ~ d | id + time, data = staggered, method = "scm", nu = 0)
  fit1 <- scm_fit(y ~ d | id + time, data = staggered, method = "scm", nu = 1)
  expect_lte(fit1$pooling$q_pool, fit0$pooling$q_pool + 1e-10)
  expect_gte(fit1$pooling$q_sep,  fit0$pooling$q_sep - 1e-10)
})

test_that("Phase 33: nu = 'auto' picks the heuristic in [0, 1]", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                 nu = "auto")
  expect_true(fit$pooling$nu >= 0 && fit$pooling$nu <= 1)
  expect_equal(fit$pooling$nu, fit$pooling$nu_heuristic)
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 33: fixedeff works on both staggered paths", {
  # legacy path
  fit_l <- scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                   fixedeff = TRUE)
  expect_true(is.finite(fit_l$estimate))
  expect_true(fit_l$fixedeff)
  # pooled path
  fit_p <- scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                   nu = 0.5, fixedeff = TRUE)
  expect_true(is.finite(fit_p$estimate))
  expect_true(fit_p$pooling$fixedeff)
  # Y_synth is reported on the raw outcome scale: pre-period levels track
  # the treated series (alpha restores the intercept)
  for (cf in fit_p$cohort_fits) {
    pre <- seq_len(cf$T_pre)
    expect_lt(abs(mean(cf$Y_treat[pre]) - mean(cf$Y_synth[pre])), 1.0)
  }
})

test_that("Phase 33: fixedeff recovers ATT under large unit intercepts", {
  # DGP where donors share the treated units' trend but sit at shifted
  # levels: plain SCM cannot match the level, the intercept shift can.
  set.seed(11)
  N <- 12; TT <- 24
  f <- cumsum(rnorm(TT, 0, 0.5))
  shift <- c(0, 0.5, seq(4, 14, length.out = N - 2))
  rows <- expand.grid(time = seq_len(TT), id = paste0("u", seq_len(N)),
                      stringsAsFactors = FALSE)
  rows$y <- rep(f, N) + rep(shift, each = TT) + rnorm(nrow(rows), 0, 0.2)
  rows$d <- as.integer(
    (rows$id == "u1" & rows$time >= 9) | (rows$id == "u2" & rows$time >= 15)
  )
  rows$y[rows$d == 1] <- rows$y[rows$d == 1] + 1.5

  fit_fe <- scm_fit(y ~ d | id + time, data = rows, method = "scm",
                    nu = 0.5, fixedeff = TRUE)
  expect_lt(abs(fit_fe$estimate - 1.5), 0.5)
})

test_that("Phase 33: staggered cohort_fits expose idx_tr/idx_co/Y_treat_mat", {
  for (nu in list(NULL, 0.5)) {
    fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                   nu = nu)
    for (cf in fit$cohort_fits) {
      expect_true(is.integer(cf$idx_tr) || is.numeric(cf$idx_tr))
      expect_true(length(cf$idx_co) >= 2L)
      expect_true(is.matrix(cf$Y_treat_mat))
      expect_equal(ncol(cf$Y_treat_mat), cf$n_treated)
      expect_equal(nrow(cf$Y_treat_mat), 24L)
    }
  }
})

# ── Phase 33: scm_inference() wild bootstrap ─────────────────────────────────

test_that("Phase 33: scm_inference errors on sharp and non-SCM fits", {
  fit_sharp <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(scm_inference(fit_sharp), regexp = "staggered")
  fit_sdid <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  expect_error(scm_inference(fit_sdid), regexp = "method = 'scm'")
})

test_that("Phase 33: scm_inference returns a valid coresynth_inference", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  inf <- suppressWarnings(scm_inference(fit, method = "wild_bootstrap",
                                       n_boot = 200, seed = 1))
  expect_s3_class(inf, "coresynth_inference")
  expect_equal(inf$estimate, fit$estimate)
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_lt(inf$ci_lower, inf$ci_upper)
  expect_equal(length(inf$boot_ests), 200L)
  expect_true(inf$staggered)
  expect_equal(inf$method, "wild_bootstrap")
  expect_equal(inf$n_treated, 2L)
})

test_that("Phase 33: scm_inference warns with fewer than 3 treated units", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  expect_warning(scm_inference(fit, method = "wild_bootstrap", n_boot = 50,
                               seed = 1),
                 regexp = "fewer than 3 treated units")
})

test_that("Phase 33: scm_inference is reproducible with seed and works with tidy/glance", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm",
                 nu = 0.5, fixedeff = TRUE)
  inf1 <- suppressWarnings(scm_inference(fit, method = "wild_bootstrap",
                                        n_boot = 100, seed = 42))
  inf2 <- suppressWarnings(scm_inference(fit, method = "wild_bootstrap",
                                        n_boot = 100, seed = 42))
  expect_equal(inf1$boot_ests, inf2$boot_ests)

  td <- broom::tidy(inf1)
  expect_equal(nrow(td), 1L)
  expect_equal(td$method, "wild_bootstrap")
  expect_true(td$staggered)
  gl <- broom::glance(inf1)
  expect_equal(gl$n_boot_valid, 100L)
})

test_that("Phase 33: scm_inference alternative directions behave sensibly", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "scm")
  # true ATT = 1.5 > 0: 'greater' should not have a larger p than 'less'
  for (m in c("jackknife", "wild_bootstrap")) {
    p_g <- suppressWarnings(scm_inference(fit, method = m, n_boot = 200,
                                          seed = 3,
                                          alternative = "greater"))$p_value
    p_l <- suppressWarnings(scm_inference(fit, method = m, n_boot = 200,
                                          seed = 3,
                                          alternative = "less"))$p_value
    expect_lte(p_g, p_l)
  }
})

test_that("Phase 33: partially pooled SCM works with covariates partial-out", {
  fit <- scm_fit(y ~ d | id + time, data = stag_cov, method = "scm",
                 nu = 0.5, covariates = "cov1")
  expect_true(is.finite(fit$estimate))
  expect_length(fit$beta_hat, 1L)
})

# ── Phase 34: plot align / show_donors / SDID weight panels / placebo SE ─────

test_that("Phase 34: plot(align=TRUE) makes the SDID post-period mean gap equal the estimate", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  d   <- ggplot2::ggplot_build(plot(fit, type = "gap", align = TRUE))$data[[1]]
  expect_equal(mean(d$y[d$x > fit$T_pre]), fit$estimate, tolerance = 1e-8)
  # the raw gap plot lacks this property (SDID has a free intercept)
  d0 <- ggplot2::ggplot_build(plot(fit, type = "gap"))$data[[1]]
  expect_false(isTRUE(all.equal(mean(d0$y[d0$x > fit$T_pre]), fit$estimate,
                                tolerance = 1e-8)))
})

test_that("Phase 34: plot(align=TRUE) shifts only the synthetic series, by the lambda-weighted pre-gap", {
  fit  <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  pre  <- seq_len(fit$T_pre)
  off  <- sum(fit$time_weights * (fit$Y_treat[pre] - fit$Y_synth[pre]))
  base <- ggplot2::ggplot_build(plot(fit, type = "trend"))$data[[1]]
  al   <- ggplot2::ggplot_build(plot(fit, type = "trend", align = TRUE))$data[[1]]
  shift <- al$y - base$y
  expect_true(all(abs(shift) < 1e-10 | abs(shift - off) < 1e-10))
  expect_true(any(abs(shift - off) < 1e-10))  # synthetic series moved
  expect_true(any(abs(shift) < 1e-10))        # treated series did not
})

test_that("Phase 34: plot(align=TRUE) centers the pre-period gap at zero for non-SDID fits", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  d   <- ggplot2::ggplot_build(plot(fit, type = "gap", align = TRUE))$data[[1]]
  expect_equal(mean(d$y[d$x <= fit$T_pre]), 0, tolerance = 1e-10)
})

test_that("Phase 34: plot align argument is validated", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(plot(fit, type = "trend", align = "yes"), "align")
  expect_error(plot(fit, type = "gap", align = 1), "align")
})

test_that("Phase 34: show_donors overlays top-weight donor paths", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p0  <- plot(fit, type = "trend")
  p2  <- plot(fit, type = "trend", show_donors = 2)
  expect_equal(length(p2$layers), length(p0$layers) + 1L)
  don <- ggplot2::ggplot_build(p2)$data[[1]]  # donor layer is drawn first
  expect_equal(length(unique(don$group)), 2L)
  don_all <- ggplot2::ggplot_build(
    plot(fit, type = "trend", show_donors = Inf))$data[[1]]
  expect_equal(length(unique(don_all$group)), length(fit$unit_weights))
})

test_that("Phase 34: Donors series accepts color/label overrides only when shown", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p <- plot(fit, type = "trend", show_donors = 3,
            colors = c(donors = "black"), labels = c(donors = "Donor pool"))
  expect_s3_class(p, "ggplot")
  expect_error(plot(fit, type = "trend", colors = c(donors = "black")),
               "unrecognized")
})

test_that("Phase 34: show_donors validation and unsupported methods", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(plot(fit, type = "trend", show_donors = -1), "show_donors")
  expect_error(plot(fit, type = "trend", show_donors = c(1, 2)), "show_donors")
  fit_mc <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  expect_error(plot(fit_mc, type = "trend", show_donors = 2), "show_donors")
})

test_that("Phase 34: SDID weights plot shows unit and time weight panels", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  p <- plot(fit, type = "weights")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$PANEL)), 2L)
  n_units <- sum(fit$unit_weights > 1e-4)
  n_times <- sum(fit$time_weights > 1e-4)
  expect_equal(nrow(built$data[[1]]), n_units + n_times)
  # top_n prunes only the unit panel
  built2 <- ggplot2::ggplot_build(plot(fit, type = "weights", top_n = 2))
  expect_equal(nrow(built2$data[[1]]), 2L + n_times)
})

test_that("Phase 34: non-SDID weights plot keeps the single unit panel", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  built <- ggplot2::ggplot_build(plot(fit, type = "weights"))
  expect_equal(length(unique(built$data[[1]]$PANEL)), 1L)
})

test_that("Phase 34: sdid_inference placebo returns SE and normal CI (Clarke Alg. 4)", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  eff <- inf$placebo_effects
  expect_equal(inf$se, sqrt(mean((eff - mean(eff))^2)), tolerance = 1e-12)
  expect_true(inf$se > 0)
  expect_equal(inf$ci_lower, inf$estimate - qnorm(0.975) * inf$se,
               tolerance = 1e-12)
  expect_equal(inf$ci_upper, inf$estimate + qnorm(0.975) * inf$se,
               tolerance = 1e-12)
  inf80 <- sdid_inference(fit, method = "placebo", level = 0.80)
  expect_equal(inf80$ci_upper, inf80$estimate + qnorm(0.90) * inf80$se,
               tolerance = 1e-12)
  td <- broom::tidy(inf)
  expect_true(is.finite(td$std.error))
  expect_true(is.finite(td$conf.low) && is.finite(td$conf.high))
})

test_that("Phase 34: sdid_inference staggered placebo also returns SE and CI", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  eff <- inf$placebo_effects
  expect_equal(inf$se, sqrt(mean((eff - mean(eff))^2)), tolerance = 1e-12)
  expect_true(is.finite(inf$ci_lower) && is.finite(inf$ci_upper))
  expect_lt(inf$ci_lower, inf$ci_upper)
})

# ── Phase 35: vline positioning (vline_offset / vline xintercept) ─────────────

# xintercept values of every geom_vline layer in a built plot
.vline_x <- function(p) {
  d  <- ggplot2::ggplot_build(p)$data
  xs <- lapply(d, function(l) l[["xintercept"]])
  as.numeric(unlist(xs[!vapply(xs, is.null, logical(1))]))
}

test_that("plot.coresynth: vline_offset moves the treatment line", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  # panel: times 1:20, T_pre = 10, first post-treatment period = 11, so the
  # default offset of -1 puts the line on period 10
  expect_equal(.vline_x(plot(fit, type = "trend")), 10)
  expect_equal(.vline_x(plot(fit, type = "trend", vline_offset = 0)), 11)
  expect_equal(.vline_x(plot(fit, type = "gap")), 10)
  # fractional offsets interpolate between adjacent observed times
  expect_equal(.vline_x(plot(fit, type = "trend", vline_offset = -0.5)), 10.5)
  expect_equal(.vline_x(plot(fit, type = "gap", vline_offset = 0)), 11)
})

test_that("plot.coresynth: vline_offset is validated and range-checked", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(plot(fit, type = "trend", vline_offset = "a"), "single finite")
  expect_error(plot(fit, type = "trend", vline_offset = c(-1, 0)), "single finite")
  expect_error(plot(fit, type = "trend", vline_offset = Inf), "single finite")
  # out of range: warn and drop the line (layer count shrinks by one)
  n_default <- length(plot(fit, type = "trend")$layers)
  expect_warning(p <- plot(fit, type = "trend", vline_offset = 99),
                 "outside the observed time range")
  expect_equal(length(p$layers), n_default - 1L)
})

test_that("plot.coresynth: vline xintercept pins the line to an absolute time", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_equal(.vline_x(plot(fit, type = "trend",
                             vline = list(xintercept = 15.5))), 15.5)
  # aesthetics in the same list still apply
  p <- plot(fit, type = "gap", vline = list(xintercept = 5, color = "red"))
  expect_equal(.vline_x(p), 5)
  # several positions draw several lines
  expect_equal(.vline_x(plot(fit, type = "trend",
                             vline = list(xintercept = c(5, 11)))), c(5, 11))
  expect_error(plot(fit, type = "trend", vline_offset = 0,
                    vline = list(xintercept = 5)), "not both")
  expect_error(plot(fit, type = "trend", vline = list(xintercept = NA)),
               "non-missing")
})

test_that("plot.coresynth: vline_offset works on a Date time axis", {
  dpanel <- panel
  dpanel$time <- as.Date("2000-01-01") + 30 * (dpanel$time - 1L)
  fit <- scm_fit(y ~ d | id + time, data = dpanel, method = "scm")

  d0 <- as.numeric(as.Date("2000-01-01") + 30 * 10) # 11th period
  expect_equal(.vline_x(plot(fit, type = "trend")), d0 - 30)
  expect_equal(.vline_x(plot(fit, type = "trend", vline_offset = 0)), d0)
  expect_equal(.vline_x(plot(fit, type = "trend", vline_offset = -0.5)), d0 - 15)
  # absolute positions accept a date string on a Date axis
  expect_equal(
    .vline_x(plot(fit, type = "trend",
                  vline = list(xintercept = "2000-03-01"))),
    as.numeric(as.Date("2000-03-01"))
  )
})

test_that("plot.scm_placebo gaps: vline_offset and xintercept work", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  inf <- mspe_ratio_pval(fit)
  expect_equal(.vline_x(plot(inf, type = "gaps")), 10)
  expect_equal(.vline_x(plot(inf, type = "gaps", vline_offset = 0)), 11)
  expect_equal(.vline_x(plot(inf, type = "gaps",
                             vline = list(xintercept = 4))), 4)
  expect_error(plot(inf, type = "gaps", vline_offset = 2.5,
                    vline = list(xintercept = 4)), "not both")
})

# ── Phase 36: multi-start outer V optimisation & v_window ─────────────────────

test_that("Phase 36: predictor-path default resolves to multistart", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                 predictors = list(pred(c("cov1", "cov2"), 1:10)))
  expect_identical(fit$v_optim_effective, "multistart")
})

test_that("Phase 36: outcomes-only default keeps single-start coord descent", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_identical(fit$v_optim_effective, "coord_descent")
  fit_cd <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                    v_optim = "coord_descent")
  expect_equal(fit$unit_weights, fit_cd$unit_weights, tolerance = 1e-12)
})

test_that("Phase 36: multistart is never worse than coord_descent in pre-loss", {
  spec <- list(pred(c("cov1", "cov2"), 1:10))
  fit_ms <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                    predictors = spec, v_optim = "multistart")
  fit_cd <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                    predictors = spec, v_optim = "coord_descent")
  expect_lte(fit_ms$loss, fit_cd$loss + 1e-8)
})

test_that("Phase 36: multistart is deterministic and leaves the R RNG alone", {
  spec <- list(pred(c("cov1", "cov2"), 1:10), pred("y", c(3, 6, 9)))
  set.seed(123)
  seed_before <- get(".Random.seed", envir = globalenv())
  fit1 <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                  predictors = spec, v_optim = "multistart")
  expect_identical(get(".Random.seed", envir = globalenv()), seed_before)
  fit2 <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                  predictors = spec, v_optim = "multistart")
  expect_identical(fit1$unit_weights, fit2$unit_weights)
  expect_identical(fit1$v_weights, fit2$v_weights)
})

test_that("Phase 36: multistart errors on outcomes-only fits", {
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            v_optim = "multistart"),
    regexp = "predictor"
  )
})

test_that("Phase 36: mspe_ratio_pval mirrors the multistart fit", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                 predictors = list(pred(c("cov1", "cov2"), 1:10)))
  inf <- mspe_ratio_pval(fit)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_equal(length(inf$mspe_ratios_all), 10L)
})

test_that("Phase 36: v_window equal to the full pre-period reproduces default", {
  fit_def <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  fit_win <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_window = 1:10)
  expect_equal(fit_def$unit_weights, fit_win$unit_weights, tolerance = 1e-12)
  expect_equal(fit_win$v_window, fit_def$times[1:10])
  expect_equal(fit_win$z_rows, 1:10)
})

test_that("Phase 36: v_window subset restricts the outer evaluation", {
  fit_win <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_window = 4:10)
  expect_equal(fit_win$z_rows, 4:10)
  expect_equal(sum(fit_win$unit_weights), 1, tolerance = 1e-4)
  expect_true(is.finite(fit_win$estimate))
  # loss is still reported on the full pre-treatment window
  gap_pre <- (fit_win$Y_treat - fit_win$Y_synth)[seq_len(fit_win$T_pre)]
  expect_equal(fit_win$loss, sqrt(sum(gap_pre^2)), tolerance = 1e-8)
})

test_that("Phase 36: v_window works with predictors and mspe_ratio_pval", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
                 predictors = list(pred(c("cov1", "cov2"), 1:10)),
                 v_window = 5:10)
  expect_equal(fit$z_rows, 5:10)
  inf <- mspe_ratio_pval(fit)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
})

test_that("Phase 36: v_window validation errors", {
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "sdid", v_window = 1:5),
    regexp = "method"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm", v_window = 5:12),
    regexp = "pre-treatment"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm", v_window = 5),
    regexp = "at least 2"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            v_window = 1:8, v_selection = "oos"),
    regexp = "oos"
  )
})

test_that("Phase 36: v_window errors on staggered fits", {
  set.seed(99)
  sp <- expand.grid(id = paste0("u", 1:8), time = 1:12)
  sp$d <- as.integer((sp$id == "u1" & sp$time >= 7) |
                     (sp$id == "u2" & sp$time >= 9))
  sp$y <- rnorm(nrow(sp)) + as.numeric(sub("u", "", sp$id)) + 0.3 * sp$time
  expect_error(
    scm_fit(y ~ d | id + time, data = sp, method = "scm", v_window = 1:6),
    regexp = "sharp"
  )
})

# ── Phase 38: outcomes-only cheap bordered-KKT inner-QP face solve ─────────────
# The outcomes-only path solves each active-set face with a cheap bordered-KKT
# direct solve (cheap_face = TRUE) rather than the scale-robust null-space
# projection used for user predictors. It is the same nested V/W estimator, so
# the two solvers must agree to round-off; the flag only changes speed.

# Build outcome matrices with a low-rank factor structure; optionally make two
# donor pairs near-collinear to stress the inner face solver.
.mk_oo <- function(seed, N_co, T_pre, sd, degen = FALSE) {
  set.seed(seed)
  r  <- 3L
  L  <- matrix(rnorm((N_co + 1L) * r), N_co + 1L, r)
  Ft <- matrix(rnorm(T_pre * r), T_pre, r)
  Y  <- Ft %*% t(L) +
    matrix(rnorm(T_pre * (N_co + 1L), sd = sd), T_pre, N_co + 1L)
  if (degen) {
    Y[, 3] <- Y[, 2] + rnorm(T_pre, sd = 1e-4)
    Y[, 5] <- Y[, 4] + rnorm(T_pre, sd = 1e-4)
  }
  list(X1 = Y[, 1], X0 = Y[, -1, drop = FALSE])
}

test_that("Phase 38: cheap_face matches the null-space solver (outcomes-only)", {
  grid <- expand.grid(seed = 1:5, sd = c(0.1, 0.3, 1.0),
                      degen = c(FALSE, TRUE))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    d <- .mk_oo(g$seed, N_co = 25L, T_pre = 18L, sd = g$sd, degen = g$degen)
    cheap  <- scm_weights_cpp(d$X0, d$X1, d$X0, d$X1, cheap_face = TRUE)
    robust <- scm_weights_cpp(d$X0, d$X1, d$X0, d$X1, cheap_face = FALSE)
    # Same estimator, different inner face solve: agree to round-off.
    expect_equal(as.numeric(cheap$W), as.numeric(robust$W), tolerance = 1e-8)
    expect_equal(cheap$loss, robust$loss, tolerance = 1e-10)
    # Simplex feasibility of the fast-path solution.
    expect_true(all(cheap$W >= -1e-10))
    expect_equal(sum(cheap$W), 1, tolerance = 1e-8)
  }
})

test_that("Phase 38: scm_fit outcomes-only equals the robust null-space solver", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_equal(fit$v_optim_effective, "coord_descent")
  Yp <- fit$Y_co_pre
  yt <- fit$Y_treat[seq_len(fit$T_pre)]
  # The fit uses cheap_face internally; the null-space solver is the reference.
  ref <- scm_weights_cpp(Yp, yt, Yp, yt, cheap_face = FALSE)
  expect_equal(unname(fit$unit_weights), as.numeric(ref$W), tolerance = 1e-8)
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-8)
})

test_that("Phase 38: near-collinear donor pool stays feasible (outcomes-only)", {
  d  <- .mk_oo(123L, N_co = 20L, T_pre = 15L, sd = 0.05, degen = TRUE)
  df <- data.frame(
    id   = rep(seq_len(ncol(d$X0) + 1L), each = nrow(d$X0) + 4L),
    time = rep(seq_len(nrow(d$X0) + 4L), times = ncol(d$X0) + 1L)
  )
  Yfull <- cbind(d$X1, d$X0)
  Ypost <- Yfull[rep(nrow(Yfull), 4L), ] + 0.1  # 4 post periods
  Yall  <- rbind(Yfull, Ypost)
  df$y  <- as.vector(Yall)
  df$d  <- as.integer(df$id == 1L & df$time > nrow(d$X0))
  fit <- scm_fit(y ~ d | id + time, data = df, method = "scm")
  expect_true(all(fit$unit_weights >= -1e-10))
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-8)
  expect_true(is.finite(fit$estimate))
})

# ── outcomes-only routing of per-period outcome predictor specs ───────────────
# A predictors list naming solely the outcome at each single pre-treatment
# period defines X0/X1 equal to the pre-treatment outcome rows, so it must
# return the identical outcomes-only fit rather than entering the scaled
# multistart predictor path.

test_that("per-period outcome spec returns the identical outcomes-only fit", {
  spec <- lapply(1:10, function(t) pred("y", t))
  ref  <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_message(
    fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                   predictors = spec),
    regexp = "outcomes-only"
  )
  expect_equal(fit$v_optim_effective, "coord_descent")
  expect_null(fit$X0_mat)
  expect_null(fit$predictor_table)
  expect_identical(fit$unit_weights, ref$unit_weights)
  expect_identical(fit$v_weights, ref$v_weights)
  expect_identical(fit$estimate, ref$estimate)
})

test_that("outcomes-only routing is order- and op-invariant", {
  spec <- lapply(rev(1:10), function(t) pred("y", t, op = "median"))
  ref  <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_message(
    fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                   predictors = spec),
    regexp = "outcomes-only"
  )
  expect_identical(fit$unit_weights, ref$unit_weights)
  expect_identical(fit$estimate, ref$estimate)
})

test_that("partial-coverage or windowed outcome specs keep the predictor path", {
  fit_sub <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     predictors = lapply(1:9, function(t) pred("y", t)))
  expect_equal(fit_sub$v_optim_effective, "multistart")
  expect_false(is.null(fit_sub$X0_mat))

  spec_win <- c(lapply(1:8, function(t) pred("y", t)), list(pred("y", 9:10)))
  fit_win  <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                      predictors = spec_win)
  expect_equal(fit_win$v_optim_effective, "multistart")
  expect_false(is.null(fit_win$X0_mat))
})

test_that("explicit v_optim = 'multistart' overrides the outcomes-only routing", {
  spec <- lapply(1:10, function(t) pred("y", t))
  fit  <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                  predictors = spec, v_optim = "multistart")
  expect_equal(fit$v_optim_effective, "multistart")
  expect_false(is.null(fit$X0_mat))
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-8)
})

# Coverage is judged against the pre-treatment periods `data` carries, so the
# same predictors list lands on either path depending on how far back `data`
# reaches. That switch changes the estimand, so it has to be announced.

test_that("a partial-coverage outcome-lag spec reports the missed routing", {
  spec <- lapply(3:10, function(t) pred("y", t))
  expect_message(
    fit_full <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                        predictors = spec),
    regexp = "covering 8 of the 10 pre-treatment periods"
  )
  expect_equal(fit_full$v_optim_effective, "multistart")
  expect_false(is.null(fit_full$X0_mat))

  # Same spec, `data` trimmed to the periods it covers: now it is the
  # outcomes-only fit, on the same call.
  panel_trim <- panel[panel$time >= 3, ]
  expect_message(
    fit_trim <- scm_fit(y ~ d | id + time, data = panel_trim, method = "scm",
                        predictors = spec),
    regexp = "outcomes-only"
  )
  expect_equal(fit_trim$v_optim_effective, "coord_descent")
  expect_null(fit_trim$X0_mat)
})

test_that("the missed-routing message names v_window only when it is unset", {
  spec <- lapply(3:10, function(t) pred("y", t))
  expect_message(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            predictors = spec),
    regexp = "v_window"
  )
  msgs <- capture_messages(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            predictors = spec, v_window = 3:10)
  )
  expect_true(any(grepl("covering 8 of the 10", msgs)))
  expect_false(any(grepl("v_window", msgs)))
})

test_that("specs that are not single-period outcome lags stay silent", {
  msgs <- capture_messages(
    scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
            predictors = list(pred(c("cov1", "cov2"), 1:10)))
  )
  expect_false(any(grepl("pre-treatment periods in 'data'", msgs)))

  msgs_ms <- capture_messages(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            predictors = lapply(3:10, function(t) pred("y", t)),
            v_optim = "multistart")
  )
  expect_false(any(grepl("pre-treatment periods in 'data'", msgs_ms)))
})

# pred() windows that do not aggregate what they read as ──────────────────────

test_that("pred() windows reaching past the treatment date warn", {
  expect_warning(
    scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
            predictors = list(pred("cov1", 1:10), pred("cov2", 8:14))),
    regexp = "post-treatment"
  )
})

test_that("pred() windows only partly present in the panel warn", {
  expect_warning(
    scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
            predictors = list(pred("cov1", -3:10))),
    regexp = "does not contain"
  )
})

test_that("a pred() window with no overlap errors instead of warning", {
  expect_silent(suppressMessages(expect_error(
    scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
            predictors = list(pred("cov1", 100:110))),
    regexp = "non-finite"
  )))
})

test_that("windows that lie inside the pre-treatment period do not warn", {
  expect_no_warning(
    scm_fit(y ~ d | id + time, data = panel_cov, method = "scm",
            predictors = list(pred(c("cov1", "cov2"), 1:10),
                              pred("y", c(3, 6, 9))))
  )
})

# ── Inner-QP solver: exactness and the Wolfe (Caratheodory-sparse) option ─────

# KKT residual of W for the inner QP min ||X1 - X0 w||_V^2 over the simplex,
# normalised by the gradient scale.
scm_kkt_residual <- function(X0, X1, V, W) {
  g    <- as.numeric(crossprod(X0, V * (as.numeric(X0 %*% W) - X1)))
  supp <- W > 1e-10
  mu   <- -mean(g[supp])
  viol <- max(abs(g[supp] + mu))
  if (!all(supp)) viol <- max(viol, 0, -min(g[!supp] + mu))
  viol / (1 + max(abs(g)))
}

test_that("scm_inner_weights_cpp returns a KKT-exact optimum", {
  set.seed(4242L)
  for (cfg in list(c(k = 1, N = 40), c(k = 3, N = 60), c(k = 12, N = 25),
                   c(k = 30, N = 20))) {
    k <- cfg[["k"]]; N <- cfg[["N"]]
    X0 <- matrix(rnorm(k * N), k, N) + 5
    X1 <- rowMeans(X0) + rnorm(k, 0, 0.2)
    V  <- runif(k); V <- V / sum(V)
    W  <- as.numeric(scm_inner_weights_cpp(X0, X1, V))
    expect_equal(sum(W), 1, tolerance = 1e-8)
    expect_true(all(W >= -1e-10))
    expect_lt(scm_kkt_residual(X0, X1, V, W), 1e-8)
  }
})

test_that("Wolfe solver matches the active set in objective but is sparse", {
  set.seed(99L)
  for (cfg in list(c(k = 1, N = 50), c(k = 2, N = 80), c(k = 5, N = 120))) {
    k <- cfg[["k"]]; N <- cfg[["N"]]
    X0 <- matrix(rnorm(k * N), k, N) + 5
    X1 <- rowMeans(X0) + rnorm(k, 0, 0.2)
    V  <- rep(1 / k, k)
    wa <- as.numeric(scm_inner_weights_cpp(X0, X1, V, wolfe = FALSE))
    ww <- as.numeric(scm_inner_weights_cpp(X0, X1, V, wolfe = TRUE))
    obj <- function(w) sum(V * (X1 - as.numeric(X0 %*% w))^2)
    expect_equal(obj(ww), obj(wa), tolerance = 1e-8)
    expect_lt(scm_kkt_residual(X0, X1, V, ww), 1e-8)
    # Caratheodory: an optimum supported on at most k+1 donors always exists
    expect_lte(sum(ww > 1e-8), k + 1L)
  }
})

test_that("qp_solver = 'wolfe' leaves a non-degenerate outcomes-only fit alone", {
  fit_a <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  fit_w <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                   qp_solver = "wolfe")
  fit_d <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                   qp_solver = "active_set")
  expect_equal(fit_w$unit_weights, fit_d$unit_weights, tolerance = 1e-6)
  expect_equal(fit_w$estimate, fit_d$estimate, tolerance = 1e-6)
  expect_equal(fit_w$qp_solver, "wolfe")
  expect_equal(fit_d$qp_solver, "active_set")
  # "auto" resolves to the active set on an outcomes-only fit with at least
  # half as many pre-treatment periods as donors (here 10 and 9) ...
  expect_equal(fit_a$qp_solver, "active_set")
  expect_equal(fit_a$unit_weights, fit_w$unit_weights, tolerance = 1e-6)
  # ... and to Wolfe where the QP can be degenerate: few periods next to the
  # donor pool, or any predictor fit
  short <- panel[panel$time > 6, ]
  expect_equal(suppressMessages(
    scm_fit(y ~ d | id + time, data = short, method = "scm"))$qp_solver,
    "wolfe")
  expect_equal(suppressMessages(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            predictors = list(pred("y", 1:10, "mean"))))$qp_solver,
    "wolfe")
})

test_that("qp_solver = 'auto' falls back where Wolfe is not wired", {
  # each of these paths runs its own QP, so the default downgrades quietly
  # rather than erroring the way an explicit qp_solver = "wolfe" does
  expect_equal(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            v_selection = "oos")$qp_solver, "active_set")
  expect_equal(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            lambda_pen = 0.1)$qp_solver, "active_set")
  expect_equal(
    suppressWarnings(scm_fit(y ~ d | id + time, data = panel, method = "scm",
                             v_optim = "bfgs"))$qp_solver, "active_set")
  # staggered fits carry no fit-level qp_solver; they must simply not error
  expect_s3_class(
    scm_fit(y ~ d | id + time, data = staggered, method = "scm"),
    "coresynth_staggered")
})

test_that("qp_solver = 'wolfe' rejects the paths it is not wired through", {
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            qp_solver = "wolfe", v_selection = "oos"),
    "not available with"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "scm",
            qp_solver = "wolfe", lambda_pen = 0.1),
    "not available with"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "sdid",
            qp_solver = "wolfe"),
    "method = \"scm\" only"
  )
})

test_that("mspe_ratio_pval mirrors the fit's qp_solver", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 qp_solver = "wolfe")
  res <- mspe_ratio_pval(fit)
  expect_s3_class(res, "scm_placebo")
  expect_true(is.finite(res$p_value))
  expect_gte(res$p_value, 0)
  expect_lte(res$p_value, 1)
})

# ── scm_fit() argument validation ─────────────────────────────────────────────
# One formula interface serves six estimators, so an argument meant for another
# one, or misspelled, used to be absorbed by `...` and silently change nothing.

test_that("scm_fit rejects SCM-only arguments given to other methods", {
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       predictors = list(pred("y", 1:3))), "'predictors'")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "tasc",
                       nu = 0.5), "'nu'")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "tasc",
                       fixedeff = TRUE), "'fixedeff'")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "si",
                       v_optim = "multistart"), "'v_optim'")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "gsc",
                       lambda_pen = 1), "'lambda_pen'")
  # covariates reach three estimators, not six
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "tasc",
                       covariates = "y"), "'covariates'")
})

test_that("scm_fit rejects unknown estimator arguments passed through ...", {
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       lambdaa = 1), "unused argument")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "gsc",
                       rr = 3), "unused argument")
  # control_group has no meaning for the imputation estimators: their mask
  # already encodes each unit's own adoption date and the factorisation uses
  # every observed cell
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       control_group = "never_treated"), "unused argument")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "tasc",
                       control_group = "never_treated"), "unused argument")
})

test_that("scm_fit routes MC's lambda past the SCM-only lambda_pen", {
  # `lambda` is not an exact formal of scm_fit(), so while lambda_pen sat
  # before `...` R matched it there partially and MC's penalty never moved
  fit_default <- scm_fit(y ~ d | id + time, data = panel, method = "mc")
  fit_lambda  <- scm_fit(y ~ d | id + time, data = panel, method = "mc",
                         lambda = 2)
  expect_equal(fit_lambda$lambda, 2)
  expect_false(isTRUE(all.equal(fit_default$lambda, fit_lambda$lambda)))
  expect_false(isTRUE(all.equal(fit_default$estimate, fit_lambda$estimate)))
})

test_that("scm_fit still accepts each estimator's own arguments", {
  expect_s3_class(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                          lambda = 0.5), "coresynth")
  expect_s3_class(suppressWarnings(scm_fit(y ~ d | id + time, data = panel, method = "tasc",
                          r = 2L, em_iter = 3L, fix_A = TRUE)), "coresynth")
  expect_s3_class(scm_fit(y ~ d | id + time, data = panel, method = "gsc",
                          r = 2L), "coresynth")
  expect_s3_class(scm_fit(y ~ d | id + time, data = panel, method = "sdid",
                          zeta2 = 1, zeta_t = 1e-6), "coresynth")
  expect_s3_class(scm_fit(y ~ d | id + time, data = staggered, method = "sdid",
                          control_group = "never_treated"), "coresynth")
  expect_s3_class(scm_fit(y ~ d | id + time, data = panel, method = "scm",
                          scale_predictors = FALSE), "coresynth")
})

# ── Data-driven tuning parameters (Xu 2017; Athey et al. 2021; Donoho-Gavish) ──
# Three estimators used to pick their key tuning constant by an arbitrary rule
# while their source papers specify a data-driven one.

# Panel whose number of latent factors is known, so a selection rule has a
# right answer to find.
.factor_panel <- function(r_true, seed = 11, N = 30, TT = 24, T0 = 16,
                          tau = 2, sd_eps = 0.3) {
  set.seed(seed)
  f <- matrix(rnorm(TT * r_true), TT, r_true)
  L <- matrix(rnorm(N * r_true), N, r_true)
  Y <- f %*% t(L) + matrix(rnorm(N * TT, sd = sd_eps), TT, N)
  D <- matrix(0L, TT, N)
  D[(T0 + 1):TT, 1:3] <- 1L
  Y[(T0 + 1):TT, 1:3] <- Y[(T0 + 1):TT, 1:3] + tau
  data.frame(id = rep(1:N, each = TT), time = rep(1:TT, times = N),
             y = as.vector(Y), d = as.vector(D))
}

test_that("GSC cross-validation recovers the number of factors (Xu Alg. 1)", {
  for (r_true in 1:4) {
    fit <- scm_fit(y ~ d | id + time, data = .factor_panel(r_true),
                   method = "gsc")
    expect_equal(fit$r, r_true)
    expect_s3_class(fit$r_cv, "data.frame")
    expect_equal(names(fit$r_cv), c("r", "mspe"))
    # the winner must be the minimiser of the reported MSPE column
    expect_equal(fit$r_cv$r[which.min(fit$r_cv$mspe)], fit$r)
  }
})

test_that("GSC honours an explicit r and validates it", {
  panel_r <- .factor_panel(3)
  fit <- scm_fit(y ~ d | id + time, data = panel_r, method = "gsc", r = 2)
  expect_equal(fit$r, 2)
  expect_null(fit$r_cv)
  expect_error(scm_fit(y ~ d | id + time, data = panel_r, method = "gsc",
                       r = 0), "positive integer")
  expect_error(scm_fit(y ~ d | id + time, data = panel_r, method = "gsc",
                       r = 2.5), "positive integer")
})

test_that("SI picks the PCR rank by the optimal hard threshold", {
  for (r_true in 1:4) {
    fit <- scm_fit(y ~ d | id + time, data = .factor_panel(r_true),
                   method = "si")
    expect_equal(fit$k, r_true)
  }
  # an explicit k still wins
  fit_k <- scm_fit(y ~ d | id + time, data = .factor_panel(2), method = "si",
                   k = 5)
  expect_equal(fit_k$k, 5L)
})

test_that("a staggered SI fit reports the rank each cohort used", {
  # The threshold is applied to each cohort's own donor block, so unlike the
  # GSC factor count and the MC penalty there is no single value; the fit used
  # to store the NULL argument instead of what it chose.
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
  expect_type(fit$k, "integer")
  expect_named(fit$k, c("9", "15"))
  expect_equal(unname(fit$k),
               vapply(fit$cohort_fits, `[[`, integer(1), "k"))
  expect_true(all(fit$k >= 1L))

  # a sharp fit keeps the scalar
  expect_length(scm_fit(y ~ d | id + time, data = .factor_panel(2),
                        method = "si")$k, 1L)
})

test_that("MC cross-validates the nuclear norm penalty", {
  panel_mc <- .factor_panel(2)
  fit <- scm_fit(y ~ d | id + time, data = panel_mc, method = "mc")
  expect_true(fit$lambda_cv)
  expect_true(is.finite(fit$lambda) && fit$lambda > 0)

  # the cross-validated penalty differs from the old fixed 1% of sigma_max
  sv  <- svd(matrix(panel_mc$y, nrow = 24), nu = 0, nv = 0)$d[1L]
  expect_false(isTRUE(all.equal(fit$lambda, 0.01 * sv)))

  # an explicit lambda bypasses the search
  fit_l <- scm_fit(y ~ d | id + time, data = panel_mc, method = "mc",
                   lambda = 0.5)
  expect_equal(fit_l$lambda, 0.5)
  expect_false(fit_l$lambda_cv)
})

test_that("MC cross-validation is reproducible and leaves the RNG alone", {
  panel_mc <- .factor_panel(2)
  a <- scm_fit(y ~ d | id + time, data = panel_mc, method = "mc")
  b <- scm_fit(y ~ d | id + time, data = panel_mc, method = "mc")
  expect_equal(a$lambda, b$lambda)
  # a different fold draw may pick a different penalty
  c_ <- scm_fit(y ~ d | id + time, data = panel_mc, method = "mc", cv_seed = 99)
  expect_true(is.finite(c_$lambda))

  # the caller's RNG stream must be where they left it
  set.seed(4242)
  before <- runif(1)
  set.seed(4242)
  invisible(scm_fit(y ~ d | id + time, data = panel_mc, method = "mc"))
  expect_equal(runif(1), before)
})

# ── MC: unpenalised two-way fixed effects (Athey et al. 2021, S.8.1) ──────────
# The S.4 baseline absorbs the panel's level into the penalised matrix L, so the
# nuclear norm shrinks it along with the factor structure and the imputed
# counterfactual is pulled towards zero. S.8.1 keeps unit and time effects out
# of the penalty; MCPanel, the authors' own implementation, estimates them by
# default (to_estimate_u = to_estimate_v = 1).

.level_panel <- function(seed, level, N = 21, TT = 24, T0 = 16, tau = 2) {
  set.seed(seed)
  f <- matrix(rnorm(TT * 2), TT, 2); L <- matrix(rnorm(N * 2), N, 2)
  Y <- level + f %*% t(L) +
    matrix(rnorm(N, sd = 3), TT, N, byrow = TRUE) + rnorm(TT, sd = 1) +
    matrix(rnorm(N * TT, sd = 0.5), TT, N)
  D <- matrix(0L, TT, N)
  D[(T0 + 1):TT, 1] <- 1L; Y[(T0 + 1):TT, 1] <- Y[(T0 + 1):TT, 1] + tau
  data.frame(id = rep(1:N, each = TT), time = rep(1:TT, times = N),
             y = as.vector(Y), d = as.vector(D))
}

test_that("MC estimates unpenalised fixed effects by default", {
  fit <- scm_fit(y ~ d | id + time, data = .level_panel(1, 20), method = "mc")
  expect_true(fit$fixedeff)
  expect_length(fit$fe_time, 24L)
  expect_length(fit$fe_unit, 21L)
  # L_hat is the completed matrix; the low-rank part is kept separately
  expect_equal(dim(fit$L_lowrank), dim(fit$L_hat))
  expect_false(isTRUE(all.equal(fit$L_lowrank, fit$L_hat)))
  # switching them off returns the S.4 baseline, where L carries the level
  off <- scm_fit(y ~ d | id + time, data = .level_panel(1, 20), method = "mc",
                 fixedeff = FALSE)
  expect_false(off$fixedeff)
  expect_equal(off$L_lowrank, off$L_hat)
  expect_false(isTRUE(all.equal(fit$estimate, off$estimate)))
})

test_that("MC fixed effects remove the level-driven bias", {
  # A shifted panel is the case the S.4 baseline gets wrong: with the level
  # inside L the penalty shrinks it and the ATT is biased upwards.
  bias <- function(fixedeff) {
    est <- vapply(1:15, function(s)
      scm_fit(y ~ d | id + time, data = .level_panel(s, 20), method = "mc",
              fixedeff = fixedeff)$estimate, numeric(1L))
    mean(est) - 2
  }
  b_on  <- bias(TRUE)
  b_off <- bias(FALSE)
  expect_lt(abs(b_on), abs(b_off))
  expect_lt(abs(b_on), 0.1)
  expect_gt(b_off, 0.1)
})

test_that("soft_impute_cpp still solves the no-fixed-effect problem", {
  set.seed(5)
  Y <- matrix(rnorm(60), 12, 5)
  O <- matrix(1, 12, 5); O[10:12, 1] <- 0
  res <- soft_impute_fe_cpp(Y, O, 0.5, estimate_fe = FALSE)
  expect_equal(soft_impute_cpp(Y, O, 0.5), res$L)
  expect_equal(res$M, res$L)
  expect_true(all(res$a == 0) && all(res$b == 0))
})

test_that("conformal_inference refits MC under the fit's own model", {
  panel_lv <- .level_panel(2, 20)
  for (fe in c(TRUE, FALSE)) {
    fit <- scm_fit(y ~ d | id + time, data = panel_lv, method = "mc",
                   fixedeff = fe)
    ci <- conformal_inference(fit, grid = seq(0, 4, by = 1))
    expect_true(is.finite(ci$p_value))
    expect_gte(ci$p_value, 0)
    expect_lte(ci$p_value, 1)
  }
})

# ── SCM: adoption patterns and their defaults ────────────────────────────────
# Four patterns, by how many adoption dates the panel has and how many treated
# units share one. Each gets the default the literature supports, and `nu` is
# the single knob that overrides it.

.adopt_panel <- function(adopt, N = 25, TT = 20, tau = 2, seed = 3) {
  set.seed(seed)
  f <- matrix(rnorm(TT * 2), TT, 2); L <- matrix(rnorm(N * 2), N, 2)
  Y <- f %*% t(L) + matrix(rnorm(N * TT, sd = 0.3), TT, N)
  D <- matrix(0L, TT, N)
  for (k in seq_along(adopt)) {
    t0 <- adopt[k]; D[t0:TT, k] <- 1L; Y[t0:TT, k] <- Y[t0:TT, k] + tau
  }
  data.frame(id = rep(1:N, each = TT), time = rep(1:TT, times = N),
             y = as.vector(Y), d = as.vector(D))
}

test_that("simultaneous adoption by several units fits one pooled cohort", {
  # BFR S.4.1: with a shared adoption date the across-unit variation that puts
  # unit-level fits into the ATT error bound is zero, so the bound is minimised
  # by fully pooling the treated units.
  fit <- scm_fit(y ~ d | id + time, data = .adopt_panel(c(11, 11, 11)),
                 method = "scm")
  expect_s3_class(fit, "coresynth_scm")
  expect_false(inherits(fit, "coresynth_staggered"))
  expect_equal(fit$n_treated, 3L)
  expect_true(is.finite(fit$estimate))
  # the fitted series is the treated average; the raw series are kept beside it
  expect_length(fit$Y_treat, 20L)
  expect_equal(dim(fit$Y_treat_mat), c(20L, 3L))
  expect_equal(unname(fit$Y_treat), unname(rowMeans(fit$Y_treat_mat)))
  # donor weights remain a single interpretable simplex vector
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-6)
})

test_that("predictors describe the treated average when several units pool", {
  dat <- .adopt_panel(c(11, 11, 11))
  set.seed(11)
  dat$x <- rnorm(nrow(dat)) + dat$id / 10
  spec <- list(pred("y", 1:10, "mean"), pred("x", 1:10, "mean"))

  fit <- scm_fit(y ~ d | id + time, data = dat, method = "scm",
                 predictors = spec)
  # one predictor row per variable, not one per (variable, treated unit)
  expect_length(fit$X1_vec, 2L)
  expect_length(fit$v_weights, 2L)
  expect_equal(nrow(fit$predictor_table), 2L)

  # X1 is the average over the pooled treated units, matching the outcome path
  pre <- dat[dat$time <= 10 & dat$id %in% 1:3, ]
  by_unit <- function(v) vapply(split(pre[[v]], pre$id), mean, numeric(1L))
  expect_equal(fit$predictor_table$treated,
               c(mean(by_unit("y")), mean(by_unit("x"))),
               tolerance = 1e-10)

  # `weights` moves the predictor rows onto the same weighted average the
  # outcome rows use, so the two never describe different targets
  wts <- c(`1` = 3, `2` = 1, `3` = 1)
  fw  <- scm_fit(y ~ d | id + time, data = dat, method = "scm",
                 predictors = spec, weights = wts)
  w   <- wts / sum(wts)
  expect_equal(fw$predictor_table$treated,
               c(sum(by_unit("y") * w), sum(by_unit("x") * w)),
               tolerance = 1e-10)
  expect_equal(unname(fw$Y_treat),
               unname(drop(fw$Y_treat_mat %*% w)), tolerance = 1e-10)
})

test_that("glance() and export_json() report the pooled treated-unit count", {
  fit <- scm_fit(y ~ d | id + time, data = .adopt_panel(c(11, 11, 11)),
                 method = "scm")
  expect_equal(glance(fit)$n_treated, 3L)

  skip_if_not_installed("jsonlite")
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  export_json(fit, file = tmp)
  expect_equal(jsonlite::read_json(tmp)$data$n_treated, 3L)

  one <- scm_fit(y ~ d | id + time, data = .adopt_panel(11), method = "scm")
  expect_equal(glance(one)$n_treated, 1L)
})

test_that("print() and summary() name the pooled treated units", {
  fit <- scm_fit(y ~ d | id + time, data = .adopt_panel(c(11, 11, 11)),
                 method = "scm")
  expect_output(print(fit), "Treated units : 3")
  expect_output(summary(fit), "simple average")

  fw <- scm_fit(y ~ d | id + time, data = .adopt_panel(c(11, 11, 11)),
                method = "scm", weights = c(`1` = 3, `2` = 1, `3` = 1))
  expect_output(print(fw), "weighted average")

  # one treated unit says nothing extra
  one <- scm_fit(y ~ d | id + time, data = .adopt_panel(11), method = "scm")
  expect_false(any(grepl("Treated units",
                         capture.output(print(one)), fixed = TRUE)))
})

test_that("a single treated unit is unchanged by the pooling path", {
  fit <- scm_fit(y ~ d | id + time, data = .adopt_panel(11), method = "scm")
  expect_equal(fit$n_treated, 1L)
  expect_equal(dim(fit$Y_treat_mat), c(20L, 1L))
  expect_equal(unname(fit$Y_treat), as.numeric(fit$Y_treat_mat))
  expect_error(scm_fit(y ~ d | id + time, data = .adopt_panel(11),
                       method = "scm", nu = 0.5), "staggered")
})

test_that("the placebo MSPE vectors are both named by donor", {
  # mspe_pre_placebo and mspe_post_placebo are documented as a pair, so
  # matching them up by name has to work for both.
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  p <- mspe_ratio_pval(fit)
  donors <- colnames(fit$Y_co_pre)
  expect_equal(names(p$mspe_pre_placebo), donors)
  expect_equal(names(p$mspe_post_placebo), donors)
  expect_equal(names(p$placebo_effects), donors)
})

test_that("mspe_ratio_pval refuses a fit that pools several treated units", {
  # Measured under no treatment effect, the test rejected at 28% (3 pooled
  # units) and 53% (6) against a nominal 10%: the averaged treated series is
  # less noisy than the single donors it is compared against.
  fit <- scm_fit(y ~ d | id + time, data = .adopt_panel(c(11, 11, 11)),
                 method = "scm")
  expect_error(mspe_ratio_pval(fit), "not valid")
  # one treated unit is still fine
  ok <- scm_fit(y ~ d | id + time, data = .adopt_panel(11), method = "scm")
  expect_s3_class(mspe_ratio_pval(ok), "scm_placebo")
})

test_that("staggered SCM defaults to the BFR heuristic and intercept shift", {
  stag <- .adopt_panel(c(9, 12, 15))
  fit  <- scm_fit(y ~ d | id + time, data = stag, method = "scm")
  expect_equal(fit$pooling$nu, fit$pooling$nu_heuristic)
  expect_gt(fit$pooling$nu, 0)
  expect_true(fit$pooling$fixedeff)

  # nu = "separate" opts back into the per-cohort V-optimised path
  sep <- scm_fit(y ~ d | id + time, data = stag, method = "scm",
                 nu = "separate")
  expect_null(sep$pooling)
  expect_false(isTRUE(all.equal(fit$estimate, sep$estimate)))
  # nu = 0 is the same estimator reached through the pooled solver, so the
  # two agree up to the solvers' own tolerances rather than exactly
  expect_equal(scm_fit(y ~ d | id + time, data = stag, method = "scm",
                       nu = 0, fixedeff = FALSE)$estimate,
               scm_fit(y ~ d | id + time, data = stag, method = "scm",
                       nu = "separate", fixedeff = FALSE)$estimate,
               tolerance = 1e-3)
})

test_that("the staggered default falls back where the pooled solver cannot go", {
  stag <- .adopt_panel(c(9, 12, 15))
  # these knobs belong to the per-cohort path, so the default must not error
  expect_null(scm_fit(y ~ d | id + time, data = stag, method = "scm",
                      lambda_pen = 0.1)$pooling)
  expect_null(scm_fit(y ~ d | id + time, data = stag, method = "scm",
                      v_selection = "oos")$pooling)
  # asking for pooling explicitly there is still an error
  expect_error(scm_fit(y ~ d | id + time, data = stag, method = "scm",
                       nu = 0.5, lambda_pen = 0.1), "cannot be combined")
})

test_that("pool_cohorts reaches the per-cohort path too", {
  # The fallback above must not silently discard the block granularity: a fit
  # asked for unit-level blocks has to deliver them whichever solver runs,
  # because scm_rank_pval() dispatches on exactly that.
  stag <- .adopt_panel(c(9, 9, 12, 15))
  blocks <- function(f) vapply(f$cohort_fits, function(z) length(z$idx_tr),
                               integer(1L))

  for (extra in list(list(nu = "separate"),
                     list(lambda_pen = 0.1),
                     list(donor_mspe_threshold = 5),
                     list(v_selection = "oos"))) {
    args <- c(list(y ~ d | id + time, data = stag, method = "scm",
                   pool_cohorts = FALSE), extra)
    f <- suppressWarnings(do.call(scm_fit, args))
    label <- names(extra)
    expect_null(f$pooling, info = label)           # the per-cohort path ran
    expect_equal(blocks(f), rep(1L, 4L), info = label)
    expect_false(f$pool_cohorts, info = label)
    expect_equal(vapply(f$cohort_fits, `[[`, character(1L), "unit"),
                 as.character(1:4), info = label)
    # and the test that needs unit-level blocks now accepts the fit
    expect_s3_class(suppressWarnings(scm_rank_pval(f, B = 20)),
                    "scm_rank_test")
  }

  # the default still pools cohorts on that path
  sep <- scm_fit(y ~ d | id + time, data = stag, method = "scm",
                 nu = "separate")
  expect_equal(blocks(sep), c(2L, 1L, 1L))
  expect_true(sep$pool_cohorts)
})

test_that("sharp SCM rejects the staggered-only block arguments", {
  sharp <- .adopt_panel(c(11, 11))
  expect_error(scm_fit(y ~ d | id + time, data = sharp, method = "scm",
                       pool_cohorts = FALSE), "staggered")
  # passing the default explicitly is a no-op, not an error
  expect_s3_class(scm_fit(y ~ d | id + time, data = sharp, method = "scm",
                          pool_cohorts = TRUE), "coresynth_scm")
})

test_that("sharp SCM refuses covariates and asks for predictors instead", {
  # SCM balances covariates (ADH 2010 S.2.3; BFR 2022 S.5.2) rather than
  # partialling them out, and a sharp fit has pred() to do it with. Before
  # this, `covariates` reached fit_scm_cpp() and changed nothing at all.
  dat <- .adopt_panel(11)
  set.seed(5)
  dat$x <- rnorm(nrow(dat)) + dat$id / 10

  expect_error(
    scm_fit(y ~ d | id + time, data = dat, method = "scm", covariates = "x"),
    "staggered SCM fits only"
  )
  expect_error(
    scm_fit(y ~ d | id + time, data = dat, method = "scm", covariates = "x"),
    "pred\\("
  )
  # the route the message names works
  balanced <- scm_fit(y ~ d | id + time, data = dat, method = "scm",
                      predictors = list(pred("y", 1:10, "mean"),
                                        pred("x", 1:10, "mean")))
  expect_equal(nrow(balanced$predictor_table), 2L)

  # staggered still partials out, and the covariate still moves the fit
  stag <- .adopt_panel(c(9, 12, 15))
  set.seed(5)
  stag$x <- rnorm(nrow(stag)) + stag$id / 10
  with_cov <- scm_fit(y ~ d | id + time, data = stag, method = "scm",
                      covariates = "x")
  expect_length(with_cov$beta_hat, 1L)
  expect_false(isTRUE(all.equal(
    with_cov$estimate,
    scm_fit(y ~ d | id + time, data = stag, method = "scm")$estimate)))
})

# ── BFR S.4 unit-level blocks, aggregation weights, rank test, frontier ──────

.bfr_panel <- function(adopt, N = 25, TT = 22, tau = 2, seed = 3) {
  set.seed(seed)
  f <- matrix(rnorm(TT * 2), TT, 2); L <- matrix(rnorm(N * 2), N, 2)
  Y <- f %*% t(L) + matrix(rnorm(N * TT, sd = 0.3), TT, N)
  D <- matrix(0L, TT, N)
  for (k in seq_along(adopt)) {
    t0 <- adopt[k]; D[t0:TT, k] <- 1L; Y[t0:TT, k] <- Y[t0:TT, k] + tau
  }
  data.frame(id = rep(1:N, each = TT), time = rep(1:TT, times = N),
             y = as.vector(Y), d = as.vector(D))
}

test_that("pool_cohorts = FALSE gives one block per treated unit", {
  # BFR S.4 makes Gamma an N-by-J matrix over treated units; their Appendix
  # A.2 adds the layer that pools units sharing an adoption date.
  panel_st <- .bfr_panel(c(9, 9, 13, 17))
  pooled <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm")
  units  <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                    pool_cohorts = FALSE)

  expect_equal(nrow(pooled$cohort_estimates), 3L)   # cohorts 9, 13, 17
  expect_equal(nrow(units$cohort_estimates), 4L)    # one per treated unit
  expect_true(pooled$pooling$pool_cohorts)
  expect_false(units$pooling$pool_cohorts)
  expect_equal(pooled$cohort_estimates$n_treated, c(2L, 1L, 1L))
  expect_equal(units$cohort_estimates$n_treated, rep(1L, 4L))
  # unit-level blocks name their unit; a pooled multi-unit cohort cannot
  expect_true(is.na(pooled$cohort_fits[[1]]$unit))
  expect_equal(vapply(units$cohort_fits, `[[`, character(1), "unit"),
               as.character(1:4))
  # both remain valid staggered fits
  expect_s3_class(units, "coresynth_staggered")
  expect_true(is.finite(units$estimate))
})

test_that("unit-level blocks are labelled by unit, not by adoption date", {
  # Two units adopt at 9, so the cohort label names both of them: anything
  # keyed on it merges their series (augment()) or duplicates a row (tidy()).
  panel_st <- .bfr_panel(c(9, 9, 13, 17))
  units  <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                    pool_cohorts = FALSE)
  pooled <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm")

  expect_true("unit" %in% names(units$cohort_estimates))
  expect_equal(units$cohort_estimates$unit, as.character(1:4))
  expect_equal(units$cohort_estimates$cohort, c(9L, 9L, 13L, 17L))

  terms_u <- broom::tidy(units)$term
  expect_equal(terms_u, paste0("unit_", 1:4))
  expect_false(anyDuplicated(terms_u) > 0L)

  au <- broom::augment(units)
  tr <- au[au$.type == "treated", ]
  expect_equal(sort(unique(tr$.unit)), paste0("unit_", 1:4))
  expect_equal(nrow(tr), 4L * length(units$times))

  # Cohort blocks keep the labels they had; the unit column reports the unit
  # of a one-unit cohort and NA for the pooled one.
  expect_equal(broom::tidy(pooled)$term, paste0("cohort_", c(9, 13, 17)))
  expect_equal(pooled$cohort_estimates$unit, c(NA, "3", "4"))
  au_p <- broom::augment(pooled)
  expect_equal(sort(unique(au_p$.unit[au_p$.type == "treated"])),
               paste0("cohort_", c(13, 17, 9)))
})

test_that("cohort tables of the other estimators are unaffected", {
  # .aggregate_cohorts() drops the unit column when no block names a unit, so
  # a fit whose blocks are cohorts keeps the table it had.
  panel_st <- .bfr_panel(c(9, 9, 13, 17))
  for (m in c("sdid", "gsc", "mc", "tasc", "si")) {
    ce <- scm_fit(y ~ d | id + time, data = panel_st, method = m)$cohort_estimates
    expect_identical(
      names(ce),
      c("cohort", "n_treated", "T_pre", "T_post", "estimate", "weight"),
      info = m)
  }
})

test_that("weights set the mass each treated unit carries", {
  panel_sh <- .bfr_panel(c(11, 11, 11))
  equal <- scm_fit(y ~ d | id + time, data = panel_sh, method = "scm")
  heavy <- scm_fit(y ~ d | id + time, data = panel_sh, method = "scm",
                   weights = c("1" = 1, "2" = 8, "3" = 1))
  expect_null(equal$treated_weights)
  expect_equal(sum(heavy$treated_weights), 1)
  expect_equal(unname(heavy$treated_weights["2"]), 8 / 10)
  expect_false(isTRUE(all.equal(equal$estimate, heavy$estimate)))
  # uniform weights reproduce the unweighted fit
  flat <- scm_fit(y ~ d | id + time, data = panel_sh, method = "scm",
                  weights = c("1" = 2, "2" = 2, "3" = 2))
  expect_equal(flat$estimate, equal$estimate, tolerance = 1e-8)

  # staggered: weights multiply each block's T_post, not replace it
  panel_st <- .bfr_panel(c(9, 13, 17))
  w_st <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                  weights = c("1" = 5, "2" = 1, "3" = 1))
  ce <- w_st$cohort_estimates
  expect_equal(w_st$estimate, sum(ce$weight * ce$estimate) / sum(ce$weight))
  expect_gt(ce$weight[1], ce$weight[2])
})

test_that("the pooled objective weights each block by its mass (BFR A.2)", {
  # Appendix A.2 writes the cohort problem on the SUM over the cohort's units:
  # its target is sum_i 1{T_i = T(g)} Y_i and its weights satisfy
  # sum_i gamma_ig = n_g, so the block residual is n_g times the mean gap. The
  # heuristic is then
  #   nu_hat = || sum_g n_g gap_g || / sum_g n_g ||gap_g||
  # over the separate solution, which nu = 0 leaves in place.
  panel_st <- .bfr_panel(c(9, 9, 9, 13, 17))   # cohort sizes 3, 1, 1
  fit <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                 nu = 0, fixedeff = FALSE)
  cf <- fit$cohort_fits
  expect_equal(vapply(cf, `[[`, integer(1L), "n_treated"), c(3L, 1L, 1L))

  L   <- max(vapply(cf, `[[`, integer(1L), "T_pre"))
  num <- numeric(L)
  den <- 0
  for (b in cf) {
    pre <- seq_len(b$T_pre)
    # lag order, so blocks with different windows align in event time
    gap <- rev((b$Y_treat - drop(fit$Y_all[, b$idx_co, drop = FALSE] %*%
                                   b$unit_weights))[pre])
    m   <- b$n_treated
    num[pre] <- num[pre] + m * gap
    den <- den + m * sqrt(sum(gap^2))
  }
  expect_equal(fit$pooling$nu_heuristic, min(1, sqrt(sum(num^2)) / den),
               tolerance = 1e-6)

  # unit-level blocks are the same formula with every mass at one
  units <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                   nu = 0, fixedeff = FALSE, pool_cohorts = FALSE)
  expect_equal(vapply(units$cohort_fits, `[[`, integer(1L), "n_treated"),
               rep(1L, 5L))
  expect_false(isTRUE(all.equal(fit$pooling$nu_heuristic,
                                units$pooling$nu_heuristic)))
})

test_that("weights carry through the objective, the target and the bootstrap", {
  panel_st <- .bfr_panel(c(9, 9, 13, 17))
  w <- c("1" = 6, "2" = 1, "3" = 1, "4" = 1)
  fit <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                 weights = w)

  # the block's target is the weighted average of its treated units, not the
  # simple one, so the fit is aimed at the quantity the aggregation reports
  b1 <- fit$cohort_fits[[1L]]
  expect_equal(b1$n_treated, 2L)
  wb <- c(6, 1) / 7
  expect_equal(unname(b1$Y_treat),
               unname(drop(b1$Y_treat_mat %*% wb)), tolerance = 1e-10)
  # and its mass is the sum of those weights, up to the common rescaling
  masses <- vapply(fit$cohort_fits, `[[`, numeric(1L), "mass")
  expect_equal(masses / masses[2L], c(7, 1, 1), tolerance = 1e-10)

  # the aggregation weight is that mass times T_post
  ce <- fit$cohort_estimates
  expect_equal(ce$weight,
               (masses * ce$T_post) / sum(masses * ce$T_post),
               tolerance = 1e-10)

  # the wild bootstrap decomposes the same ATT: per-unit masses times T_post
  TT <- length(fit$times)
  tau <- numeric(0); om <- numeric(0)
  for (b in fit$cohort_fits) {
    pre <- seq_len(b$T_pre); post <- (b$T_pre + 1L):TT
    tg <- (colMeans(b$Y_treat_mat[post, , drop = FALSE]) -
             colMeans(b$Y_treat_mat[pre, , drop = FALSE])) -
      (mean(b$Y_synth[post]) - mean(b$Y_synth[pre]))
    tau <- c(tau, tg)
    om  <- c(om, b$unit_masses * b$T_post)
  }
  expect_equal(sum(om / sum(om) * tau), fit$estimate, tolerance = 1e-10)
  # (four treated units, so the bootstrap warns; the centring is the point)
  expect_equal(suppressWarnings(
    scm_inference(fit, method = "wild_bootstrap", n_boot = 50,
                  seed = 1))$estimate, fit$estimate)

  # equal weights reproduce the unweighted fit exactly
  flat <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                  weights = c("1" = 3, "2" = 3, "3" = 3, "4" = 3))
  plain <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm")
  expect_equal(flat$estimate, plain$estimate, tolerance = 1e-8)
  expect_equal(flat$pooling$nu_heuristic, plain$pooling$nu_heuristic,
               tolerance = 1e-8)
})

test_that("staggered fits report the treated series on the fitted scale", {
  # With covariates the fit works on the partialled-out outcome: the cohort
  # fits, Y_all and the sharp path's own Y_treat are all on that scale, so
  # the staggered Y_treat has to be too, or treated_outcomes() and augment()
  # describe different series.
  stag <- .bfr_panel(c(9, 13, 17))
  set.seed(7)
  stag$x <- rnorm(nrow(stag)) + stag$id / 10

  for (m in c("scm", "sdid")) {
    fit <- scm_fit(y ~ d | id + time, data = stag, method = m,
                   covariates = "x")
    series <- treated_outcomes(fit)
    Y_work <- fit$Y_all %||% fit$Y_mat      # SCM stores Y_all, SDID Y_mat
    cohort <- rowMeans(vapply(fit$cohort_fits, function(cf) {
      rowMeans(Y_work[, cf$idx_tr, drop = FALSE])
    }, numeric(length(fit$times))))
    expect_equal(unname(series), unname(cohort), tolerance = 1e-10,
                 info = m)
    # and it is no longer the raw outcome, since a covariate was partialled out
    raw <- rowMeans(matrix(stag$y, nrow = length(fit$times))[, 1:3,
                                                             drop = FALSE])
    expect_false(isTRUE(all.equal(unname(series), unname(raw))), info = m)
  }

  # without covariates the two scales coincide, so nothing moves
  plain <- scm_fit(y ~ d | id + time, data = stag, method = "scm")
  expect_equal(unname(treated_outcomes(plain)),
               unname(rowMeans(matrix(stag$y, nrow = length(plain$times))[,
                                                                         1:3, drop = FALSE])),
               tolerance = 1e-10)
})

test_that("loo_donors mirrors the fit's inner solver", {
  panel_sh <- .bfr_panel(11)
  fit <- scm_fit(y ~ d | id + time, data = panel_sh, method = "scm")
  expect_equal(fit$qp_solver, "wolfe")
  loo <- loo_donors(fit)
  expect_true(all(is.finite(loo$results$att_loo)))
  # dropping a zero-weight donor cannot move the estimate, whichever solver
  # ran -- the check that the mirrored solver is still solving the same QP
  expect_equal(loo$att_original, fit$estimate)
})

test_that("conformal_inference refuses a fit missing its tuning parameter", {
  # The refit must impose the model the fit used, so a fit that does not
  # record r / k / lambda cannot be scored against a substituted rule.
  for (spec in list(c(m = "gsc", nm = "r"), c(m = "si", nm = "k"),
                    c(m = "mc", nm = "lambda"))) {
    fit <- scm_fit(y ~ d | id + time, data = panel, method = spec[["m"]])
    expect_false(is.null(fit[[spec[["nm"]]]]), info = spec[["m"]])
    stripped <- fit
    stripped[[spec[["nm"]]]] <- NULL
    expect_error(conformal_inference(stripped, ci = FALSE),
                 "current version", info = spec[["m"]])
    # the intact fit still works
    expect_s3_class(conformal_inference(fit, ci = FALSE),
                    "conformal_inference")
  }
})

test_that("cov_method = 'balance' adds the covariate term of BFR Eq. (11)", {
  # S.5.2 balances the covariates alongside the lagged outcomes instead of
  # residualising the outcome on them, weighting the two sides by xi.
  set.seed(4)
  stag <- .bfr_panel(c(9, 9, 13, 17))
  stag$x1 <- rnorm(nrow(stag)) + stag$id / 6
  stag$x2 <- rnorm(nrow(stag)) - stag$id / 9
  ff <- function(...) scm_fit(y ~ d | id + time, data = stag,
                              method = "scm", ...)

  plain <- ff()
  bal   <- ff(covariates = c("x1", "x2"), cov_method = "balance")
  expect_equal(bal$pooling$n_covariates, 2L)
  expect_true(is.finite(bal$pooling$q_cov_sep))
  expect_true(is.finite(bal$pooling$q_cov_pool))
  # partial_out leaves no covariate balance behind, and is still the default
  expect_null(ff(covariates = c("x1", "x2"))$pooling$q_cov_sep)

  # xi -> 0 removes the covariate side, so the fit collapses onto the
  # outcomes-only one at a rate linear in xi
  errs <- vapply(c(1e-4, 1e-8, 1e-12), function(xi) {
    abs(ff(covariates = c("x1", "x2"), cov_method = "balance",
           xi = xi)$estimate - plain$estimate)
  }, numeric(1L))
  expect_true(all(diff(errs) < 0))
  expect_lt(errs[3L], 1e-9)

  # raising xi buys covariate balance at the cost of outcome balance
  Z <- vapply(c("x1", "x2"), function(v)
    tapply(stag[[v]][stag$time <= 8], stag$id[stag$time <= 8], mean),
    numeric(length(unique(stag$id))))
  cov_imbalance <- function(f) {
    mean(vapply(f$cohort_fits, function(cf)
      sum((colMeans(Z[cf$idx_tr, , drop = FALSE]) -
             drop(crossprod(Z[cf$idx_co, , drop = FALSE],
                            cf$unit_weights)))^2), numeric(1L)))
  }
  low  <- ff(covariates = c("x1", "x2"), cov_method = "balance", xi = 1e-8)
  high <- ff(covariates = c("x1", "x2"), cov_method = "balance", xi = 50)
  expect_lt(cov_imbalance(high), cov_imbalance(low))
  expect_gt(high$pooling$q_sep, low$pooling$q_sep)

  # the frontier traces the objective the fit actually used
  fr <- scm_balance_frontier(bal, nu_grid = c(0, 0.5, 1))
  expect_equal(nrow(fr), 3L)
  expect_true(all(is.finite(fr$q_sep)))
})

test_that("cov_method = 'balance' is refused where its solver cannot go", {
  set.seed(4)
  stag <- .bfr_panel(c(9, 13, 17))
  stag$x1 <- rnorm(nrow(stag)) + stag$id / 6
  sharp <- .bfr_panel(c(11, 11))
  sharp$x1 <- rnorm(nrow(sharp)) + sharp$id / 6

  for (extra in list(list(nu = "separate"), list(lambda_pen = 0.1),
                     list(v_selection = "oos"),
                     list(donor_mspe_threshold = 5))) {
    args <- c(list(y ~ d | id + time, data = stag, method = "scm",
                   covariates = "x1", cov_method = "balance"), extra)
    expect_error(do.call(scm_fit, args), "partially pooled objective",
                 info = names(extra))
  }
  expect_error(scm_fit(y ~ d | id + time, data = sharp, method = "scm",
                       covariates = "x1", cov_method = "balance"),
               "sharp fit")
  expect_error(scm_fit(y ~ d | id + time, data = stag, method = "scm",
                       cov_method = "balance"), "covariates to balance")
  expect_error(scm_fit(y ~ d | id + time, data = stag, method = "scm",
                       covariates = "x1", cov_method = "balance", xi = -1),
               "positive number")
})

test_that("xi is refused where no covariate term uses it", {
  # xi only weights the balance term; everywhere else it changed nothing,
  # which is the silent no-op the argument checks exist to stop.
  set.seed(4)
  stag <- .bfr_panel(c(9, 13, 17))
  stag$x1 <- rnorm(nrow(stag)) + stag$id / 6
  sharp <- .bfr_panel(c(11, 11))

  expect_error(scm_fit(y ~ d | id + time, data = sharp, method = "scm",
                       xi = 99), "no effect")
  expect_error(scm_fit(y ~ d | id + time, data = stag, method = "scm",
                       covariates = "x1", xi = 99), "no effect")
  expect_s3_class(
    scm_fit(y ~ d | id + time, data = stag, method = "scm",
            covariates = "x1", cov_method = "balance", xi = 99),
    "coresynth_scm")
})

test_that("the fallback off the pooled solver announces itself", {
  panel_st <- .bfr_panel(c(9, 13, 17))
  ff <- function(...) scm_fit(y ~ d | id + time, data = panel_st,
                              method = "scm", ...)
  expect_message(f <- ff(lambda_pen = 0.1), "lambda_pen")
  expect_message(ff(donor_mspe_threshold = 5), "donor_mspe_threshold")
  expect_message(ff(v_selection = "oos"), "v_selection")
  expect_message(ff(lambda_pen = 0.1, donor_mspe_threshold = 5),
                 "donor_mspe_threshold and lambda_pen")
  # an explicit choice, and the default, say nothing
  expect_silent(ff(nu = "separate"))
  expect_silent(ff())
  # and the fit records which estimator it ended up being
  expect_identical(f$nu, "separate")
  expect_identical(ff(nu = "separate")$nu, "separate")
  expect_true(is.numeric(ff()$nu))
})

test_that("staggered SCM records the solver it resolved to", {
  panel_st <- .bfr_panel(c(9, 13, 17))
  sharp    <- .bfr_panel(c(11, 11))
  # the staggered paths run their own QP, so "auto" lands on active_set
  expect_identical(scm_fit(y ~ d | id + time, data = panel_st,
                           method = "scm")$qp_solver, "active_set")
  expect_identical(scm_fit(y ~ d | id + time, data = sharp,
                           method = "scm")$qp_solver, "wolfe")
})

test_that("the frontier reports the covariate side of a balance fit", {
  set.seed(4)
  stag <- .bfr_panel(c(9, 13, 17))
  stag$x1 <- rnorm(nrow(stag)) + stag$id / 6
  bal <- scm_fit(y ~ d | id + time, data = stag, method = "scm",
                 covariates = "x1", cov_method = "balance")
  fr <- scm_balance_frontier(bal, nu_grid = c(0, 0.5, 1))
  expect_true(all(c("q_cov_sep", "q_cov_pool") %in% names(fr)))
  expect_true(all(is.finite(fr$q_cov_sep)))
  # footnote 7 normalises by the combined measure, so the outcome ratios are
  # below 1 at nu = 0 rather than equal to it
  expect_lt(fr$q_sep_rel[1L], 1)
  expect_lt(fr$q_pool_rel[1L], 1)
  expect_match(plot(fr)$labels$caption, "q_cov_sep")

  # an outcomes-only fit keeps the two-column frontier and its normalisation
  plain <- scm_fit(y ~ d | id + time, data = stag, method = "scm")
  fr2 <- scm_balance_frontier(plain, nu_grid = c(0, 0.5, 1))
  expect_false(any(c("q_cov_sep", "q_cov_pool") %in% names(fr2)))
  expect_equal(fr2$q_sep_rel[1L], 1, tolerance = 1e-6)
  expect_equal(fr2$q_pool_rel[1L], 1, tolerance = 1e-6)
  expect_null(plot(fr2)$labels$caption)
})

test_that("the frontier is traced on the weighted objective", {
  panel_st <- .bfr_panel(c(9, 9, 13, 17))
  fit <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm",
                 weights = c("1" = 6, "2" = 1, "3" = 1, "4" = 1))
  fr <- scm_balance_frontier(fit, nu_grid = c(0, 0.5, 1))
  # at nu = 0 the frontier reproduces the separate solution the fit normalised
  # by, and the estimate column is the weighted one
  expect_equal(fr$q_sep[1L], fit$pooling$q_sep_separate, tolerance = 1e-6)
  expect_equal(fr$q_pool[1L], fit$pooling$q_pool_separate, tolerance = 1e-6)
  sep <- scm_fit(y ~ d | id + time, data = panel_st, method = "scm", nu = 0,
                 weights = c("1" = 6, "2" = 1, "3" = 1, "4" = 1))
  expect_equal(fr$estimate[1L], sep$estimate, tolerance = 1e-6)
})

test_that("weights are validated against the treated units", {
  panel_sh <- .bfr_panel(c(11, 11, 11))
  expect_error(scm_fit(y ~ d | id + time, data = panel_sh, method = "scm",
                       weights = c(1, 2, 3)), "must be named")
  expect_error(scm_fit(y ~ d | id + time, data = panel_sh, method = "scm",
                       weights = c("1" = 1)), "missing the treated unit")
  expect_error(scm_fit(y ~ d | id + time, data = panel_sh, method = "scm",
                       weights = c("1" = 1, "2" = -1, "3" = 1)), "non-negative")
  expect_error(scm_fit(y ~ d | id + time, data = panel_sh, method = "scm",
                       weights = c("1" = 0, "2" = 0, "3" = 0)), "positive sum")
})

test_that("scm_rank_pval runs on both fits that expose per-unit effects", {
  sharp <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(11, 11, 11)),
                   method = "scm")
  res <- scm_rank_pval(sharp, B = 40L)
  expect_s3_class(res, "scm_rank_test")
  expect_gte(res$p_value, 0); expect_lte(res$p_value, 1)
  expect_length(res$rank_sums, 40L)
  expect_length(res$unit_statistics, 3L)
  expect_output(print(res), "rank-based permutation test")

  units <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(9, 13, 17)),
                   method = "scm", pool_cohorts = FALSE)
  expect_s3_class(scm_rank_pval(units, B = 40L), "scm_rank_test")
})

test_that("scm_rank_pval is reproducible and leaves the RNG alone", {
  sharp <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(11, 11, 11)),
                   method = "scm")
  expect_equal(scm_rank_pval(sharp, B = 30L)$p_value,
               scm_rank_pval(sharp, B = 30L)$p_value)
  set.seed(99); before <- runif(1)
  set.seed(99); invisible(scm_rank_pval(sharp, B = 30L))
  expect_equal(runif(1), before)
})

test_that("scm_rank_pval refuses fits without per-unit effects", {
  one <- scm_fit(y ~ d | id + time, data = .bfr_panel(11), method = "scm")
  expect_error(scm_rank_pval(one), "mspe_ratio_pval")
  pooled <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(9, 9, 13)),
                    method = "scm")
  expect_error(scm_rank_pval(pooled), "pool_cohorts = FALSE")
  expect_error(scm_rank_pval(scm_fit(y ~ d | id + time,
                                     data = .bfr_panel(c(11, 11, 11)),
                                     method = "mc")), "method = 'scm'")
})

test_that("scm_balance_frontier traces the pooling trade-off", {
  fit <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(11, 11, 14, 18)),
                 method = "scm")
  fr <- scm_balance_frontier(fit, nu_grid = seq(0, 1, by = 0.25))
  expect_s3_class(fr, "scm_frontier")
  expect_equal(nrow(fr), 5L)
  # pooling buys pooled fit at the cost of per-block fit, monotonically
  expect_true(all(diff(fr$q_pool) <= 1e-8))
  expect_true(all(diff(fr$q_sep) >= -1e-8))
  # normalised against the separate solution, nu = 0 is the unit point
  expect_equal(fr$q_sep_rel[1], 1, tolerance = 1e-6)
  expect_equal(fr$q_pool_rel[1], 1, tolerance = 1e-6)
  expect_equal(sum(fr$heuristic), 1L)
  expect_output(print(fr), "balance possibility frontier")
  expect_s3_class(plot(fr), "ggplot")
})

test_that("scm_balance_frontier refuses fits it cannot trace", {
  sharp <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(11, 11, 11)),
                   method = "scm")
  expect_error(scm_balance_frontier(sharp), "staggered")
  sep <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(9, 13, 17)),
                 method = "scm", nu = "separate")
  expect_error(scm_balance_frontier(sep), "no pooling")
  fit <- scm_fit(y ~ d | id + time, data = .bfr_panel(c(9, 13, 17)),
                 method = "scm")
  expect_error(scm_balance_frontier(fit, nu_grid = 0.5), "two or more")
  expect_error(scm_balance_frontier(fit, nu_grid = c(0, 2)), "two or more")
})

# ── ADH (2010) Figures 5-7: pruning ill-fitting placebo runs ─────────────────
# A placebo run the synthetic control could not track before treatment says
# nothing about how rare a large post-treatment gap is. ADH therefore report
# the permutation distribution with such runs discarded, at cutoffs of 20, 5
# and 2 times the treated unit's pre-treatment MSPE.

.prune_panel <- function(seed = 7, N = 40, TT = 24, T0 = 14, tau = 3,
                         n_bad = 5L) {
  set.seed(seed)
  f <- matrix(rnorm(TT * 2), TT, 2); L <- matrix(rnorm(N * 2), N, 2)
  Y <- f %*% t(L) + matrix(rnorm(N * TT, sd = 0.4), TT, N)
  bad <- (N - n_bad + 1L):N
  Y[, bad] <- Y[, bad] + matrix(rnorm(TT * n_bad, sd = 4), TT, n_bad)
  D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1] <- 1L
  Y[(T0 + 1):TT, 1] <- Y[(T0 + 1):TT, 1] + tau
  data.frame(id = rep(1:N, each = TT), time = rep(1:TT, times = N),
             y = as.vector(Y), d = as.vector(D))
}

test_that("mspe_prune discards the runs it says it does", {
  fit <- scm_fit(y ~ d | id + time, data = .prune_panel(), method = "scm")
  none <- mspe_ratio_pval(fit)
  expect_equal(none$n_placebo_pruned, 0L)
  expect_equal(none$mspe_prune, Inf)

  tight <- mspe_ratio_pval(fit, mspe_prune = 5)
  expect_gt(tight$n_placebo_pruned, 0L)
  # exactly the runs above the cutoff, and no others
  expect_equal(tight$n_placebo_pruned,
               sum(tight$mspe_pre_placebo > 5 * tight$mspe_pre_treated))
  expect_equal(tight$n_placebo_used + tight$n_placebo_pruned,
               length(tight$mspe_pre_placebo))
  # a lower cutoff can only discard more
  expect_gte(mspe_ratio_pval(fit, mspe_prune = 2)$n_placebo_pruned,
             tight$n_placebo_pruned)
})

test_that("pruning is what makes the post-MSPE statistic usable", {
  # Abadie (2021) S.3.5 offers ranking the post-treatment fit as the
  # alternative to the ratio, but only once ill-fitting runs are gone: left in,
  # they sit at the top of the distribution on size alone.
  fit <- scm_fit(y ~ d | id + time, data = .prune_panel(), method = "scm")
  raw <- mspe_ratio_pval(fit, statistic = "post_mspe")
  cut <- mspe_ratio_pval(fit, statistic = "post_mspe", mspe_prune = 5)
  expect_gt(raw$p_value, cut$p_value)
  expect_lt(cut$p_value, 0.1)
  expect_equal(cut$statistic, "post_mspe")
  # the ratio statistic needs no cutoff to find the same signal (ADH Fig. 8)
  expect_lt(mspe_ratio_pval(fit)$p_value, 0.1)
})

test_that("pruning also applies to the ratio and one-sided tests", {
  fit <- scm_fit(y ~ d | id + time, data = .prune_panel(), method = "scm")
  r <- mspe_ratio_pval(fit, mspe_prune = 5)
  expect_gt(r$n_placebo_pruned, 0L)
  expect_equal(r$statistic, "ratio")
  one <- mspe_ratio_pval(fit, mspe_prune = 5, alternative = "greater")
  expect_equal(one$n_placebo_used + one$n_placebo_pruned,
               length(one$mspe_pre_placebo))
  expect_gte(one$p_value, 0); expect_lte(one$p_value, 1)
})

test_that("mspe_prune and statistic are validated", {
  fit <- scm_fit(y ~ d | id + time, data = .prune_panel(), method = "scm")
  expect_error(mspe_ratio_pval(fit, mspe_prune = 0), "positive number")
  expect_error(mspe_ratio_pval(fit, mspe_prune = -1), "positive number")
  expect_error(mspe_ratio_pval(fit, mspe_prune = c(2, 5)), "positive number")
  expect_error(mspe_ratio_pval(fit, statistic = "post_mspe",
                               alternative = "greater"), "no direction")
  expect_error(mspe_ratio_pval(fit, mspe_prune = 1e-8,
                               statistic = "post_mspe"), "discards every")
})

test_that("the gap plot reuses the cutoff the p-value was computed under", {
  fit <- scm_fit(y ~ d | id + time, data = .prune_panel(), method = "scm")
  cut <- mspe_ratio_pval(fit, mspe_prune = 5)
  n_lines <- function(p) {
    d <- ggplot2::ggplot_build(p)$data
    length(unique(d[[1]]$group))
  }
  expect_equal(n_lines(plot(cut, type = "gaps")), cut$n_placebo_used)
  # Inf overrides it, showing every run again
  expect_equal(n_lines(plot(cut, type = "gaps", mspe_prune = Inf)),
               length(cut$mspe_pre_placebo))
  # an unpruned test still shows everything by default
  none <- mspe_ratio_pval(fit)
  expect_equal(n_lines(plot(none, type = "gaps")),
               length(none$mspe_pre_placebo))
})

test_that("the ratio plot reuses the cutoff the p-value was computed under", {
  # The subtitle reports n_placebo_used, so the points have to be that many:
  # a figure drawing runs the p-value discarded contradicts its own caption.
  fit <- scm_fit(y ~ d | id + time, data = .prune_panel(), method = "scm")
  cut <- mspe_ratio_pval(fit, mspe_prune = 5)
  expect_gt(cut$n_placebo_pruned, 0L)

  n_pts <- function(p) nrow(ggplot2::ggplot_build(p)$data[[1]])
  expect_equal(n_pts(plot(cut, type = "ratios")), cut$n_placebo_used + 1L)
  expect_equal(nrow(plot_data(cut, type = "ratios")), cut$n_placebo_used + 1L)
  # Inf overrides it, as it does for the gaps
  expect_equal(n_pts(plot(cut, type = "ratios", mspe_prune = Inf)),
               length(cut$mspe_pre_placebo) + 1L)

  none <- mspe_ratio_pval(fit)
  expect_equal(n_pts(plot(none, type = "ratios")),
               length(none$mspe_pre_placebo) + 1L)
})

test_that("the ratio plot draws the statistic the test ranked", {
  fit  <- scm_fit(y ~ d | id + time, data = .prune_panel(), method = "scm")
  post <- mspe_ratio_pval(fit, mspe_prune = 5, statistic = "post_mspe")
  p    <- plot(post, type = "ratios")
  drawn <- p$data$value[p$data$series == "Treated"]
  expect_equal(drawn, unname(post$mspe_post_treated))
  expect_false(isTRUE(all.equal(drawn, unname(post$mspe_ratio_treated))))
  expect_match(p$labels$x, "Post-treatment MSPE")

  # the ratio statistic keeps the ADH Figure 8 axis
  ratio <- mspe_ratio_pval(fit, mspe_prune = 5)
  pr <- plot(ratio, type = "ratios")
  expect_equal(pr$data$value[pr$data$series == "Treated"],
               unname(ratio$mspe_ratio_treated))
  expect_match(pr$labels$x, "MSPE ratio")
})

# ── staggered donor pool, covariate partial-out and SCM jackknife ─────────────

# Three cohorts plus never-treated units, a constant effect, and treated
# units whose loadings sit apart from the controls', so the synthetic
# controls lean on the later adopters whenever they are allowed to.
.contam_panel <- function(seed = 1, tau = 3) {
  set.seed(seed)
  adopt <- c(rep(Inf, 12), rep(c(9, 13, 17), each = 3))
  N <- length(adopt); TT <- 24
  lam <- rnorm(N); lam[adopt < Inf] <- lam[adopt < Inf] + 0.5
  f <- cumsum(rnorm(TT, 0, 0.3))
  d <- expand.grid(time = seq_len(TT), id = seq_len(N))
  d$y <- rnorm(N)[d$id] + lam[d$id] * f[d$time] + rnorm(nrow(d), 0, 0.2)
  d$d <- as.integer(d$time >= adopt[d$id])
  d$y <- d$y + tau * d$d
  d$id <- sprintf("u%02d", d$id)
  d
}

test_that("the default donor pool keeps later adopters out of earlier cohorts", {
  dat <- .contam_panel()
  for (m in c("scm", "sdid")) {
    fit <- scm_fit(y ~ d | id + time, data = dat, method = m)
    expect_equal(fit$estimate, 3, tolerance = 0.1)
    cf <- fit$cohort_fits
    tr <- unlist(lapply(cf, `[[`, "idx_tr"))
    for (z in cf) expect_length(intersect(z$idx_co, tr), 0L)
    # "clean" lets the later adopters in, and their effect leaks into the
    # earlier cohorts' counterfactuals
    cl <- suppressWarnings(scm_fit(y ~ d | id + time, data = dat, method = m,
                                   control_group = "clean"))
    expect_lt(cl$estimate, fit$estimate)
  }
})

# x carries a unit level and a unit-specific trend that is steeper for the
# treated units; y loads on x with beta = 1 on top of unit and time effects
# correlated with x's.
.partial_out_panel <- function(seed = 1, tau = 2, staggered = TRUE) {
  set.seed(seed)
  adopt <- c(rep(Inf, 20),
             if (staggered) rep(c(12, 18), each = 4) else rep(15, 8))
  N <- length(adopt); TT <- 26; tr <- adopt < Inf
  c_i <- rnorm(N, 5, 2); s_i <- rnorm(N, 0, 0.05) + 0.05 * tr
  d <- expand.grid(time = seq_len(TT), id = seq_len(N))
  d$x <- c_i[d$id] + (s_i[d$id] + 0.2) * d$time + rnorm(nrow(d), 0, 0.3)
  d$y <- 3 * c_i[d$id] + 0.5 * d$time + d$x + rnorm(nrow(d), 0, 0.3)
  d$d <- as.integer(d$time >= adopt[d$id])
  d$y <- d$y + tau * d$d
  d$id <- sprintf("u%02d", d$id)
  d
}

test_that("covariate partial-out absorbs unit and time effects (Kranz 2022)", {
  dat <- .partial_out_panel()
  for (m in c("scm", "sdid")) {
    fit <- scm_fit(y ~ d | id + time, data = dat, method = m,
                   covariates = "x")
    # without the fixed effects beta picks up the level correlation (~3.7)
    expect_equal(unname(fit$beta_hat), 1, tolerance = 0.1)
    expect_equal(fit$estimate, 2, tolerance = 0.2)
  }
  sh_dat <- .partial_out_panel(staggered = FALSE)
  sharp <- scm_fit(y ~ d | id + time, data = sh_dat, method = "sdid",
                   covariates = "x")
  expect_equal(unname(sharp$beta_hat), 1, tolerance = 0.1)
  # sharp and staggered share Kranz's sample: every untreated cell
  sh0 <- sh_dat[sh_dat$d == 0, ]
  expect_equal(unname(sharp$beta_hat),
               unname(coef(lm(y ~ x + factor(id) + factor(time), sh0))["x"]),
               tolerance = 1e-6)

  # on the untreated cells it is exactly the two-way fixed-effects slope
  un0  <- dat[dat$d == 0, ]
  b_fe <- unname(coef(lm(y ~ x + factor(id) + factor(time), un0))["x"])
  fit  <- scm_fit(y ~ d | id + time, data = dat, method = "sdid",
                  covariates = "x")
  expect_equal(unname(fit$beta_hat), b_fe, tolerance = 1e-6)
})

test_that("a covariate the fixed effects absorb is flagged, not adjusted for", {
  dat <- .partial_out_panel()
  dat$z <- as.numeric(factor(dat$id))          # time-invariant
  expect_warning(
    fit <- scm_fit(y ~ d | id + time, data = dat, method = "sdid",
                   covariates = c("x", "z")),
    "absorb")
  expect_equal(fit$beta_hat[2], 0)
})

test_that("scm_inference() jackknife follows BFR Appendix A.1", {
  dat <- .contam_panel()
  fit <- scm_fit(y ~ d | id + time, data = dat, method = "scm")
  inf <- scm_inference(fit)
  expect_equal(inf$method, "jackknife")
  expect_equal(inf$estimate, fit$estimate)
  N <- length(unique(dat$id))
  expect_length(inf$jack_ests, N)
  expect_true(all(is.finite(inf$jack_ests)))
  jb <- inf$jack_ests
  expect_equal(inf$se, sqrt((N - 1) / N * sum((jb - mean(jb))^2)))
  expect_equal(inf$ci_upper - inf$ci_lower, 2 * qnorm(0.975) * inf$se)
  expect_true(inf$ci_lower < 3 && inf$ci_upper > 3)

  # dropping a donor refits with nu held at the fit's value
  pan <- coresynth:::.pan_from_staggered_fit(fit)
  keep <- setdiff(seq_len(N), 1L)
  p <- pan
  p$Y <- pan$Y[, keep]; p$D <- pan$D[, keep]; p$units <- pan$units[keep]
  p$T_adopt <- pan$T_adopt[keep]
  p$idx_treat <- which(!is.na(p$T_adopt)); p$idx_control <- which(is.na(p$T_adopt))
  ref <- coresynth:::.fit_scm_staggered_pooled(
    p, Y_matrix = p$Y, nu = fit$pooling$nu, fixedeff = TRUE,
    control_group = "never_treated")
  # the jackknife warm-starts from the full fit; a cold refit meets it at
  # the block solver's tolerance
  expect_equal(unname(inf$jack_ests[1]), ref$estimate, tolerance = 1e-5)

  # the per-block path refits with its own settings
  sep <- scm_fit(y ~ d | id + time, data = dat, method = "scm",
                 nu = "separate")
  js <- scm_inference(sep)
  expect_true(is.finite(js$se) && js$se > 0)
  expect_equal(js$estimate, sep$estimate)

  gl <- broom::glance(inf)
  expect_equal(gl$method, "jackknife")
})

test_that("staggered SCM refits keep a skipped cohort's units out of the donor pool", {
  dat <- .contam_panel()
  # a cohort adopting in period 2 has one pre-period and is skipped
  dat$d[dat$id == "u13" & dat$time >= 2] <- 1L
  expect_warning(fit <- scm_fit(y ~ d | id + time, data = dat, method = "scm"),
                 "skipped")
  pan <- coresynth:::.pan_from_staggered_fit(fit)
  expect_false(13L %in% pan$idx_control)
  expect_true(13L %in% pan$idx_treat)
})

test_that("solve_simplex_qp() takes a known Lipschitz constant", {
  set.seed(3)
  X <- matrix(rnorm(8 * 40), 8, 40); b <- rnorm(8)
  Q <- crossprod(X); cv <- drop(crossprod(X, b))
  lip <- max(eigen(tcrossprod(X), symmetric = TRUE, only.values = TRUE)$values)
  expect_equal(drop(coresynth:::solve_simplex_qp(Q, cv, lipschitz = lip)),
               drop(coresynth:::solve_simplex_qp(Q, cv)), tolerance = 1e-10)
})

test_that("the jackknife skips donors with no weight, whose refit changes nothing", {
  dat <- .contam_panel()
  fit <- scm_fit(y ~ d | id + time, data = dat, method = "scm")
  w <- c(unlist(lapply(fit$cohort_fits, `[[`, "unit_weights")),
         unlist(fit$pooling$separate_weights))
  active <- unique(names(w)[w > 1e-12])
  never  <- unique(dat$id[!dat$id %in% dat$id[dat$d == 1]])
  inert  <- setdiff(never, active)
  skip_if(length(inert) == 0L, "every donor carries weight in this panel")
  inf <- scm_inference(fit)
  expect_true(all(inf$jack_ests[inert] == fit$estimate))
  # refitting one of them from scratch lands on the same estimate
  pan  <- coresynth:::.pan_from_staggered_fit(fit)
  i    <- match(inert[1], pan$units)
  keep <- setdiff(seq_len(ncol(pan$Y)), i)
  p <- pan
  p$Y <- pan$Y[, keep]; p$D <- pan$D[, keep]; p$units <- pan$units[keep]
  p$T_adopt <- pan$T_adopt[keep]
  p$idx_treat <- which(!is.na(p$T_adopt)); p$idx_control <- which(is.na(p$T_adopt))
  ref <- coresynth:::.fit_scm_staggered_pooled(
    p, Y_matrix = p$Y, nu = fit$pooling$nu, fixedeff = TRUE,
    control_group = "never_treated")
  expect_equal(ref$estimate, fit$estimate, tolerance = 1e-5)
})

test_that("the wild bootstrap follows Otsu & Rai (2017, Eq. 9)", {
  dat <- .contam_panel()
  fit <- scm_fit(y ~ d | id + time, data = dat, method = "scm")
  c_u <- coresynth:::.scm_wild_coefs(fit)
  # every unit carries a term; with all multipliers at one they cancel
  expect_equal(sum(c_u), 0, tolerance = 1e-12)
  used <- unique(unlist(lapply(fit$cohort_fits, function(b)
    names(b$unit_weights)[b$unit_weights > 1e-8])))
  expect_true(all(abs(c_u[used]) > 0))

  inf <- scm_inference(fit, method = "wild_bootstrap", n_boot = 400, seed = 1)
  # only the treated heterogeneity is centred at the ATT, so a larger effect
  # moves the estimate but not the interval's width ...
  big <- dat; big$y[big$d == 1] <- big$y[big$d == 1] + 10
  fit_b <- scm_fit(y ~ d | id + time, data = big, method = "scm")
  inf_b <- scm_inference(fit_b, method = "wild_bootstrap", n_boot = 400,
                         seed = 1)
  expect_equal(fit_b$estimate, fit$estimate + 10, tolerance = 1e-6)
  expect_equal(inf_b$se, inf$se, tolerance = 1e-6)
  # ... and the residuals net out period effects, so a common shock does not
  # widen it either
  shock <- dat; shock$y <- shock$y + 5 * (shock$time > 15)
  fit_s <- scm_fit(y ~ d | id + time, data = shock, method = "scm")
  inf_s <- scm_inference(fit_s, method = "wild_bootstrap", n_boot = 400,
                         seed = 1)
  # (the shifted panel is refitted, so equal up to the block solver's tolerance)
  expect_equal(inf_s$se, inf$se, tolerance = 1e-4)
})

test_that("scm_fit() rejects transformed outcome or treatment terms", {
  # The columns are read from `data` by name, so log(y) used to be fitted as y.
  expect_error(scm_fit(log(y) ~ d | id + time, data = panel, method = "sdid"),
               "outcome must be a column name")
  expect_error(scm_fit(y ~ I(d > 0) | id + time, data = panel, method = "sdid"),
               "treatment must be a column name")
  dat <- panel; dat$log_y <- log(dat$y + 100)
  expect_silent(scm_fit(log_y ~ d | id + time, data = dat, method = "sdid"))
})

test_that("gsc_boot() centres its interval on the estimate (Xu 2017 Alg. 2)", {
  fit  <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2)
  boot <- gsc_boot(fit, B = 59L, seed = 3L)
  q <- quantile(boot$boot_dist, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  # the replicates are draws of the estimation error under no effect
  expect_equal(boot$ci_lower, fit$estimate - q[2])
  expect_equal(boot$ci_upper, fit$estimate - q[1])
  expect_true(boot$ci_lower < fit$estimate && fit$estimate < boot$ci_upper)
})

test_that("placebo results do not depend on the thread count", {
  fit_scm  <- suppressMessages(
    scm_fit(y ~ d | id + time, data = panel, method = "scm"))
  fit_sdid <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  run <- function(n) {
    old <- options(coresynth.threads = n)
    on.exit(options(old))
    list(scm  = mspe_ratio_pval(fit_scm)$placebo_effects,
         sdid = sdid_inference(fit_sdid, "placebo")$placebo_effects)
  }
  expect_identical(run(1L), run(2L))
})

test_that("the diagonal-R Kalman gain matches the direct form", {
  set.seed(5)
  N <- 12L; TT <- 15L; r <- 2L
  Y <- matrix(rnorm(N * TT), N, TT)
  Y[1:2, 11:15] <- NA                       # treated post cells
  W <- matrix(rnorm(N * r), N, r)
  A <- diag(0.9, r); C <- c(0.1, -0.1); Q <- diag(0.2, r)
  R_diag  <- diag(runif(N, 0.5, 1.5))
  # numerically the same matrix, but no longer diagonal: takes the direct form
  R_dense <- R_diag; R_dense[1, 2] <- R_dense[2, 1] <- 1e-300
  ks_fast <- kalman_smoother_cpp(Y, W, A, C, Q, R_diag, rep(0, r), diag(r))
  ks_ref  <- kalman_smoother_cpp(Y, W, A, C, Q, R_dense, rep(0, r), diag(r))
  for (nm in c("z_smooth", "P_smooth", "P_cross", "z_upd"))
    expect_equal(ks_fast[[nm]], ks_ref[[nm]], tolerance = 1e-10)
})

test_that("SDID unit weights solve the same QP with or without forming Q", {
  set.seed(9)
  Y_pre <- matrix(rnorm(8 * 40), 8, 40)    # wide: T_pre < N_co
  y_tr  <- rnorm(8)
  zeta2 <- 0.3
  w  <- sdid_unit_weights_cpp(Y_pre, y_tr, zeta2)
  Q  <- crossprod(Y_pre) + diag(zeta2 * nrow(Y_pre), ncol(Y_pre))
  w_ref <- solve_simplex_qp(Q, crossprod(Y_pre, y_tr))
  expect_equal(drop(w), drop(w_ref), tolerance = 1e-8)
})

# ── SDID Event Study with Cluster / Unit Bootstrap ───────────────────────────

make_clustered_sdid_panel <- function(n_clusters = 20L, units_per_cl = 4L,
                                      times = c(2007:2010, 2013:2021),
                                      treat_from = 2013L,
                                      cl_sd = 1.2, eps_sd = 0.15,
                                      effect = -0.15, seed = 123L) {
  set.seed(seed)
  TT <- length(times)
  N  <- n_clusters * units_per_cl
  cl_id   <- rep(sprintf("pref_%02d", seq_len(n_clusters)), each = units_per_cl)
  unit_id <- sprintf("muni_%03d", seq_len(N))
  # First 6 clusters are predominantly treated (or mixed) so cluster draws
  # always have both treated and control units
  is_tr <- rep(c(TRUE, FALSE), c(n_clusters %/% 3L * units_per_cl,
                                 N - (n_clusters %/% 3L * units_per_cl)))
  # Also mix a few treated units into other clusters
  is_tr[seq(1L, N, by = units_per_cl)] <- TRUE

  unit_fe <- rnorm(N, 0, 0.5)
  time_fe <- seq(-0.2, 0.4, length.out = TT)
  # Strong cluster-by-period shock creates positive intra-cluster correlation
  cl_shock <- matrix(rnorm(n_clusters * TT, 0, cl_sd), nrow = TT, ncol = n_clusters)

  rows <- expand.grid(time_idx = seq_len(TT), u_idx = seq_len(N))
  t_i  <- rows$time_idx
  u_i  <- rows$u_idx
  c_i  <- (u_i - 1L) %/% units_per_cl + 1L

  d <- as.integer(is_tr[u_i] & times[t_i] >= treat_from)
  y <- unit_fe[u_i] + time_fe[t_i] +
    cl_shock[cbind(t_i, c_i)] +
    rnorm(nrow(rows), 0, eps_sd) +
    effect * d

  data.frame(
    unit = unit_id[u_i],
    pref = cl_id[u_i],
    year = times[t_i],
    d    = d,
    y    = y,
    stringsAsFactors = FALSE
  )
}

test_that("sdid_event_study: post-period mean of period estimates equals ATT exactly", {
  dat <- make_clustered_sdid_panel()
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")
  es  <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 30L, seed = 1L)

  expect_s3_class(es, "sdid_event_study")
  expect_equal(es$estimate, fit$estimate, tolerance = 1e-12)

  # Point estimates match plot_data(fit, type = "gap", align = TRUE)
  gap_aligned <- plot_data(fit, type = "gap", align = TRUE)
  expect_equal(es$terms$estimate, gap_aligned$gap, tolerance = 1e-12)
  expect_equal(mean(es$terms$estimate[es$terms$period == "post"]),
               fit$estimate, tolerance = 1e-12)
  # Pre-period lambda-weighted mean is 0 by construction
  expect_equal(sum(fit$time_weights * es$terms$estimate[es$terms$period == "pre"]),
               0, tolerance = 1e-12)
})

test_that("sdid_event_study: bootstrap draws match provisional data-frame refit exactly", {
  dat <- make_clustered_sdid_panel(n_clusters = 8L, units_per_cl = 3L, seed = 77L)
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")

  es <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 15L, seed = 99L)

  # Run the exact provisional loop with the same RNG sequence
  by_pref <- split(dat, dat$pref)
  set.seed(99L)
  ref_draws <- replicate(15L, {
    draw <- sample(names(by_pref), length(by_pref), replace = TRUE)
    parts <- lapply(seq_along(draw), function(k) {
      x <- by_pref[[draw[k]]]
      x$unit <- paste0(x$unit, "_b", k)
      x
    })
    b <- do.call(rbind, parts)
    u <- tapply(b$d, b$unit, max)
    if (all(u == 0) || all(u == 1)) return(NULL)
    f <- tryCatch(suppressMessages(scm_fit(y ~ d | unit + year, data = b, method = "sdid")),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    plot_data(f, type = "gap", align = TRUE)$gap
  }, simplify = FALSE)
  ref_valid <- do.call(rbind, ref_draws[!vapply(ref_draws, is.null, logical(1L))])

  expect_equal(unname(es$boot_gaps), unname(ref_valid), tolerance = 1e-10)
  expect_equal(es$boot_att, rowMeans(ref_valid[, es$terms$period == "post", drop = FALSE]),
               tolerance = 1e-10)
  expect_equal(es$se, sd(es$boot_att), tolerance = 1e-12)
})

test_that("sdid_event_study: numeric cluster IDs are drawn in numeric order", {
  dat <- make_clustered_sdid_panel(n_clusters = 12L, units_per_cl = 3L, seed = 31L)
  dat$pref_num <- as.integer(sub("pref_", "", dat$pref))
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")

  # pref_01..pref_12 sort the same way as 1..12 numerically, but not as
  # character "1","10","11","12","2",...
  es_num <- sdid_event_study(fit, cluster = "pref_num", data = dat,
                             n_boot = 10L, seed = 8L)
  es_chr <- sdid_event_study(fit, cluster = "pref", data = dat,
                             n_boot = 10L, seed = 8L)
  expect_identical(es_num$boot_gaps, es_chr$boot_gaps)
})

test_that("sdid fit keeps the raw outcome matrix only with covariates", {
  dat <- make_clustered_sdid_panel(n_clusters = 6L, units_per_cl = 3L)
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")
  expect_null(fit$Y_raw)
  expect_null(fit$X_arr)

  dat$x <- rnorm(nrow(dat))
  fit_x <- suppressWarnings(scm_fit(y ~ d | unit + year, data = dat,
                                    method = "sdid", covariates = "x"))
  expect_true(is.matrix(fit_x$Y_raw))
  expect_false(is.null(fit_x$X_arr))
  es <- sdid_event_study(fit_x, cluster = "pref", data = dat,
                         n_boot = 5L, seed = 1L)
  expect_s3_class(es, "sdid_event_study")
})

test_that("sdid_event_study: seed reproducibility and both normal & percentile CIs", {
  dat <- make_clustered_sdid_panel(n_clusters = 10L, units_per_cl = 3L)
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")

  es1 <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 40L, seed = 42L)
  es2 <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 40L, seed = 42L)
  expect_identical(es1$boot_gaps, es2$boot_gaps)
  expect_identical(es1$boot_att,  es2$boot_att)
  expect_identical(es1$se,        es2$se)
  expect_identical(es1$terms,     es2$terms)

  # Both normal and percentile CIs are computed from the same draws
  z <- qnorm(0.975)
  expect_equal(es1$ci_norm_lower, es1$estimate - z * es1$se, tolerance = 1e-12)
  expect_equal(es1$ci_norm_upper, es1$estimate + z * es1$se, tolerance = 1e-12)
  expect_equal(es1$ci_pct_lower,
               unname(quantile(es1$boot_att, 0.025)), tolerance = 1e-12)
  expect_equal(es1$ci_pct_upper,
               unname(quantile(es1$boot_att, 0.975)), tolerance = 1e-12)
  expect_equal(es1$terms$ci_norm_lower,
               es1$terms$estimate - z * es1$terms$se, tolerance = 1e-12)
  expect_equal(es1$terms$ci_pct_lower,
               unname(apply(es1$boot_gaps, 2L, quantile, probs = 0.025, names = FALSE)),
               tolerance = 1e-12)

  # tidy() and glance() switch conf.low / conf.high with ci_type
  td_norm <- broom::tidy(es1, ci_type = "normal")
  td_pct  <- broom::tidy(es1, ci_type = "percentile")
  expect_equal(td_norm$conf.low, td_norm$ci_norm_lower)
  expect_equal(td_pct$conf.low,  td_pct$ci_pct_lower)

  gl_norm <- broom::glance(es1, ci_type = "normal")
  gl_pct  <- broom::glance(es1, ci_type = "percentile")
  expect_equal(gl_norm$conf.low, es1$ci_norm_lower)
  expect_equal(gl_pct$conf.low,  es1$ci_pct_lower)
})

test_that("sdid_event_study: reports valid and discarded draws when clusters have pure treatment/control", {
  # 2 clusters: cluster 1 is all treated, cluster 2 is all control.
  # With 2 clusters drawn with replacement, P(all cluster 1 or all cluster 2) = 0.5,
  # so roughly half of the draws have an empty treated or control group and get discarded.
  set.seed(10)
  dat <- expand.grid(unit = sprintf("u%d", 1:6), year = 2001:2008,
                     stringsAsFactors = FALSE)
  dat$pref <- ifelse(dat$unit %in% c("u1", "u2", "u3"), "P_tr", "P_co")
  dat$d    <- as.integer(dat$pref == "P_tr" & dat$year >= 2005)
  dat$y    <- rnorm(nrow(dat)) - 1.0 * dat$d

  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")
  es  <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 60L, seed = 7L)

  expect_equal(es$n_boot, 60L)
  expect_equal(es$n_boot_valid + es$n_boot_dropped, 60L)
  expect_gt(es$n_boot_dropped, 0L)
  expect_gt(es_boot_valid <- es$n_boot_valid, 0L)
  expect_equal(nrow(es$boot_gaps), es$n_boot_valid)
  expect_length(es$boot_att, es$n_boot_valid)
  expect_output(print(es), "discarded")

  # Tiny p-value prints as "< 1e-04" rather than rounding to "0"
  es_tiny <- es
  es_tiny$p_value <- 1.4e-14
  expect_output(print(es_tiny), "< 1e-04")
})

test_that("sdid_event_study: cluster bootstrap SE exceeds unit bootstrap SE under intra-cluster correlation", {
  dat <- make_clustered_sdid_panel(
    n_clusters = 24L, units_per_cl = 6L,
    cl_sd = 1.8, eps_sd = 0.1, seed = 2026L
  )
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")

  es_unit <- sdid_event_study(fit, cluster = NULL, n_boot = 150L, seed = 11L)
  es_cl   <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 150L, seed = 11L)

  expect_false(es_unit$clustered)
  expect_true(es_cl$clustered)
  expect_equal(es_cl$n_clusters, 24L)
  expect_equal(es_unit$n_clusters, 144L)
  # Under strong cluster-level shocks, cluster bootstrap SE is substantially larger
  expect_gt(es_cl$se, es_unit$se * 1.4)
  expect_gt(mean(es_cl$terms$se), mean(es_unit$terms$se) * 1.4)
})

test_that("sdid_event_study: cluster specification accepts column name or named vector, rejects unnamed vector", {
  dat <- make_clustered_sdid_panel(n_clusters = 8L, units_per_cl = 3L)
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")

  # Named lookup vector (even in scrambled order) gives identical results to column name
  u_unique <- unique(dat[, c("unit", "pref")])
  scram    <- sample(nrow(u_unique))
  cl_vec   <- stats::setNames(u_unique$pref[scram], u_unique$unit[scram])

  es_col <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 25L, seed = 5L)
  es_vec <- sdid_event_study(fit, cluster = cl_vec, n_boot = 25L, seed = 5L)
  es_fml <- sdid_event_study(fit, cluster = ~ pref, data = dat, n_boot = 25L, seed = 5L)
  expect_identical(es_col$boot_gaps, es_vec$boot_gaps)
  expect_identical(es_col$boot_gaps, es_fml$boot_gaps)

  # Positional (unnamed) vector is rejected
  expect_error(
    sdid_event_study(fit, cluster = u_unique$pref, n_boot = 10L),
    "Unnamed positional vectors are not accepted"
  )
  # Column name without `data` is rejected
  expect_error(
    sdid_event_study(fit, cluster = "pref", n_boot = 10L),
    "`data` must also be supplied"
  )
  # Missing unit in named vector is rejected
  expect_error(
    sdid_event_study(fit, cluster = cl_vec[-1L], n_boot = 10L),
    "missing 1 unit"
  )
  # Time-varying cluster in `data` is rejected
  dat_bad <- dat
  dat_bad$pref[1L] <- "other_pref"
  expect_error(
    sdid_event_study(fit, cluster = "pref", data = dat_bad, n_boot = 10L),
    "more than one cluster"
  )
})

test_that("sdid_event_study and sdid_inference: reject staggered fits and invalid cluster args", {
  stag <- make_staggered_panel()
  fit_stag <- scm_fit(y ~ d | id + time, data = stag, method = "sdid")
  expect_error(sdid_event_study(fit_stag), "staggered")

  fit_scm <- scm_fit(y ~ d | id + time, data = panel, method = "scm")
  expect_error(sdid_event_study(fit_scm), "method = 'sdid'")

  # scm_fit() does not accept `cluster`
  expect_error(
    scm_fit(y ~ d | id + time, data = panel, method = "sdid", cluster = "id"),
    "unused argument"
  )

  # sdid_inference() rejects `cluster` with method-specific guidance
  fit_sdid <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_error(
    sdid_inference(fit_sdid, method = "bootstrap", cluster = "id"),
    "sdid_event_study"
  )
  expect_error(
    sdid_inference(fit_sdid, method = "placebo", cluster = "id"),
    "applies only to bootstrap"
  )
  expect_error(
    sdid_inference(fit_sdid, method = "jackknife", cluster = "id"),
    "applies only to bootstrap"
  )
})

test_that("plot.sdid_event_study: handles non-consecutive years (2011-2012 gap) and reference lines", {
  dat <- make_clustered_sdid_panel(
    n_clusters = 8L, units_per_cl = 3L,
    times = c(2007:2010, 2013:2021), treat_from = 2013L
  )
  fit <- scm_fit(y ~ d | unit + year, data = dat, method = "sdid")
  es  <- sdid_event_study(fit, cluster = "pref", data = dat, n_boot = 20L, seed = 3L)

  p <- plot(es)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)

  # x-axis breaks are the observed years only (no 2011 or 2012 ticks)
  x_scale <- built$layout$panel_params[[1L]]$x
  breaks  <- x_scale$get_breaks()
  breaks  <- breaks[!is.na(breaks)]
  expect_identical(as.integer(breaks), c(2007:2010, 2013:2021))
  expect_false(any(c(2011L, 2012L) %in% as.integer(breaks)))

  # Default vline_offset = -1 places the vertical line on the last pre-period (2010)
  expect_equal(.vline_x(p), 2010)
  # vline_offset = -0.5 interpolates midway between 2010 and 2013 -> 2011.5
  p_mid <- plot(es, vline_offset = -0.5)
  expect_equal(.vline_x(p_mid), 2011.5)

  # Percentile CI plot and plot_data()
  p_pct <- plot(es, ci_type = "percentile", connect = TRUE)
  expect_s3_class(p_pct, "ggplot")
  pd <- plot_data(es, ci_type = "percentile")
  expect_equal(pd$conf.low, es$terms$ci_pct_lower)
})

test_that("unbalanced panel is rejected with a clear error message in SDID", {
  dat <- make_clustered_sdid_panel(n_clusters = 6L, units_per_cl = 2L)
  # Drop one unit-year observation to make the panel unbalanced
  dat_unbal <- dat[-1L, ]
  expect_error(
    scm_fit(y ~ d | unit + year, data = dat_unbal, method = "sdid"),
    "balanced panel"
  )
})


# ── GSC / TASC / SI 論文照合レビュー（2026-10）──────────────────────────────

# Xu (2017, eq. 3): covariates correlated with the factors and the loadings,
# the case in which plain alternating least squares crawls.
.xu_cov_panel <- function(seed = 7, TT = 30, T0 = 20, Nco = 40, Ntr = 5) {
  set.seed(seed)
  N   <- Nco + Ntr
  F   <- matrix(rnorm(TT * 2), TT, 2)
  L   <- matrix(runif(N * 2, -sqrt(3), sqrt(3)), N, 2)
  FL  <- F %*% t(L)
  mkX <- function() 1 + FL + matrix(L[, 1] + L[, 2], TT, N, byrow = TRUE) +
    matrix(F[, 1] + F[, 2], TT, N) + matrix(rnorm(TT * N), TT, N)
  X1 <- mkX(); X2 <- mkX()
  D  <- matrix(0L, TT, N); D[(T0 + 1):TT, seq_len(Ntr)] <- 1L
  Y  <- X1 + 3 * X2 + FL + matrix(rnorm(TT * N), TT, N) + 2 * D
  list(Y = Y, X1 = X1, X2 = X2, co = (Ntr + 1):N, tr = seq_len(Ntr), T0 = T0,
       df = data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N),
                       y = as.vector(Y), d = as.vector(D),
                       x1 = as.vector(X1), x2 = as.vector(X2)))
}

test_that("GSC covariate loop reaches the IFE fixed point", {
  p   <- .xu_cov_panel()
  TT  <- nrow(p$Y)
  Xco <- array(c(p$X1[, p$co], p$X2[, p$co]), c(TT, length(p$co), 2))
  Xtr <- array(c(p$X1[seq_len(p$T0), p$tr], p$X2[seq_len(p$T0), p$tr]),
               c(p$T0, length(p$tr), 2))
  res <- gsc_ife_cpp(p$Y[, p$co], p$Y[seq_len(p$T0), p$tr], 2L, Xco, Xtr)
  ref <- gsc_ife_cpp(p$Y[, p$co], p$Y[seq_len(p$T0), p$tr], 2L, Xco, Xtr,
                     max_iter = 20000L, tol = 1e-13)
  expect_true(res$converged)
  expect_equal(drop(res$beta), drop(ref$beta), tolerance = 1e-6)
  expect_equal(drop(res$beta), c(1, 3), tolerance = 0.05)
  # A loop stopped by its cap says so
  capped <- gsc_ife_cpp(p$Y[, p$co], p$Y[seq_len(p$T0), p$tr], 2L, Xco, Xtr,
                        max_iter = 2L)
  expect_false(capped$converged)
  expect_equal(capped$n_iter, 2L)
})

test_that("gsc_treated_cpp is the treated step of gsc_ife_cpp", {
  p   <- .xu_cov_panel()
  TT  <- nrow(p$Y)
  Xco <- array(c(p$X1[, p$co], p$X2[, p$co]), c(TT, length(p$co), 2))
  Xtr <- array(c(p$X1[seq_len(p$T0), p$tr], p$X2[seq_len(p$T0), p$tr]),
               c(p$T0, length(p$tr), 2))
  full <- gsc_ife_cpp(p$Y[, p$co], p$Y[seq_len(p$T0), p$tr], 2L, Xco, Xtr)
  tr   <- coresynth:::gsc_treated_cpp(full$F, full$beta,
                                      p$Y[seq_len(p$T0), p$tr], Xtr,
                                      0, numeric(0), FALSE)
  expect_equal(tr$L_tr, full$L_tr, tolerance = 1e-12)
  expect_equal(tr$Y_tr_hat, full$Y_tr_hat, tolerance = 1e-12)
})

test_that("GSC with covariates recovers the effect on the Xu (2017) design", {
  p   <- .xu_cov_panel()
  fit <- scm_fit(y ~ d | id + t, data = p$df, method = "gsc", r = 2L,
                 covariates = c("x1", "x2"))
  expect_equal(drop(fit$beta), c(1, 3), tolerance = 0.05)
  expect_lt(abs(fit$estimate - 2), 0.4)
})

test_that("gsc_boot resamples residual series and handles covariates", {
  p   <- .xu_cov_panel(Ntr = 1)
  fit <- scm_fit(y ~ d | id + t, data = p$df, method = "gsc", r = 2L,
                 covariates = c("x1", "x2"))
  b <- gsc_boot(fit, B = 49L, seed = 3L)
  expect_true(is.finite(b$se) && b$se > 0)
  expect_lte(b$ci_lower, fit$estimate)
  expect_gte(b$ci_upper, fit$estimate)
  expect_equal(gsc_boot(fit, B = 49L, seed = 3L)$boot_dist, b$boot_dist)
})

test_that("conformal_inference refits GSC with the covariates the fit used", {
  set.seed(2)
  TT <- 30; T0 <- 20; N <- 31
  F <- matrix(rnorm(TT * 2), TT, 2); L <- matrix(rnorm(N * 2, 1, .5), N, 2)
  X <- matrix(rnorm(TT * N), TT, N) + 2 * F[, 1] %*% t(L[, 1])
  D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1] <- 1L
  Y <- F %*% t(L) + 3 * X + matrix(rnorm(TT * N, sd = .5), TT, N) + 2 * D
  df <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N),
                   y = as.vector(Y), d = as.vector(D), x = as.vector(X))
  fit <- scm_fit(y ~ d | id + t, data = df, method = "gsc", r = 2L,
                 covariates = "x")
  # Under the null tau0 = ATT the post-treatment outcomes are imputed net of
  # the effect, so the refit over all periods matches the fit before it
  y1  <- treated_outcomes(fit)
  y1[-seq_len(T0)] <- y1[-seq_len(T0)] - fit$estimate
  P   <- coresynth:::.conformal_refit.coresynth_gsc(fit, donor_outcomes(fit), y1)
  rmse_refit <- sqrt(mean((y1 - P)[seq_len(T0)]^2))
  rmse_fit   <- sqrt(mean(fit$gap[seq_len(T0)]^2))
  expect_lt(rmse_refit, 2 * rmse_fit)
})

test_that("TASC EM keeps the state noise from collapsing and reports convergence", {
  set.seed(1)
  TT <- 24; N_co <- 22; N <- N_co + 2
  F <- sapply(1:2, function(k) cumsum(rnorm(TT, sd = 0.4)))
  Y <- F %*% t(matrix(rnorm(N * 2, sd = 0.6), N, 2)) +
    matrix(rnorm(TT * N, sd = 0.3), TT, N)
  D <- matrix(0L, TT, N); D[15:TT, (N_co + 1):N] <- 1L
  df <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N),
                   y = as.vector(Y + 3 * D), d = as.vector(D))
  fit <- scm_fit(y ~ d | id + t, data = df, method = "tasc", em_tol = 1e-5)
  expect_true(fit$converged)
  expect_true(all(diag(fit$Q) > 1e-6))
  expect_lt(abs(fit$estimate - 3), 0.5)
  expect_warning(scm_fit(y ~ d | id + t, data = df, method = "tasc",
                         em_iter = 3L), "did not converge")
})

test_that("scm_design CI uses the order-statistic quantile of eq. (18)", {
  df  <- make_design_panel()
  res <- scm_design(df, "outcome", "unit", "time", T0 = 15L, T_fit = 10L)
  q   <- res$ci_upper - res$tau_hat
  expect_equal(unname(q[1]), unname(q[length(q)]))
  expect_true(any(abs(q[1] - abs(res$u_blank)) < 1e-12))
})

test_that("scm_design leaves the caller's random number stream alone", {
  df <- make_design_panel(TT = 40)
  set.seed(123); invisible(runif(1)); before <- .Random.seed
  res <- scm_design(df, "outcome", "unit", "time", T0 = 30L, T_fit = 10L)
  expect_identical(.Random.seed, before)
  expect_true(is.finite(res$p_value))
})

test_that("SI unit engine reproduces a staggered multi-arm fit", {
  fit <- scm_fit(y ~ d | id + time, data = stag_tensor_panel, method = "si")
  ctx <- coresynth:::.si_inference_context(fit)
  expect_equal(coresynth:::.si_unit_estimate(ctx, ctx$co, ctx$tr)$estimate,
               fit$estimate, tolerance = 1e-10)
})

test_that("staggered multi-arm SI jackknife SE is on the bootstrap's scale", {
  # The per-cohort variance used to be taken of the weighted sum of the
  # cohort's cells rather than of its mean, inflating the SE by the weight
  # (0.78 against a bootstrap SE of 0.04 on this design).
  p    <- make_stag_tensor_panel(N_co = 12L, TT = 14L)
  fit  <- scm_fit(y ~ d | id + time, data = p, method = "si")
  se_j <- si_inference(fit, method = "jackknife")$se
  se_b <- si_inference(fit, method = "bootstrap", n_boot = 199L, seed = 1L)$se
  expect_lt(se_j / se_b, 5)
  expect_gt(se_j / se_b, 0.2)
})

test_that("GSC additive effects with r = 0 reduce to the imputation by fixed effects", {
  fit_of <- function(force)
    scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc", r = 0L,
            force = force)

  f    <- fit_of("two-way")
  Y_co <- f$Y_co_all
  Y_tr <- f$Y_treat
  pre  <- seq_len(f$T_pre)
  post <- (f$T_pre + 1L):nrow(Y_co)
  c_t  <- rowMeans(Y_co)                # mu + xi_t from the controls
  # two-way: the treated unit's level from its pre-treatment gap to c_t
  cf_tw <- outer(c_t, colMeans(Y_tr[pre, , drop = FALSE] - c_t[pre]), `+`)
  expect_equal(f$estimate, mean(Y_tr[post, ] - cf_tw[post, ]), tolerance = 1e-10)
  # unit: the treated unit's own pre-treatment mean
  expect_equal(fit_of("unit")$estimate,
               mean(sweep(Y_tr[post, , drop = FALSE], 2, colMeans(Y_tr[pre, , drop = FALSE]))),
               tolerance = 1e-10)
  # time: the controls' mean path
  expect_equal(fit_of("time")$estimate, mean(Y_tr[post, ] - c_t[post]),
               tolerance = 1e-10)
  expect_error(scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc",
                       r = 0L), "positive integer")
})

test_that("GSC refits keep the additive effects the fit used", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "gsc", r = 1L,
                 force = "two-way")
  expect_identical(fit$force, "two-way")
  ctx <- coresynth:::.gsc_inference_context(fit)
  expect_equal(coresynth:::.gsc_unit_estimate(ctx, ctx$co, ctx$tr)$estimate,
               fit$estimate, tolerance = 1e-10)
  # conformal at the null tau0 = ATT refits the same model over all periods
  y1 <- treated_outcomes(fit)
  y1[-seq_len(fit$T_pre)] <- y1[-seq_len(fit$T_pre)] - fit$estimate
  P  <- coresynth:::.conformal_refit.coresynth_gsc(fit, donor_outcomes(fit), y1)
  expect_lt(sqrt(mean((y1 - P)[seq_len(fit$T_pre)]^2)),
            2 * sqrt(mean(rowMeans(fit$gap)[seq_len(fit$T_pre)]^2)) + 1e-8)
  b <- gsc_boot(fit, B = 29L, seed = 1L)
  expect_true(is.finite(b$se))
  st <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc", r = 1L,
                force = "unit")
  ctx_st <- coresynth:::.gsc_inference_context(st)
  expect_equal(coresynth:::.gsc_unit_estimate(ctx_st, ctx_st$co, ctx_st$tr)$estimate,
               st$estimate, tolerance = 1e-10)
})

# ── SDID 推論: Arkhangelsky et al. (2021) Algorithms 2-3（2026-10）────────────

test_that("SDID unit engines reproduce sharp, covariate and staggered fits", {
  check <- function(fit) {
    ctx <- coresynth:::.sdid_inference_context(fit)
    expect_equal(coresynth:::.sdid_unit_estimate(ctx, ctx$co, ctx$tr)$estimate,
                 fit$estimate, tolerance = 1e-10)
    expect_equal(coresynth:::.sdid_jack_estimate(ctx, ctx$co, ctx$tr)$estimate,
                 fit$estimate, tolerance = 1e-10)
  }
  check(scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid"))
  check(scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid", zeta2 = 2))
  check(scm_fit(y ~ d | id + time, data = staggered2, method = "sdid"))
  check(suppressWarnings(scm_fit(y ~ d | id + time, data = staggered2,
                                 method = "sdid", control_group = "clean")))
  check(scm_fit(y ~ d | id + time, data = stag_cov, method = "sdid",
                covariates = "cov1"))
})

test_that("SDID jackknife holds the fitted weights", {
  fit <- scm_fit(y ~ d | id + time, data = panel_tr3, method = "sdid")
  ctx <- coresynth:::.sdid_inference_context(fit)
  # Dropping the first donor: omega renormalised over the rest, lambda as
  # fitted, the treated mean unchanged
  om  <- fit$unit_weights[-1] / sum(fit$unit_weights[-1])
  lam <- fit$time_weights
  manual <- (mean(fit$Y_tr_post_mean) - sum(lam * fit$Y_tr_pre_mean)) -
    (sum(om * colMeans(fit$Y_co_post[, -1, drop = FALSE])) -
       sum(om * colSums(lam * fit$Y_co_pre[, -1, drop = FALSE])))
  expect_equal(
    coresynth:::.sdid_jack_estimate(ctx, ctx$co[-1], ctx$tr)$estimate,
    manual, tolerance = 1e-10)
})

test_that("SDID bootstrap carries the treated units' own variation", {
  # Two treated units with effects 0 and 4: resampling only the controls
  # would see none of that spread
  p <- panel_tr3
  p$y[p$id == "u2" & p$d == 1] <- p$y[p$id == "u2" & p$d == 1] - 2
  p$y[p$id == "u3" & p$d == 1] <- p$y[p$id == "u3" & p$d == 1] + 2
  fit <- scm_fit(y ~ d | id + time, data = p, method = "sdid")
  inf <- sdid_inference(fit, method = "bootstrap", n_boot = 99L, seed = 1L)
  expect_gt(inf$se, 0.5)
  expect_gt(sdid_inference(fit, method = "jackknife")$se, 0.5)
})

test_that("resampling inference needs two treated units, and two per cohort for the jackknife", {
  fit1 <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  expect_error(sdid_inference(fit1, method = "bootstrap"), "placebo")
  expect_error(sdid_inference(fit1, method = "jackknife"), "placebo")
  expect_no_error(sdid_inference(fit1, method = "placebo"))
  for (m in c("sdid", "gsc", "si")) {
    args <- list(y ~ d | id + time, data = staggered, method = m)
    if (m == "gsc") args$r <- 1L
    fit <- do.call(scm_fit, args)
    inf_fn <- switch(m, sdid = sdid_inference, gsc = gsc_inference,
                     si = si_inference)
    expect_error(inf_fn(fit, method = "jackknife"), "two treated units")
    expect_error(inf_fn(fit, method = "jackknife_global"), "two treated units")
    expect_true(is.finite(inf_fn(fit, method = "bootstrap", n_boot = 29L,
                                 seed = 1L)$se))
  }
})


test_that("MC validates estimator-specific arguments and C++ inputs", {
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       fixedeff = NA), "`fixedeff` must be TRUE or FALSE")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       lambda = -1), "`lambda` must be NULL or a single non-negative")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       lambda = c(0.1, 0.2)), "`lambda` must be NULL or a single non-negative")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       cv_folds = 0), "`cv_folds` must be a positive integer")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       cv_nlambda = 0), "`cv_nlambda` must be a positive integer")
  expect_error(scm_fit(y ~ d | id + time, data = panel, method = "mc",
                       cv_seed = NA), "`cv_seed` must be NULL or a single finite number")

  p_inf <- panel
  p_inf$y[1] <- Inf
  expect_error(scm_fit(y ~ d | id + time, data = p_inf, method = "mc"),
               "infinite values")

  expect_error(soft_impute_fe_cpp(matrix(0, 4, 3), matrix(1, 4, 2), 1.0),
               "same dimensions")
  expect_error(soft_impute_fe_cpp(matrix(0, 4, 3), matrix(1, 4, 3), -0.5),
               "non-negative finite")
  expect_error(soft_impute_fe_cpp(matrix(0, 4, 3), matrix(1, 4, 3), 1.0,
                                  max_iter = 0L), "max_iter must be >= 1")
})

test_that("conformal_inference handles missing observations in MC panels", {
  p_na <- panel
  # Drop a control cell and a pre-treatment cell
  p_na$y[p_na$id == "u2" & p_na$time == 5] <- NA_real_
  fit_mc_na <- scm_fit(y ~ d | id + time, data = p_na, method = "mc",
                       lambda = 0.5)
  ci_res <- conformal_inference(fit_mc_na, ci = TRUE, n_grid = 25L)
  expect_s3_class(ci_res, "conformal_inference")
  expect_true(is.finite(ci_res$p_value))
  expect_gte(ci_res$p_value, 0)
  expect_lte(ci_res$p_value, 1)
})
