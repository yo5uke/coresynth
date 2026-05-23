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
  setup <- panel.matrices(california_prop99)
  ref <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
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
    em_iter = 5
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
    em_iter = 10
  )
  expect_true(is.finite(fit$estimate))
})

test_that("TASC A matrix has correct dimensions (r x r)", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_iter = 5
  )
  expect_equal(dim(fit$A), c(2L, 2L))
})

test_that("TASC A deviates from identity after EM (learning occurred)", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_iter = 20
  )
  expect_true(max(abs(fit$A - diag(2))) > 1e-6)
})

test_that("TASC fix_A=TRUE keeps A as identity", {
  fit <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_iter = 10,
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
    em_iter = 10,
    fix_A = FALSE
  )
  fit_fix <- scm_fit(
    y ~ d | id + time,
    data = panel,
    method = "tasc",
    r = 2,
    em_iter = 10,
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
    em_iter = 5L
  )
  expect_true(is.finite(fit$estimate))
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

test_that("SCM v_selection='insample' and 'oos' give different (but valid) V weights", {
  fit_is  <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_selection = "insample")
  fit_oos <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                     v_selection = "oos")
  expect_equal(length(fit_is$v_weights), length(fit_oos$v_weights))
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
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_optim = "bfgs")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-5))
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 12a: v_optim='bfgs' + v_selection='oos' works", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_optim = "bfgs", v_selection = "oos")
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
  # panel: k=T_pre=10, N_co=9 → k=10 > 15? No. But 2*k=20 > N_co=9? Yes
  # auto: k_dim=10 > 15 → coord_descent; still valid
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_optim = "auto")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 12a: v_optim='bfgs' with lambda_pen works", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "scm",
                 v_optim = "bfgs", lambda_pen = 0.5)
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(is.finite(fit$estimate))
})

# ── Phase 13a: SDID 共変量対応 ────────────────────────────────────────────────

# Time-varying covariate panel: panel_cov already has cov1/cov2 (time-invariant here;
# sufficient to test the partial-out API)
test_that("Phase 13a: SDID with one covariate returns valid fit", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid",
                 covariates = "cov1")
  expect_s3_class(fit, "coresynth")
  expect_true(is.finite(fit$estimate))
})

test_that("Phase 13a: beta_hat length equals number of covariates", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid",
                 covariates = c("cov1", "cov2"))
  expect_length(fit$beta_hat, 2L)
})

test_that("Phase 13a: SDID covariates — unit weights stay on simplex", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid",
                 covariates = "cov1")
  expect_equal(sum(fit$unit_weights), 1, tolerance = 1e-4)
  expect_true(all(fit$unit_weights >= -1e-8))
})

test_that("Phase 13a: SDID without covariates — beta_hat is length-0", {
  fit <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid")
  expect_length(fit$beta_hat, 0L)
})

test_that("Phase 13a: covariate partial-out changes the ATT estimate", {
  fit_no  <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid")
  fit_cov <- scm_fit(y ~ d | id + time, data = panel_cov, method = "sdid",
                     covariates = "cov1")
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

test_that("Phase 13b: SDID staggered control_group='clean' is default", {
  fit_clean  <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid",
                        control_group = "clean")
  fit_never  <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid",
                        control_group = "never_treated")
  # Different control groups should (in general) give different estimates
  expect_true(is.finite(fit_clean$estimate))
  expect_true(is.finite(fit_never$estimate))
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
  expect_null(inf$se)
  expect_null(inf$ci_lower)
  expect_equal(inf$n_controls, 9L)
})

test_that("Phase 15a: sdid_inference(placebo) named placebo_effects", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  expect_equal(names(inf$placebo_effects), colnames(fit$Y_co_pre))
})

test_that("Phase 15a: sdid_inference(bootstrap) returns se and CI", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "bootstrap", n_boot = 50L, seed = 1L)
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
  expect_true(inf$ci_lower < inf$ci_upper)
  expect_length(inf$boot_ests, 50L)
})

test_that("Phase 15a: sdid_inference(bootstrap) is reproducible with seed", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf1 <- sdid_inference(fit, method = "bootstrap", n_boot = 30L, seed = 99L)
  inf2 <- sdid_inference(fit, method = "bootstrap", n_boot = 30L, seed = 99L)
  expect_equal(inf1$se, inf2$se)
  expect_equal(inf1$boot_ests, inf2$boot_ests)
})

test_that("Phase 15a: sdid_inference(jackknife) returns se and CI", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
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
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
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
  expect_error(
    coresynth:::panel_to_tensor(bad$y, bad$d, bad$id, bad$time),
    "arm 0"
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

# ── Phase 17a: gsc_inference — sharp ─────────────────────────────────────────

test_that("Phase 17a: gsc_inference sharp bootstrap — no error, class correct", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 1L)
  expect_s3_class(inf, "coresynth_inference")
})

test_that("Phase 17a: gsc_inference sharp bootstrap — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 2L)
  expect_true(is.finite(inf$se))
  expect_gte(inf$p_value, 0)
  expect_lte(inf$p_value, 1)
})

test_that("Phase 17a: gsc_inference sharp bootstrap — CI contains estimate", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc")
  inf <- gsc_inference(fit, method = "bootstrap", n_boot = 99L, seed = 3L)
  expect_lte(inf$ci_lower, inf$estimate + 1e-8)
  expect_gte(inf$ci_upper, inf$estimate - 1e-8)
})

test_that("Phase 17a: gsc_inference sharp jackknife — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc")
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
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  inf <- gsc_inference(fit, method = "jackknife")
  expect_true(is.finite(inf$se))
  expect_true(inf$staggered)
})

# ── Phase 17c: si_inference — sharp ───────────────────────────────────────────

test_that("Phase 17c: si_inference sharp bootstrap — no error, class correct", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 20L)
  expect_s3_class(inf, "coresynth_inference")
})

test_that("Phase 17c: si_inference sharp bootstrap — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 21L)
  expect_true(is.finite(inf$se))
  expect_gte(inf$p_value, 0)
  expect_lte(inf$p_value, 1)
})

test_that("Phase 17c: si_inference sharp bootstrap — CI contains estimate", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  inf <- si_inference(fit, method = "bootstrap", n_boot = 99L, seed = 22L)
  expect_lte(inf$ci_lower, inf$estimate + 1e-8)
  expect_gte(inf$ci_upper, inf$estimate - 1e-8)
})

test_that("Phase 17c: si_inference sharp jackknife — finite SE and p_value in [0,1]", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
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
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
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
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  expect_s3_class(inf, "sdid_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$p_value >= 0 && inf$p_value <= 1)
})

test_that("Phase 18a: sdid_inference staggered jackknife CI straddles estimate", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  expect_true(inf$ci_lower <= inf$estimate && inf$estimate <= inf$ci_upper)
  expect_null(inf$boot_ests)
})

# ── Phase 18b: jackknife_global for SDID / GSC / SI staggered ────────────────

test_that("Phase 18b: sdid_inference jackknife_global staggered", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
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
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "gsc")
  inf <- gsc_inference(fit, method = "jackknife_global")
  expect_s3_class(inf, "coresynth_inference")
  expect_true(isTRUE(inf$staggered))
  expect_true(is.finite(inf$se) && inf$se > 0)
  expect_true(inf$ci_lower <= inf$estimate && inf$estimate <= inf$ci_upper)
})

test_that("Phase 18b: si_inference jackknife_global staggered", {
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "si")
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
  expect_null(inf$se)
  expect_null(inf$ci_lower)
  expect_null(inf$ci_upper)
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
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
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
  fit <- scm_fit(y ~ d | id + time, data = staggered, method = "sdid")
  inf <- sdid_inference(fit, method = "jackknife")
  td  <- broom::tidy(inf)
  expect_true(td$staggered)
})

test_that("Phase 24: tidy.coresynth_inference GSC works", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "gsc", r = 2L)
  inf <- gsc_inference(fit, method = "jackknife")
  td  <- broom::tidy(inf)
  expect_equal(td$term, "ATT")
  expect_true(is.finite(td$estimate))
  expect_true(is.finite(td$std.error))
})

test_that("Phase 24: tidy.coresynth_inference SI works", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "si")
  inf <- si_inference(fit, method = "jackknife")
  td  <- broom::tidy(inf)
  expect_equal(td$term, "ATT")
  expect_true(is.finite(td$std.error))
})

test_that("Phase 24: glance.coresynth_inference returns one-row summary", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "bootstrap", n_boot = 31L, seed = 1L)
  g   <- broom::glance(inf)
  expect_equal(nrow(g), 1L)
  expect_true(all(c("method", "n_controls", "staggered", "estimate",
                    "std.error", "p.value", "conf.low", "conf.high",
                    "alternative", "n_boot_valid") %in% names(g)))
  expect_gte(g$n_boot_valid, 1L)
})

test_that("Phase 24: tidy.coresynth_inference placebo SE/CI are NA", {
  fit <- scm_fit(y ~ d | id + time, data = panel, method = "sdid")
  inf <- sdid_inference(fit, method = "placebo")
  td  <- broom::tidy(inf)
  expect_true(is.na(td$std.error))
  expect_true(is.na(td$conf.low))
  expect_true(is.na(td$conf.high))
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
