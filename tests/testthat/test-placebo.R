## In-space placebo (permutation) tests: predictor-spec consistency between
## treated and placebo refits, and agreement of the fast outcomes-only solver
## with the reference nested optimiser.

library(testthat)

# Balanced panel: 1 treated unit (u1), rest control, treatment from T_pre+1.
make_placebo_panel <- function(N = 10, T = 20, T_pre = 10, effect = 2.0,
                               seed = 42) {
  set.seed(seed)
  times <- 1:T
  units <- paste0("u", 1:N)
  f <- cumsum(rnorm(T, 0, 0.5))
  lam <- rnorm(N, 1, 0.3)
  rows <- expand.grid(time = times, id = units, stringsAsFactors = FALSE)
  rows$y <- as.vector(outer(f, lam)) + rnorm(nrow(rows), 0, 0.3)
  rows$d <- as.integer(rows$id == "u1" & rows$time > T_pre)
  rows$y[rows$d == 1] <- rows$y[rows$d == 1] + effect
  rows[order(rows$id, rows$time), ]
}

panel <- make_placebo_panel()

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
