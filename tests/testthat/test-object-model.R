# Object model: class tags, accessor generics, and S3 dispatch --------------

make_sharp_panel <- function(seed = 42L) {
  set.seed(seed)
  panel <- expand.grid(unit = 1:12, year = 1:24)
  panel$treated <- as.integer(panel$unit == 1 & panel$year > 18)
  panel$gdp <- panel$unit + 0.5 * panel$year +
    rnorm(nrow(panel)) + 3 * panel$treated
  panel
}

make_staggered_panel <- function(seed = 42L) {
  set.seed(seed)
  stag <- expand.grid(unit = 1:12, year = 1:24)
  adopt <- c(16, 19, 21, rep(NA, 9))
  stag$treated <- as.integer(!is.na(adopt[stag$unit]) &
                               stag$year >= adopt[stag$unit])
  stag$gdp <- stag$unit + 0.5 * stag$year +
    rnorm(nrow(stag)) + 2.5 * stag$treated
  stag
}

make_multiarm_panel <- function(seed = 42L) {
  set.seed(seed)
  ma <- expand.grid(unit = 1:14, year = 1:24)
  arm_of <- c(1, 1, 2, 2, rep(0, 10))
  ma$treated <- ifelse(arm_of[ma$unit] > 0 & ma$year > 18,
                       arm_of[ma$unit], 0L)
  ma$gdp <- ma$unit + 0.4 * ma$year + rnorm(nrow(ma)) +
    ifelse(ma$treated == 1, 2, ifelse(ma$treated == 2, -1.5, 0))
  ma
}

test_that("sharp fits carry method subclass but no structural subclass", {
  panel <- make_sharp_panel()
  for (m in c("scm", "sdid", "gsc", "mc", "tasc", "si")) {
    fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = m)
    expect_identical(class(fit), c(paste0("coresynth_", m), "coresynth"))
  }
})

test_that("staggered fits gain the coresynth_staggered subclass", {
  stag <- make_staggered_panel()
  for (m in c("scm", "sdid", "gsc")) {
    fit <- scm_fit(gdp ~ treated | unit + year, data = stag, method = m)
    expect_s3_class(fit, "coresynth_staggered")
    expect_s3_class(fit, paste0("coresynth_", m))
    expect_s3_class(fit, "coresynth")
    # staggered tag must precede the method tag so it wins dispatch
    expect_lt(match("coresynth_staggered", class(fit)),
              match(paste0("coresynth_", m), class(fit)))
  }
})

test_that("multi-arm SI fits gain the coresynth_multiarm subclass", {
  fit <- scm_fit(gdp ~ treated | unit + year, data = make_multiarm_panel(),
                 method = "si")
  expect_s3_class(fit, "coresynth_multiarm")
  expect_s3_class(fit, "coresynth_si")
  expect_false(inherits(fit, "coresynth_staggered"))
})

test_that("accessors return consistent series for every sharp method", {
  panel <- make_sharp_panel()
  for (m in c("scm", "sdid", "gsc", "mc", "si")) {
    fit <- scm_fit(gdp ~ treated | unit + year, data = panel, method = m)
    TT  <- length(fit$times)

    y1 <- treated_outcomes(fit)
    y0 <- synthetic_outcomes(fit)
    expect_length(y1, TT)
    expect_length(y0, TT)
    expect_true(all(is.finite(y1)))
    expect_true(all(is.finite(y0)))

    Yco <- donor_outcomes(fit)
    expect_true(is.matrix(Yco))
    expect_identical(nrow(Yco), TT)
    expect_identical(ncol(Yco), 11L)
  }
})

test_that("treated - synthetic reproduces the stored gap (SCM)", {
  fit <- scm_fit(gdp ~ treated | unit + year, data = make_sharp_panel(),
                 method = "scm")
  expect_equal(treated_outcomes(fit) - synthetic_outcomes(fit),
               unname(fit$gap))
})

test_that("donor_outcomes returns NULL where donors are not stored", {
  panel <- make_sharp_panel()
  fit_tasc <- scm_fit(gdp ~ treated | unit + year, data = panel,
                      method = "tasc")
  expect_null(donor_outcomes(fit_tasc))

  stag <- make_staggered_panel()
  fit_stag <- scm_fit(gdp ~ treated | unit + year, data = stag,
                      method = "scm")
  expect_null(donor_outcomes(fit_stag))  # per-cohort data live in cohort_fits
})

test_that("donor_outcomes falls back to field sniffing for legacy objects", {
  fit <- scm_fit(gdp ~ treated | unit + year, data = make_sharp_panel(),
                 method = "scm")
  legacy <- fit
  class(legacy) <- "coresynth"  # object saved by an older package version
  expect_identical(donor_outcomes(legacy), donor_outcomes(fit))
})

test_that("tidy dispatches on structural subclasses", {
  stag <- make_staggered_panel()
  fit_stag <- scm_fit(gdp ~ treated | unit + year, data = stag,
                      method = "sdid")
  td <- broom::tidy(fit_stag)
  expect_true(all(td$type == "cohort_estimate"))
  expect_identical(nrow(td), nrow(fit_stag$cohort_estimates))

  fit_ma <- scm_fit(gdp ~ treated | unit + year, data = make_multiarm_panel(),
                    method = "si")
  td_ma <- broom::tidy(fit_ma)
  # Sharp multi-arm fits fall through to the unit-weight/empty representation
  expect_s3_class(td_ma, "data.frame")

  fit_sharp <- scm_fit(gdp ~ treated | unit + year, data = make_sharp_panel(),
                       method = "scm")
  td_sharp <- broom::tidy(fit_sharp)
  expect_true(all(td_sharp$type == "unit_weight"))
})

test_that("augment dispatches to the staggered method", {
  stag <- make_staggered_panel()
  fit <- scm_fit(gdp ~ treated | unit + year, data = stag, method = "sdid")
  au <- broom::augment(fit)
  expect_true(".cohort" %in% names(au))
  expect_gt(nrow(au), 0L)
})

test_that("augment(include_donors) returns donor rows for SI fits", {
  fit <- scm_fit(gdp ~ treated | unit + year, data = make_sharp_panel(),
                 method = "si")
  au <- broom::augment(fit, include_donors = TRUE)
  expect_setequal(unique(au$.type), c("treated", "control"))
  expect_identical(sum(au$.type == "control"), 11L * 24L)
})

test_that("print and summary methods dispatch for multi-arm fits", {
  fit <- scm_fit(gdp ~ treated | unit + year, data = make_multiarm_panel(),
                 method = "si")
  expect_output(print(fit), "Multi-arm SI \\(K = 2 arms\\)")
  expect_output(summary(fit), "Per-arm ATT:")

  fit_sharp <- scm_fit(gdp ~ treated | unit + year, data = make_sharp_panel(),
                       method = "scm")
  out <- capture.output(print(fit_sharp))
  expect_false(any(grepl("Multi-arm", out)))
})

test_that("conformal inference rejects staggered fits via class", {
  stag <- make_staggered_panel()
  fit <- scm_fit(gdp ~ treated | unit + year, data = stag, method = "sdid")
  expect_error(conformal_inference(fit), "sharp")
})

test_that("plot gives a clear error for staggered fits", {
  stag <- make_staggered_panel()
  fit <- scm_fit(gdp ~ treated | unit + year, data = stag, method = "sdid")
  expect_error(plot(fit, type = "trend"), "per cohort")
})
