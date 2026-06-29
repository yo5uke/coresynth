# Staggered adoption

``` r

library(coresynth)
```

When units adopt treatment at different times, coresynth fits each
adoption **cohort** separately and aggregates the cohort ATTs with
weights proportional to `N_treated × T_post` (Clarke et al. 2023):

``` math
\hat\tau = \frac{\sum_g N_{tr,g}\, T_{post,g}\, \hat\tau_g}{\sum_g N_{tr,g}\, T_{post,g}}.
```

All six estimators support this; staggered timing is **detected
automatically** from the treatment column — no extra flag is needed.

## A staggered panel

Here `u1` is treated from period 11 and `u2` from period 16.

``` r

set.seed(42)
N <- 12; TT <- 20
f   <- cumsum(rnorm(TT, 0, 0.5))
lam <- rnorm(N, 1, 0.3)
dat <- expand.grid(time = seq_len(TT), id = paste0("u", seq_len(N)))
dat$y <- as.vector(outer(f, lam)) + rnorm(nrow(dat), 0, 0.3)

dat$d <- 0L
dat$d[dat$id == "u1" & dat$time > 10] <- 1L
dat$d[dat$id == "u2" & dat$time > 15] <- 1L
dat$y[dat$d == 1] <- dat$y[dat$d == 1] + 2.0   # true ATT = 2.0
```

## Fitting

``` r

fit <- scm_fit(y ~ d | id + time, data = dat, method = "sdid")
fit$estimate          # aggregate ATT
#> [1] 1.872811
fit$staggered         # TRUE
#> [1] TRUE
```

Per-cohort detail is in `cohort_estimates`:

``` r

fit$cohort_estimates
#>   cohort n_treated T_pre T_post estimate    weight
#> 1     11         1    10     10 1.839366 0.6666667
#> 2     16         1    15      5 1.939702 0.3333333
```

## Choosing the control group

`control_group` controls which units serve as donors for each cohort:

- `"clean"` (default) — never-treated units **plus** not-yet-treated
  units (those adopting later than the current cohort).
- `"never_treated"` — never-treated units only.

``` r

fit_clean <- scm_fit(y ~ d | id + time, data = dat, method = "sdid",
                     control_group = "clean")
fit_nt    <- scm_fit(y ~ d | id + time, data = dat, method = "sdid",
                     control_group = "never_treated")
c(clean = fit_clean$estimate, never_treated = fit_nt$estimate)
#>         clean never_treated 
#>      1.872811      1.962067
```

## Across estimators

The same call works for every method:

``` r

methods <- c("scm", "sdid", "gsc", "mc", "tasc", "si")
sapply(methods, function(m)
  scm_fit(y ~ d | id + time, data = dat, method = m)$estimate)
#>      scm     sdid      gsc       mc     tasc       si 
#> 1.835085 1.872811 1.760712 2.292777 3.153267 1.934303
```

## Inference under staggered adoption

[`sdid_inference()`](https://yo5uke.com/coresynth/reference/sdid_inference.md),
[`gsc_inference()`](https://yo5uke.com/coresynth/reference/gsc_inference.md),
and
[`si_inference()`](https://yo5uke.com/coresynth/reference/si_inference.md)
all extend to staggered fits. `jackknife_global` is the
staggered-specific variant that removes each unique control unit across
**all** cohorts at once, correctly capturing cross-cohort correlation.

``` r

library(broom)
tidy(sdid_inference(fit, method = "bootstrap", n_boot = 100, seed = 1))
#>   term estimate  std.error statistic     p.value conf.low conf.high    method
#> 1  ATT 1.872811 0.09200243  20.35611 4.09877e-92 1.711835  2.058967 bootstrap
#>   alternative n_controls staggered
#> 1   two.sided       10.5      TRUE
tidy(sdid_inference(fit, method = "jackknife_global"))
#>   term estimate std.error statistic      p.value conf.low conf.high
#> 1  ATT 1.872811 0.1185776  15.79397 3.423562e-56 1.640404  2.105219
#>             method alternative n_controls staggered
#> 1 jackknife_global   two.sided         11      TRUE
```

## Notes

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html) for
  staggered fits operate per cohort; the aggregate synthetic series is
  not defined (`Y_synth = NULL`).
- SCM staggered adoption supports `covariates =` (partial-out) but not
  the `predictors = pred(...)` interface.
- SI additionally supports staggered **and** multi-arm simultaneously —
  see
  [Estimators](https://yo5uke.com/coresynth/articles/estimators.html#multi-arm-si).
