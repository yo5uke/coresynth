# Inference

``` r

library(coresynth)
library(broom)
```

coresynth ships method-appropriate inference for each estimator, plus a
cross-cutting **conformal** procedure. This article surveys them on a
single synthetic panel.

``` r

set.seed(123)
N <- 12; TT <- 20; T_pre <- 12
f   <- cumsum(rnorm(TT, 0, 0.5))
lam <- rnorm(N, 1, 0.3)
dat <- expand.grid(time = seq_len(TT), id = paste0("u", seq_len(N)))
dat$y <- as.vector(outer(f, lam)) + rnorm(nrow(dat), 0, 0.3)
dat$d <- as.integer(dat$id == "u1" & dat$time > T_pre)
dat$y[dat$d == 1] <- dat$y[dat$d == 1] + 2.0   # true ATT = 2.0

fit_scm  <- scm_fit(y ~ d | id + time, data = dat, method = "scm")
fit_sdid <- scm_fit(y ~ d | id + time, data = dat, method = "sdid")
fit_gsc  <- scm_fit(y ~ d | id + time, data = dat, method = "gsc", r = 2)
fit_si   <- scm_fit(y ~ d | id + time, data = dat, method = "si")
```

## SCM — MSPE ratio permutation test

*Abadie et al. (2010).*
[`mspe_ratio_pval()`](https://yo5uke.com/coresynth/reference/mspe_ratio_pval.md)
compares the treated unit’s post/pre MSPE ratio against placebo ratios
from every donor. It returns a plain list.

``` r

mspe <- mspe_ratio_pval(fit_scm, alternative = "two.sided")
c(p_value = mspe$p_value, ratio_obs = mspe$ratio_obs)
#>    p_value 
#> 0.08333333
```

## SDID — four inference methods

*Arkhangelsky et al. (2021); Clarke et al. (2023).*
[`sdid_inference()`](https://yo5uke.com/coresynth/reference/sdid_inference.md)
supports `"placebo"`, `"bootstrap"`, `"jackknife"`, and
`"jackknife_global"`. The result is `broom`-friendly.

``` r

inf_plac <- sdid_inference(fit_sdid, method = "placebo")
inf_boot <- sdid_inference(fit_sdid, method = "bootstrap", n_boot = 100, seed = 1)

tidy(inf_plac)
#>   term estimate std.error statistic    p.value conf.low conf.high  method
#> 1  ATT 1.821899        NA        NA 0.08333333       NA        NA placebo
#>   alternative n_controls staggered
#> 1   two.sided         11     FALSE
tidy(inf_boot)
#>   term estimate std.error statistic      p.value conf.low conf.high    method
#> 1  ATT 1.821899 0.1287232  14.15362 1.773922e-45  1.60038  2.027952 bootstrap
#>   alternative n_controls staggered
#> 1   two.sided         11     FALSE
```

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns a
one-row data frame with `estimate`, `std.error`, `p.value`, and a CI —
ready to [`rbind()`](https://rdrr.io/r/base/cbind.html) into a results
table. [`glance()`](https://generics.r-lib.org/reference/glance.html)
gives a compact summary:

``` r

glance(inf_boot)
#>      method n_controls staggered estimate std.error      p.value conf.low
#> 1 bootstrap         11     FALSE 1.821899 0.1287232 1.773922e-45  1.60038
#>   conf.high alternative n_boot_valid
#> 1  2.027952   two.sided          100
```

## GSC — parametric and non-parametric

*Xu (2017).*
[`gsc_boot()`](https://yo5uke.com/coresynth/reference/gsc_boot.md) is a
parametric bootstrap under H0 (sharp fits only);
[`gsc_inference()`](https://yo5uke.com/coresynth/reference/gsc_inference.md)
offers non-parametric bootstrap / jackknife.

``` r

gb <- gsc_boot(fit_gsc, B = 100, alpha = 0.05)
c(ci_lower = gb$ci_lower, ci_upper = gb$ci_upper, p_value = gb$p_value)
#>   ci_lower   ci_upper    p_value 
#> -0.3934672  0.2557330  0.0000000

tidy(gsc_inference(fit_gsc, method = "jackknife"))
#>   term estimate  std.error statistic p.value conf.low conf.high    method
#> 1  ATT 1.966522 0.04177117  47.07845       0 1.884652  2.048392 jackknife
#>   alternative n_controls staggered
#> 1   two.sided         11     FALSE
```

## SI — bootstrap / jackknife

*Agarwal et al. (2025).*
[`si_inference()`](https://yo5uke.com/coresynth/reference/si_inference.md)
mirrors the GSC non-parametric API and also handles staggered and
multi-arm fits.

``` r

tidy(si_inference(fit_si, method = "bootstrap", n_boot = 100, seed = 1))
#>   term estimate  std.error statistic      p.value conf.low conf.high    method
#> 1  ATT 2.009438 0.05357169  37.50932 6.49154e-308 1.881553  2.080903 bootstrap
#>   alternative n_controls staggered
#> 1   two.sided         11     FALSE
```

## Conformal inference (any sharp fit)

*Chernozhukov, Wüthrich & Zhu (2021).*
[`conformal_inference()`](https://yo5uke.com/coresynth/reference/conformal_inference.md)
works across `scm` / `sdid` / `gsc` / `mc` / `si` sharp fits. Under H0:
τ = τ0 it re-imputes the treated post-period as `Y1 - τ0`,
**re-estimates the counterfactual on all T periods**, and inverts a
moving-block permutation test of the residuals to get a p-value and a
confidence interval.

``` r

conf <- conformal_inference(fit_scm, tau0 = 0, level = 0.95)
tidy(conf)
#>   term estimate std.error statistic p.value   conf.low conf.high    method
#> 1  ATT 1.698439        NA        NA     0.1 -0.1072756  3.162243 conformal
#>   alternative n_controls staggered
#> 1   two.sided         11     FALSE
```

The p-value tests the sharp null τ0 = 0; the CI is obtained by test
inversion over a grid. Because the permutation uses \|Π\| = T cyclic
shifts, the smallest attainable p-value is 1/T.

## Assembling a results table

Since every inference object tidies to the same columns, comparing
methods is a single [`rbind()`](https://rdrr.io/r/base/cbind.html):

``` r

do.call(rbind, list(
  cbind(method = "sdid (placebo)",   tidy(inf_plac)[c("estimate","p.value","conf.low","conf.high")]),
  cbind(method = "sdid (bootstrap)", tidy(inf_boot)[c("estimate","p.value","conf.low","conf.high")]),
  cbind(method = "scm (conformal)",  tidy(conf)[c("estimate","p.value","conf.low","conf.high")])
))
#>             method estimate      p.value   conf.low conf.high
#> 1   sdid (placebo) 1.821899 8.333333e-02         NA        NA
#> 2 sdid (bootstrap) 1.821899 1.773922e-45  1.6003801  2.027952
#> 3  scm (conformal) 1.698439 1.000000e-01 -0.1072756  3.162243
```
