# Estimators

``` r

library(coresynth)
```

coresynth exposes six estimators through
[`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) plus an
experimental-design variant
([`scm_design()`](https://yo5uke.com/coresynth/reference/scm_design.md)).
This article covers the method-specific options for each. For the basics
of the shared interface, see [Get
started](https://yo5uke.com/coresynth/articles/coresynth.md). For
inference, see
[Inference](https://yo5uke.com/coresynth/articles/inference.md); for
multi-period adoption, see [Staggered
adoption](https://yo5uke.com/coresynth/articles/staggered.md).

We use one synthetic panel throughout, with two auxiliary covariates so
we can demonstrate predictor- and covariate-based matching.

``` r

set.seed(1)
N <- 12; TT <- 20; T_pre <- 12
f   <- cumsum(rnorm(TT, 0, 0.5))
lam <- rnorm(N, 1, 0.3)
dat <- expand.grid(time = seq_len(TT), id = paste0("u", seq_len(N)))
dat$y      <- as.vector(outer(f, lam)) + rnorm(nrow(dat), 0, 0.3)
dat$income <- as.vector(outer(f * 0.5, lam)) + rnorm(nrow(dat), 0, 0.3)  # predictor
dat$x      <- rnorm(nrow(dat))                                           # time-varying cov
dat$d      <- as.integer(dat$id == "u1" & dat$time > T_pre)
dat$y[dat$d == 1] <- dat$y[dat$d == 1] + 2.0   # true ATT = 2.0
```

## SCM — Synthetic Control Method

*Abadie, Diamond & Hainmueller (2010).* Donor weights on the simplex
match the treated unit’s pre-treatment trajectory.

### Predictors via `pred()`

By default (`predictors = NULL`) all pre-treatment outcomes form the
predictor matrix. To match on covariates and specific outcome lags
instead, build a list of `pred(vars, times, op)` specs:

``` r

fit_scm <- scm_fit(
  y ~ d | id + time, data = dat, method = "scm",
  predictors = list(
    pred("income", 1:T_pre),          # mean income over the pre-period
    pred("y", T_pre),                  # outcome in the last pre-period
    pred("y", 1:4, op = "mean")        # outcome averaged over early pre-period
  )
)
fit_scm$estimate
#> [1] 2.038808
```

### Out-of-sample V selection, donor filtering, penalisation

``` r

# Abadie (2021) S.3.2: split the pre-period to choose V out of sample
fit_oos <- scm_fit(y ~ d | id + time, data = dat, method = "scm",
                   v_selection = "oos")

# Abadie (2021) S.4: drop poorly-fitting donors by individual MSPE ratio
fit_filt <- scm_fit(y ~ d | id + time, data = dat, method = "scm",
                    donor_mspe_threshold = 5)

# Abadie & L'Hour (2021): penalised SCM (auto-selects the penalty out of sample)
fit_pen <- scm_fit(y ~ d | id + time, data = dat, method = "scm",
                   lambda_pen = "auto")

c(oos = fit_oos$estimate, filtered = fit_filt$estimate, penalised = fit_pen$estimate)
#>       oos  filtered penalised 
#>  2.583619  2.073836  2.073836
```

### Augmented SCM

[`augment_scm()`](https://yo5uke.com/coresynth/reference/augment_scm.md)
applies the ridge bias-correction of Ben-Michael, Feller & Rothstein
(2021):

``` r

aug <- augment_scm(fit_scm)
c(scm = aug$att_scm, augmented = aug$att_aug)
#>       scm augmented 
#>  2.038808  2.115741
```

## SDID — Synthetic Difference-in-Differences

*Arkhangelsky et al. (2021).* Combines unit and time weights with a
DiD-style double differencing. Time-varying covariates are partialled
out first (Clarke et al. 2023):

``` r

fit_sdid <- scm_fit(y ~ d | id + time, data = dat, method = "sdid",
                    covariates = "x")
fit_sdid$estimate
#> [1] 1.975303
```

## GSC — Generalised Synthetic Control

*Xu (2017).* Interactive fixed effects estimated by SVD; `r` sets the
number of latent factors. With covariates the full EM algorithm is used
(E-step SVD, M-step ridge OLS), and `fit$beta` holds the covariate
coefficients.

``` r

fit_gsc <- scm_fit(y ~ d | id + time, data = dat, method = "gsc",
                   r = 2, covariates = "x")
c(estimate = fit_gsc$estimate, beta = unname(fit_gsc$beta))
#>   estimate       beta 
#> 1.89883014 0.01877872
```

`r` sensitivity is easy to sweep:

``` r

sapply(1:3, function(rr)
  scm_fit(y ~ d | id + time, data = dat, method = "gsc", r = rr)$estimate)
#> [1] 1.912170 1.905871 1.875384
```

## MC — Matrix Completion

*Athey et al. (2021).* Nuclear-norm-regularised completion
(Soft-Impute). The penalty `lambda` defaults to `0.01 * sigma_max(Y)`;
pass a number to override.

``` r

fit_mc_auto <- scm_fit(y ~ d | id + time, data = dat, method = "mc")
fit_mc_man  <- scm_fit(y ~ d | id + time, data = dat, method = "mc", lambda = 1.0)
c(auto = fit_mc_auto$estimate, manual = fit_mc_man$estimate)
#>     auto   manual 
#> 2.165341 2.165341
```

> On small donor pools MC’s nuclear-norm shrinkage often pushes the ATT
> above the truth — a known property (Mazumder et al. 2010), not an
> implementation bug.

## TASC — Time-Aware Synthetic Control

*Rho et al. (2026).* A state-space model fitted by Kalman EM. `r` sets
the latent state dimension, `em_iter` the EM iterations, and
`fix_A = TRUE` constrains the transition matrix to be constant.

``` r

fit_tasc      <- scm_fit(y ~ d | id + time, data = dat, method = "tasc",
                         r = 2, em_iter = 10)
fit_tasc_fixA <- scm_fit(y ~ d | id + time, data = dat, method = "tasc",
                         r = 2, em_iter = 10, fix_A = TRUE)
c(free_A = fit_tasc$estimate, fixed_A = fit_tasc_fixA$estimate)
#>   free_A  fixed_A 
#> 2.045796 2.049410
```

## SI — Synthetic Interventions

*Agarwal et al. (2025).* SI-PCR uses a truncated SVD of the donor
pre-period. `k` sets the rank (default `floor(sqrt(min(T_pre, N_co)))`).

``` r

fit_si <- scm_fit(y ~ d | id + time, data = dat, method = "si")
fit_si$estimate
#> [1] 1.925724
```

### Multi-arm SI

SI uniquely supports multiple treatment arms (`K > 1`). When the
treatment column takes values `0, 1, ..., K`,
[`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md)
automatically routes to the multi-arm path: the control-arm SVD basis is
shared across all arms (the unit factors are arm-invariant). Use
`panel_to_tensor()` to prepare/inspect the tensor structure.

``` r

set.seed(7)
Nc <- 10; T2 <- 16; Tp <- 10
U <- matrix(rnorm(T2 * 2), T2, 2)
V <- matrix(rnorm((Nc + 4) * 2), Nc + 4, 2)
Lam <- rbind(c(1, 1), c(1.4, .6), c(.6, 1.4))   # arm 0, 1, 2
arm_unit <- c(rep(0L, Nc), 1L, 1L, 2L, 2L)
g <- expand.grid(time = seq_len(T2), id = seq_len(Nc + 4))
g$arm <- arm_unit[g$id]
g$d   <- ifelse(g$arm == 0L, 0L, ifelse(g$time > Tp, g$arm, 0L))
g$y   <- mapply(function(t, i, a)
  sum(U[t, ] * Lam[a + 1L, ] * V[i, ]) + 50 + rnorm(1, sd = .3),
  g$time, g$id, g$arm)

fit_multi <- scm_fit(y ~ d | id + time, data = g, method = "si")
fit_multi$arm_estimates       # per-arm ATTs
#>          1          2 
#> -0.2139400  0.1224433
fit_multi$estimate            # weighted aggregate
#> [1] -0.04574832
```

## SCM-Design — experimental design

*Abadie & Zhao (2026).*
[`scm_design()`](https://yo5uke.com/coresynth/reference/scm_design.md)
has a different, design-oriented interface (it selects which units to
expose to a planned intervention). It offers `"base"`,
`"weakly_targeted"`, and `"unit_level"` variants with blank-period
permutation tests and split-conformal CIs.

``` r

des <- scm_design(
  data = dat, outcome = "y", unit = "id", time = "time",
  T0 = T_pre, design = "base"
)
des
#> === scm_design (Abadie & Zhao 2026) ===
#> Design variant : base 
#> Treated units  : u3 
#> ATT estimate   : 0.0083 
#> p-value        : NA (no blank periods)
```
