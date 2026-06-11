# Parametric Bootstrap Inference for GSC (Xu 2017 S.3)

Generates the null distribution of the ATT under H0 (no treatment
effect) by parametric resampling from the estimated IFE factor model.
Under H0, both the control panel and treated unit are generated from the
fitted factor model with homoskedastic noise. When the fit includes
covariate adjustment (beta), the covariate contribution is included in
the simulated DGP and re-estimated in each bootstrap replicate.

## Usage

``` r
gsc_boot(fit, B = 499L, alpha = 0.05, seed = NULL)
```

## Arguments

- fit:

  A `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) with
  `method = "gsc"`.

- B:

  Bootstrap replications (default 499L).

- alpha:

  Significance level for the confidence interval (default 0.05).

- seed:

  RNG seed for reproducibility (default NULL).

## Value

A list with:

- `p_value`: Two-sided p-value: mean(\|ATT\*\| \>= \|ATT_obs\|)

- `ci_lower`: Lower bound of (1-alpha)\*100% bootstrap CI

- `ci_upper`: Upper bound of (1-alpha)\*100% bootstrap CI

- `se`: Bootstrap standard error

- `boot_dist`: Numeric vector of length B (bootstrap ATT\* values)

- `att_obs`: Observed ATT from the original fit
