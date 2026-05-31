# Non-parametric Inference for GSC (Xu 2017)

Estimates SE and confidence intervals for the ATT via non-parametric
cluster bootstrap or jackknife over control units. Works for both sharp
and staggered GSC fits. For staggered fits, bootstrap resamples each
cohort's control pool independently, and jackknife uses a per-cohort LOO
with delta-method variance aggregation.

## Usage

``` r
gsc_inference(
  fit,
  method = c("bootstrap", "jackknife", "jackknife_global"),
  n_boot = 499L,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  seed = NULL
)
```

## Arguments

- fit:

  A `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) with
  `method = "gsc"`.

- method:

  `"bootstrap"` (default) or `"jackknife"`.

- n_boot:

  Number of bootstrap replications (default 499L; ignored for
  jackknife).

- level:

  Confidence level (default 0.95).

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`.

- seed:

  RNG seed for reproducibility (default NULL).

## Value

A list of class `coresynth_inference`.

## Details

Note: [`gsc_boot()`](https://yo5uke.com/coresynth/reference/gsc_boot.md)
performs a *parametric* bootstrap under H0 (hypothesis testing).
`gsc_inference()` provides *non-parametric* SE and CIs suitable for
inference about the ATT magnitude.
