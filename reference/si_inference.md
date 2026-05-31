# Non-parametric Inference for SI (Agarwal et al. 2025)

Estimates SE and confidence intervals for the ATT via non-parametric
cluster bootstrap or jackknife over control units. Works for both sharp
and staggered SI fits. For staggered fits, bootstrap resamples each
cohort's control pool independently, and jackknife uses a per-cohort LOO
with delta-method variance aggregation.

## Usage

``` r
si_inference(
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
  `method = "si"`.

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
