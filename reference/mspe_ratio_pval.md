# Permutation Inference via MSPE Ratio for SCM

Computes the Abadie et al. (2010) / Abadie (2021) permutation p-value.
For each control unit, a leave-one-out synthetic control is fitted.

## Usage

``` r
mspe_ratio_pval(
  fit,
  mspe_threshold = 0,
  max_iter = 100L,
  tol = 1e-04,
  use_covariates = FALSE,
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- fit:

  A `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) with
  `method = "scm"`.

- mspe_threshold:

  Minimum pre-treatment MSPE for including a control unit in the
  two-sided test. Ignored for one-sided tests. Default: 0 (no
  filtering).

- max_iter:

  Passed to
  [`scm_placebo_cpp()`](https://yo5uke.com/coresynth/reference/scm_placebo_cpp.md).
  Default 100L.

- tol:

  Passed to
  [`scm_placebo_cpp()`](https://yo5uke.com/coresynth/reference/scm_placebo_cpp.md).
  Default 1e-4.

- use_covariates:

  If `TRUE` and the fit used predictor covariates, applies the same
  covariate spec to each placebo unit (R-level loop). Default `FALSE`
  (faster C++ outcomes-only placebos).

- alternative:

  Direction of the alternative hypothesis: `"two.sided"` (default) uses
  the MSPE ratio statistic; `"greater"` tests whether the treatment
  increased the outcome; `"less"` tests whether the treatment decreased
  the outcome. One-sided tests use the signed ATT as the test statistic.

## Value

A list with:

- `p_value`: Permutation p-value between 0 and 1

- `mspe_ratio_treated`: MSPE_post / MSPE_pre for the treated unit
  (two.sided only)

- `mspe_ratios_all`: Named numeric vector (treated first, then
  controls); two.sided only

- `placebo_effects`: Named N_co-vector of placebo ATT estimates

- `treated_effect`: ATT estimate for the treated unit

- `n_placebo_used`: Number of control units used

## Details

When `alternative = "two.sided"` (default), the test statistic is the
post/pre MSPE ratio, following Abadie et al. (2010). When
`alternative = "greater"` or `"less"`, the test statistic is the signed
average post-treatment gap (ATT), giving a one-sided permutation test as
recommended by Abadie (2021) S.3.5 for improved power when the direction
of the treatment effect is known.
