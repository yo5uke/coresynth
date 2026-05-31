# Conformal Inference for Synthetic Control Estimators

Implements the permutation-based conformal inference procedure of
Chernozhukov, Wuthrich & Zhu (2021, JASA). The test inverts a sharp null
\\H_0: \tau = \tau_0\\ by imputing the treated post-treatment
counterfactual as \\Y\_{1t} - \tau_0\\, re-estimating the counterfactual
proxy on **all** \\T\\ periods (imposing the null), and computing a
moving-block permutation p-value from the estimated residuals. A
confidence interval is obtained by test inversion over a grid of
candidate \\\tau_0\\.

## Usage

``` r
conformal_inference(
  fit,
  tau0 = 0,
  q = 1,
  alternative = c("two.sided", "greater", "less"),
  ci = TRUE,
  level = 0.95,
  grid = NULL,
  n_grid = 200L,
  grid_mult = 4,
  ...
)
```

## Arguments

- fit:

  A `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md).

- tau0:

  Null value of the ATT for the reported p-value (default 0).

- q:

  Exponent of the \\S_q\\ test statistic
  (`S_q = (T_post^{-1} \sum |u_t|^q)^{1/q}`). Default 1, robust to
  heavy-tailed data (CWZ 2021). Used only for
  `alternative = "two.sided"`; one-sided tests use the signed mean
  post-treatment residual.

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`.

- ci:

  Logical; construct a confidence interval by test inversion (default
  `TRUE`).

- level:

  Confidence level for the interval (default 0.95).

- grid:

  Optional numeric vector of candidate \\\tau_0\\ values for test
  inversion. When `NULL` (default), a grid of `n_grid` points is centred
  on the point estimate with half-width `grid_mult` times the
  pre-treatment residual standard deviation.

- n_grid:

  Number of grid points when `grid = NULL` (default 200).

- grid_mult:

  Half-width multiplier when `grid = NULL` (default 4).

- ...:

  Unused.

## Value

A list of class `c("conformal_inference", "coresynth_inference")` with
`estimate`, `se` (`NA`; conformal has no SE), `p_value` (at `tau0`),
`ci_lower`, `ci_upper`, `method` (`"conformal"`), `n_controls`,
`alternative`, `staggered` (`FALSE`), plus `tau0`, `q`, `grid`, and
`p_grid` (p-values along the grid). Compatible with
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
[`glance()`](https://generics.r-lib.org/reference/glance.html).

## Details

Supported for **sharp** (single-cohort) fits with `method` in
`c("scm", "sdid", "gsc", "mc", "si")`. Staggered, multi-arm, and `tasc`
fits are not supported (use
[`sdid_inference()`](https://yo5uke.com/coresynth/reference/sdid_inference.md),
[`gsc_inference()`](https://yo5uke.com/coresynth/reference/gsc_inference.md),
or
[`si_inference()`](https://yo5uke.com/coresynth/reference/si_inference.md)
instead).

## References

Chernozhukov, V., Wuthrich, K., & Zhu, Y. (2021). An Exact and Robust
Conformal Inference Method for Counterfactual and Synthetic Controls.
*Journal of the American Statistical Association*, 116(536), 1849-1864.
