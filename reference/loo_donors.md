# Leave-One-Out Donor Robustness for SCM

Iteratively re-estimates the synthetic control excluding one
contributing donor at a time, holding the predictor weights V fixed at
their baseline values (Abadie, Diamond & Hainmueller 2015, footnote 20).
The spread of the leave-one-out ATT estimates shows how much the result
hinges on any single donor.

## Usage

``` r
loo_donors(fit, weight_threshold = 1e-06)
```

## Arguments

- fit:

  A sharp `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) with
  `method = "scm"`.

- weight_threshold:

  Only donors whose baseline weight exceeds this value are dropped
  (removing a zero-weight donor cannot change the fit). Default `1e-6`.

## Value

A list with:

- `att_original`: baseline ATT

- `results`: data.frame with one row per excluded donor (`donor`,
  `weight`, `att_loo`)

- `att_range`: range of the leave-one-out ATTs

## Details

For penalised fits (`lambda_pen` used), the same penalty is re-applied
in each leave-one-out QP.

## See also

[`placebo_in_time()`](https://yo5uke.com/coresynth/reference/placebo_in_time.md),
[`mspe_ratio_pval()`](https://yo5uke.com/coresynth/reference/mspe_ratio_pval.md)
