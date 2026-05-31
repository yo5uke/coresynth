# Tidy an inference result

Coerces a `coresynth_inference` or `sdid_inference` object to a one-row
tidy `data.frame` with broom-style column names so it can be combined
with regression output for paper tables.

## Usage

``` r
# S3 method for class 'coresynth_inference'
tidy(x, conf.int = TRUE, ...)
```

## Arguments

- x:

  A `coresynth_inference` (or `sdid_inference`) object returned by
  [`sdid_inference()`](https://yo5uke.com/coresynth/reference/sdid_inference.md),
  [`gsc_inference()`](https://yo5uke.com/coresynth/reference/gsc_inference.md),
  or
  [`si_inference()`](https://yo5uke.com/coresynth/reference/si_inference.md).

- conf.int:

  Logical. Include `conf.low`/`conf.high` columns when CI is available
  (default `TRUE`). Permutation placebo SE/CI are `NA`.

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `term`, `estimate`, `std.error`,
`statistic`, `p.value`, `conf.low`, `conf.high`, `method`,
`alternative`, `n_controls`, `staggered`.
