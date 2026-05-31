# Export coresynth Results to JSON

Generates a comprehensive, standardized JSON record covering all six
estimators. Suitable for reproducibility workflows (Xu & Yang 2026) and
downstream tooling. Pass the result of
[`mspe_ratio_pval()`](https://yo5uke.com/coresynth/reference/mspe_ratio_pval.md)
or [`gsc_boot()`](https://yo5uke.com/coresynth/reference/gsc_boot.md)
via the `inference` argument to include inference results.

## Usage

``` r
export_json(x, file = "coresynth_results.json", inference = NULL, digits = 6L)
```

## Arguments

- x:

  A `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md).

- file:

  Output file path. Default `"coresynth_results.json"`. Pass `NULL` to
  skip writing and return the R list invisibly.

- inference:

  Optional list from
  [`mspe_ratio_pval()`](https://yo5uke.com/coresynth/reference/mspe_ratio_pval.md)
  or [`gsc_boot()`](https://yo5uke.com/coresynth/reference/gsc_boot.md).
  When provided, populates the `inference` section and updates
  `estimate` with `p_value`, `se`, `ci_lower`, `ci_upper`.

- digits:

  Number of significant digits applied to numeric values (default 6L).

## Value

Invisibly, the R list that was (or would be) serialized.
