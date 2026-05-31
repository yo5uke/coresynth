# Predictor Specification for SCM

Creates a single predictor specification for use in
[`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) with
`method = "scm"`. Pass a [`list()`](https://rdrr.io/r/base/list.html) of
`pred()` calls as the `predictors` argument to define the full covariate
matrix.

## Usage

``` r
pred(vars, times, op = "mean")
```

## Arguments

- vars:

  Character vector of variable names. All variables share the same
  `times` window and `op` operator. Use separate `pred()` calls for
  variables with different time windows.

- times:

  Numeric/integer vector of time values to aggregate over.

- op:

  Aggregation operator applied to each variable over `times`. One of
  `"mean"` (default), `"median"`, or `"sum"`.

## Value

A `pred_spec` object (a named list with class `"pred_spec"`).

## See also

[`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) for the
`predictors` argument that consumes a
[`list()`](https://rdrr.io/r/base/list.html) of `pred_spec` objects.

## Examples

``` r
# Three variables averaged over the same window
pred(c("lnincome", "retprice", "age15to24"), 1980:1988)
#> pred(lnincome, retprice, age15to24, 1980:1988, op = "mean")

# Single variable at a specific year
pred("cigsale", 1975)
#> pred(cigsale, 1975, op = "mean")

# Single variable averaged over a range
pred("beer", 1984:1988)
#> pred(beer, 1984:1988, op = "mean")

# Abadie, Diamond & Hainmueller (2010) California Prop 99 style:
# combine several covariates aggregated over different windows plus
# three outcome lags at specific years. Pass the list to scm_fit(..., predictors = ...).
if (FALSE) { # \dontrun{
predictors <- list(
  pred(c("lnincome", "retprice", "age15to24"), 1980:1988),
  pred("beer",    1984:1988),
  pred("cigsale", 1988),
  pred("cigsale", 1980),
  pred("cigsale", 1975)
)
fit <- scm_fit(cigsale ~ treated | state + year,
               data = prop99, method = "scm",
               predictors = predictors)
} # }
```
