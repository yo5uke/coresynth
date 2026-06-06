# Plot an scm_design object

Plot an scm_design object

## Usage

``` r
# S3 method for class 'scm_design'
plot(x, type = c("outcome", "gap"), ...)
```

## Arguments

- x:

  An `scm_design` object.

- type:

  `"outcome"` (default): synthetic treated vs synthetic control outcome
  series over all periods. `"gap"`: estimated treatment effect in the
  experimental periods, with split-conformal confidence intervals.

- ...:

  Currently ignored.

## Value

A `ggplot` object: for `type = "outcome"`, the synthetic treated and
synthetic control outcome series; for `type = "gap"`, the estimated
treatment effect over the experimental periods with split-conformal
confidence intervals. The object is returned for printing or further
customisation.
