# Plot a coresynth model

Plot a coresynth model

## Usage

``` r
# S3 method for class 'coresynth'
plot(x, type = c("trend", "gap", "weights"), ...)
```

## Arguments

- x:

  A `coresynth` object.

- type:

  One of `"trend"` (observed vs synthetic), `"gap"` (ATT over time), or
  `"weights"` (donor unit weight bar chart).

- ...:

  Ignored.

## Value

A `ggplot2` plot object.
