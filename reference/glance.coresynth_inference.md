# Glance at an inference result

One-row summary of a `coresynth_inference` (or `sdid_inference`) object.

## Usage

``` r
# S3 method for class 'coresynth_inference'
glance(x, ...)
```

## Arguments

- x:

  An inference object.

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `method`, `n_controls`, `staggered`,
`estimate`, `std.error`, `p.value`, `conf.low`, `conf.high`,
`alternative`, `n_boot_valid`.
