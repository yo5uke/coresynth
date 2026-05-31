# Inference for Synthetic Difference-in-Differences

Computes standard errors and p-values for a SDID estimate using one of
three methods: permutation placebo test (Algorithm 4), cluster bootstrap
(Algorithm 2), or leave-one-out jackknife (Algorithm 3), following
Clarke et al. (2023).

## Usage

``` r
sdid_inference(
  fit,
  method = c("placebo", "bootstrap", "jackknife", "jackknife_global"),
  n_boot = 200L,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  seed = NULL
)
```

## Arguments

- fit:

  A `coresynth` object with `method = "sdid"` (sharp adoption only).

- method:

  Inference method: `"placebo"` (permutation), `"bootstrap"`, or
  `"jackknife"`.

- n_boot:

  Number of bootstrap replications (only for `method = "bootstrap"`).

- level:

  Confidence level for the interval (only for `method = "bootstrap"` or
  `"jackknife"`).

- alternative:

  Direction of the alternative hypothesis: `"two.sided"`, `"greater"`,
  or `"less"`.

- seed:

  Integer seed for reproducibility (only for `method = "bootstrap"`).

## Value

A list with:

- `estimate`: The SDID point estimate.

- `se`: Standard error (bootstrap / jackknife only).

- `p_value`: Permutation or normal-approximation p-value.

- `ci_lower`, `ci_upper`: Confidence interval bounds (bootstrap /
  jackknife).

- `method`: The inference method used.

- `n_controls`: Number of control units.

- `alternative`: The alternative hypothesis direction.

- `placebo_effects`: Named vector of LOO placebo effects (placebo only).

- `boot_ests`: Bootstrap estimate distribution (bootstrap only).
