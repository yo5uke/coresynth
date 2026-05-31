# Fast Leave-One-Out Placebo Test for SCM (Abadie et al. 2010)

For each control unit, treats it as pseudo-treated and fits SCM weights
from the remaining N_co-1 donors. Returns MSPE components for
constructing MSPE-ratio permutation p-values in R.

## Usage

``` r
scm_placebo_cpp(Y_pre, Y_post, max_iter = 100L, tol = 1e-04)
```

## Arguments

- Y_pre:

  Control pre-treatment outcomes (T_pre x N_co)

- Y_post:

  Control post-treatment outcomes (T_post x N_co)

- max_iter:

  Outer coordinate-descent iterations (default 100)

- tol:

  Convergence tolerance for V updates (default 1e-4)

## Value

A list with:

- `mspe_pre`: N_co-vector of pre-treatment MSPE per placebo unit

- `mspe_post`: N_co-vector of post-treatment MSPE per placebo unit

- `effects`: N_co-vector of mean post-period gap per placebo unit
