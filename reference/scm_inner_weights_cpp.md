# SCM Inner Weights (QP Given V)

Solves the inner-loop QP for SCM: given a fixed diagonal metric matrix
V, finds donor weights W on the simplex minimising the V-weighted
covariate loss.

## Usage

``` r
scm_inner_weights_cpp(X0, X1, V_diag)
```

## Arguments

- X0:

  Covariate matrix for control units (k x N_co)

- X1:

  Covariate vector for the treated unit (k x 1)

- V_diag:

  Diagonal of the metric matrix V (k x 1, non-negative, need not sum to
  1)

## Value

Donor weight vector W (N_co x 1) on the unit simplex
