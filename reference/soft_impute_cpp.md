# Fast Matrix Completion using Soft-Impute Algorithm

Solves: min_L (1/2) \|\|O o (Y - L)\|\|*F^2 + lambda \* \|\|L\|\|*\* via
iterative SVD soft-thresholding (Mazumder, Hastie, Tibshirani 2010).
Note: lambda is NOT normalized by \|O\|. Default lambda = 0.01 \*
sigma_max(Y).

## Usage

``` r
soft_impute_cpp(Y, O, lambda, max_iter = 1000L, tol = 1e-05)
```

## Arguments

- Y:

  Observed outcome matrix (N x T). Unobserved entries should be 0.

- O:

  Binary mask matrix (N x T): 1 = observed, 0 = missing (treated post).

- lambda:

  Nuclear norm penalty (soft-threshold on singular values).

- max_iter:

  Maximum iterations.

- tol:

  Convergence tolerance (relative Frobenius norm change).

## Value

A numeric matrix of the same dimension as `Y` (N x T): the completed
low-rank matrix `L` that minimises the soft-thresholded nuclear-norm
objective.
