# Fast Interactive Fixed Effects (IFE) for Generalized Synthetic Control

Implements Xu (2017) IFE model with optional covariate adjustment. When
X_co has p \> 0 slices, runs an EM loop alternating between: E-step:
truncated SVD of Y_tilde = Y_co - X_co \* beta M-step: panel OLS to
update beta given current factors When X_co has 0 slices (default),
falls back to the plain 3-step estimator.

## Usage

``` r
gsc_ife_cpp(Y_co, Y_tr_pre, r, X_co, X_tr_pre, max_iter = 50L, tol = 1e-06)
```

## Arguments

- Y_co:

  Control units outcome matrix (T x N_co)

- Y_tr_pre:

  Treated units pre-treatment outcomes (T_pre x N_tr)

- r:

  Number of latent factors (must be \<= min(T, N_co))

- X_co:

  Time-varying covariate cube (T x N_co x p). Pass an empty cube (0
  slices) for the covariate-free estimator.

- X_tr_pre:

  Time-varying covariate cube for treated units in the pre-treatment
  window (T_pre x N_tr x p). Required for correct Step 2 loading
  estimation per Xu (2017): lambda_hat is estimated from Y_tr_pre -
  X_tr_pre \* beta (covariate- demeaned). Pass an empty cube (0 slices)
  to skip demeaning (backward-compatible, but biased when beta != 0).

- max_iter:

  Maximum EM iterations (default 50)

- tol:

  Convergence tolerance on relative beta change (default 1e-6)
