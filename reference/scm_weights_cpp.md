# SCM Outer Weights (Joint Optimization of W and V)

Jointly optimises donor weights W (on the simplex) and the diagonal
metric matrix V via coordinate descent on the pre-treatment prediction
MSPE, following Abadie, Diamond & Hainmueller (2010).

## Usage

``` r
scm_weights_cpp(X0, X1, Z0, Z1, max_iter = 100L, tol = 1e-04, t_train = -1L)
```

## Arguments

- X0:

  Covariate matrix for control units (k x N_co, typically pre-treatment
  outcomes)

- X1:

  Covariate vector for the treated unit (k x 1)

- Z0:

  Outcome matrix for control units in the pre-treatment window (T_pre x
  N_co)

- Z1:

  Outcome vector for the treated unit in the pre-treatment window (T_pre
  x 1)

- max_iter:

  Maximum coordinate-descent iterations (default 100)

- tol:

  Convergence tolerance on MSPE improvement (default 1e-4)

- t_train:

  Training window length for out-of-sample V selection. -1 (default):
  in-sample V selection (original behaviour). Positive: use rows
  0..(t_train-1) of Z for fitting W, rows t_train..(T_pre-1) as the
  validation window for V selection, then refit W on full data.

## Value

A list with:

- `W`: Donor weight vector (N_co x 1) on the unit simplex

- `V`: Optimal metric diagonal (k x 1, normalised to sum to 1)

- `loss`: Final pre-treatment prediction loss (full pre-treatment
  window)

## Details

When `t_train > 0`, uses out-of-sample V selection per Abadie (2021)
§3.2: V is selected by minimising MSPE on a validation window (rows
t_train..T_pre-1 of Z), while W is fitted on the training window (rows
0..t_train-1 of X when X and Z have the same row count, i.e. the
outcomes-only case). After selecting V\*, W is refit on the full data.
