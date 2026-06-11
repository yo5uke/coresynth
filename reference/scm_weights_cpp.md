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

  Validation-window split for V selection. -1 (default): V selected on
  the full Z window (in-sample). Positive: rows t_train..(T_pre-1) of Z
  form the validation window used to select V (W is fitted on the full X
  throughout); after selecting V\*, W is refit and the reported loss
  uses the full Z window.

## Value

A list with:

- `W`: Donor weight vector (N_co x 1) on the unit simplex

- `V`: Optimal metric diagonal (k x 1, normalised to sum to 1)

- `loss`: Final pre-treatment prediction loss (full pre-treatment
  window)

## Details

When `t_train > 0`, V is selected by minimising MSPE on a validation
window (rows t_train..T_pre-1 of Z) while W is fitted on the full
predictor matrix X. This is appropriate when X is a fixed predictor
matrix that contains no validation-period outcome information (the
user-supplied predictors case). For the outcomes-only case the proper
Abadie (2021) S.3.2 train/validation split is implemented in R
(`.scm_oos_outcomes()`): candidate W(V) are fitted on training-half
outcomes only, by passing the training rows as X and the validation rows
as Z to this function with `t_train = -1`.
