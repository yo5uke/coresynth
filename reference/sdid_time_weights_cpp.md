# Calculate SDID Time Weights (lambda)

Solves the time-weight QP (with implicit intercept lambda_0 concentrated
out):

## Usage

``` r
sdid_time_weights_cpp(Y_pre_co, Y_post_target, zeta_t)
```

## Arguments

- Y_pre_co:

  Pre-treatment outcomes for control units, row-demeaned (T_pre x N_co)

- Y_post_target:

  Post-treatment mean per control unit, demeaned (N_co x 1)

- zeta_t:

  Ridge penalty for time weights (paper: 1e-6 \* sigma_hat)

## Value

A numeric vector of length `T_pre` holding the SDID time weights
`lambda` (non-negative and summing to one).

## Details

min over lambda in Delta_pre: \|\|Y_post_target - Y_pre_co^T
lambda\|\|^2 + zeta_t^2 \* N_co \* \|\|lambda\|\|^2

The caller is responsible for pre-demeaning Y_pre_co (row-wise) and
Y_post_target (subtract the cross-unit mean) to concentrate out
lambda_0, as described in Arkhangelsky et al. (2021) Algorithm 1, Eq.
(2.3).
