# Calculate SDID Unit Weights (omega)

Solves the regularized QP: min over omega in Delta: sum_t (sum_i omega_i
Y_it - Y_tr_t)^2 + zeta^2 \* T_pre \* \|\|omega\|\|^2

## Usage

``` r
sdid_unit_weights_cpp(Y_pre, Y_tr_pre, zeta2)
```

## Arguments

- Y_pre:

  Pre-treatment outcome matrix for control units (T_pre x N_co)

- Y_tr_pre:

  Pre-treatment outcome vector for treated unit (T_pre x 1), averaged if
  multiple

- zeta2:

  Ridge penalty parameter (zeta^2). The code internally multiplies by
  T_pre per the paper.

## Details

This corresponds to equation (5) in Arkhangelsky et al. (2021).
