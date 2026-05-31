# SI-PCR: Synthetic Interventions via Principal Component Regression

Implements the SI-PCR estimator of Agarwal et al. (2025). Uses the top-k
SVD of pre-treatment control outcomes to find donor weights that predict
each treated unit's pre-treatment trajectory, then applies those weights
to post-treatment control outcomes.

## Usage

``` r
si_pcr_cpp(Y_pre_co, Y_post_co, Y_pre_tr, k)
```

## Arguments

- Y_pre_co:

  Pre-treatment control outcomes (T_pre x N_co)

- Y_post_co:

  Post-treatment control outcomes (T_post x N_co)

- Y_pre_tr:

  Pre-treatment treated outcomes (T_pre x N_tr)

- k:

  Number of SVD components to retain

## Value

A list with:

- `W`: Donor weight matrix (N_co x N_tr)

- `Y_hat`: Counterfactual post-treatment outcomes (T_post x N_tr)
