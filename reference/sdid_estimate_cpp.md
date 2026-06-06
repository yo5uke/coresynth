# Calculate SDID Estimate (tau_sdid)

Given unit weights omega and time weights lambda, computes the SDID
estimator as a weighted two-way difference:

## Usage

``` r
sdid_estimate_cpp(Y_pre_co, Y_post_co, Y_pre_tr, Y_post_tr, omega, lambda)
```

## Arguments

- Y_pre_co:

  Control pre-treatment outcomes (T_pre x N_co)

- Y_post_co:

  Control post-treatment outcomes (T_post x N_co)

- Y_pre_tr:

  Treated pre-treatment outcomes (T_pre x 1)

- Y_post_tr:

  Treated post-treatment outcomes (T_post x 1)

- omega:

  Unit weights (N_co x 1)

- lambda:

  Time weights (T_pre x 1)

## Value

A single numeric value: the SDID treatment-effect estimate `tau_sdid`.

## Details

tau_sdid = (Y_tr_post_mean - Y_tr_pre_wt) - (Y_co_post_wt - Y_co_pre_wt)
