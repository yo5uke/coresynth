# Fast Placebo Test for SDID

For each control unit, treats it as the "pseudo-treated" unit and
estimates the leave-one-out SDID effect. The distribution of these
placebo effects provides a permutation-based null distribution for
inference.

## Usage

``` r
sdid_placebo_cpp(Y_pre, Y_post, time_weights, zeta2)
```

## Arguments

- Y_pre:

  Control units pre-treatment outcomes (T_pre x N_co)

- Y_post:

  Control units post-treatment outcomes (T_post x N_co)

- time_weights:

  Lambda weights for pre-treatment periods (T_pre x 1)

- zeta2:

  Ridge penalty (same as used in the main estimate)

## Value

A numeric vector of length `N_co`. Each element is the leave-one-out
placebo SDID effect obtained by treating that control unit as the
pseudo-treated unit; the vector serves as a permutation-based null
distribution for inference.
