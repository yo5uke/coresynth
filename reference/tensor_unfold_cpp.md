# Tensor Unfolding (Matricization) for Synthetic Interventions

Tensor Unfolding (Matricization) for Synthetic Interventions

## Usage

``` r
tensor_unfold_cpp(T_cube, mode)
```

## Arguments

- T_cube:

  A 3D array (cube) of dimensions (n1, n2, n3)

- mode:

  The mode to unfold along (1, 2, or 3)

## Value

A numeric matrix: the mode-`mode` unfolding (matricization) of `T_cube`,
with dimensions `n1 x (n2 * n3)`, `n2 x (n1 * n3)`, or `n3 x (n1 * n2)`
for `mode` 1, 2, or 3 respectively.
