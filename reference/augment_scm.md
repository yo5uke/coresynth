# Augmented Synthetic Control Method (Ridge ASCM)

Applies a ridge-regression-based bias correction to a fitted SCM object,
following Ben-Michael, Feller & Rothstein (2021, JASA). The corrected
estimator is:

## Usage

``` r
augment_scm(fit, lambda_ridge = NULL)
```

## Arguments

- fit:

  A `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) with
  `method = "scm"`.

- lambda_ridge:

  Ridge penalty (non-negative). `NULL` (default) selects the penalty by
  leave-one-out cross-validation on the control units.

## Value

A list with:

- `att_aug`: Augmented ATT estimate

- `delta`: Bias correction term (m_tr_post - sum_j W_j m_j_post)

- `att_scm`: Original SCM ATT for comparison

- `lambda_ridge`: Ridge penalty used

- `beta_hat`: Ridge regression coefficients (length T_pre)

## Details

tau_aug = tau_SCM + (m_tr_post - sum_j W_j \* m_j_post)

where m_i_post = Y_pre_i' beta_hat is the ridge outcome model prediction
for unit i's mean post-treatment outcome, and beta_hat is estimated by
ridge regression across control units.
