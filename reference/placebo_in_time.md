# In-Time Placebo (Backdating) Test for SCM

Re-estimates the synthetic control after artificially backdating the
treatment to a pre-treatment period, following Abadie, Diamond &
Hainmueller (2015) and Abadie & Vives-i-Bastida (2022, principle 7:
"out-of-sample validation is key"). Only pre-treatment data enter the
exercise, so the placebo gap after `t0_placebo` is uncontaminated by the
actual intervention. A credible design shows no sizable divergence at
the backdated treatment time.

## Usage

``` r
placebo_in_time(fit, t0_placebo = NULL)
```

## Arguments

- fit:

  A sharp `coresynth` object from
  [`scm_fit()`](https://yo5uke.com/coresynth/reference/scm_fit.md) with
  `method = "scm"`.

- t0_placebo:

  Backdated treatment period as a 1-based position in `fit$times`; must
  satisfy `2 <= t0_placebo < T_pre`. Default `floor(T_pre / 2)`.

## Value

A list with:

- `t0_placebo`: the backdated treatment period used

- `times`: time values of the pre-treatment window

- `unit_weights`: placebo donor weights

- `Y_treat`, `Y_synth`, `gap`: series over the pre-treatment window

- `placebo_att`: mean placebo gap over `(t0_placebo, T_pre]`

- `fit_rmspe`: RMSPE over the placebo fitting window `1..t0_placebo`

- `eval_rmspe`: RMSPE over the placebo post window `(t0_placebo, T_pre]`

## Details

The refit uses the outcomes of periods `1..t0_placebo` as predictors
(the `predictors = NULL` default), regardless of how the original fit
was specified, because user-supplied
[`pred()`](https://yo5uke.com/coresynth/reference/pred.md) windows
cannot be lagged automatically (ADH 2015 lag their predictors by hand).

## See also

[`mspe_ratio_pval()`](https://yo5uke.com/coresynth/reference/mspe_ratio_pval.md)
for in-space placebos,
[`loo_donors()`](https://yo5uke.com/coresynth/reference/loo_donors.md)
for donor-robustness checks.
