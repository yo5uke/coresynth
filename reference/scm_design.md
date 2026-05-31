# Experimental Synthetic Control Design

Selects which units to assign to the treatment arm (and which to the
control arm) in a planned experiment, following Abadie and Zhao (2026).
Both sets of units are chosen by minimising the distance between their
weighted-average predictor vectors and the population-average predictor
vector \\\bar{X}\\, so the resulting estimates are less susceptible to
post-randomisation bias than pure random assignment.

## Usage

``` r
scm_design(
  data,
  outcome,
  unit,
  time,
  T0,
  T_fit = NULL,
  m_min = 1L,
  m_max = 1L,
  f = NULL,
  predictors = NULL,
  design = c("base", "weakly_targeted", "unit_level"),
  beta = 1,
  xi = 1,
  alpha = 0.05,
  normalize = TRUE,
  max_subsets = 100000L
)
```

## Arguments

- data:

  Long-format data frame (one row per unit–time).

- outcome:

  Name of the outcome column.

- unit:

  Name of the unit identifier column.

- time:

  Name of the time identifier column.

- T0:

  Last pre-experimental period (a value present in the time column).
  Periods after `T0` are the experimental periods.

- T_fit:

  Number of fitting periods, counted from the **start** of the
  pre-experimental phase. Defaults to `NULL`, which uses all
  pre-experimental periods for fitting (no blank periods; inference
  disabled). When `T_fit` is smaller than the total number of
  pre-experimental periods, the remaining periods become blank periods
  used for inference.

- m_min:

  Minimum number of units assigned to treatment (default 1).

- m_max:

  Maximum number of units assigned to treatment (default 1).

- f:

  Named numeric vector of population weights \\f_j\\. Defaults to
  uniform weights \\1/J\\. Will be normalised to sum to 1.

- predictors:

  A [`list()`](https://rdrr.io/r/base/list.html) of
  [`pred()`](https://yo5uke.com/coresynth/reference/pred.md)
  specifications that define the predictor matrix \\X_j\\. Defaults to
  `NULL`, which uses all fitting-period outcome values as predictors.

- design:

  Design formulation: `"base"` (default), `"weakly_targeted"`, or
  `"unit_level"`.

- beta:

  Trade-off parameter \\\beta \> 0\\ for the Weakly targeted design
  (default 1).

- xi:

  Trade-off parameter \\\xi \> 0\\ for the Unit-level design (default
  1).

- alpha:

  Significance level for confidence intervals (default 0.05).

- normalize:

  If `TRUE` (default), each row of the predictor matrix is divided by
  its cross-unit standard deviation before optimisation, so predictors
  measured on different scales contribute equally.

- max_subsets:

  Maximum number of treatment-set candidates to evaluate before
  switching to random sampling (default 100 000).

## Value

An object of class `"scm_design"` with components:

- `treated_units`: unit identifiers selected for treatment

- `control_units`: unit identifiers in the control pool

- `w`: J-length weight vector for the synthetic treated unit (sums to 1)

- `v`: J-length weight vector for the synthetic control unit (sums to 1)

- `tau_hat`: estimated treatment effects for each experimental period

- `p_value`: permutation p-value (NA when blank periods are unavailable)

- `ci_lower`, `ci_upper`: per-period split-conformal confidence interval

- `Y_synth_tr`, `Y_synth_co`: synthetic treated/control series (all
  periods)

- `estimate`: ATT (mean of `tau_hat`)

## Details

Three design formulations are available:

- **`"base"`** (eq. 7): both the synthetic treated and the synthetic
  control independently target the population average \\\bar{X}\\.

- **`"weakly_targeted"`** (eq. 9): the synthetic treated targets
  \\\bar{X}\\; the synthetic control targets the synthetic treated
  predictor vector (controlled by `beta`).

- **`"unit_level"`** (eq. 10): each treated unit gets its own synthetic
  control; the aggregate control weight is a convex combination
  (controlled by `xi`).

Inference uses "blank periods" — pre-experimental periods whose outcomes
were *not* used to estimate the weights. Set `T_fit` strictly smaller
than the number of pre-experimental periods to enable the permutation
test and split- conformal confidence intervals from Section 3 of Abadie
and Zhao (2026).

## References

Abadie, A. and Zhao, J. (2026). "Synthetic Controls for Experimental
Design." MIT Working Paper.
