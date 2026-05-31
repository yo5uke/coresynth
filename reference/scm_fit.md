# Fit a Synthetic Control Method Model

Unified formula interface for Synthetic Control and related causal
inference methods. The formula syntax is:

## Usage

``` r
scm_fit(
  formula,
  data,
  method = c("scm", "sdid", "gsc", "mc", "tasc", "si"),
  predictors = NULL,
  covariates = NULL,
  v_selection = c("insample", "oos"),
  donor_mspe_threshold = Inf,
  lambda_pen = NULL,
  v_optim = c("coord_descent", "auto", "bfgs"),
  ...
)
```

## Arguments

- formula:

  A `Formula` object, e.g. `y ~ D | unit + time`.

- data:

  A `data.frame` in **long** format (one row per unit-time).

- method:

  One of `"scm"`, `"sdid"`, `"gsc"`, `"mc"`, `"tasc"`, `"si"`.

- predictors:

  A [`list()`](https://rdrr.io/r/base/list.html) of
  [`pred()`](https://yo5uke.com/coresynth/reference/pred.md)
  specifications that define the predictor matrix for SCM (see Abadie et
  al. 2010, S.2.3). Each
  [`pred()`](https://yo5uke.com/coresynth/reference/pred.md) entry
  aggregates one or more variables over a time window. Pass `NULL`
  (default) to use all pre-treatment outcome periods as predictors.
  Applies to `method = "scm"` only.

- covariates:

  An optional named `list` of additional time-varying covariates to
  partial out before estimation. Each element is a character string
  naming a column in `data`. Supported for `method = "sdid"`, `"scm"`,
  and `"gsc"`.

- v_selection:

  V matrix selection method for `method = "scm"`. `"insample"` (default)
  follows Abadie et al. (2010): V is chosen by minimising in-sample
  pre-treatment MSPE. `"oos"` follows Abadie (2021) S.3.2: the
  pre-treatment window is split in half; V is selected to minimise MSPE
  on the validation half, then W is refit on the full window.

- donor_mspe_threshold:

  Donor pool filtering threshold (Abadie 2021 S.4). For `method = "scm"`
  only. Each donor's individual pre-treatment MSPE (using that donor
  alone as the counterfactual) is divided by the minimum such MSPE
  across all donors. Donors whose ratio exceeds this threshold are
  excluded from estimation. `Inf` (default) disables filtering.

- lambda_pen:

  Penalised SCM parameter (Abadie & L'Hour 2021, JASA). For
  `method = "scm"` only. `NULL` (default) runs standard unpenalised SCM.
  `"auto"` selects the penalty via out-of-sample pre-treatment MSPE on
  the same validation window as `v_selection = "oos"`. A non-negative
  number uses that value directly.

- v_optim:

  Outer V-optimisation method for `method = "scm"`. `"coord_descent"`
  (default) uses the existing C++ coordinate descent with 11-point grid
  search – fastest when `k = T_pre` is large (outcomes-only). `"bfgs"`
  uses R's L-BFGS-B, which requires only O(k^2) inner QP calls and is
  faster when `k` is small (e.g. external predictors with k \<= 15).
  `"auto"` selects `"bfgs"` when `k <= 15`, otherwise `"coord_descent"`.

- ...:

  Additional arguments forwarded to the specific method (e.g. `r`,
  `lambda`, `zeta2`).

## Value

An object of classes `c("coresynth_<method>", "coresynth")`. All methods
return at minimum:

- `method`: estimator name

- `estimate`: average treatment effect (ATT)

- `times`: time index vector

- `T_pre`: number of pre-treatment periods

- `Y_treat`: treated unit outcome series

- `gap`: treatment effect series (Y_treat - counterfactual)

## Details

`outcome ~ treatment | unit_id + time_id`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- scm_fit(gdp ~ treated | country + year, data = panel_data, method = "sdid")
summary(fit)
plot(fit, type = "gap")
} # }
```
