#' @section Parallel computation:
#' The in-space placebo refits of [mspe_ratio_pval()] and [sdid_inference()],
#' and the multi-start search of [scm_fit()], run in parallel with OpenMP.
#' Set `options(coresynth.threads = n)` to cap the number of threads. Left
#' unset, the OpenMP default applies (all cores, or `OMP_NUM_THREADS` /
#' `OMP_THREAD_LIMIT` when set), except under `R CMD check --as-cran`, where at
#' most two threads are used. The results do not depend on the thread count.
#'
#' @keywords internal
#' @useDynLib coresynth, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases median optim pnorm qnorm quantile reorder rnorm sd setNames var
#' @importFrom utils combn packageVersion
"_PACKAGE"

utils::globalVariables(c(
  "time", "value", "series", "weight",
  "group", "tau_hat", "ci_lower", "ci_upper",
  "sd", "gap", "ratio", "label", "panel", "predictor",
  "xv", "yv", "nu", "estimate", "conf.low", "conf.high"
))
