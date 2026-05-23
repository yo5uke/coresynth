#' @keywords internal
#' @useDynLib coresynth, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases median optim pnorm qnorm quantile reorder rnorm sd setNames var
#' @importFrom utils combn packageVersion
"_PACKAGE"

utils::globalVariables(c(
  "time", "value", "series", "weight",
  "group", "tau_hat", "ci_lower", "ci_upper",
  "sd"
))
