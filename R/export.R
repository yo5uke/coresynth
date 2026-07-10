#' Export coresynth Results to JSON
#'
#' Generates a comprehensive, standardized JSON record covering all six
#' estimators. Suitable for reproducibility workflows (Xu & Yang 2026) and
#' downstream tooling. Pass the result of [mspe_ratio_pval()] or [gsc_boot()]
#' via the `inference` argument to include inference results.
#'
#' @param x A `coresynth` object from [scm_fit()].
#' @param file Output file path. Default `"coresynth_results.json"`. Pass
#'   `NULL` to skip writing and return the R list invisibly.
#' @param inference Optional list from [mspe_ratio_pval()] or [gsc_boot()].
#'   When provided, populates the `inference` section and updates `estimate`
#'   with `p_value`, `se`, `ci_lower`, `ci_upper`.
#' @param digits Number of significant digits applied to numeric values
#'   (default 6L).
#' @return Invisibly, the R list that was (or would be) serialized.
#' @export
export_json <- function(x, file = "coresynth_results.json",
                        inference = NULL, digits = 6L) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for JSON export. Please install it.")
  if (!inherits(x, "coresynth"))
    stop("x must be a coresynth object.")

  rd <- function(v) if (is.numeric(v)) signif(v, digits) else v

  # ── 1. meta ──────────────────────────────────────────────────────────────────
  meta <- list(
    package    = "coresynth",
    version    = as.character(packageVersion("coresynth")),
    method     = x$method,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  # ── 2. data ───────────────────────────────────────────────────────────────────
  T_total <- length(x$times)
  T_pre   <- x$T_pre
  T_post  <- T_total - T_pre

  unit_names_co <- if (!is.null(x$unit_weights)) names(x$unit_weights)
                   else if (!is.null(x$Y_co_all)) colnames(x$Y_co_all)
                   else NULL
  n_controls <- if (!is.null(x$unit_weights)) length(x$unit_weights)
                else if (!is.null(x$L_co))    nrow(x$L_co)
                else NULL
  n_treated <- if (!is.null(x$L_tr))      nrow(x$L_tr)
               else if (is.matrix(x$Y_treat)) ncol(x$Y_treat)
               else 1L

  data_section <- list(
    n_controls = n_controls,
    n_treated  = n_treated,
    T_total    = T_total,
    T_pre      = T_pre,
    T_post     = T_post,
    times      = as.numeric(x$times),
    unit_names = unit_names_co
  )

  # ── 3. estimate ───────────────────────────────────────────────────────────────
  estimate_section <- list(
    att      = rd(x$estimate),
    se       = if (!is.null(x$se))       rd(x$se)       else NULL,
    p_value  = if (!is.null(x$p_value))  rd(x$p_value)  else NULL,
    ci_lower = if (!is.null(x$ci_lower)) rd(x$ci_lower) else NULL,
    ci_upper = if (!is.null(x$ci_upper)) rd(x$ci_upper) else NULL
  )
  if (!is.null(inference)) {
    if (!is.null(inference$p_value))  estimate_section$p_value  <- rd(inference$p_value)
    if (!is.null(inference$se))       estimate_section$se       <- rd(inference$se)
    if (!is.null(inference$ci_lower)) estimate_section$ci_lower <- rd(inference$ci_lower)
    if (!is.null(inference$ci_upper)) estimate_section$ci_upper <- rd(inference$ci_upper)
  }

  # ── 4. weights ────────────────────────────────────────────────────────────────
  weights_section <- list()
  if (!is.null(x$unit_weights))
    weights_section$unit <- rd(as.list(x$unit_weights))
  if (!is.null(x$time_weights))
    weights_section$time <- rd(as.numeric(x$time_weights))
  if (!is.null(x$v_weights))
    weights_section$covariate <- rd(as.list(x$v_weights))

  # ── 5. time_series ────────────────────────────────────────────────────────────
  flatten_ts <- function(m) {
    if (is.matrix(m)) {
      if (ncol(m) == 1L) return(rd(drop(m)))
      return(lapply(seq_len(ncol(m)), function(j) rd(m[, j])))
    }
    rd(as.numeric(m))
  }
  Y_synth_out <- if (!is.null(x$Y_synth))   x$Y_synth
                 else if (!is.null(x$Y_tr_hat)) x$Y_tr_hat
                 else if (!is.null(x$Y_hat)) {
                   # TASC's Y_hat spans all N units; keep the treated columns
                   if (identical(x$method, "tasc") && !is.null(x$idx_tr))
                     x$Y_hat[, x$idx_tr, drop = FALSE]
                   else x$Y_hat
                 }
                 else NULL

  ts_section <- list(
    times   = as.numeric(x$times),
    Y_treat = flatten_ts(x$Y_treat),
    Y_synth = if (!is.null(Y_synth_out)) flatten_ts(Y_synth_out) else NULL,
    gap     = flatten_ts(x$gap)
  )

  # ── 6. method_specific ───────────────────────────────────────────────────────
  ms_section <- switch(x$method,
    "scm"  = list(
      loss      = rd(x$loss),
      v_weights = rd(as.list(x$v_weights))
    ),
    "sdid" = list(
      zeta2      = rd(x$zeta2),
      sigma2_hat = rd(x$sigma2_hat),
      omega0     = rd(x$omega0),
      lambda0    = rd(x$lambda0)
    ),
    "gsc"  = list(
      r               = x$r,
      singular_values = rd(as.numeric(x$singular_values)),
      F               = lapply(seq_len(ncol(x$F)),   function(j) rd(x$F[, j])),
      L_co            = lapply(seq_len(ncol(x$L_co)), function(j) rd(x$L_co[, j])),
      L_tr            = lapply(seq_len(ncol(x$L_tr)), function(j) rd(x$L_tr[, j]))
    ),
    "mc"   = list(lambda = rd(x$lambda)),
    "tasc" = list(
      r = x$r,
      A = lapply(seq_len(nrow(x$A)), function(i) rd(x$A[i, ]))
    ),
    "si"   = list(
      k       = x$k,
      weights = rd(as.numeric(x$unit_weights))
    ),
    list()
  )

  # ── 7. inference ─────────────────────────────────────────────────────────────
  inf_section <- list()
  if (!is.null(inference)) {
    if (!is.null(inference$mspe_ratios_all)) {
      inf_section <- list(
        type               = "mspe_permutation",
        p_value            = rd(inference$p_value),
        mspe_ratio_treated = rd(inference$mspe_ratio_treated),
        mspe_ratios_all    = rd(as.numeric(inference$mspe_ratios_all)),
        n_placebo_used     = inference$n_placebo_used,
        placebo_effects    = rd(as.numeric(inference$placebo_effects))
      )
    } else if (!is.null(inference$boot_dist)) {
      inf_section <- list(
        type      = "parametric_bootstrap",
        B         = length(inference$boot_dist),
        p_value   = rd(inference$p_value),
        se        = rd(inference$se),
        ci_lower  = rd(inference$ci_lower),
        ci_upper  = rd(inference$ci_upper),
        boot_dist = rd(inference$boot_dist)
      )
    }
  }

  # ── Assemble ──────────────────────────────────────────────────────────────────
  result <- list(
    meta            = meta,
    data            = data_section,
    estimate        = estimate_section,
    weights         = if (length(weights_section) > 0L) weights_section else NULL,
    time_series     = ts_section,
    method_specific = if (length(ms_section) > 0L) ms_section else NULL,
    inference       = if (length(inf_section) > 0L) inf_section else NULL
  )

  if (!is.null(file))
    jsonlite::write_json(result, path = file, auto_unbox = TRUE,
                         pretty = TRUE, null = "null")
  invisible(result)
}
