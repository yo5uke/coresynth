## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.6.0 (x86_64-w64-mingw32)
* win-builder: R-devel (via `devtools::check_win_devel()`)

## Submission notes

* This is a resubmission. Version 0.2.0 was submitted but has not yet been
  published. Version 0.2.1 fixes a data leak in `v_selection = "oos"`
  (Abadie 2021 S.3.2) and adds Synth-style predictor scaling
  (`scale_predictors`); both change numerical results for affected fits and
  are documented in NEWS.md.
* The package contains compiled C++ code (via Rcpp / RcppArmadillo) and uses
  OpenMP for parallelism where available; this is declared through
  `SHLIB_OPENMP_*FLAGS` in `src/Makevars`.
* The "Possibly misspelled words in DESCRIPTION" reported by the incoming
  feasibility check are author surnames (e.g. Abadie, Arkhangelsky, Hainmueller)
  and standard method acronyms (SCM, SDID, GSC, TASC). These are correct.

## Downstream dependencies

There are currently no downstream dependencies, as this is a new package.
