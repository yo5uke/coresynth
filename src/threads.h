#ifndef CORESYNTH_THREADS_H
#define CORESYNTH_THREADS_H

#include <Rinternals.h>
#include <cstdlib>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

// Number of OpenMP threads for the parallel loops. Must be called on the main
// thread, before entering a parallel region: it reads R state.
//
//  * options(coresynth.threads = n) wins when set to a positive number;
//  * otherwise, under R CMD check --as-cran (which sets _R_CHECK_LIMIT_CORES_),
//    at most two threads, as the CRAN Repository Policy requires;
//  * otherwise the OpenMP default, which honours OMP_NUM_THREADS and
//    OMP_THREAD_LIMIT.
inline int cs_num_threads() {
#ifdef _OPENMP
  int n = omp_get_max_threads();
  SEXP opt = Rf_GetOption1(Rf_install("coresynth.threads"));
  if (opt != R_NilValue && Rf_length(opt) == 1 &&
      (Rf_isInteger(opt) || Rf_isReal(opt))) {
    double v = Rf_asReal(opt);
    if (R_FINITE(v) && v >= 1.0) return v < n ? (int)v : n;
  }
  const char* lim = std::getenv("_R_CHECK_LIMIT_CORES_");
  if (lim != nullptr && *lim != '\0' && std::strcmp(lim, "false") != 0 &&
      std::strcmp(lim, "FALSE") != 0) {
    return n < 2 ? n : 2;
  }
  return n;
#else
  return 1;
#endif
}

#endif
