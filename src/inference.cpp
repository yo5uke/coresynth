#include <RcppArmadillo.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Declaration of SDID unit weight solver
arma::vec sdid_unit_weights_cpp(const arma::mat& Y_pre, const arma::vec& Y_tr_pre, double zeta2);

// Forward declarations of internal SCM weight solvers (defined in scm.cpp)
arma::vec scm_weights_vec_internal(const arma::mat& X0, const arma::vec& X1,
                                    const arma::mat& Z0, const arma::vec& Z1,
                                    int max_iter, double tol, bool multistart);
void scm_multistart_batch(const std::vector<arma::mat>& X0d,
                          const std::vector<arma::vec>& X1d,
                          const std::vector<arma::mat>& Z0d,
                          const std::vector<arma::vec>& Z1d,
                          int max_iter, double tol,
                          std::vector<arma::vec>& W_out,
                          std::vector<char>& ok_out);

//' Fast Placebo Test for SDID
//'
//' For each control unit, treats it as the "pseudo-treated" unit and
//' estimates the leave-one-out SDID effect. The distribution of these
//' placebo effects provides a permutation-based null distribution for inference.
//'
//' @param Y_pre     Control units pre-treatment outcomes (T_pre x N_co)
//' @param Y_post    Control units post-treatment outcomes (T_post x N_co)
//' @param time_weights Lambda weights for pre-treatment periods (T_pre x 1)
//' @param zeta2     Ridge penalty (same as used in the main estimate)
//' @return A numeric vector of length `N_co`. Each element is the
//'   leave-one-out placebo SDID effect obtained by treating that control unit
//'   as the pseudo-treated unit; the vector serves as a permutation-based null
//'   distribution for inference.
//' @export
// [[Rcpp::export]]
arma::vec sdid_placebo_cpp(const arma::mat& Y_pre, const arma::mat& Y_post,
                           const arma::vec& time_weights, double zeta2) {
  int N_co  = Y_pre.n_cols;
  arma::vec placebo_effects(N_co);

  #pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < N_co; i++) {
    // Pseudo-treated unit
    arma::vec y_tr_pre  = Y_pre.col(i);
    arma::vec y_tr_post = Y_post.col(i);

    // Donor pool (all other control units)
    // Use arma::regspace for safe integer index generation
    arma::uvec all_idx   = arma::regspace<arma::uvec>(0, N_co - 1);
    arma::uvec donor_idx = all_idx(arma::find(all_idx != (arma::uword)i));

    arma::mat Y_donor_pre  = Y_pre.cols(donor_idx);
    arma::mat Y_donor_post = Y_post.cols(donor_idx);

    // Compute unit weights for this leave-one-out placebo
    arma::vec unit_w = sdid_unit_weights_cpp(Y_donor_pre, y_tr_pre, zeta2);

    // SDID-style DiD using the provided time weights
    double tr_pre_wt    = arma::dot(time_weights, y_tr_pre);
    double tr_post_mean = arma::mean(y_tr_post);

    arma::vec co_pre_wt   = Y_donor_pre.t() * time_weights;        // N_co-1
    arma::vec co_post_mean = arma::mean(Y_donor_post, 0).t();       // N_co-1

    double synth_pre_wt   = arma::dot(unit_w, co_pre_wt);
    double synth_post_mean = arma::dot(unit_w, co_post_mean);

    placebo_effects(i) = (tr_post_mean - tr_pre_wt) - (synth_post_mean - synth_pre_wt);
  }

  return placebo_effects;
}

//' Fast Leave-One-Out Placebo Test for SCM (Abadie et al. 2010)
//'
//' For each control unit, treats it as pseudo-treated and fits SCM weights
//' from the remaining N_co-1 donors. Returns MSPE components for constructing
//' MSPE-ratio permutation p-values in R.
//'
//' @param Y_pre   Control pre-treatment outcomes (T_pre x N_co)
//' @param Y_post  Control post-treatment outcomes (T_post x N_co)
//' @param max_iter Outer coordinate-descent iterations (default 100)
//' @param tol      Convergence tolerance for V updates (default 1e-4)
//' @param z_rows  Optional 1-based pre-period row indices of the outer
//'   evaluation window (the `v_window` of the treated fit), so each placebo
//'   refit optimises V on the same window. `NULL` (default) uses all rows.
//'   MSPE components are always computed on the full pre/post windows.
//' @return A list with:
//'   * `mspe_pre`:  N_co-vector of pre-treatment MSPE per placebo unit
//'   * `mspe_post`: N_co-vector of post-treatment MSPE per placebo unit
//'   * `effects`:   N_co-vector of mean post-period gap per placebo unit
//'   * `gaps`:      (T_pre + T_post) x N_co matrix of placebo gap paths
//' @export
// [[Rcpp::export]]
Rcpp::List scm_placebo_cpp(const arma::mat& Y_pre, const arma::mat& Y_post,
                            int max_iter = 100, double tol = 1e-4,
                            Rcpp::Nullable<Rcpp::IntegerVector> z_rows = R_NilValue) {
  int N_co = Y_pre.n_cols;
  arma::vec mspe_pre(N_co), mspe_post(N_co), effects(N_co);
  arma::mat gaps(Y_pre.n_rows + Y_post.n_rows, N_co);
  arma::uvec all_idx = arma::regspace<arma::uvec>(0, N_co - 1);

  // Convert the window outside the parallel region (no Rcpp API in threads).
  bool has_window = z_rows.isNotNull();
  arma::uvec zr;
  if (has_window) {
    Rcpp::IntegerVector zri(z_rows);
    zr.set_size(zri.size());
    for (int i = 0; i < zri.size(); i++) zr(i) = (arma::uword)(zri[i] - 1);
  }

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < N_co; i++) {
    arma::vec y_pre_i  = Y_pre.col(i);
    arma::vec y_post_i = Y_post.col(i);
    arma::uvec donors  = all_idx(arma::find(all_idx != (arma::uword)i));

    // T_pre x (N_co-1): same layout as Y_co_pre in fit_scm_cpp (k=T_pre rows, N_donors cols)
    arma::mat Y_d_pre  = Y_pre.cols(donors);
    arma::mat Y_d_post = Y_post.cols(donors);

    // Predictors (X) stay the full pre-period outcomes; only the outer
    // evaluation window is restricted -- mirroring the treated fit.
    arma::mat Z0_fit;
    arma::vec Z1_fit;
    if (has_window) {
      Z0_fit = Y_d_pre.rows(zr);
      Z1_fit = y_pre_i(zr);
    } else {
      Z0_fit = Y_d_pre;
      Z1_fit = y_pre_i;
    }

    arma::vec w = scm_weights_vec_internal(Y_d_pre, y_pre_i,
                                           Z0_fit, Z1_fit,
                                           max_iter, tol, false);

    arma::vec synth_pre  = Y_d_pre  * w;
    arma::vec synth_post = Y_d_post * w;

    mspe_pre(i)  = arma::mean(arma::square(y_pre_i  - synth_pre));
    mspe_post(i) = arma::mean(arma::square(y_post_i - synth_post));
    effects(i)   = arma::mean(y_post_i - synth_post);
    gaps.col(i)  = arma::join_cols(y_pre_i - synth_pre, y_post_i - synth_post);
  }

  return Rcpp::List::create(
    Rcpp::Named("mspe_pre")  = mspe_pre,
    Rcpp::Named("mspe_post") = mspe_post,
    Rcpp::Named("effects")   = effects,
    Rcpp::Named("gaps")      = gaps
  );
}

//' Fast Leave-One-Out Placebo Test for SCM with a Predictor Specification
//'
//' Covariate-spec counterpart of [scm_placebo_cpp()]: for each control unit,
//' treats it as pseudo-treated with its own predictor column `X0[, i]` and
//' fits the nested V/W optimisation against the remaining donors' predictors
//' `X0[, -i]`, evaluating the prediction loss on pre-treatment outcomes.
//' Each leave-one-out problem is identical to a [scm_weights_cpp()] call on
//' the same submatrices; iterations are independent and run in parallel
//' under OpenMP.
//'
//' @param X0      Predictor matrix for control units (k x N_co), on the same
//'   scale as the treated fit (SD-scaled when `scale_predictors = TRUE`)
//' @param Y_pre   Control pre-treatment outcomes (T_pre x N_co)
//' @param Y_post  Control post-treatment outcomes (T_post x N_co)
//' @param max_iter Outer coordinate-descent iterations (default 100)
//' @param tol      Convergence tolerance for V updates (default 1e-4)
//' @param z_rows  Optional 1-based pre-period row indices of the outer
//'   evaluation window (the `v_window` of the treated fit), so each placebo
//'   refit optimises V on the same window. `NULL` (default) uses all rows.
//'   MSPE components are always computed on the full pre/post windows.
//' @param multistart If `TRUE`, each placebo refit uses the same
//'   deterministic multi-start outer search as the treated fit, keeping the
//'   permutation test symmetric.
//' @return A list with:
//'   * `mspe_pre`:  N_co-vector of pre-treatment MSPE per placebo unit
//'   * `mspe_post`: N_co-vector of post-treatment MSPE per placebo unit
//'   * `effects`:   N_co-vector of mean post-period gap per placebo unit
//'   * `gaps`:      (T_pre + T_post) x N_co matrix of placebo gap paths
//'   A placebo unit whose solver fails yields NaN entries.
//' @export
// [[Rcpp::export]]
Rcpp::List scm_placebo_x_cpp(const arma::mat& X0,
                              const arma::mat& Y_pre, const arma::mat& Y_post,
                              int max_iter = 100, double tol = 1e-4,
                              Rcpp::Nullable<Rcpp::IntegerVector> z_rows = R_NilValue,
                              bool multistart = false) {
  int N_co = Y_pre.n_cols;
  arma::vec mspe_pre(N_co), mspe_post(N_co), effects(N_co);
  arma::mat gaps(Y_pre.n_rows + Y_post.n_rows, N_co);
  arma::uvec all_idx = arma::regspace<arma::uvec>(0, N_co - 1);

  // Convert the window outside the parallel region (no Rcpp API in threads).
  bool has_window = z_rows.isNotNull();
  arma::uvec zr;
  if (has_window) {
    Rcpp::IntegerVector zri(z_rows);
    zr.set_size(zri.size());
    for (int i = 0; i < zri.size(); i++) zr(i) = (arma::uword)(zri[i] - 1);
  }

  if (multistart) {
    // Flat (donor x candidate) batch: materialise the per-donor matrices
    // up front (cheap copies), let scm_multistart_batch schedule every
    // candidate pipeline of every donor across one thread pool, then
    // compute the MSPE components. Per-donor refit cost spans two orders
    // of magnitude on real data, so donor-level parallelism alone leaves
    // most threads idle behind the slowest donor.
    std::vector<arma::mat> X0d(N_co), Z0f(N_co);
    std::vector<arma::vec> x1d(N_co), Z1f(N_co);
    for (int i = 0; i < N_co; i++) {
      arma::uvec donors = all_idx(arma::find(all_idx != (arma::uword)i));
      X0d[i] = X0.cols(donors);
      x1d[i] = X0.col(i);
      arma::mat Y_d_pre = Y_pre.cols(donors);
      arma::vec y_pre_i = Y_pre.col(i);
      if (has_window) {
        Z0f[i] = Y_d_pre.rows(zr);
        Z1f[i] = y_pre_i(zr);
      } else {
        Z0f[i] = Y_d_pre;
        Z1f[i] = y_pre_i;
      }
    }

    std::vector<arma::vec> W(N_co);
    std::vector<char> ok;
    scm_multistart_batch(X0d, x1d, Z0f, Z1f, max_iter, tol, W, ok);

    for (int i = 0; i < N_co; i++) {
      // A donor whose refit failed entirely carries NaN through all
      // fields (same contract as the solver-failure catch below).
      if (!ok[i] || !W[i].is_finite()) {
        mspe_pre(i)  = arma::datum::nan;
        mspe_post(i) = arma::datum::nan;
        effects(i)   = arma::datum::nan;
        gaps.col(i).fill(arma::datum::nan);
        continue;
      }
      arma::uvec donors = all_idx(arma::find(all_idx != (arma::uword)i));
      arma::vec y_pre_i  = Y_pre.col(i);
      arma::vec y_post_i = Y_post.col(i);
      arma::mat Y_d_pre  = Y_pre.cols(donors);
      arma::mat Y_d_post = Y_post.cols(donors);

      arma::vec synth_pre  = Y_d_pre  * W[i];
      arma::vec synth_post = Y_d_post * W[i];

      mspe_pre(i)  = arma::mean(arma::square(y_pre_i  - synth_pre));
      mspe_post(i) = arma::mean(arma::square(y_post_i - synth_post));
      effects(i)   = arma::mean(y_post_i - synth_post);
      gaps.col(i)  = arma::join_cols(y_pre_i - synth_pre,
                                     y_post_i - synth_post);
    }

    return Rcpp::List::create(
      Rcpp::Named("mspe_pre")  = mspe_pre,
      Rcpp::Named("mspe_post") = mspe_post,
      Rcpp::Named("effects")   = effects,
      Rcpp::Named("gaps")      = gaps
    );
  }

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < N_co; i++) {
    arma::uvec donors = all_idx(arma::find(all_idx != (arma::uword)i));
    // Exceptions must not escape the OpenMP region; a failed solver marks
    // this unit NaN (same contract as the previous R-level tryCatch).
    try {
      arma::vec x1_i     = X0.col(i);
      arma::mat X0_d     = X0.cols(donors);
      arma::vec y_pre_i  = Y_pre.col(i);
      arma::vec y_post_i = Y_post.col(i);
      arma::mat Y_d_pre  = Y_pre.cols(donors);
      arma::mat Y_d_post = Y_post.cols(donors);

      arma::mat Z0_fit;
      arma::vec Z1_fit;
      if (has_window) {
        Z0_fit = Y_d_pre.rows(zr);
        Z1_fit = y_pre_i(zr);
      } else {
        Z0_fit = Y_d_pre;
        Z1_fit = y_pre_i;
      }

      arma::vec w = scm_weights_vec_internal(X0_d, x1_i,
                                             Z0_fit, Z1_fit,
                                             max_iter, tol, false);

      arma::vec synth_pre  = Y_d_pre  * w;
      arma::vec synth_post = Y_d_post * w;

      mspe_pre(i)  = arma::mean(arma::square(y_pre_i  - synth_pre));
      mspe_post(i) = arma::mean(arma::square(y_post_i - synth_post));
      effects(i)   = arma::mean(y_post_i - synth_post);
      gaps.col(i)  = arma::join_cols(y_pre_i - synth_pre, y_post_i - synth_post);
    } catch (...) {
      mspe_pre(i)  = arma::datum::nan;
      mspe_post(i) = arma::datum::nan;
      effects(i)   = arma::datum::nan;
      gaps.col(i).fill(arma::datum::nan);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mspe_pre")  = mspe_pre,
    Rcpp::Named("mspe_post") = mspe_post,
    Rcpp::Named("effects")   = effects,
    Rcpp::Named("gaps")      = gaps
  );
}
