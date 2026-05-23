#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Declaration of SDID unit weight solver
arma::vec sdid_unit_weights_cpp(const arma::mat& Y_pre, const arma::vec& Y_tr_pre, double zeta2);

// Forward declaration of internal SCM weight solver (defined in scm.cpp)
arma::vec scm_weights_vec_internal(const arma::mat& X0, const arma::vec& X1,
                                    const arma::mat& Z0, const arma::vec& Z1,
                                    int max_iter, double tol);

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
//' @return A list with:
//'   * `mspe_pre`:  N_co-vector of pre-treatment MSPE per placebo unit
//'   * `mspe_post`: N_co-vector of post-treatment MSPE per placebo unit
//'   * `effects`:   N_co-vector of mean post-period gap per placebo unit
//' @export
// [[Rcpp::export]]
Rcpp::List scm_placebo_cpp(const arma::mat& Y_pre, const arma::mat& Y_post,
                            int max_iter = 100, double tol = 1e-4) {
  int N_co = Y_pre.n_cols;
  arma::vec mspe_pre(N_co), mspe_post(N_co), effects(N_co);
  arma::uvec all_idx = arma::regspace<arma::uvec>(0, N_co - 1);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < N_co; i++) {
    arma::vec y_pre_i  = Y_pre.col(i);
    arma::vec y_post_i = Y_post.col(i);
    arma::uvec donors  = all_idx(arma::find(all_idx != (arma::uword)i));

    // T_pre x (N_co-1): same layout as Y_co_pre in fit_scm_cpp (k=T_pre rows, N_donors cols)
    arma::mat Y_d_pre  = Y_pre.cols(donors);
    arma::mat Y_d_post = Y_post.cols(donors);

    arma::vec w = scm_weights_vec_internal(Y_d_pre, y_pre_i,
                                           Y_d_pre, y_pre_i,
                                           max_iter, tol);

    arma::vec synth_pre  = Y_d_pre  * w;
    arma::vec synth_post = Y_d_post * w;

    mspe_pre(i)  = arma::mean(arma::square(y_pre_i  - synth_pre));
    mspe_post(i) = arma::mean(arma::square(y_post_i - synth_post));
    effects(i)   = arma::mean(y_post_i - synth_post);
  }

  return Rcpp::List::create(
    Rcpp::Named("mspe_pre")  = mspe_pre,
    Rcpp::Named("mspe_post") = mspe_post,
    Rcpp::Named("effects")   = effects
  );
}
