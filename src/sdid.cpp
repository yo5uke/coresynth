#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Declaration of the solver from optim.cpp
arma::vec solve_simplex_qp(const arma::mat& Q, const arma::vec& c, int max_iter = 10000, double tol = 1e-6);

//' Calculate SDID Unit Weights (omega)
//'
//' Solves the regularized QP:
//' min over omega in Delta: sum_t (sum_i omega_i Y_it - Y_tr_t)^2 + zeta^2 * T_pre * ||omega||^2
//'
//' This corresponds to equation (5) in Arkhangelsky et al. (2021).
//'
//' @param Y_pre Pre-treatment outcome matrix for control units (T_pre x N_co)
//' @param Y_tr_pre Pre-treatment outcome vector for treated unit (T_pre x 1), averaged if multiple
//' @param zeta2 Ridge penalty parameter (zeta^2). The code internally multiplies by T_pre per the paper.
//' @export
// [[Rcpp::export]]
arma::vec sdid_unit_weights_cpp(const arma::mat& Y_pre, const arma::vec& Y_tr_pre, double zeta2) {
  int T_pre = Y_pre.n_rows;

  // Q = Y_pre^T * Y_pre + zeta^2 * T_pre * I   (Eq. 5 in Arkhangelsky et al. 2021)
  arma::mat Q = Y_pre.t() * Y_pre;
  Q.diag() += zeta2 * T_pre;   // <-- paper-correct scaling

  // c = Y_pre^T * Y_tr
  arma::vec c = Y_pre.t() * Y_tr_pre;

  return solve_simplex_qp(Q, c);
}

//' Calculate SDID Time Weights (lambda)
//'
//' Solves the time-weight QP (with implicit intercept lambda_0 concentrated out):
//'
//' min over lambda in Delta_pre: ||Y_post_target - Y_pre_co^T lambda||^2 + zeta_t^2 * N_co * ||lambda||^2
//'
//' The caller is responsible for pre-demeaning Y_pre_co (row-wise) and
//' Y_post_target (subtract the cross-unit mean) to concentrate out lambda_0,
//' as described in Arkhangelsky et al. (2021) Algorithm 1, Eq. (2.3).
//'
//' @param Y_pre_co  Pre-treatment outcomes for control units, row-demeaned (T_pre x N_co)
//' @param Y_post_target Post-treatment mean per control unit, demeaned (N_co x 1)
//' @param zeta_t    Ridge penalty for time weights (paper: 1e-6 * sigma_hat)
//' @export
// [[Rcpp::export]]
arma::vec sdid_time_weights_cpp(const arma::mat& Y_pre_co,
                                const arma::vec& Y_post_target,
                                double zeta_t) {
  int N_co  = Y_pre_co.n_cols;

  // QP: min ||Y_pre_co^T lambda - Y_post_target||^2 + zeta_t^2 * N_co * ||lambda||^2
  arma::mat Q = Y_pre_co * Y_pre_co.t();   // T_pre x T_pre
  Q.diag() += zeta_t * zeta_t * N_co;

  arma::vec c = Y_pre_co * Y_post_target;  // T_pre x 1

  return solve_simplex_qp(Q, c);
}


//' Calculate SDID Estimate (tau_sdid)
//'
//' Given unit weights omega and time weights lambda, computes the SDID
//' estimator as a weighted two-way difference:
//'
//' tau_sdid = (Y_tr_post_mean - Y_tr_pre_wt) - (Y_co_post_wt - Y_co_pre_wt)
//'
//' @param Y_pre_co   Control pre-treatment outcomes (T_pre x N_co)
//' @param Y_post_co  Control post-treatment outcomes (T_post x N_co)
//' @param Y_pre_tr   Treated pre-treatment outcomes (T_pre x 1)
//' @param Y_post_tr  Treated post-treatment outcomes (T_post x 1)
//' @param omega      Unit weights (N_co x 1)
//' @param lambda     Time weights (T_pre x 1)
//' @export
// [[Rcpp::export]]
double sdid_estimate_cpp(
    const arma::mat& Y_pre_co,  const arma::mat& Y_post_co,
    const arma::vec& Y_pre_tr,  const arma::vec& Y_post_tr,
    const arma::vec& omega, const arma::vec& lambda) {

  // Treated: lambda-weighted pre mean, simple post mean
  double tr_pre_wt   = arma::dot(lambda, Y_pre_tr);
  double tr_post_mean = arma::mean(Y_post_tr);

  // Synthetic control: lambda-weighted pre mean, omega-weighted post mean
  arma::vec co_pre_wt   = Y_pre_co.t() * lambda;   // N_co x 1
  arma::rowvec co_post_mean_row = arma::mean(Y_post_co, 0); // 1 x N_co
  arma::vec co_post_mean = co_post_mean_row.t();

  double synth_pre_wt   = arma::dot(omega, co_pre_wt);
  double synth_post_mean = arma::dot(omega, co_post_mean);

  // SDID DiD estimate
  return (tr_post_mean - tr_pre_wt) - (synth_post_mean - synth_pre_wt);
}
