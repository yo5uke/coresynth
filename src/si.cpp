#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Tensor Unfolding (Matricization) for Synthetic Interventions
//'
//' @param T_cube A 3D array (cube) of dimensions (n1, n2, n3)
//' @param mode The mode to unfold along (1, 2, or 3)
//' @export
// [[Rcpp::export]]
arma::mat tensor_unfold_cpp(const arma::cube& T_cube, int mode) {
  int n1 = T_cube.n_rows;
  int n2 = T_cube.n_cols;
  int n3 = T_cube.n_slices;
  
  if (mode == 1) {
    // Mode-1 unfolding: size (n1) x (n2 * n3)
    arma::mat res(n1, n2 * n3);
    for(int k = 0; k < n3; k++) {
      res.cols(k * n2, (k + 1) * n2 - 1) = T_cube.slice(k);
    }
    return res;
  } 
  else if (mode == 2) {
    // Mode-2 unfolding: size (n2) x (n1 * n3)
    arma::mat res(n2, n1 * n3);
    for(int k = 0; k < n3; k++) {
      res.cols(k * n1, (k + 1) * n1 - 1) = T_cube.slice(k).t();
    }
    return res;
  }
  else if (mode == 3) {
    // Mode-3 unfolding: size (n3) x (n1 * n2)
    arma::mat res(n3, n1 * n2);
    for(int k = 0; k < n3; k++) {
      // vec() vectorizes column-by-column
      res.row(k) = arma::vectorise(T_cube.slice(k)).t();
    }
    return res;
  }
  
  Rcpp::stop("Mode must be 1, 2, or 3");
}

//' SI-PCR: Synthetic Interventions via Principal Component Regression
//'
//' Implements the SI-PCR estimator of Agarwal et al. (2025).
//' Uses the top-k SVD of pre-treatment control outcomes to find donor
//' weights that predict each treated unit's pre-treatment trajectory,
//' then applies those weights to post-treatment control outcomes.
//'
//' @param Y_pre_co  Pre-treatment control outcomes (T_pre x N_co)
//' @param Y_post_co Post-treatment control outcomes (T_post x N_co)
//' @param Y_pre_tr  Pre-treatment treated outcomes (T_pre x N_tr)
//' @param k         Number of SVD components to retain
//' @return A list with:
//'   * `W`: Donor weight matrix (N_co x N_tr)
//'   * `Y_hat`: Counterfactual post-treatment outcomes (T_post x N_tr)
//' @export
// [[Rcpp::export]]
Rcpp::List si_pcr_cpp(const arma::mat& Y_pre_co,
                      const arma::mat& Y_post_co,
                      const arma::mat& Y_pre_tr,
                      int k) {
  int T_pre = Y_pre_co.n_rows;
  int N_co  = Y_pre_co.n_cols;
  int T_post = Y_post_co.n_rows;
  int N_tr  = Y_pre_tr.n_cols;

  if (k < 1)
    Rcpp::stop("k must be >= 1");
  if (k > std::min(T_pre, N_co))
    Rcpp::stop("k (%d) must be <= min(T_pre=%d, N_co=%d)", k, T_pre, N_co);
  if (Y_pre_co.n_rows != Y_pre_tr.n_rows)
    Rcpp::stop("T_pre mismatch between control and treated matrices");
  if (Y_pre_co.n_cols != Y_post_co.n_cols)
    Rcpp::stop("N_co mismatch between pre and post control matrices");

  // Economy SVD of Y_pre_co^T (N_co x T_pre)
  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, Y_pre_co.t());

  // Truncate to rank k
  arma::mat U_k = U.cols(0, k - 1);  // N_co x k
  arma::vec S_k = S.subvec(0, k - 1);
  arma::mat V_k = V.cols(0, k - 1);  // T_pre x k

  // Pseudo-inverse of rank-k truncation: Pinv = U_k * diag(1/S_k) * V_k^T
  // Guard near-zero singular values (same threshold as MATLAB pinv)
  double tol_si = arma::datum::eps *
                  static_cast<double>(std::max(T_pre, N_co)) * S_k(0);
  arma::vec S_k_inv(k);
  for (int l = 0; l < k; l++)
    S_k_inv(l) = (S_k(l) > tol_si) ? (1.0 / S_k(l)) : 0.0;
  arma::mat Pinv = U_k * arma::diagmat(S_k_inv) * V_k.t();

  // Compute weights and counterfactuals for each treated unit
  arma::mat W(N_co, N_tr);
  arma::mat Y_hat(T_post, N_tr);

  for (int i = 0; i < N_tr; i++) {
    arma::vec w_i    = Pinv * Y_pre_tr.col(i);  // N_co x 1
    W.col(i)         = w_i;
    Y_hat.col(i)     = Y_post_co * w_i;          // T_post x 1
  }

  return Rcpp::List::create(
    Rcpp::Named("W")     = W,
    Rcpp::Named("Y_hat") = Y_hat
  );
}
