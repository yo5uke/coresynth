#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Fast Matrix Completion using Soft-Impute Algorithm
//'
//' Solves: min_L (1/2) ||O o (Y - L)||_F^2 + lambda * ||L||_*
//' via iterative SVD soft-thresholding (Mazumder, Hastie, Tibshirani 2010).
//' Note: lambda is NOT normalized by |O|. Default lambda = 0.01 * sigma_max(Y).
//'
//' @param Y       Observed outcome matrix (N x T). Unobserved entries should be 0.
//' @param O       Binary mask matrix (N x T): 1 = observed, 0 = missing (treated post).
//' @param lambda  Nuclear norm penalty (soft-threshold on singular values).
//' @param max_iter Maximum iterations.
//' @param tol     Convergence tolerance (relative Frobenius norm change).
//' @return A numeric matrix of the same dimension as `Y` (N x T): the
//'   completed low-rank matrix `L` that minimises the soft-thresholded
//'   nuclear-norm objective.
//' @export
// [[Rcpp::export]]
arma::mat soft_impute_cpp(const arma::mat& Y, const arma::mat& O,
                          double lambda, int max_iter = 1000, double tol = 1e-5) {
  arma::mat L = arma::zeros(Y.n_rows, Y.n_cols);

  for(int iter = 0; iter < max_iter; iter++) {
    // Fill-in: use observed Y where observed, current L estimate where missing
    arma::mat Z = O % Y + (1.0 - O) % L;

    // Economy SVD of Z
    arma::mat U;
    arma::vec D;
    arma::mat V;
    bool ok = arma::svd_econ(U, D, V, Z);
    if(!ok) {
      Rcpp::warning("SVD failed to converge in soft_impute; returning current estimate.");
      break;
    }

    // Soft-threshold singular values
    arma::vec D_thresh = arma::max(D - lambda, arma::zeros<arma::vec>(D.n_elem));

    arma::mat L_new = U * arma::diagmat(D_thresh) * V.t();

    // Relative Frobenius change (guard against zero denominator)
    double denom = std::max(arma::norm(L, "fro"), 1e-12);
    double rel_change = arma::norm(L_new - L, "fro") / denom;
    L = L_new;

    if(iter > 0 && rel_change < tol) break;
  }

  return L;
}
