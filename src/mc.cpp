#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Matrix Completion with Optional Two-Way Fixed Effects
//'
//' Solves
//'   `min over L, a, b of (1/2) * squaredFrobenius(O o (Y - L - a1' - 1b')) + lambda * nuclearNorm(L)`
//' by alternating a back-fitting update of the fixed effects with one
//' soft-thresholded SVD step for `L` (Mazumder, Hastie & Tibshirani 2010).
//'
//' The row and column effects are the unpenalised fixed effects of Athey et
//' al. (2021, S.8.1): "we regularize the estimates of L*, but do not wish to
//' regularize the estimates of the fixed effects". Leaving them inside `L`
//' shrinks the panel's level towards zero along with its factor structure,
//' which biases the imputed counterfactual whenever the level is far from
//' zero. `estimate_fe = false` reproduces the S.4 baseline model `Y = L + e`.
//'
//' Both effects absorb an additive constant, so `a` and `b` are identified
//' only up to a shift between them; their sum, which is what enters the
//' completed matrix, is unaffected.
//'
//' @param Y       Observed outcome matrix (T x N). Unobserved entries should be 0.
//' @param O       Binary mask matrix (T x N): 1 = observed, 0 = missing (treated post).
//' @param lambda  Nuclear norm penalty (soft-threshold on singular values).
//' @param estimate_fe Estimate unpenalised row and column effects.
//' @param max_iter Maximum iterations.
//' @param tol     Convergence tolerance (relative Frobenius norm change).
//' @param fe_sweeps Back-fitting sweeps per iteration. The mask makes the two
//'   effects non-orthogonal, so one closed-form pass does not solve for both.
//' @return A list with `L` (the penalised low-rank part), `a` (row effects),
//'   `b` (column effects), and `M`, the completed matrix `L + a 1' + 1 b'`.
//' @export
// [[Rcpp::export]]
Rcpp::List soft_impute_fe_cpp(const arma::mat& Y, const arma::mat& O,
                              double lambda, bool estimate_fe = true,
                              int max_iter = 1000, double tol = 1e-5,
                              int fe_sweeps = 5) {
  const arma::uword T = Y.n_rows, N = Y.n_cols;
  arma::mat L = arma::zeros(T, N);
  arma::vec a = arma::zeros(T);
  arma::vec b = arma::zeros(N);

  // Observation counts per row and column; guarded so an all-missing row or
  // column leaves its effect at zero instead of dividing by zero.
  arma::vec  n_row = arma::sum(O, 1);
  arma::rowvec n_col = arma::sum(O, 0);
  arma::vec  n_row_s = arma::clamp(n_row, 1.0, arma::datum::inf);
  arma::rowvec n_col_s = arma::clamp(n_col, 1.0, arma::datum::inf);

  for (int iter = 0; iter < max_iter; iter++) {
    if (estimate_fe) {
      arma::mat R = O % (Y - L);
      for (int s = 0; s < fe_sweeps; s++) {
        a = arma::sum(R - (O.each_row() % b.t()), 1) / n_row_s;
        b = (arma::sum(R - (O.each_col() % a), 0) / n_col_s).t();
      }
    }
    arma::mat FE = arma::repmat(a, 1, N) + arma::repmat(b.t(), T, 1);

    // Fill-in: observed cells carry Y net of the fixed effects, missing cells
    // carry the current L
    arma::mat Z = O % (Y - FE) + (1.0 - O) % L;

    arma::mat U;
    arma::vec D;
    arma::mat V;
    if (!arma::svd_econ(U, D, V, Z)) {
      Rcpp::warning("SVD failed to converge in soft_impute; returning current estimate.");
      break;
    }

    arma::vec D_thresh = arma::max(D - lambda, arma::zeros<arma::vec>(D.n_elem));
    arma::mat L_new = U * arma::diagmat(D_thresh) * V.t();

    double denom = std::max(arma::norm(L, "fro"), 1e-12);
    double rel_change = arma::norm(L_new - L, "fro") / denom;
    L = L_new;

    if (iter > 0 && rel_change < tol) break;
  }

  arma::mat FE = arma::repmat(a, 1, N) + arma::repmat(b.t(), T, 1);
  return Rcpp::List::create(
    Rcpp::Named("L") = L,
    Rcpp::Named("a") = a,
    Rcpp::Named("b") = b,
    Rcpp::Named("M") = L + FE
  );
}

//' Fast Matrix Completion using Soft-Impute Algorithm
//'
//' Solves: min_L (1/2) ||O o (Y - L)||_F^2 + lambda * ||L||_*
//' via iterative SVD soft-thresholding (Mazumder, Hastie, Tibshirani 2010).
//' Note: lambda is NOT normalized by |O|.
//'
//' This is [soft_impute_fe_cpp()] with the fixed effects switched off.
//'
//' @param Y       Observed outcome matrix (T x N). Unobserved entries should be 0.
//' @param O       Binary mask matrix (T x N): 1 = observed, 0 = missing (treated post).
//' @param lambda  Nuclear norm penalty (soft-threshold on singular values).
//' @param max_iter Maximum iterations.
//' @param tol     Convergence tolerance (relative Frobenius norm change).
//' @return A numeric matrix of the same dimension as `Y` (T x N): the
//'   completed low-rank matrix `L` that minimises the soft-thresholded
//'   nuclear-norm objective.
//' @export
// [[Rcpp::export]]
arma::mat soft_impute_cpp(const arma::mat& Y, const arma::mat& O,
                          double lambda, int max_iter = 1000, double tol = 1e-5) {
  Rcpp::List res = soft_impute_fe_cpp(Y, O, lambda, false, max_iter, tol, 0);
  return Rcpp::as<arma::mat>(res["L"]);
}
