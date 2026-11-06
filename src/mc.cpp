#include <RcppArmadillo.h>
#include <vector>
#include "threads.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Pure-Armadillo Soft-Impute core (OpenMP-safe). Returns false if SVD fails.
static bool soft_impute_fe_core(const arma::mat& Y, const arma::mat& O,
                                double lambda, bool estimate_fe,
                                int max_iter, double tol, int fe_sweeps,
                                arma::mat& L, arma::vec& a, arma::vec& b,
                                arma::mat& M) {
  const arma::uword T = Y.n_rows, N = Y.n_cols;
  L.zeros(T, N);
  a.zeros(T);
  b.zeros(N);

  // Observation counts per row and column; guarded so an all-missing row or
  // column leaves its effect at zero instead of dividing by zero.
  arma::vec n_row_s = arma::clamp(arma::sum(O, 1), 1.0, arma::datum::inf);
  arma::vec n_col_s = arma::clamp(arma::sum(O, 0).t(), 1.0, arma::datum::inf);
  arma::mat Ot;
  if (estimate_fe && fe_sweeps > 0) {
    Ot = O.t();
  }

  arma::mat OY = O % Y;
  arma::mat one_minus_O = 1.0 - O;
  arma::mat Z(T, N);
  arma::mat L_new(T, N);
  arma::mat U;
  arma::vec D;
  arma::mat V;
  bool svd_ok = true;
  double norm_L = 0.0;

  for (int iter = 0; iter < max_iter; iter++) {
    if (estimate_fe) {
      arma::mat R = O % (Y - L);
      arma::vec r_row = arma::sum(R, 1);
      arma::vec r_col = arma::sum(R, 0).t();
      for (int s = 0; s < fe_sweeps; s++) {
        a = (r_row - O * b) / n_row_s;
        b = (r_col - Ot * a) / n_col_s;
      }
      Z = Y;
      Z.each_col() -= a;
      Z.each_row() -= b.t();
      Z = O % Z + one_minus_O % L;
    } else {
      Z = OY + one_minus_O % L;
    }

    if (!arma::svd_econ(U, D, V, Z)) {
      svd_ok = false;
      break;
    }

    arma::vec D_thresh = arma::clamp(D - lambda, 0.0, arma::datum::inf);
    arma::uword r_pos = arma::accu(D_thresh > 0.0);
    if (r_pos == 0) {
      L_new.zeros();
    } else {
      arma::mat U_scaled = U.head_cols(r_pos);
      U_scaled.each_row() %= D_thresh.head(r_pos).t();
      L_new = U_scaled * V.head_cols(r_pos).t();
    }

    double denom = std::max(norm_L, 1e-12);
    double rel_change = arma::norm(L_new - L, "fro") / denom;
    L = L_new;
    norm_L = (r_pos == 0) ? 0.0 : arma::norm(L, "fro");

    if (iter > 0 && rel_change < tol) break;
  }

  M = L;
  if (estimate_fe) {
    M.each_col() += a;
    M.each_row() += b.t();
  }
  return svd_ok;
}

//' Matrix Completion with Optional Two-Way Fixed Effects
//'
//' Solves
//'   `min over L, a, b of (1/2) * squaredFrobenius(O o (Y - L - a1' - 1b')) + lambda * nuclearNorm(L)`
//' by alternating a back-fitting update of the fixed effects with one
//' soft-thresholded SVD step for `L` (Mazumder, Hastie & Tibshirani 2010;
//' Athey et al. 2021, S.4 Eq. 4.3).
//'
//' The row and column effects are the unpenalised time and unit fixed effects
//' of Athey et al. (2021, S.4 Eq. 4.3; see also S.8.1): "we regularize the
//' estimates of L*, but do not wish to regularize the estimates of the fixed
//' effects". Leaving them inside `L` shrinks the panel's level towards zero
//' along with its factor structure, which biases the imputed counterfactual
//' whenever the level is far from zero. `estimate_fe = false` reproduces the
//' pure low-rank baseline `Y = L + e` (Mazumder et al. 2010, Eq. 2).
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
  if (Y.n_rows == 0 || Y.n_cols == 0) {
    Rcpp::stop("Y must have at least 1 row and 1 column.");
  }
  if (O.n_rows != Y.n_rows || O.n_cols != Y.n_cols) {
    Rcpp::stop("O must have the same dimensions as Y.");
  }
  if (!std::isfinite(lambda) || lambda < 0.0) {
    Rcpp::stop("lambda must be a non-negative finite number.");
  }
  if (max_iter < 1) {
    Rcpp::stop("max_iter must be >= 1.");
  }
  if (!std::isfinite(tol) || tol <= 0.0) {
    Rcpp::stop("tol must be a positive finite number.");
  }
  if (fe_sweeps < 0) {
    Rcpp::stop("fe_sweeps must be >= 0.");
  }

  arma::mat L, M;
  arma::vec a, b;
  if (!soft_impute_fe_core(Y, O, lambda, estimate_fe, max_iter, tol, fe_sweeps,
                           L, a, b, M)) {
    Rcpp::warning("SVD failed to converge in soft_impute; returning current estimate.");
  }

  return Rcpp::List::create(
    Rcpp::Named("L") = L,
    Rcpp::Named("a") = a,
    Rcpp::Named("b") = b,
    Rcpp::Named("M") = M
  );
}

// Parallel evaluation of the (fold x lambda) grid in .mc_cv_lambda().
// Fold masks (keep_list, held_list) are drawn in R first so the RNG state and
// fold assignments are identical to the sequential loop.
// [[Rcpp::export]]
arma::mat mc_cv_cpp(const arma::mat& Y_obs, const Rcpp::List& keep_list,
                    const Rcpp::List& held_list, const arma::vec& lambdas,
                    bool estimate_fe = true) {
  const int n_folds = keep_list.size();
  const int n_lam   = lambdas.n_elem;
  const arma::uword T = Y_obs.n_rows, N = Y_obs.n_cols;

  std::vector<arma::mat> O_folds(n_folds, arma::zeros(T, N));
  std::vector<arma::uvec> held_folds(n_folds);
  for (int k = 0; k < n_folds; k++) {
    Rcpp::IntegerVector keep_k = keep_list[k];
    Rcpp::IntegerVector held_k = held_list[k];
    for (int j = 0; j < keep_k.size(); j++) {
      O_folds[k]((arma::uword)(keep_k[j] - 1)) = 1.0;
    }
    held_folds[k].set_size(held_k.size());
    for (int j = 0; j < held_k.size(); j++) {
      held_folds[k](j) = (arma::uword)(held_k[j] - 1);
    }
  }

  arma::mat err(n_folds, n_lam);
  err.fill(arma::datum::nan);
  const int total = n_folds * n_lam;
  const int n_threads = cs_num_threads();

  #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
  for (int idx = 0; idx < total; idx++) {
    try {
      int k = idx / n_lam;
      int i = idx % n_lam;
      if (held_folds[k].n_elem == 0) continue;
      arma::mat L, M;
      arma::vec a, b;
      soft_impute_fe_core(Y_obs, O_folds[k], lambdas(i), estimate_fe,
                          1000, 1e-5, 5, L, a, b, M);
      arma::vec diff = Y_obs.elem(held_folds[k]) - M.elem(held_folds[k]);
      err(k, i) = arma::mean(arma::square(diff));
    } catch (...) {
      // Leave err(k, i) as NaN on failure
    }
  }

  return err;
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
