#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Declarations of solvers from optim.cpp
arma::vec solve_simplex_qp(const arma::mat& Q, const arma::vec& c, int max_iter = 10000, double tol = 1e-6);
arma::vec solve_simplex_qp_lr(const arma::mat& B, const arma::vec& b, int max_iter = 10000, double tol = 1e-6);

// Inner Optimization for SCM: find W given V
//' SCM Inner Weights (QP Given V)
//'
//' Solves the inner-loop QP for SCM: given a fixed diagonal metric matrix V,
//' finds donor weights W on the simplex minimising the V-weighted covariate loss.
//'
//' @param X0     Covariate matrix for control units (k x N_co)
//' @param X1     Covariate vector for the treated unit (k x 1)
//' @param V_diag Diagonal of the metric matrix V (k x 1, non-negative, need not sum to 1)
//' @return Donor weight vector W (N_co x 1) on the unit simplex
//' @export
// [[Rcpp::export]]
arma::vec scm_inner_weights_cpp(const arma::mat& X0, const arma::vec& X1, const arma::vec& V_diag) {
  int k    = (int)X0.n_rows;
  int N_co = (int)X0.n_cols;

  if (2 * k < N_co) {
    // Low-rank path (k < N_co/2): form B = diag(sqrt(V)) * X0, b = sqrt(V) % X1.
    // ||X1 - X0 W||_V^2 = ||b - B W||^2  (algebraically equivalent).
    // Avoids N_co x N_co matrix; uses eig_sym(k x k) for Lipschitz constant.
    arma::vec sqV = arma::sqrt(arma::clamp(V_diag, 0.0, arma::datum::inf));
    arma::mat B   = arma::diagmat(sqV) * X0;   // k x N_co
    arma::vec bv  = sqV % X1;                  // k x 1
    return solve_simplex_qp_lr(B, bv);
  }

  // Full path (k >= N_co/2): original implementation.
  arma::mat V = arma::diagmat(V_diag);
  arma::mat Q = X0.t() * V * X0;
  arma::vec c = X0.t() * V * X1;
  return solve_simplex_qp(Q, c);
}

// Internal helper: returns only the weight vector (no Rcpp::List overhead).
// Called by scm_placebo_cpp in inference.cpp via forward declaration.
arma::vec scm_weights_vec_internal(const arma::mat& X0, const arma::vec& X1,
                                    const arma::mat& Z0, const arma::vec& Z1,
                                    int max_iter, double tol) {
  int k = X0.n_rows;
  arma::vec V_diag = arma::ones<arma::vec>(k) / k;
  arma::vec best_W = scm_inner_weights_cpp(X0, X1, V_diag);
  double best_loss = arma::norm(Z1 - Z0 * best_W, 2);

  for (int iter = 0; iter < max_iter; iter++) {
    double max_change = 0.0;
    for (int j = 0; j < k; j++) {
      double orig_v         = V_diag(j);
      double local_best_v   = orig_v;
      double local_best_loss = best_loss;
      arma::vec local_best_W = best_W;
      arma::vec grid = arma::linspace<arma::vec>(0.0, 1.0, 11);
      for (arma::uword g = 0; g < grid.n_elem; g++) {
        V_diag(j) = grid(g);
        double vsum = arma::sum(V_diag);
        if (vsum < 1e-14) { V_diag(j) = orig_v; continue; }  // restore before skipping
        arma::vec V_norm = V_diag / vsum;
        arma::vec cand_W = scm_inner_weights_cpp(X0, X1, V_norm);
        double cand_loss = arma::norm(Z1 - Z0 * cand_W, 2);
        if (cand_loss < local_best_loss) {
          local_best_loss = cand_loss;
          local_best_v    = grid(g);
          local_best_W    = cand_W;
        }
      }
      V_diag(j) = local_best_v;
      double vsum_after = arma::sum(V_diag);
      if (vsum_after > 1e-14) V_diag = V_diag / vsum_after;
      if (local_best_loss < best_loss) {
        max_change = std::max(max_change, best_loss - local_best_loss);
        best_loss  = local_best_loss;
        best_W     = local_best_W;
      }
    }
    if (max_change < tol) break;
  }
  return best_W;
}

// Outer Optimization for SCM using simple Coordinate Descent
//' SCM Outer Weights (Joint Optimization of W and V)
//'
//' Jointly optimises donor weights W (on the simplex) and the diagonal
//' metric matrix V via coordinate descent on the pre-treatment prediction
//' MSPE, following Abadie, Diamond & Hainmueller (2010).
//'
//' When `t_train > 0`, uses out-of-sample V selection per Abadie (2021)
//' §3.2: V is selected by minimising MSPE on a validation window
//' (rows t_train..T_pre-1 of Z), while W is fitted on the training window
//' (rows 0..t_train-1 of X when X and Z have the same row count, i.e. the
//' outcomes-only case). After selecting V*, W is refit on the full data.
//'
//' @param X0      Covariate matrix for control units (k x N_co, typically pre-treatment outcomes)
//' @param X1      Covariate vector for the treated unit (k x 1)
//' @param Z0      Outcome matrix for control units in the pre-treatment window (T_pre x N_co)
//' @param Z1      Outcome vector for the treated unit in the pre-treatment window (T_pre x 1)
//' @param max_iter Maximum coordinate-descent iterations (default 100)
//' @param tol     Convergence tolerance on MSPE improvement (default 1e-4)
//' @param t_train Training window length for out-of-sample V selection.
//'   -1 (default): in-sample V selection (original behaviour).
//'   >0: use rows 0..(t_train-1) of Z for fitting W, rows t_train..(T_pre-1)
//'   as the validation window for V selection, then refit W on full data.
//' @return A list with:
//'   * `W`: Donor weight vector (N_co x 1) on the unit simplex
//'   * `V`: Optimal metric diagonal (k x 1, normalised to sum to 1)
//'   * `loss`: Final pre-treatment prediction loss (full pre-treatment window)
//' @export
// [[Rcpp::export]]
Rcpp::List scm_weights_cpp(const arma::mat& X0, const arma::vec& X1,
                            const arma::mat& Z0, const arma::vec& Z1,
                            int max_iter = 100, double tol = 1e-4,
                            int t_train = -1) {
  int k     = X0.n_rows;
  int T_pre = (int)Z0.n_rows;

  bool do_oos = (t_train > 0 && t_train < T_pre);

  // In OOS mode: W is always fitted on the full predictor matrix X; only the
  // MSPE evaluation window is restricted to the validation period of Z.
  // This avoids the dimension mismatch that arises when X and Z have the same
  // number of rows (outcomes-only case) and V has k=T_pre entries.
  arma::mat Z0_eval;
  arma::vec Z1_eval;

  if (do_oos) {
    // Validation window: rows t_train..(T_pre-1) of Z
    Z0_eval = Z0.rows(t_train, T_pre - 1);
    Z1_eval = Z1.subvec(t_train, T_pre - 1);
  } else {
    Z0_eval = Z0;
    Z1_eval = Z1;
  }

  arma::vec V_diag = arma::ones(k) / k;

  // W is always fitted on the full X (inner QP uses full predictor matrix)
  arma::vec best_W = scm_inner_weights_cpp(X0, X1, V_diag);
  double best_loss = arma::norm(Z1_eval - Z0_eval * best_W, 2);

  // Coordinate Descent over V
  for(int iter = 0; iter < max_iter; iter++) {
    double max_change = 0.0;

    for(int j = 0; j < k; j++) {
      double orig_v          = V_diag(j);
      double local_best_v    = orig_v;
      double local_best_loss = best_loss;
      arma::vec local_best_W = best_W;

      // Grid search over coordinate j — hold all other V elements fixed
      arma::vec grid = arma::linspace<arma::vec>(0.0, 1.0, 11);
      for(int g = 0; g < (int)grid.n_elem; g++) {
        V_diag(j) = grid(g);
        double vsum = arma::sum(V_diag);
        if (vsum < 1e-14) { V_diag(j) = orig_v; continue; }
        arma::vec V_norm = V_diag / vsum;  // normalised copy — V_diag unchanged

        // W fitted on full X; MSPE evaluated on validation window (OOS) or full Z (in-sample)
        arma::vec cand_W    = scm_inner_weights_cpp(X0, X1, V_norm);
        double    cand_loss = arma::norm(Z1_eval - Z0_eval * cand_W, 2);

        if(cand_loss < local_best_loss) {
          local_best_loss = cand_loss;
          local_best_v    = grid(g);
          local_best_W    = cand_W;
        }
      }
      // Restore j-th element to best found, then re-normalise
      V_diag(j) = local_best_v;
      double vsum_after = arma::sum(V_diag);
      if (vsum_after > 1e-14) V_diag = V_diag / vsum_after;

      if(local_best_loss < best_loss) {
        max_change = std::max(max_change, std::abs(local_best_loss - best_loss));
        best_loss  = local_best_loss;
        best_W     = local_best_W;
      }
    }

    if(max_change < tol) break;
  }

  // OOS mode: refit W on full pre-treatment data with selected V*
  arma::vec final_W;
  double    final_loss;
  if (do_oos) {
    final_W    = scm_inner_weights_cpp(X0, X1, V_diag);
    final_loss = arma::norm(Z1 - Z0 * final_W, 2);
  } else {
    final_W    = best_W;
    final_loss = best_loss;
  }

  return Rcpp::List::create(
    Rcpp::Named("W")    = final_W,
    Rcpp::Named("V")    = V_diag,
    Rcpp::Named("loss") = final_loss
  );
}
