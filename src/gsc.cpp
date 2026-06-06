#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Fast Interactive Fixed Effects (IFE) for Generalized Synthetic Control
//'
//' Implements Xu (2017) IFE model with optional covariate adjustment.
//' When X_co has p > 0 slices, runs an EM loop alternating between:
//'   E-step: truncated SVD of Y_tilde = Y_co - X_co * beta
//'   M-step: panel OLS to update beta given current factors
//' When X_co has 0 slices (default), falls back to the plain 3-step estimator.
//'
//' @param Y_co     Control units outcome matrix (T x N_co)
//' @param Y_tr_pre Treated units pre-treatment outcomes (T_pre x N_tr)
//' @param r        Number of latent factors (must be <= min(T, N_co))
//' @param X_co     Time-varying covariate cube (T x N_co x p). Pass an empty
//'                 cube (0 slices) for the covariate-free estimator.
//' @param X_tr_pre Time-varying covariate cube for treated units in the
//'                 pre-treatment window (T_pre x N_tr x p). Required for
//'                 correct Step 2 loading estimation per Xu (2017): lambda_hat
//'                 is estimated from Y_tr_pre - X_tr_pre * beta (covariate-
//'                 demeaned). Pass an empty cube (0 slices) to skip demeaning
//'                 (backward-compatible, but biased when beta != 0).
//' @param max_iter Maximum EM iterations (default 50)
//' @param tol      Convergence tolerance on relative beta change (default 1e-6)
//' @return A list with components:
//'   * `F`: estimated time factors (T x r).
//'   * `L_co`: control-unit factor loadings (N_co x r).
//'   * `L_tr`: treated-unit factor loadings (N_tr x r).
//'   * `Y_tr_hat`: estimated treated-unit counterfactual outcomes (T x N_tr).
//'   * `singular_values`: singular values from the final truncated SVD.
//'   * `beta`: estimated covariate coefficients (p x 1), empty when no
//'     covariates are supplied.
//' @export
// [[Rcpp::export]]
Rcpp::List gsc_ife_cpp(const arma::mat& Y_co,
                        const arma::mat& Y_tr_pre,
                        int r,
                        const arma::cube& X_co,
                        const arma::cube& X_tr_pre,
                        int max_iter = 50,
                        double tol   = 1e-6) {
  int T     = (int)Y_co.n_rows;
  int N_co  = (int)Y_co.n_cols;
  int T_pre = (int)Y_tr_pre.n_rows;
  int p     = (int)X_co.n_slices;

  if (r > std::min(T, N_co))
    Rcpp::stop("r must be <= min(T, N_co)");
  if (T_pre < r)
    Rcpp::stop("T_pre must be >= r to identify treated loadings");

  arma::vec beta = arma::zeros<arma::vec>(p);
  arma::mat U, V, F, L_co;
  arma::vec D;

  if (p == 0) {
    // ── Plain 3-step estimator (no covariates, backward-compatible) ──────────
    arma::svd_econ(U, D, V, Y_co);
    F    = U.cols(0, r - 1) * arma::diagmat(D.subvec(0, r - 1));
    L_co = V.cols(0, r - 1);

  } else {
    // ── EM loop: alternating SVD (E-step) and panel OLS (M-step) ─────────────
    // Initial E-step with beta=0
    arma::svd_econ(U, D, V, Y_co);
    F    = U.cols(0, r - 1) * arma::diagmat(D.subvec(0, r - 1));
    L_co = V.cols(0, r - 1);

    int n_obs = T * N_co;

    for (int iter = 0; iter < max_iter; iter++) {
      arma::vec beta_old = beta;

      // M-step: OLS regression of residuals R_it on covariates x_{it}
      // R_it = Y_co(t,i) - F(t,:) . L_co(i,:)
      // X_stack(t*N_co + i, j) = X_co(t, i, j)
      arma::mat X_stack(n_obs, p);
      arma::vec R_stack(n_obs);

      for (int t = 0; t < T; t++) {
        for (int i = 0; i < N_co; i++) {
          int idx = t * N_co + i;
          R_stack(idx) = Y_co(t, i) - arma::dot(F.row(t), L_co.row(i));
          for (int j = 0; j < p; j++) {
            X_stack(idx, j) = X_co(t, i, j);
          }
        }
      }

      // Ridge-regularised OLS: beta = (X'X + ridge*I)^{-1} X'R
      double ridge = 1e-8;
      arma::mat XtX = X_stack.t() * X_stack + ridge * arma::eye(p, p);
      arma::vec XtR = X_stack.t() * R_stack;
      arma::vec beta_new = arma::solve(XtX, XtR);

      // E-step: update F, L_co from demeaned Y_tilde = Y_co - X_co * beta
      arma::mat Y_tilde = Y_co;
      for (int j = 0; j < p; j++)
        Y_tilde -= X_co.slice(j) * beta_new(j);

      arma::svd_econ(U, D, V, Y_tilde);
      F    = U.cols(0, r - 1) * arma::diagmat(D.subvec(0, r - 1));
      L_co = V.cols(0, r - 1);

      // Convergence: relative change in beta
      double delta = arma::norm(beta_new - beta_old, 2) /
                     (1.0 + arma::norm(beta_old, 2));
      beta = beta_new;
      if (delta < tol) break;
    }
  }

  // ── Treated loadings: Xu (2017) Step 2 ─────────────────────────────────────
  // lambda_hat_i = (F0' F0)^{-1} F0' (Y0_i - X0_i beta)
  // Demean Y_tr_pre by X_tr_pre * beta when covariates are provided.
  arma::mat F_pre   = F.rows(0, T_pre - 1);
  arma::mat Y_tr_dm = Y_tr_pre;
  int p_tr = (int)X_tr_pre.n_slices;
  if (p_tr > 0) {
    for (int j = 0; j < p_tr; j++)
      Y_tr_dm -= X_tr_pre.slice(j) * beta(j);
  }
  arma::mat L_tr_t  = arma::pinv(F_pre) * Y_tr_dm;
  arma::mat L_tr    = L_tr_t.t();
  arma::mat Y_tr_hat = F * L_tr_t;

  return Rcpp::List::create(
    Rcpp::Named("F")               = F,
    Rcpp::Named("L_co")            = L_co,
    Rcpp::Named("L_tr")            = L_tr,
    Rcpp::Named("Y_tr_hat")        = Y_tr_hat,
    Rcpp::Named("singular_values") = D,
    Rcpp::Named("beta")            = beta
  );
}
