#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Kalman Filter and RTS Smoother (TASC)
//'
//' Implements the Kalman filter (forward pass) and Rauch-Tung-Striebel smoother
//' (backward pass) for the state-space model in Rho et al. (2026):
//'
//'   State:       z_{t+1} = A z_t + C + eta_t,  eta_t ~ N(0, Q)
//'   Observation: y_t     = W z_t + eps_t,       eps_t ~ N(0, R)
//'
//' Observation rows with NA (treated post-intervention) are automatically
//' dropped at each time step so only control-unit rows update the filter.
//'
//' The P update uses the numerically stable Joseph form:
//'   P_{t|t} = (I - K W_obs) P_{t|t-1} (I - K W_obs)^T + K R_obs K^T
//'
//' @param Y  Observed data matrix (N x T). Use NA for unobserved entries.
//' @param W  Observation / loading matrix (N x r)
//' @param A  State transition matrix (r x r). Pass diag(r) for random-walk dynamics.
//' @param C  State drift vector (r x 1)
//' @param Q  State noise covariance (r x r)
//' @param R  Observation noise covariance (N x N, diagonal in practice)
//' @param z0 Initial state mean (r x 1)
//' @param P0 Initial state covariance (r x r)
//' @return A list with z_smooth, P_smooth, P_cross, z_pred, z_upd.
//'   P_cross is an r x r x (T-1) cube. Slice t (C++ 0-indexed, t=0,...,T-2)
//'   stores P_{t+1, t | T} (0-indexed), i.e. P_{t+2, t+1 | T} in 1-indexed
//'   Shumway-Stoffer notation. Formula: P_{t+1|T} * J_t^T (eq. 6.68-6.69).
//' @export
// [[Rcpp::export]]
Rcpp::List kalman_smoother_cpp(const arma::mat& Y,
                               const arma::mat& W,
                               const arma::mat& A,
                               const arma::vec& C,
                               const arma::mat& Q,
                               const arma::mat& R,
                               const arma::vec& z0,
                               const arma::mat& P0) {
  int T = Y.n_cols;
  int r = z0.n_elem;
  int N = Y.n_rows;

  // Storage
  arma::mat  z_pred(r, T),  z_upd(r, T);
  arma::cube P_pred(r, r, T), P_upd(r, r, T);

  arma::vec  z_curr = z0;
  arma::mat  P_curr = P0;
  arma::mat  I_r = arma::eye(r, r);

  // ── Forward Kalman Filter ──
  for(int t = 0; t < T; t++) {
    // Predict
    arma::vec z_p = A * z_curr + C;
    arma::mat P_p = A * P_curr * A.t() + Q;
    z_pred.col(t) = z_p;
    P_pred.slice(t) = P_p;

    // Find finite (observed) rows
    arma::vec y_t = Y.col(t);
    arma::uvec obs_idx = arma::find_finite(y_t);

    if(obs_idx.n_elem > 0) {
      arma::mat W_obs = W.rows(obs_idx);
      arma::vec y_obs = y_t.elem(obs_idx);
      arma::mat R_obs = R.submat(obs_idx, obs_idx);

      // Innovation covariance S and Kalman gain K
      arma::mat S = W_obs * P_p * W_obs.t() + R_obs;
      arma::mat K = P_p * W_obs.t() * arma::inv_sympd(S);

      // State update
      z_curr = z_p + K * (y_obs - W_obs * z_p);

      // Covariance update — Joseph form for numerical stability
      arma::mat IKW = I_r - K * W_obs;
      P_curr = IKW * P_p * IKW.t() + K * R_obs * K.t();

      // Enforce symmetry
      P_curr = 0.5 * (P_curr + P_curr.t());
    } else {
      z_curr = z_p;
      P_curr = P_p;
    }

    z_upd.col(t)    = z_curr;
    P_upd.slice(t)  = P_curr;
  }

  // ── Backward RTS Smoother ──
  arma::mat  z_smooth(r, T);
  arma::cube P_smooth(r, r, T);
  arma::cube P_cross(r, r, T - 1);  // P_{t+1,t|T}: lag-one cross-covariances

  z_smooth.col(T - 1)    = z_upd.col(T - 1);
  P_smooth.slice(T - 1)  = P_upd.slice(T - 1);

  for(int t = T - 2; t >= 0; t--) {
    arma::mat P_pred_next = P_pred.slice(t + 1);
    // Smoother gain J_t = P_{t|t} * A^T * P_{t+1|t}^{-1}
    arma::mat J = P_upd.slice(t) * A.t() * arma::inv_sympd(P_pred_next);

    z_smooth.col(t) = z_upd.col(t)
                    + J * (z_smooth.col(t + 1) - z_pred.col(t + 1));

    P_smooth.slice(t) = P_upd.slice(t)
                      + J * (P_smooth.slice(t + 1) - P_pred_next) * J.t();

    // Enforce symmetry
    P_smooth.slice(t) = 0.5 * (P_smooth.slice(t) + P_smooth.slice(t).t());

    // Lag-one cross-covariance: P_{t+1, t | T} = P_{t+1|T} * J_t^T
    // (Shumway & Stoffer, eq. 6.68-6.69)
    P_cross.slice(t) = P_smooth.slice(t + 1) * J.t();
  }

  return Rcpp::List::create(
    Rcpp::Named("z_smooth") = z_smooth,
    Rcpp::Named("P_smooth") = P_smooth,
    Rcpp::Named("P_cross")  = P_cross,
    Rcpp::Named("z_pred")   = z_pred,
    Rcpp::Named("z_upd")    = z_upd
  );
}
