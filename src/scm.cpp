#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Declarations of solvers from optim.cpp
arma::vec solve_simplex_qp(const arma::mat& Q, const arma::vec& c, int max_iter = 10000, double tol = 1e-6, Rcpp::Nullable<Rcpp::NumericVector> x0 = R_NilValue);
arma::vec solve_simplex_qp_lr(const arma::mat& B, const arma::vec& b, int max_iter = 10000, double tol = 1e-6);
arma::vec proj_simplex(arma::vec y);

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

// Active-set solver (Lawson-Hanson style) for min 0.5 w'Q_g w - c'w on the
// simplex, where Q_g = Q + dv * r r' (rank-1 term applied implicitly).
// Starting from a warm-start support A, alternates:
//  * solve the equality-constrained QP on A (stationarity + sum-to-one);
//  * primal infeasible (negative weight) -> drop the most negative coordinate;
//  * dual infeasible (an inactive coordinate with negative multiplier) ->
//    add the worst violator.
// Terminates only at a KKT-verified exact optimum (up to linear-solve
// precision). Returns false on solve failure or pivot-cap overrun so the
// caller can fall back to FISTA; x_out is untouched in that case.
// SCM solutions are sparse, so |A| stays small and each pivot costs
// O(|A|^3 + N^2) -- orders of magnitude below FISTA's iteration count on
// the ill-conditioned Q_g typical of collinear donor outcomes.
static bool active_set_simplex_rank1(const arma::mat& Q, const arma::vec& r,
                                     double dv, const arma::vec& c,
                                     arma::uvec A, arma::vec& x_out,
                                     int max_pivots = 100) {
  const arma::uword N = Q.n_rows;
  if (A.n_elem == 0) return false;

  for (int pivot = 0; pivot < max_pivots; pivot++) {
    arma::uword m = A.n_elem;
    // Stationarity on the active set: Q_AA w_A + mu * 1 = c_A, 1'w_A = 1
    arma::vec rA = r(A);
    arma::mat K(m + 1, m + 1, arma::fill::zeros);
    K.submat(0, 0, m - 1, m - 1) = Q.submat(A, A) + dv * (rA * rA.t());
    K.col(m).head(m).ones();
    K.row(m).head(m).ones();
    arma::vec rhs(m + 1);
    rhs.head(m) = c(A);
    rhs(m)      = 1.0;

    arma::vec sol;
    if (!arma::solve(sol, K, rhs, arma::solve_opts::no_approx)) return false;
    arma::vec wA = sol.head(m);
    double    mu = sol(m);

    if (wA.min() < -1e-12) {
      if (m <= 1) return false;
      A.shed_row(wA.index_min());
      continue;
    }

    arma::vec w(N, arma::fill::zeros);
    w(A) = arma::clamp(wA, 0.0, arma::datum::inf);
    double s = arma::accu(w);
    if (s <= 0.0) return false;
    w /= s;  // remove clamp dust (<= 1e-12 per coordinate)

    // Dual feasibility on the inactive set: lambda_i = g_i + mu >= 0
    arma::vec g   = Q * w + (dv * arma::dot(r, w)) * r - c;
    arma::vec lam = g + mu;
    lam.elem(A).zeros();  // active coordinates are stationary by construction
    double eps = 1e-9 * (1.0 + arma::norm(g, "inf"));
    arma::uword worst = lam.index_min();
    if (lam(worst) >= -eps) {
      x_out = w;
      return true;
    }
    A.insert_rows(A.n_elem, arma::uvec{worst});
  }
  return false;
}

// Simplex QP for Q_g = Q + dv * r r' (rank-1 term applied implicitly), warm
// started from x. Strategy: try the exact active-set solve seeded with the
// warm start's support; if it cannot verify KKT (singular subproblem, pivot
// cap), fall back to warm-started FISTA with periodic active-set retries.
// Pure Armadillo (no Rcpp types): safe to call inside OpenMP threads.
// L must be an upper bound on lambda_max(Q_g); a loose bound only shrinks
// the FISTA step size, it does not change the minimiser.
static arma::vec fista_simplex_rank1(const arma::mat& Q, const arma::vec& r,
                                     double dv, const arma::vec& c,
                                     double L, arma::vec x,
                                     int max_iter = 10000, double tol = 1e-6) {
  x = proj_simplex(x);

  arma::vec x_as;
  if (active_set_simplex_rank1(Q, r, dv, c, arma::find(x > 1e-10), x_as)) {
    return x_as;
  }

  if (L < 1e-14) L = 1.0;
  double t = 1.0 / L;
  arma::vec y = x;
  double t_acc = 1.0;

  for (int iter = 0; iter < max_iter; iter++) {
    arma::vec x_prev = x;
    arma::vec grad = Q * y + (dv * arma::dot(r, y)) * r - c;
    x = proj_simplex(y - t * grad);
    // Adaptive restart (gradient scheme) -- see solve_simplex_qp.
    if (arma::dot(grad, x - x_prev) > 0.0) {
      t_acc = 1.0;
    }
    double t_acc_next = (1.0 + std::sqrt(1.0 + 4.0 * t_acc * t_acc)) / 2.0;
    y = x + ((t_acc - 1.0) / t_acc_next) * (x - x_prev);
    t_acc = t_acc_next;
    if (arma::norm(x - x_prev, 2) < tol) break;
    // Retry the exact solve once the support has moved on a little.
    if (iter > 0 && iter % 50 == 0) {
      if (active_set_simplex_rank1(Q, r, dv, c, arma::find(x > 1e-10), x_as)) {
        return x_as;
      }
    }
  }
  return x;
}

// Coordinate-descent core for the nested V/W SCM optimisation. Shared by
// scm_weights_cpp (R-facing, arbitrary evaluation window) and
// scm_weights_vec_internal (placebo loop). Semantics are the classic ADH
// coordinate descent (grid {0, 0.1, ..., 1} per V coordinate, strict-improve
// accept rule, per-coordinate renormalisation of V), but each inner QP
// exploits structure:
//  * scale invariance: the QP argmin under V/sum(V) equals the argmin under
//    the unnormalised V, so per-grid-point renormalisation of Q is skipped;
//  * rank-1 structure: changing one V coordinate shifts Q by dv * r_j r_j',
//    applied implicitly in the QP instead of rebuilding Q = X0' V X0
//    (O(k N^2)) for every grid point;
//  * warm starts: each QP starts from the previous grid point's solution
//    (active-set exact solve, FISTA fallback -- see fista_simplex_rank1);
//  * Lipschitz bound: lambda_max is computed by eig_sym once per sweep and
//    tracked through rank-1 updates via lambda_max(Q + dv r r') <=
//    lambda_max(Q) + max(dv, 0) ||r||^2 (a valid upper bound, so FISTA
//    convergence is unaffected).
// V_diag must come in normalised (sum 1); it is updated in place and leaves
// normalised. Returns the best weight vector; best_loss is updated in place.
static arma::vec scm_coord_descent_core(const arma::mat& X0, const arma::vec& X1,
                                        const arma::mat& Z0_eval,
                                        const arma::vec& Z1_eval,
                                        int max_iter, double tol,
                                        arma::vec& V_diag, double& best_loss) {
  int k = X0.n_rows;
  arma::vec best_W = scm_inner_weights_cpp(X0, X1, V_diag);
  best_loss = arma::norm(Z1_eval - Z0_eval * best_W, 2);

  // Gram caches for the current (normalised) V_diag
  arma::mat Q_base = X0.t() * arma::diagmat(V_diag) * X0;  // N x N
  arma::vec c_base = X0.t() * (V_diag % X1);               // N
  arma::vec row_n2 = arma::sum(arma::square(X0), 1);       // k: ||r_j||^2

  arma::vec grid = arma::linspace<arma::vec>(0.0, 1.0, 11);

  for (int iter = 0; iter < max_iter; iter++) {
    // Exact Lipschitz constant once per sweep; rank-1 bounds accumulate
    // slack within the sweep only.
    arma::mat Q_sym = (Q_base + Q_base.t()) / 2.0;
    double L_base = arma::max(arma::eig_sym(Q_sym));
    if (L_base < 1e-14) L_base = 1.0;

    double max_change = 0.0;
    for (int j = 0; j < k; j++) {
      double orig_v          = V_diag(j);
      double local_best_v    = orig_v;
      double local_best_loss = best_loss;
      arma::vec local_best_W = best_W;
      arma::vec r_j  = X0.row(j).t();
      arma::vec warm = best_W;

      for (arma::uword g = 0; g < grid.n_elem; g++) {
        if (std::abs(grid(g) - orig_v) < 1e-15) continue;  // current point: ties only
        double vsum = arma::sum(V_diag) - orig_v + grid(g);
        if (vsum < 1e-14) continue;
        double dv = grid(g) - orig_v;
        arma::vec c_g = c_base + (dv * X1(j)) * r_j;
        double    L_g = L_base + std::max(dv, 0.0) * row_n2(j);
        arma::vec cand_W = fista_simplex_rank1(Q_base, r_j, dv, c_g, L_g, warm);
        warm = cand_W;
        double cand_loss = arma::norm(Z1_eval - Z0_eval * cand_W, 2);
        if (cand_loss < local_best_loss) {
          local_best_loss = cand_loss;
          local_best_v    = grid(g);
          local_best_W    = cand_W;
        }
      }

      if (local_best_loss < best_loss) {
        double dv_acc = local_best_v - orig_v;
        Q_base += dv_acc * (r_j * r_j.t());
        c_base += (dv_acc * X1(j)) * r_j;
        L_base += std::max(dv_acc, 0.0) * row_n2(j);
        V_diag(j) = local_best_v;
        double vsum_after = arma::sum(V_diag);
        if (vsum_after > 1e-14) {
          V_diag /= vsum_after;
          Q_base /= vsum_after;
          c_base /= vsum_after;
          L_base /= vsum_after;
        }
        max_change = std::max(max_change, best_loss - local_best_loss);
        best_loss  = local_best_loss;
        best_W     = local_best_W;
      }
    }
    if (max_change < tol) break;
  }
  return best_W;
}

// Internal helper: returns only the weight vector (no Rcpp::List overhead).
// Called by scm_placebo_cpp in inference.cpp via forward declaration.
arma::vec scm_weights_vec_internal(const arma::mat& X0, const arma::vec& X1,
                                    const arma::mat& Z0, const arma::vec& Z1,
                                    int max_iter, double tol) {
  int k = X0.n_rows;
  arma::vec V_diag = arma::ones<arma::vec>(k) / k;
  double best_loss = 0.0;
  return scm_coord_descent_core(X0, X1, Z0, Z1, max_iter, tol,
                                V_diag, best_loss);
}

// Outer Optimization for SCM using simple Coordinate Descent
//' SCM Outer Weights (Joint Optimization of W and V)
//'
//' Jointly optimises donor weights W (on the simplex) and the diagonal
//' metric matrix V via coordinate descent on the pre-treatment prediction
//' MSPE, following Abadie, Diamond & Hainmueller (2010).
//'
//' When `t_train > 0`, V is selected by minimising MSPE on a validation
//' window (rows t_train..T_pre-1 of Z) while W is fitted on the full
//' predictor matrix X. This is appropriate when X is a fixed predictor
//' matrix that contains no validation-period outcome information (the
//' user-supplied predictors case). For the outcomes-only case the proper
//' Abadie (2021) S.3.2 train/validation split is implemented in R
//' (`.scm_oos_outcomes()`): candidate W(V) are fitted on training-half
//' outcomes only, by passing the training rows as X and the validation
//' rows as Z to this function with `t_train = -1`.
//'
//' @param X0      Covariate matrix for control units (k x N_co, typically pre-treatment outcomes)
//' @param X1      Covariate vector for the treated unit (k x 1)
//' @param Z0      Outcome matrix for control units in the pre-treatment window (T_pre x N_co)
//' @param Z1      Outcome vector for the treated unit in the pre-treatment window (T_pre x 1)
//' @param max_iter Maximum coordinate-descent iterations (default 100)
//' @param tol     Convergence tolerance on MSPE improvement (default 1e-4)
//' @param t_train Validation-window split for V selection.
//'   -1 (default): V selected on the full Z window (in-sample).
//'   Positive: rows t_train..(T_pre-1) of Z form the validation window used
//'   to select V (W is fitted on the full X throughout); after selecting V*,
//'   W is refit and the reported loss uses the full Z window.
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

  // Coordinate descent over V; W always fitted on the full X (see
  // scm_coord_descent_core for the fast inner-QP machinery).
  double best_loss = 0.0;
  arma::vec best_W = scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                                            max_iter, tol, V_diag, best_loss);

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
