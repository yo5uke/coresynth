#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Projection onto the unit simplex
// min ||x - y||^2 s.t. sum(x) = 1, x >= 0
// Using algorithm by Duchi et al. (2008)
// [[Rcpp::export]]
arma::vec proj_simplex(arma::vec y) {
  int n = y.n_elem;
  arma::vec u = arma::sort(y, "descend");
  
  double cssv = 0.0;
  double rho = 0.0;
  
  for(int i = 0; i < n; i++) {
    cssv += u(i);
    if(u(i) > (cssv - 1.0) / (i + 1.0)) {
      rho = (cssv - 1.0) / (i + 1.0);
    } else {
      break; // Since u is sorted descending, the condition will fail for the rest
    }
  }
  
  arma::vec x = arma::max(y - rho, arma::zeros(n));
  return x;
}

// Projected Gradient Descent for Quadratic Programming on the unit simplex
// min 0.5 * x^T Q x - c^T x
// s.t. sum(x) = 1, x >= 0
// x0: optional warm start (projected onto the simplex). Block coordinate
// descent callers (partially pooled staggered SCM) resolve near-identical
// QPs every sweep; restarting FISTA from the previous block solution cuts
// the iteration count by orders of magnitude on ill-conditioned Q.
// [[Rcpp::export]]
arma::vec solve_simplex_qp(const arma::mat& Q, const arma::vec& c, int max_iter = 10000, double tol = 1e-6,
                           Rcpp::Nullable<Rcpp::NumericVector> x0 = R_NilValue) {
  int n = c.n_elem;
  arma::vec x;
  if (x0.isNotNull()) {
    x = proj_simplex(Rcpp::as<arma::vec>(x0.get()));
  } else {
    x = arma::ones(n) / n; // initialize uniformly
  }

  // Lipschitz constant for step size: largest eigenvalue of Q (optimal for FISTA).
  // eig_sym is O(n^3) but gives the tightest step size, minimising FISTA iterations.
  // norm_inf is cheaper to compute but yields a smaller step → more iterations;
  // benchmarking shows eig_sym wins overall for the n <= 100 range typical in SCM.
  // Symmetrise before eig_sym: when X0 has fewer rows than columns (k < N_co),
  // Q = X0'VX0 is rank-deficient; floating-point errors in the near-zero
  // eigenspace can make Q appear slightly asymmetric to Armadillo's strict check.
  arma::mat Q_sym = (Q + Q.t()) / 2.0;
  arma::vec eigval = arma::eig_sym(Q_sym);
  double L = arma::max(eigval);
  if(L < 1e-14) L = 1.0;  // guard against degenerate Q
  double t = 1.0 / L;
  
  arma::vec y = x;
  double t_acc = 1.0;

  for(int iter = 0; iter < max_iter; iter++) {
    arma::vec x_prev = x;

    // Gradient step
    arma::vec grad = Q * y - c;

    // Projection step
    x = proj_simplex(y - t * grad);

    // Adaptive restart (O'Donoghue & Candes 2015, gradient scheme): when the
    // momentum direction disagrees with the descent direction, the
    // acceleration is overshooting and slowing convergence on ill-conditioned
    // Q (the dominant cost in SCM, where Q = X0'VX0 is highly collinear).
    // Resetting t_acc removes the stale momentum. This only changes the path
    // to the optimum, not the optimum itself.
    if (arma::dot(grad, x - x_prev) > 0.0) {
      t_acc = 1.0;
    }

    // FISTA acceleration
    double t_acc_next = (1.0 + std::sqrt(1.0 + 4.0 * t_acc * t_acc)) / 2.0;
    y = x + ((t_acc - 1.0) / t_acc_next) * (x - x_prev);
    t_acc = t_acc_next;

    if(arma::norm(x - x_prev, 2) < tol) {
      break;
    }
  }

  return x;
}

// Low-rank FISTA for inner SCM QP when k < N_co/2.
// Solves min_{w in Delta} ||B*w - b||^2  (B = sqrt(V)*X0: k x N_co)
// without forming the N_co x N_co matrix Q = B'B explicitly.
// Lipschitz constant: lambda_max(B*B') via eig_sym(k x k) [O(k^3)] vs O(N_co^3).
// Gradient per iter: B'*(B*y - b) [O(k*N_co)] vs Q*y [O(N_co^2)].
// Only activated when 2*k < N_co (guarantees per-iteration speedup).
arma::vec solve_simplex_qp_lr(const arma::mat& B, const arma::vec& b,
                               int max_iter, double tol) {
  int N_co = B.n_cols;
  arma::vec x = arma::ones<arma::vec>(N_co) / N_co;

  arma::mat BBt     = B * B.t();                    // k x k
  arma::mat BBt_sym = (BBt + BBt.t()) / 2.0;        // symmetry guard
  arma::vec eigval  = arma::eig_sym(BBt_sym);
  double L = arma::max(eigval);
  if (L < 1e-14) L = 1.0;
  double t = 1.0 / L;

  arma::vec y = x;
  double t_acc = 1.0;

  for (int iter = 0; iter < max_iter; iter++) {
    arma::vec x_prev = x;
    arma::vec grad   = B.t() * (B * y - b);         // O(k * N_co)
    x = proj_simplex(y - t * grad);
    // Adaptive restart (gradient scheme) -- see solve_simplex_qp for rationale.
    if (arma::dot(grad, x - x_prev) > 0.0) {
      t_acc = 1.0;
    }
    double t_acc_next = (1.0 + std::sqrt(1.0 + 4.0 * t_acc * t_acc)) / 2.0;
    y = x + ((t_acc - 1.0) / t_acc_next) * (x - x_prev);
    t_acc = t_acc_next;
    if (arma::norm(x - x_prev, 2) < tol) break;
  }
  return x;
}
