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

// Accelerated projected gradient (FISTA) on the unit simplex, shared by the
// dense and factored QP entry points below. `grad(y)` returns Q y - c; `L` is
// an upper bound on lambda_max(Q). Pure Armadillo, so it is safe inside
// OpenMP threads.
template <typename Grad>
static arma::vec fista_simplex_core(Grad grad_fn, arma::vec x, double L,
                                    int max_iter, double tol) {
  if (L < 1e-14) L = 1.0;  // guard against degenerate Q
  double t = 1.0 / L;

  arma::vec y = x;
  double t_acc = 1.0;

  for (int iter = 0; iter < max_iter; iter++) {
    arma::vec x_prev = x;

    // Gradient step
    arma::vec grad = grad_fn(y);

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

    if (arma::norm(x - x_prev, 2) < tol) {
      break;
    }
  }

  return x;
}

// Largest eigenvalue of a symmetric matrix. Symmetrise before eig_sym: when
// the matrix is a rank-deficient Gram (e.g. Q = X0'VX0 with fewer rows than
// columns), floating-point errors in the near-zero eigenspace can make it
// look slightly asymmetric to Armadillo's strict check.
static double lambda_max_sym(const arma::mat& M) {
  arma::mat M_sym = (M + M.t()) / 2.0;
  arma::vec eigval = arma::eig_sym(M_sym);
  return arma::max(eigval);
}

// Projected Gradient Descent for Quadratic Programming on the unit simplex
// min 0.5 * x^T Q x - c^T x
// s.t. sum(x) = 1, x >= 0
// x0: optional warm start (projected onto the simplex). Block coordinate
// descent callers (partially pooled staggered SCM) resolve near-identical
// QPs every sweep; restarting FISTA from the previous block solution cuts
// the iteration count by orders of magnitude on ill-conditioned Q.
// lipschitz: the largest eigenvalue of Q when the caller already knows it,
// typically from a factor Q = F'F with few rows, where lambda_max(FF') is the
// same number at a fraction of the O(n^3) eig_sym below. Non-positive (the
// default) computes it here.
// [[Rcpp::export]]
arma::vec solve_simplex_qp(const arma::mat& Q, const arma::vec& c, int max_iter = 10000, double tol = 1e-6,
                           Rcpp::Nullable<Rcpp::NumericVector> x0 = R_NilValue,
                           double lipschitz = -1.0) {
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
  double L = lipschitz;
  if (L <= 0.0) L = lambda_max_sym(Q);

  return fista_simplex_core(
    [&](const arma::vec& y) -> arma::vec { return Q * y - c; },
    x, L, max_iter, tol);
}

// Ridge-penalised least squares on the unit simplex:
//   min_w 0.5 ||A w||^2 + 0.5 * ridge * ||w||^2 - (A'b)'w,  w in simplex,
// i.e. solve_simplex_qp with Q = A'A + ridge * I and c = A'b, without forming
// Q when A is wide. The SDID weight problems have this shape, and the unit
// weights are wide (T_pre rows, N_co columns): there Q is N_co x N_co of rank
// at most T_pre, so the factor gives the gradient A'(A y) + ridge * y in
// O(T_pre N_co) instead of O(N_co^2), and the Lipschitz constant
// lambda_max(AA') + ridge from a T_pre x T_pre eigenproblem instead of an
// O(N_co^3) one. The FISTA iterates are the same as for the dense form up to
// round-off. Pure Armadillo (OpenMP-safe); starts from the uniform vector.
arma::vec solve_simplex_ridge_ls(const arma::mat& A, const arma::vec& b,
                                 double ridge, int max_iter, double tol) {
  const arma::uword n = A.n_cols;
  arma::vec x = arma::ones(n) / n;
  arma::vec c = A.t() * b;

  if (A.n_rows >= n) {
    arma::mat Q = A.t() * A;
    Q.diag() += ridge;
    return fista_simplex_core(
      [&](const arma::vec& y) -> arma::vec { return Q * y - c; },
      x, lambda_max_sym(Q), max_iter, tol);
  }

  double L = lambda_max_sym(A * A.t()) + ridge;
  return fista_simplex_core(
    [&](const arma::vec& y) -> arma::vec {
      return A.t() * (A * y) + ridge * y - c;
    },
    x, L, max_iter, tol);
}
