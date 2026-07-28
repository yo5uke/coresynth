#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <random>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Declarations of solvers from optim.cpp
arma::vec proj_simplex(arma::vec y);

// Helmert-style orthonormal basis of the null space of the sum constraint
// {z : 1'z = 0}: column j has 1/sqrt(j(j+1)) in its first j rows,
// -j/sqrt(j(j+1)) at row j+1, zeros below. Analytic (no SVD), deterministic.
static arma::mat sum_null_basis(arma::uword m) {
  arma::mat Z(m, m - 1, arma::fill::zeros);
  for (arma::uword j = 1; j < m; j++) {
    double s = 1.0 / std::sqrt((double)j * (double)(j + 1));
    for (arma::uword i = 0; i < j; i++) Z(i, j - 1) = s;
    Z(j, j - 1) = -(double)j * s;
  }
  return Z;
}

// Build the factored inner-QP data (B, b) for a given predictor weighting.
// V need not be normalised: the QP argmin is invariant to a positive rescaling
// of V, which is what lets the coordinate sweep skip renormalising per grid
// point.
static inline void scm_make_factor(const arma::mat& X0, const arma::vec& X1,
                                   const arma::vec& V, arma::mat& B,
                                   arma::vec& b) {
  arma::vec sqV = arma::sqrt(arma::clamp(V, 0.0, arma::datum::inf));
  B = arma::diagmat(sqV) * X0;   // k x N
  b = sqV % X1;                  // k
}

// ---- Factored inner QP -----------------------------------------------------
// Every inner QP of the nested V/W optimisation is
//   min_{w in simplex} 0.5 w' (X0' diag(V) X0) w - (X0' (V % X1))' w
//     == 0.5 ||b - B w||^2 + const,   B = diag(sqrt(V)) X0 (k x N),
//                                     b = sqrt(V) % X1     (k).
// The solvers below carry the k x N factor (B, b) instead of the N x N Gram
// Q = B'B. Q has rank at most k, so with few predictors and many donors it is
// maximally rank deficient -- exactly the regime where forming it is both the
// most expensive and the least informative thing to do. From the factor the
// gradient Q w - c = B'(B w - b) costs O(kN) instead of O(N^2), and the face
// solve below costs O(k^2 m) instead of O(m^3). A one-coordinate change of V
// is expressed by rebuilding B in O(kN), which replaces the rank-1 update of
// Q the callers used to carry.

// Equality-constrained face problem on the active set A:
//   min ||b - B_A w_A||^2  s.t.  1'w_A = 1,
// returning the face minimiser w_A and the multiplier mu of the sum
// constraint (stationarity: G w_A + mu 1 = c_A with G = B_A'B_A, c_A = B_A'b).
// Three interchangeable routes, selected on cost alone -- all three return the
// same minimiser up to round-off:
//  * cheap_face: bordered-KKT direct solve of the (m+1)x(m+1) saddle system,
//    tried first in the outcomes-only regime and falling through to the exact
//    routes below when it fails (see the note at the branch).
//  * m > k: null-space projection carried in the k-dimensional factor.
//    w_A = w0 + pinv(B_A P) (b - B_A w0) with w0 = 1/m and P = I - 11'/m.
//    For any orthonormal basis Z of {z : 1'z = 0} we have B_A P = (B_A Z) Z',
//    and Z' is a co-isometry, so pinv(B_A P) = Z pinv(B_A Z): this is exactly
//    the Helmert null-space solution below, including its rank truncation,
//    at O(k^2 m) instead of O(m^3). This is the path that carries predictor
//    fits, where m runs to the size of the donor pool.
//  * otherwise: the explicit Helmert null-space route on the m x m face Gram,
//    which is the cheaper of the two once m <= k.
// Rank truncation is at 1e-9 of the leading eigenvalue of the face Gram in
// the null-space routes (equivalently sqrt(1e-9) of the leading singular
// value of the factor): curvature below that carries no statistical
// information. The caller's dual-feasibility check gates every result.
static bool face_solve(const arma::mat& BA, const arma::vec& b,
                       bool cheap_face, arma::vec& wA, double& mu) {
  const arma::uword m = BA.n_cols;
  const arma::uword k = BA.n_rows;

  if (m == 1) {
    wA = arma::vec{1.0};
    mu = arma::dot(BA.col(0), b) - arma::dot(BA.col(0), BA.col(0));
    return true;
  }

  if (cheap_face) {
    // Fast attempt: bordered-KKT direct solve of the (m+1)x(m+1) saddle
    // system. In the outcomes-only regime V is dense and small faces are well
    // conditioned, so this succeeds on the faces that carry the fit and is
    // markedly cheaper than the null-space routes. It cannot be the only
    // route: once m exceeds the rank of the face Gram the saddle system is
    // singular and the solve fails outright, and it is a knife edge that a
    // rescaling of the predictor units can flip. On failure we fall through
    // to the exact null-space routes below rather than giving up -- the
    // caller would otherwise fall back to FISTA, which stops at its 1e-6
    // tolerance and leaves the outer V search reading an inexact objective.
    arma::mat Gc  = BA.t() * BA;
    arma::vec cAc = BA.t() * b;
    arma::mat K(m + 1, m + 1, arma::fill::zeros);
    K.submat(0, 0, m - 1, m - 1) = Gc;
    K.col(m).head(m).ones();
    K.row(m).head(m).ones();
    arma::vec rhs(m + 1);
    rhs.head(m) = cAc;
    rhs(m)      = 1.0;
    arma::vec sol;
    if (arma::solve(sol, K, rhs, arma::solve_opts::no_approx)) {
      wA = sol.head(m);
      mu = sol(m);
      return true;
    }
  }

  if (m > k) {
    arma::vec w0(m);
    w0.fill(1.0 / m);
    arma::vec rowmean = arma::mean(BA, 1);      // == BA * w0
    arma::vec r = b - rowmean;
    arma::mat M = BA;
    M.each_col() -= rowmean;                    // BA * P, k x m

    arma::mat U, Vs;
    arma::vec sv;
    if (!arma::svd_econ(U, sv, Vs, M)) return false;
    arma::vec z(m, arma::fill::zeros);
    if (sv.n_elem > 0) {
      double smax = sv.max();
      double scut = 3.16227766016837933e-5 * smax;   // sqrt(1e-9) * smax
      arma::vec Ur = U.t() * r;
      for (arma::uword i = 0; i < sv.n_elem; i++) {
        Ur(i) = (sv(i) > scut) ? Ur(i) / sv(i) : 0.0;
      }
      z = Vs * Ur;
    }
    wA = w0 + z;
    if (!wA.is_finite()) return false;
    mu = -arma::mean(BA.t() * (BA * wA - b));
    return true;
  }

  arma::mat G  = BA.t() * BA;
  arma::vec cA = BA.t() * b;

  arma::mat Z_ns = sum_null_basis(m);
  arma::vec w0(m);
  w0.fill(1.0 / m);
  arma::mat H = Z_ns.t() * G * Z_ns;
  arma::vec bz = Z_ns.t() * (cA - G * w0);

  // Fast path: Cholesky solve when H is comfortably well conditioned, judged
  // by the diagonal ratio of the factor (a scale-robust rcond proxy -- ulp
  // perturbations cannot flip it, unlike the internal failure threshold of a
  // plain linear solve). Degenerate path: eigendecomposition with the rank
  // truncation described above. At the branch boundary the two routes agree,
  // so the branch itself is benign.
  arma::vec u;
  arma::mat R;
  bool fast = arma::chol(R, arma::symmatu(H));
  if (fast) {
    double dmin = R.diag().min();
    double dmax = R.diag().max();
    fast = (dmin > 0.0) && (dmin >= 3.2e-5 * dmax);  // rcond(H) >~ 1e-9
    if (fast) {
      u = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), bz));
    }
  }
  if (!fast) {
    arma::vec ev;
    arma::mat evec;
    if (!arma::eig_sym(ev, evec, arma::symmatu(H))) return false;
    double ev_max = ev.max();
    if (ev_max <= 0.0) {
      u = arma::zeros<arma::vec>(m - 1);  // flat face: w0 is stationary
    } else {
      double cutoff = 1e-9 * ev_max;
      arma::vec bt = evec.t() * bz;
      for (arma::uword i = 0; i < ev.n_elem; i++) {
        bt(i) = (ev(i) > cutoff) ? bt(i) / ev(i) : 0.0;
      }
      u = evec * bt;
    }
  }

  wA = w0 + Z_ns * u;
  if (!wA.is_finite()) return false;
  mu = -arma::mean(G * wA - cA);
  return true;
}

// Active-set solver (Lawson-Hanson style) for min ||b - B w||^2 on the
// simplex. Starting from a warm-start support A, alternates:
//  * solve the equality-constrained face problem on A (see face_solve);
//  * primal infeasible (negative weight) -> drop the most negative coordinate;
//  * dual infeasible (an inactive coordinate with negative multiplier) ->
//    add the worst violator.
// Terminates only at a KKT-verified exact optimum (up to linear-solve
// precision). Returns false on solve failure or pivot-cap overrun so the
// caller can fall back to FISTA; x_out is untouched in that case.
static bool active_set_simplex(const arma::mat& B, const arma::vec& b,
                               arma::uvec A, arma::vec& x_out,
                               int max_pivots = 100,
                               bool cheap_face = false) {
  const arma::uword N = B.n_cols;
  if (A.n_elem == 0) return false;

  for (int pivot = 0; pivot < max_pivots; pivot++) {
    arma::uword m = A.n_elem;
    arma::vec wA;
    double    mu;
    if (!face_solve(B.cols(A), b, cheap_face, wA, mu)) return false;

    // Anti-cycling: after 30 pivots assume the heuristic rules are cycling
    // on a degenerate face and switch to Bland's rule (always pick the
    // lowest-index violator). Pivot paths for problems that already
    // terminated within 30 pivots are untouched.
    bool bland = (pivot >= 30);

    if (wA.min() < -1e-12) {
      if (m <= 1) return false;
      arma::uword drop = wA.index_min();
      if (bland) {
        for (arma::uword j = 0; j < m; j++) {
          if (wA(j) < -1e-12) { drop = j; break; }
        }
      }
      A.shed_row(drop);
      continue;
    }

    arma::vec w(N, arma::fill::zeros);
    w(A) = arma::clamp(wA, 0.0, arma::datum::inf);
    double s = arma::accu(w);
    if (s <= 0.0) return false;
    w /= s;  // remove clamp dust (<= 1e-12 per coordinate)

    // Dual feasibility on the inactive set: lambda_i = g_i + mu >= 0
    arma::vec g   = B.t() * (B * w - b);   // == Q w - c, at O(kN)
    arma::vec lam = g + mu;
    lam.elem(A).zeros();  // active coordinates are stationary by construction
    double eps = 1e-9 * (1.0 + arma::norm(g, "inf"));
    arma::uword worst = lam.index_min();
    if (lam(worst) >= -eps) {
      x_out = w;
      return true;
    }
    if (bland) {
      for (arma::uword j = 0; j < N; j++) {
        if (lam(j) < -eps) { worst = j; break; }
      }
    }
    A.insert_rows(A.n_elem, arma::uvec{worst});
  }
  return false;
}

// Lawson-Hanson NNLS for a tiny least-squares system: min ||rhs - A x||
// subject to x >= 0, with A of size (effective QP rank + 1) x N. Finite
// termination is the entire point -- it is the classic, proven algorithm,
// so it cannot churn the way an active-set QP method can on degenerate
// faces. Used purely as a support oracle (see nnls_rescue); the
// pseudo-inverse handles rank-deficient passive blocks and the outer cap
// is a belt-and-braces bound.
static arma::vec nnls_small(const arma::mat& A, const arma::vec& rhs) {
  const arma::uword n = A.n_cols;
  arma::vec x(n, arma::fill::zeros);
  std::vector<char> inP(n, 0);
  arma::vec grad = A.t() * rhs;
  const double tol = 1e-10 * (1.0 + arma::norm(grad, "inf"));
  const int max_outer = 3 * (int)n;

  for (int it = 0; it < max_outer; it++) {
    double best = tol;
    arma::sword j_add = -1;
    for (arma::uword j = 0; j < n; j++) {
      if (!inP[j] && grad(j) > best) { best = grad(j); j_add = (arma::sword)j; }
    }
    if (j_add < 0) break;
    inP[(arma::uword)j_add] = 1;

    for (int inner = 0; inner <= (int)n; inner++) {
      std::vector<arma::uword> Pv;
      for (arma::uword j = 0; j < n; j++) if (inP[j]) Pv.push_back(j);
      if (Pv.empty()) break;
      arma::uvec P = arma::conv_to<arma::uvec>::from(Pv);
      arma::mat AP = A.cols(P);
      arma::vec z;
      if (!arma::solve(z, AP, rhs)) {
        arma::mat APp;
        if (!arma::pinv(APp, AP)) { z = arma::zeros<arma::vec>(P.n_elem); }
        else                      { z = APp * rhs; }
      }
      if (z.n_elem == P.n_elem && z.is_finite() && z.min() > 0) {
        x.zeros();
        x(P) = z;
        break;
      }
      // Standard LH backtrack: move toward z until the first passive
      // coordinate hits zero, then drop it from the passive set.
      double alpha = arma::datum::inf;
      for (arma::uword q = 0; q < P.n_elem; q++) {
        if (z(q) <= 0) {
          double xi = x(P(q));
          double d  = xi - z(q);
          if (d > 0) alpha = std::min(alpha, xi / d);
        }
      }
      if (!std::isfinite(alpha)) alpha = 0.0;
      for (arma::uword q = 0; q < P.n_elem; q++) {
        double nx = x(P(q)) + alpha * (z(q) - x(P(q)));
        if (nx <= 1e-12) { nx = 0.0; inP[P(q)] = 0; }
        x(P(q)) = nx;
      }
    }
    grad = A.t() * (rhs - A * x);
  }
  return x;
}

// NNLS-seeded rescue for simplex QPs where the warm active-set attempt
// cannot verify KKT within its pivot cap (degenerate faces make its
// pivoting churn). The factored form is already the least-squares form the
// NNLS wants, so no eigendecomposition is needed to recover it: run the
// finite-termination NNLS on (B, b) with a soft sum-to-one row appended,
// purely as a *support oracle*, and hand that support to the exact face
// solver. Only a KKT-verified solution may be returned, so exactness
// guarantees are unchanged; this merely replaces a milliseconds-scale cold
// fallback with a microseconds-scale one on the pathological problems that
// used to dominate placebo wall time.
static bool nnls_rescue(const arma::mat& B, const arma::vec& b,
                        arma::vec& x_out) {
  const arma::uword N = B.n_cols;
  // tau = largest singular value of B, from the k x k Gram.
  arma::mat BBt = B * B.t();
  arma::vec ev;
  if (!arma::eig_sym(ev, arma::symmatu(BBt))) return false;
  double ev_max = ev.max();
  if (!(ev_max > 0)) return false;
  double tau = std::sqrt(ev_max);

  arma::mat A(B.n_rows + 1, N);
  A.head_rows(B.n_rows) = B;
  A.row(B.n_rows).fill(tau);
  arma::vec rhs(B.n_rows + 1);
  rhs.head(B.n_rows) = b;
  rhs(B.n_rows)      = tau;

  arma::vec x = nnls_small(A, rhs);
  if (!x.is_finite() || !(arma::accu(x) > 0)) return false;
  arma::uvec S = arma::find(x > 1e-8 * x.max());
  if (S.n_elem == 0) return false;
  return active_set_simplex(B, b, S, x_out, 300);
}

// Simplex QP min ||b - B w||^2, warm started from x. Strategy: try the exact
// active-set solve seeded with the warm start's support; if it cannot verify
// KKT (singular subproblem, pivot cap), fall back to warm-started FISTA with
// periodic active-set retries. When `use_nnls` is set (multi-start internals
// only), a failed warm attempt first tries the NNLS support oracle before
// resorting to FISTA.
// Pure Armadillo (no Rcpp types): safe to call inside OpenMP threads.
// L must be an upper bound on lambda_max(B'B); a loose bound only shrinks
// the FISTA step size, it does not change the minimiser.
// When `wolfe` is set the Wolfe min-norm-point method replaces the warm
// active-set attempt (see wolfe_simplex); its failures fall through to the
// same NNLS/FISTA ladder, so the solver can only ever return a KKT-verified
// optimum or the same FISTA iterate the default path would return.
static bool wolfe_simplex(const arma::mat& B, const arma::vec& b,
                          arma::vec& w_out);

static arma::vec fista_simplex(const arma::mat& B, const arma::vec& b,
                               double L, arma::vec x,
                               int max_iter = 10000, double tol = 1e-6,
                               bool use_nnls = false,
                               bool cheap_face = false,
                               bool wolfe = false) {
  x = proj_simplex(x);

  arma::vec x_as;
  if (wolfe) {
    if (wolfe_simplex(B, b, x_as)) return x_as;
  } else if (active_set_simplex(B, b, arma::find(x > 1e-10), x_as, 100,
                                cheap_face)) {
    return x_as;
  }
  if (use_nnls && nnls_rescue(B, b, x_as)) {
    return x_as;
  }

  if (L < 1e-14) L = 1.0;
  double t = 1.0 / L;
  arma::vec y = x;
  double t_acc = 1.0;

  for (int iter = 0; iter < max_iter; iter++) {
    arma::vec x_prev = x;
    arma::vec grad = B.t() * (B * y - b);
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
      if (active_set_simplex(B, b, arma::find(x > 1e-10), x_as, 100,
                             cheap_face)) {
        return x_as;
      }
    }
  }
  return x;
}

// Wolfe (1976) min-norm-point method for min ||b - B w||^2 on the simplex.
//
// Bw traces the convex hull of the donor columns as w ranges over the simplex,
// so the problem is a projection onto a polytope in the k-dimensional
// predictor space -- not an N-dimensional one. By Caratheodory's theorem some
// optimal point of that hull is a convex combination of at most k+1 donors, so
// the method keeps a small affinely independent "corral" S and grows it one
// vertex at a time:
//   major cycle: the reduced costs are g = B'(B w - b); the entering donor is
//     argmin_j g_j, and min_j g_j >= g'w certifies optimality on the simplex.
//   minor cycle: solve the affine-hull problem on S exactly (face_solve, at
//     most (k+1)-dimensional); while any coefficient is non-positive, move to
//     the boundary of the simplex along the segment and drop the blocking
//     vertex.
// Termination is finite and the returned point is KKT-verified, so this is an
// exact solver, not a heuristic. Compared with the default active set -- which
// starts from the warm start's support (all N donors on a cold start) and
// sheds one coordinate per pivot -- it grows instead of shrinks, which is the
// direction the geometry favours when k is small: every face solve stays
// O(k^3) regardless of the donor pool, and the answer is Caratheodory-sparse
// rather than one arbitrary dense point of a degenerate optimal face.
// Returns false without touching w_out if the caps are hit, so the caller can
// fall back; a false return never costs correctness, only speed.
static bool wolfe_simplex(const arma::mat& B, const arma::vec& b,
                          arma::vec& w_out) {
  const arma::uword N = B.n_cols;
  const arma::uword k = B.n_rows;
  if (N == 0) return false;

  // Start at the donor closest to the target in predictor space.
  arma::vec d2 = arma::sum(arma::square(B.each_col() - b), 0).t();
  if (!d2.is_finite()) return false;
  std::vector<arma::uword> S{d2.index_min()};
  arma::vec alpha = arma::vec{1.0};

  // The corral cannot exceed k+1 affinely independent points, so the major
  // cycle is bounded by the number of vertices it can swap through; the caps
  // are generous multiples that only guard against degenerate cycling.
  const int max_major = (int)(20 * (k + 2)) + 200;
  const int max_minor = (int)k + 3;

  for (int major = 0; major < max_major; major++) {
    for (int minor = 0; minor < max_minor; minor++) {
      arma::uvec Su = arma::conv_to<arma::uvec>::from(S);
      arma::vec wA;
      double mu;
      if (!face_solve(B.cols(Su), b, false, wA, mu)) return false;
      if (!wA.is_finite()) return false;
      if (wA.min() > 0.0) {          // interior of the face: corral is valid
        alpha = wA / arma::accu(wA);
        break;
      }
      // Move from alpha toward wA until the first coordinate hits zero.
      double theta = 1.0;
      for (arma::uword i = 0; i < wA.n_elem; i++) {
        if (wA(i) < alpha(i)) {
          double denom = alpha(i) - wA(i);
          if (denom > 0.0) theta = std::min(theta, alpha(i) / denom);
        }
      }
      if (!(theta >= 0.0) || !std::isfinite(theta)) return false;
      arma::vec anew = (1.0 - theta) * alpha + theta * wA;
      std::vector<arma::uword> S2;
      std::vector<double> a2;
      for (arma::uword i = 0; i < anew.n_elem; i++) {
        if (anew(i) > 1e-14) { S2.push_back(S[i]); a2.push_back(anew(i)); }
      }
      if (S2.empty() || S2.size() == S.size()) return false;  // no progress
      S = S2;
      alpha = arma::vec(a2);
      alpha /= arma::accu(alpha);
      if (minor == max_minor - 1) return false;
    }

    arma::vec w(N, arma::fill::zeros);
    w(arma::conv_to<arma::uvec>::from(S)) = alpha;
    arma::vec g = B.t() * (B * w - b);
    double gw = arma::dot(g, w);
    arma::uword j = g.index_min();
    double eps = 1e-10 * (1.0 + arma::norm(g, "inf"));
    if (g(j) >= gw - eps) {          // KKT: no donor improves the objective
      w_out = w;
      return true;
    }
    if (std::find(S.begin(), S.end(), j) != S.end()) return false;  // cycling
    S.push_back(j);
    alpha.resize(alpha.n_elem + 1);
    alpha(alpha.n_elem - 1) = 0.0;
  }
  return false;
}

// Inner Optimization for SCM: find W given V
//' SCM Inner Weights (QP Given V)
//'
//' Solves the inner-loop QP for SCM: given a fixed diagonal metric matrix V,
//' finds donor weights W on the simplex minimising the V-weighted covariate
//' loss. The returned weights are a KKT-verified exact optimum whenever the
//' active-set solver converges, with accelerated projected gradient as a
//' fallback.
//'
//' @param X0     Covariate matrix for control units (k x N_co)
//' @param X1     Covariate vector for the treated unit (k x 1)
//' @param V_diag Diagonal of the metric matrix V (k x 1, non-negative, need not sum to 1)
//' @param wolfe  If `TRUE`, solve with the Wolfe min-norm-point method, which
//'   returns a Caratheodory-sparse optimum (at most k+1 donors carry weight)
//'   instead of one arbitrary point of a degenerate optimal face. `FALSE`
//'   (default) uses the warm-started active-set solver.
//' @return Donor weight vector W (N_co x 1) on the unit simplex
//' @export
// [[Rcpp::export]]
arma::vec scm_inner_weights_cpp(const arma::mat& X0, const arma::vec& X1,
                                const arma::vec& V_diag, bool wolfe = false) {
  arma::mat B;
  arma::vec b;
  scm_make_factor(X0, X1, V_diag, B, b);
  arma::mat BBt = B * B.t();
  double L = 1.0;
  arma::vec ev;
  if (arma::eig_sym(ev, arma::symmatu(BBt)) && ev.max() > 1e-14) L = ev.max();
  arma::vec x0 = arma::ones<arma::vec>(X0.n_cols) / X0.n_cols;
  return fista_simplex(B, b, L, x0, 10000, 1e-6, /*use_nnls=*/true,
                       /*cheap_face=*/false, wolfe);
}

// Coordinate-descent core for the nested V/W SCM optimisation. Shared by
// scm_weights_cpp (R-facing, arbitrary evaluation window) and
// scm_weights_vec_internal (placebo loop). Semantics are the classic ADH
// coordinate descent (grid {0, 0.1, ..., 1} per V coordinate, strict-improve
// accept rule, per-coordinate renormalisation of V), but each inner QP
// exploits structure:
//  * scale invariance: the QP argmin under V/sum(V) equals the argmin under
//    the unnormalised V, so per-grid-point renormalisation is skipped;
//  * factored form: changing one V coordinate rebuilds the k x N factor B in
//    O(kN) instead of touching an N x N Gram (see face_solve);
//  * warm starts: each QP starts from the previous grid point's solution
//    (active-set exact solve, FISTA fallback -- see fista_simplex);
//  * Lipschitz bound: lambda_max(B'B) = lambda_max(BB') is computed by
//    eig_sym on the k x k Gram once per sweep and tracked through the
//    single-coordinate change via lambda_max(Q + dv r r') <= lambda_max(Q) +
//    max(dv, 0) ||r||^2 (a valid upper bound, so FISTA convergence is
//    unaffected).
// V_diag must come in normalised (sum 1); it is updated in place and leaves
// normalised. Returns the best weight vector; best_loss is updated in place.
static arma::vec scm_coord_descent_core(const arma::mat& X0, const arma::vec& X1,
                                        const arma::mat& Z0_eval,
                                        const arma::vec& Z1_eval,
                                        int max_iter, double tol,
                                        arma::vec& V_diag, double& best_loss,
                                        bool use_nnls = false,
                                        bool cheap_face = false,
                                        bool wolfe = false) {
  int k = X0.n_rows;
  arma::vec best_W = scm_inner_weights_cpp(X0, X1, V_diag, wolfe);
  best_loss = arma::norm(Z1_eval - Z0_eval * best_W, 2);

  // With a single predictor the normalised V is the single point V = 1: every
  // grid candidate rescales the inner QP by a positive constant and leaves its
  // argmin unchanged, so the whole sweep evaluates one and the same QP. What
  // the sweep can still do is improve on the cold solve above, whose FISTA
  // tolerance leaves it short of the exact optimum -- so run that one QP
  // through the exact solver (same strict-improvement accept rule as the
  // sweep) and skip the ten redundant candidates.
  if (k == 1) {
    arma::mat B1;
    arma::vec b1;
    scm_make_factor(X0, X1, V_diag, B1, b1);
    arma::mat BBt1 = B1 * B1.t();
    double L1 = arma::max(arma::eig_sym(arma::symmatu(BBt1)));
    if (L1 < 1e-14) L1 = 1.0;
    arma::vec W1 = fista_simplex(B1, b1, L1, best_W, 10000, 1e-6,
                                 use_nnls, cheap_face, wolfe);
    double loss1 = arma::norm(Z1_eval - Z0_eval * W1, 2);
    if (loss1 < best_loss) {
      best_loss = loss1;
      best_W    = W1;
    }
    return best_W;
  }

  arma::vec row_n2 = arma::sum(arma::square(X0), 1);       // k: ||r_j||^2
  arma::vec grid   = arma::linspace<arma::vec>(0.0, 1.0, 11);

  arma::mat B_base;
  arma::vec b_base;

  for (int iter = 0; iter < max_iter; iter++) {
    // Exact Lipschitz constant once per sweep, from the k x k Gram; the
    // single-coordinate bounds accumulate slack within the sweep only.
    scm_make_factor(X0, X1, V_diag, B_base, b_base);
    arma::mat BBt = B_base * B_base.t();
    double L_base = arma::max(arma::eig_sym(arma::symmatu(BBt)));
    if (L_base < 1e-14) L_base = 1.0;

    double max_change = 0.0;
    for (int j = 0; j < k; j++) {
      double orig_v          = V_diag(j);
      double local_best_v    = orig_v;
      double local_best_loss = best_loss;
      arma::vec local_best_W = best_W;
      arma::vec warm = best_W;

      for (arma::uword g = 0; g < grid.n_elem; g++) {
        if (std::abs(grid(g) - orig_v) < 1e-15) continue;  // current point: ties only
        double vsum = arma::sum(V_diag) - orig_v + grid(g);
        if (vsum < 1e-14) continue;
        double dv = grid(g) - orig_v;
        arma::vec V_g = V_diag;
        V_g(j) = grid(g);
        arma::mat B_g;
        arma::vec b_g;
        scm_make_factor(X0, X1, V_g, B_g, b_g);
        double L_g = L_base + std::max(dv, 0.0) * row_n2(j);
        arma::vec cand_W = fista_simplex(B_g, b_g, L_g, warm, 10000, 1e-6,
                                         use_nnls, cheap_face, wolfe);
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
        L_base += std::max(dv_acc, 0.0) * row_n2(j);
        V_diag(j) = local_best_v;
        double vsum_after = arma::sum(V_diag);
        if (vsum_after > 1e-14) {
          V_diag /= vsum_after;
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

// Warm-started outer-loss evaluator: one inner QP per candidate V, seeded
// with the previous evaluation's solution. Consecutive candidates in a
// screen or Nelder-Mead run have similar V, so the active-set solve usually
// verifies KKT within a few pivots -- orders of magnitude cheaper than a
// cold FISTA solve, which is what makes a thorough multi-start search
// affordable. Pure Armadillo state (OpenMP-safe).
struct scm_outer_eval {
  const arma::mat& X0;
  const arma::vec& X1;
  const arma::mat& Z0_eval;
  const arma::vec& Z1_eval;
  arma::vec warm;    // persistent warm start across evaluations
  bool use_nnls;     // multi-start internals only (see nnls_rescue)
  bool wolfe;        // Caratheodory-sparse inner solver (opt-in)

  scm_outer_eval(const arma::mat& X0_, const arma::vec& X1_,
                 const arma::mat& Z0_, const arma::vec& Z1_,
                 bool use_nnls_ = false, bool wolfe_ = false)
    : X0(X0_), X1(X1_), Z0_eval(Z0_), Z1_eval(Z1_), use_nnls(use_nnls_),
      wolfe(wolfe_) {
    warm = arma::ones<arma::vec>(X0.n_cols) / X0.n_cols;
  }

  double operator()(const arma::vec& V_raw, arma::vec* W_out = nullptr) {
    arma::vec V = arma::clamp(V_raw, 0.0, arma::datum::inf);
    double s = arma::accu(V);
    if (s < 1e-14 || !V.is_finite()) return arma::datum::inf;
    V /= s;
    arma::mat B;
    arma::vec b;
    scm_make_factor(X0, X1, V, B, b);
    // Only a KKT-verified exact solve may report a loss: an under-converged
    // W would make the outer objective look flat (every V returns roughly
    // the warm-start's loss) and send the refiner drifting. The active-set
    // solve seeded with the previous support is exact when it succeeds;
    // when it cannot verify KKT (degenerate faces), the NNLS support oracle
    // is a finite-termination rescue that lands the exact face solve
    // without the milliseconds-scale cold solve. The cold reference solver
    // is the last resort.
    arma::vec W;
    bool got = wolfe
      ? wolfe_simplex(B, b, W)
      : active_set_simplex(B, b, arma::find(warm > 1e-10), W, 300);
    if (!got && !(use_nnls && nnls_rescue(B, b, W))) {
      W = scm_inner_weights_cpp(X0, X1, V, wolfe);
    }
    if (!W.is_finite()) return arma::datum::inf;
    warm = W;
    double loss = arma::norm(Z1_eval - Z0_eval * W, 2);
    if (!std::isfinite(loss)) return arma::datum::inf;
    if (W_out) *W_out = W;
    return loss;
  }
};

// Nelder-Mead refinement of the outer V objective on the raw non-negative
// parametrisation (the evaluator normalises internally, so the search is
// scale-free and simplex constraints need no explicit handling). Standard
// reflection/expansion/contraction/shrink; fully deterministic. Each
// function evaluation costs one warm-started inner QP.
static arma::vec nm_refine_v(scm_outer_eval& eval,
                             const arma::vec& v_init, double& f_out,
                             int max_eval = 350, double reltol = 1e-8) {
  const int k = v_init.n_elem;
  const double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;

  std::vector<arma::vec> pts(k + 1);
  arma::vec fv(k + 1);
  pts[0] = v_init;
  fv(0)  = eval(pts[0]);
  for (int i = 0; i < k; i++) {
    arma::vec p = v_init;
    p(i) += 0.1;  // 10% of the (normalised) total mass
    pts[i + 1] = p;
    fv(i + 1)  = eval(p);
  }
  int n_eval = k + 1;

  while (n_eval < max_eval) {
    arma::uvec ord = arma::sort_index(fv);
    if (fv(ord(k)) - fv(ord(0)) <=
        reltol * (std::abs(fv(ord(0))) + reltol)) break;

    arma::uword i_best = ord(0), i_worst = ord(k), i_second = ord(k - 1);
    arma::vec centroid(k, arma::fill::zeros);
    for (int i = 0; i <= k; i++) {
      if ((arma::uword)i != i_worst) centroid += pts[i];
    }
    centroid /= k;

    arma::vec x_r = centroid + alpha * (centroid - pts[i_worst]);
    double f_r = eval(x_r);
    n_eval++;

    if (f_r < fv(i_best)) {
      arma::vec x_e = centroid + gamma * (x_r - centroid);
      double f_e = eval(x_e);
      n_eval++;
      if (f_e < f_r) { pts[i_worst] = x_e; fv(i_worst) = f_e; }
      else           { pts[i_worst] = x_r; fv(i_worst) = f_r; }
    } else if (f_r < fv(i_second)) {
      pts[i_worst] = x_r; fv(i_worst) = f_r;
    } else {
      arma::vec x_c = centroid + rho * (pts[i_worst] - centroid);
      double f_c = eval(x_c);
      n_eval++;
      if (f_c < fv(i_worst)) {
        pts[i_worst] = x_c; fv(i_worst) = f_c;
      } else {
        for (int i = 0; i <= k; i++) {
          if ((arma::uword)i == i_best) continue;
          pts[i] = pts[i_best] + sigma * (pts[i] - pts[i_best]);
          fv(i)  = eval(pts[i]);
          n_eval++;
        }
      }
    }
  }

  arma::uword i_min = fv.index_min();
  f_out = fv(i_min);
  return pts[i_min];
}

// ---- Multi-start machinery -------------------------------------------------
// The nested V/W objective is non-convex and a single start can settle in a
// poor basin (2x worse pre-period SSR than reachable on the Prop99 predictor
// spec). The multi-start search below is fully deterministic: no R RNG state
// is touched, and everything is pure Armadillo (safe inside OpenMP threads).

struct scm_ms_result {
  arma::vec V, W;
  double loss = arma::datum::inf;
};

// Track the best (V, W, loss) triple. Full triples are required: coordinate
// descent's reported W is the warm-started QP solution recorded during its
// sweep, which is not exactly reproducible from V alone, so re-deriving W
// from the winning V would break the never-worse guarantee.
static void scm_ms_consider(scm_ms_result& best, const arma::vec& V_raw,
                            const arma::vec& W, double loss) {
  if (!std::isfinite(loss) || loss >= best.loss) return;
  arma::vec V = arma::clamp(V_raw, 0.0, arma::datum::inf);
  double s = arma::accu(V);
  if (s < 1e-14 || !W.is_finite()) return;
  best.loss = loss;
  best.V    = V / s;
  best.W    = W;
}

// Fixed start set: uniform, k smoothed one-hots, and 100 Dirichlet(1) draws
// from a fixed-seed mt19937 (uniforms transformed via -log(u), so the stream
// is identical across platforms). Depends only on k, so callers fitting many
// problems with the same predictor count can share one matrix.
static arma::mat scm_ms_start_matrix(int k) {
  const int n_rand = 100;
  const int n_starts = 1 + k + n_rand;
  arma::mat starts(k, n_starts);
  starts.col(0).fill(1.0 / k);
  for (int j = 0; j < k; j++) {
    if (k > 1) {
      starts.col(1 + j).fill(0.1 / (k - 1));
      starts(j, 1 + j) = 0.9;
    } else {
      starts(0, 1 + j) = 1.0;
    }
  }
  std::mt19937 rng(20260719u);
  const double inv_max = 1.0 / (double(std::mt19937::max()) + 1.0);
  for (int j = 0; j < n_rand; j++) {
    arma::vec g(k);
    for (int i = 0; i < k; i++) {
      double u = (double(rng()) + 0.5) * inv_max;
      g(i) = -std::log(u);
    }
    starts.col(1 + k + j) = g / arma::accu(g);
  }
  return starts;
}

// Screen every start at one warm-started inner QP each and pick the
// candidates to refine: the uniform start (never-worse-than-single-start
// guarantee) plus the n_top best screened starts.
static std::vector<arma::uword> scm_ms_select(const arma::mat& starts,
                                              scm_outer_eval& eval,
                                              int n_top) {
  const int n_starts = (int)starts.n_cols;
  arma::vec screen(n_starts);
  for (int j = 0; j < n_starts; j++) {
    screen(j) = eval(starts.col(j));
  }
  arma::uvec ord = arma::sort_index(screen);
  std::vector<arma::uword> cand;
  cand.push_back(0);
  for (arma::uword q = 0; q < ord.n_elem; q++) {
    if ((int)cand.size() >= n_top + 1) break;
    if (ord(q) != 0 && std::isfinite(screen(ord(q)))) cand.push_back(ord(q));
  }
  return cand;
}

// One candidate's refinement: the exact single-start path (coordinate
// descent from V = 1/k) when the candidate is the uniform start, then the
// Nelder-Mead -> coordinate-descent -> Nelder-Mead pipeline for all: NM
// makes continuous moves between the 0.1-spaced grid points, the grid
// sweep escapes flat NM valleys via large single-coordinate jumps. The
// best triple seen anywhere in the pipeline is returned.
static scm_ms_result scm_ms_pipeline(const arma::mat& X0, const arma::vec& X1,
                                     const arma::mat& Z0_eval,
                                     const arma::vec& Z1_eval,
                                     int max_iter, double tol,
                                     const arma::vec& start,
                                     bool is_uniform, bool wolfe = false) {
  scm_ms_result best;
  scm_outer_eval eval_q(X0, X1, Z0_eval, Z1_eval, /*use_nnls=*/true, wolfe);

  if (is_uniform) {
    // The uniform candidate's plain coordinate descent is the never-worse
    // reference: it must reproduce v_optim = "coord_descent" exactly, so it
    // keeps the historical solver path (use_nnls stays false).
    arma::vec V = start;
    double loss = 0.0;
    arma::vec W = scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                                         max_iter, tol, V, loss,
                                         /*use_nnls=*/false,
                                         /*cheap_face=*/false, wolfe);
    scm_ms_consider(best, V, W, loss);
  }

  arma::vec W_stage;
  double f1 = 0.0;
  arma::vec V1 = nm_refine_v(eval_q, start, f1);
  double f1x = eval_q(V1, &W_stage);
  scm_ms_consider(best, V1, W_stage, f1x);

  arma::vec V2 = arma::clamp(V1, 0.0, arma::datum::inf);
  double s2 = arma::accu(V2);
  if (s2 >= 1e-14) {
    V2 /= s2;
    double f2 = 0.0;
    arma::vec W2 = scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                                          max_iter, tol, V2, f2,
                                          /*use_nnls=*/true,
                                          /*cheap_face=*/false, wolfe);
    scm_ms_consider(best, V2, W2, f2);

    double f3 = 0.0;
    arma::vec V3 = nm_refine_v(eval_q, V2, f3);
    double f3x = eval_q(V3, &W_stage);
    scm_ms_consider(best, V3, W_stage, f3x);
  }
  return best;
}

// Deterministic multi-start driver for a single fit: screen, refine the
// uniform + top-4 candidates in parallel (suppressed inside an enclosing
// parallel region -- nested teams oversubscribe), reduce serially in fixed
// candidate order with a strict '<' so the result is deterministic
// regardless of thread scheduling.
static arma::vec scm_multistart_core(const arma::mat& X0, const arma::vec& X1,
                                     const arma::mat& Z0_eval,
                                     const arma::vec& Z1_eval,
                                     int max_iter, double tol,
                                     arma::vec& V_out, double& best_loss,
                                     bool wolfe = false) {
  const int k = X0.n_rows;
  const int n_top = 4;

  // With a single predictor there is only one normalised V, so every start in
  // the set is the same point and no refinement can move it: the search
  // degenerates to the single inner QP that coordinate descent already
  // short-circuits to.
  if (k == 1) {
    V_out = arma::ones<arma::vec>(1);
    return scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval, max_iter, tol,
                                  V_out, best_loss, /*use_nnls=*/false,
                                  /*cheap_face=*/false, wolfe);
  }

  arma::mat starts = scm_ms_start_matrix(k);

  scm_outer_eval eval(X0, X1, Z0_eval, Z1_eval, /*use_nnls=*/true, wolfe);
  std::vector<arma::uword> cand = scm_ms_select(starts, eval, n_top);

  const int n_cand = (int)cand.size();
  std::vector<scm_ms_result> slots(n_cand);

#ifdef _OPENMP
  const bool spawn_team = !omp_in_parallel();
#else
  const bool spawn_team = false;
#endif

  #pragma omp parallel for schedule(dynamic, 1) if(spawn_team)
  for (int qi = 0; qi < n_cand; qi++) {
    try {
      slots[qi] = scm_ms_pipeline(X0, X1, Z0_eval, Z1_eval, max_iter, tol,
                                  starts.col(cand[qi]), cand[qi] == 0, wolfe);
    } catch (...) {
      // Exceptions must not escape the OpenMP region; a failed pipeline
      // simply contributes nothing (its loss stays Inf).
    }
  }

  best_loss = arma::datum::inf;
  arma::vec best_W, best_V;
  for (int qi = 0; qi < n_cand; qi++) {
    if (std::isfinite(slots[qi].loss) && slots[qi].loss < best_loss &&
        !slots[qi].V.is_empty() && slots[qi].W.is_finite()) {
      best_loss = slots[qi].loss;
      best_V    = slots[qi].V;
      best_W    = slots[qi].W;
    }
  }

  if (best_V.is_empty()) {
    // All pipeline stages failed (pathological inputs): fall back to the
    // plain single-start path.
    best_V = starts.col(0);
    best_W = scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                                    max_iter, tol, best_V, best_loss,
                                    /*use_nnls=*/false, /*cheap_face=*/false,
                                    wolfe);
  }
  V_out = best_V;
  return best_W;
}

// Flat (donor x candidate) multi-start batch for the covariate placebo
// loop; called from scm_placebo_x_cpp via forward declaration. Donor-level
// parallelism alone lets one slow donor pin a single thread while the
// others go idle after finishing (per-donor refit cost spans two orders of
// magnitude on real data); flattening the task list lets every candidate
// pipeline of every donor be scheduled independently, with no nested teams
// anywhere. Per-donor results are identical to running
// scm_weights_vec_internal(multistart = TRUE) donor by donor: same screen,
// same pipelines, same fixed-order reduction. ok_out[d] = 0 marks a donor
// whose refit failed entirely (the caller maps it to NaN, matching the
// existing solver-failure contract).
void scm_multistart_batch(const std::vector<arma::mat>& X0d,
                          const std::vector<arma::vec>& X1d,
                          const std::vector<arma::mat>& Z0d,
                          const std::vector<arma::vec>& Z1d,
                          int max_iter, double tol,
                          std::vector<arma::vec>& W_out,
                          std::vector<char>& ok_out, bool wolfe) {
  const int N = (int)X0d.size();
  const int n_top = 4;
  W_out.assign(N, arma::vec());
  ok_out.assign(N, 0);
  if (N == 0) return;

  // Single predictor: the multi-start search is vacuous (see
  // scm_multistart_core), so every donor reduces to one inner QP.
  if (X0d[0].n_rows == 1) {
    #pragma omp parallel for schedule(dynamic, 1)
    for (int d = 0; d < N; d++) {
      try {
        arma::vec V = arma::ones<arma::vec>(1);
        double loss = 0.0;
        W_out[d] = scm_coord_descent_core(X0d[d], X1d[d], Z0d[d], Z1d[d],
                                          max_iter, tol, V, loss,
                                          /*use_nnls=*/false,
                                          /*cheap_face=*/false, wolfe);
        ok_out[d] = 1;
      } catch (...) {}
    }
    return;
  }

  arma::mat starts = scm_ms_start_matrix((int)X0d[0].n_rows);

  // Phase A: per-donor screening (cheap: one warm QP per start). Screening
  // wall time is also recorded as a cost proxy -- donors whose evaluations
  // are expensive (degenerate faces, cold fallbacks) are expensive in the
  // refinement pipelines too.
  std::vector<std::vector<arma::uword>> cand(N);
  std::vector<char>   screened(N, 0);
  std::vector<double> screen_cost(N, 0.0);
  #pragma omp parallel for schedule(dynamic, 1)
  for (int d = 0; d < N; d++) {
    try {
      auto t0 = std::chrono::steady_clock::now();
      scm_outer_eval ev(X0d[d], X1d[d], Z0d[d], Z1d[d], /*use_nnls=*/true,
                        wolfe);
      cand[d] = scm_ms_select(starts, ev, n_top);
      screened[d] = 1;
      screen_cost[d] = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    } catch (...) {}
  }

  // Phase B: one flat task per (donor, candidate), dispatched
  // heaviest-donor-first (LPT heuristic) so a slow donor's pipelines start
  // immediately instead of extending the tail of the schedule. Execution
  // order only affects packing: results are reduced by (donor, candidate)
  // index in Phase C, so this is exactly result-neutral.
  std::vector<int> donor_order;
  for (int d = 0; d < N; d++) {
    if (screened[d]) donor_order.push_back(d);
  }
  std::stable_sort(donor_order.begin(), donor_order.end(),
                   [&](int a, int b) { return screen_cost[a] > screen_cost[b]; });

  std::vector<int>         task_d;
  std::vector<arma::uword> task_start;
  for (size_t o = 0; o < donor_order.size(); o++) {
    const int d = donor_order[o];
    for (size_t q = 0; q < cand[d].size(); q++) {
      task_d.push_back(d);
      task_start.push_back(cand[d][q]);
    }
  }
  const int n_task = (int)task_d.size();
  std::vector<scm_ms_result> res(n_task);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int t = 0; t < n_task; t++) {
    try {
      const int d = task_d[t];
      res[t] = scm_ms_pipeline(X0d[d], X1d[d], Z0d[d], Z1d[d],
                               max_iter, tol,
                               starts.col(task_start[t]),
                               task_start[t] == 0, wolfe);
    } catch (...) {}
  }

  // Phase C: per-donor reduction in task order (== candidate order), with
  // the same strict '<' and the same single-start fallback as the
  // single-fit driver.
  std::vector<scm_ms_result> best(N);
  for (int t = 0; t < n_task; t++) {
    const int d = task_d[t];
    if (std::isfinite(res[t].loss) && res[t].loss < best[d].loss &&
        !res[t].V.is_empty() && res[t].W.is_finite()) {
      best[d] = res[t];
    }
  }
  for (int d = 0; d < N; d++) {
    if (!best[d].V.is_empty()) {
      W_out[d]  = best[d].W;
      ok_out[d] = 1;
    } else if (screened[d]) {
      try {
        arma::vec V = starts.col(0);
        double loss = 0.0;
        W_out[d] = scm_coord_descent_core(X0d[d], X1d[d], Z0d[d], Z1d[d],
                                          max_iter, tol, V, loss,
                                          /*use_nnls=*/false,
                                          /*cheap_face=*/false, wolfe);
        ok_out[d] = 1;
      } catch (...) {}
    }
  }
}

// Internal helper: returns only the weight vector (no Rcpp::List overhead).
// Called by scm_placebo_cpp / scm_placebo_x_cpp in inference.cpp via forward
// declaration. `multistart` selects the deterministic multi-start driver so
// placebo refits stay symmetric with a multi-start treated fit.
arma::vec scm_weights_vec_internal(const arma::mat& X0, const arma::vec& X1,
                                    const arma::mat& Z0, const arma::vec& Z1,
                                    int max_iter, double tol, bool multistart,
                                    bool cheap_face, bool wolfe) {
  int k = X0.n_rows;
  arma::vec V_diag = arma::ones<arma::vec>(k) / k;
  double best_loss = 0.0;
  if (multistart) {
    return scm_multistart_core(X0, X1, Z0, Z1, max_iter, tol,
                               V_diag, best_loss, wolfe);
  }
  return scm_coord_descent_core(X0, X1, Z0, Z1, max_iter, tol,
                                V_diag, best_loss, /*use_nnls=*/false,
                                cheap_face, wolfe);
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
//' @param z_rows Optional 1-based row indices of Z defining the evaluation
//'   window for the outer V optimisation (the `v_window` argument of
//'   [scm_fit()]). `NULL` (default) evaluates on the full Z window. Takes
//'   precedence over `t_train`; the reported loss always uses the full Z
//'   window.
//' @param multistart If `TRUE`, the outer V optimisation runs a
//'   deterministic multi-start search (screened start set, coordinate-descent
//'   polish, Nelder-Mead refinement) instead of a single coordinate-descent
//'   pass from the uniform V. The result is never worse (in outer loss) than
//'   the single-start path.
//' @param cheap_face If `TRUE`, the inner simplex QP first tries a cheap
//'   bordered-KKT direct solve on each active-set face, falling back to the
//'   scale-robust null-space solve when that system is singular. Worth trying
//'   only in the outcomes-only regime (no user predictors, no predictor
//'   rescaling), where V is dense and the faces that carry the fit are well
//'   conditioned; it reproduces the null-space solution to round-off.
//'   `FALSE` (default) goes straight to the null-space solve, which is what
//'   predictor fits need for scale invariance and rank-deficient faces.
//' @param wolfe If `TRUE`, every inner QP is solved with the Wolfe
//'   min-norm-point method, which returns a Caratheodory-sparse optimum (at
//'   most k+1 donors carry weight) instead of one arbitrary point of a
//'   degenerate optimal face. `FALSE` (default) uses the warm-started
//'   active-set solver.
//' @return A list with:
//'   * `W`: Donor weight vector (N_co x 1) on the unit simplex
//'   * `V`: Optimal metric diagonal (k x 1, normalised to sum to 1)
//'   * `loss`: Final pre-treatment prediction loss (full pre-treatment window)
//' @export
// [[Rcpp::export]]
Rcpp::List scm_weights_cpp(const arma::mat& X0, const arma::vec& X1,
                            const arma::mat& Z0, const arma::vec& Z1,
                            int max_iter = 100, double tol = 1e-4,
                            int t_train = -1,
                            Rcpp::Nullable<Rcpp::IntegerVector> z_rows = R_NilValue,
                            bool multistart = false,
                            bool cheap_face = false,
                            bool wolfe = false) {
  int k     = X0.n_rows;
  int T_pre = (int)Z0.n_rows;

  bool has_window = z_rows.isNotNull();
  bool do_oos     = !has_window && (t_train > 0 && t_train < T_pre);

  // The evaluation window restricts only the outer V-selection loss; W is
  // always fitted on the full predictor matrix X. This avoids the dimension
  // mismatch that arises when X and Z have the same number of rows
  // (outcomes-only case) and V has k=T_pre entries.
  arma::mat Z0_eval;
  arma::vec Z1_eval;

  if (has_window) {
    Rcpp::IntegerVector zr(z_rows);
    arma::uvec rows(zr.size());
    for (int i = 0; i < zr.size(); i++) rows(i) = (arma::uword)(zr[i] - 1);
    Z0_eval = Z0.rows(rows);
    Z1_eval = Z1(rows);
  } else if (do_oos) {
    // Validation window: rows t_train..(T_pre-1) of Z
    Z0_eval = Z0.rows(t_train, T_pre - 1);
    Z1_eval = Z1.subvec(t_train, T_pre - 1);
  } else {
    Z0_eval = Z0;
    Z1_eval = Z1;
  }

  arma::vec V_diag = arma::ones(k) / k;

  // Outer optimisation over V; W always fitted on the full X (see
  // scm_coord_descent_core / scm_multistart_core for the inner-QP machinery).
  double best_loss = 0.0;
  arma::vec best_W = multistart
    ? scm_multistart_core(X0, X1, Z0_eval, Z1_eval,
                          max_iter, tol, V_diag, best_loss, wolfe)
    : scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                             max_iter, tol, V_diag, best_loss,
                             /*use_nnls=*/false, cheap_face, wolfe);

  // OOS mode: refit W on full pre-treatment data with selected V*
  arma::vec final_W = do_oos ? scm_inner_weights_cpp(X0, X1, V_diag, wolfe)
                             : best_W;
  // Reported loss always covers the full Z window so pre-fit quality stays
  // comparable across evaluation-window choices.
  double final_loss = arma::norm(Z1 - Z0 * final_W, 2);

  return Rcpp::List::create(
    Rcpp::Named("W")    = final_W,
    Rcpp::Named("V")    = V_diag,
    Rcpp::Named("loss") = final_loss
  );
}
