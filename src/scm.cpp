#include <RcppArmadillo.h>
#include <random>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

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
    // Stationarity on the active set: G w_A + mu * 1 = c_A, 1'w_A = 1,
    // with G = Q_AA + dv * rA rA'. Solved in the null space of the sum
    // constraint (w = w0 + Z u, Z an orthonormal basis of {z : 1'z = 0})
    // rather than via the bordered KKT system: G is rank-deficient
    // whenever V has few effective predictors (k < |A|, e.g. one-hot V
    // gives rank-1 Q), and a raw KKT solve then sits on a success/failure
    // knife edge that ulp-scale input perturbations (rescaled predictor
    // units) can flip, breaking the scale-invariance of the fit. The
    // reduced problem's H = Z'GZ is PSD, so a rank-truncated
    // pseudo-inverse gives a true face minimiser whenever one exists --
    // identical to the exact solve when H is well conditioned. The
    // truncation cutoff sits well above machine noise (curvature below
    // 1e-9 of the matrix scale carries no statistical information), and
    // the dual-feasibility check below still gates every result.
    arma::vec rA = r(A);
    arma::vec cA = c(A);
    arma::mat G  = Q.submat(A, A) + dv * (rA * rA.t());

    arma::vec wA;
    double    mu;
    if (m == 1) {
      wA = arma::vec{1.0};
      mu = cA(0) - G(0, 0);
    } else {
      arma::mat Z_ns = sum_null_basis(m);
      arma::vec w0(m);
      w0.fill(1.0 / m);
      arma::mat H = Z_ns.t() * G * Z_ns;
      arma::vec b = Z_ns.t() * (cA - G * w0);

      // Fast path: Cholesky solve when H is comfortably well conditioned,
      // judged by the diagonal ratio of the factor (a scale-robust rcond
      // proxy -- ulp perturbations cannot flip it, unlike the internal
      // failure threshold of a plain linear solve). Degenerate path:
      // eigendecomposition with rank truncation well above machine noise
      // (curvature below 1e-9 of the leading eigenvalue carries no
      // statistical information); the truncated inverse is a true face
      // minimiser whenever one exists. At the branch boundary the two
      // routes agree, so the branch itself is benign.
      arma::vec u;
      arma::mat R;
      bool fast = arma::chol(R, arma::symmatu(H));
      if (fast) {
        double dmin = R.diag().min();
        double dmax = R.diag().max();
        fast = (dmin > 0.0) && (dmin >= 3.2e-5 * dmax);  // rcond(H) >~ 1e-9
        if (fast) {
          u = arma::solve(arma::trimatu(R),
                          arma::solve(arma::trimatl(R.t()), b));
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
          arma::vec bt = evec.t() * b;
          for (arma::uword i = 0; i < ev.n_elem; i++) {
            bt(i) = (ev(i) > cutoff) ? bt(i) / ev(i) : 0.0;
          }
          u = evec * bt;
        }
      }

      wA = w0 + Z_ns * u;
      if (!wA.is_finite()) return false;
      mu = -arma::mean(G * wA - cA);
    }

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
    arma::vec g   = Q * w + (dv * arma::dot(r, w)) * r - c;
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
  arma::vec row_n2;  // ||r_j||^2 per predictor row (Lipschitz bound pieces)
  arma::vec r_zero;  // zero rank-1 term: plain QP through fista_simplex_rank1
  arma::vec warm;    // persistent warm start across evaluations

  scm_outer_eval(const arma::mat& X0_, const arma::vec& X1_,
                 const arma::mat& Z0_, const arma::vec& Z1_)
    : X0(X0_), X1(X1_), Z0_eval(Z0_), Z1_eval(Z1_) {
    row_n2 = arma::sum(arma::square(X0), 1);
    r_zero = arma::zeros<arma::vec>(X0.n_cols);
    warm   = arma::ones<arma::vec>(X0.n_cols) / X0.n_cols;
  }

  double operator()(const arma::vec& V_raw, arma::vec* W_out = nullptr) {
    arma::vec V = arma::clamp(V_raw, 0.0, arma::datum::inf);
    double s = arma::accu(V);
    if (s < 1e-14 || !V.is_finite()) return arma::datum::inf;
    V /= s;
    arma::mat Q = X0.t() * arma::diagmat(V) * X0;
    arma::vec c = X0.t() * (V % X1);
    // Only a KKT-verified exact solve may report a loss: an under-converged
    // W would make the outer objective look flat (every V returns roughly
    // the warm-start's loss) and send the refiner drifting. The active-set
    // solve seeded with the previous support is exact when it succeeds;
    // otherwise fall back to the cold reference solver.
    arma::vec W;
    if (!active_set_simplex_rank1(Q, r_zero, 0.0, c,
                                  arma::find(warm > 1e-10), W, 300)) {
      W = scm_inner_weights_cpp(X0, X1, V);
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

// Deterministic multi-start driver for the outer V problem. The nested
// V/W objective is non-convex and a single start can settle in a poor basin
// (2x worse pre-period SSR than reachable on the Prop99 predictor spec), so:
//  1. screen a fixed start set -- uniform, k smoothed one-hots, and 100
//     Dirichlet(1) draws from a fixed-seed mt19937 (uniforms transformed via
//     -log(u), so the stream is identical across platforms) -- at one
//     warm-started inner QP each;
//  2. take the uniform start plus the three best screened starts as
//     candidates (keeping uniform guarantees the result is never worse than
//     the single-start path), and run each through a
//     Nelder-Mead -> coordinate-descent -> Nelder-Mead pipeline: NM makes
//     continuous moves between the 0.1-spaced grid points, the grid sweep
//     escapes flat NM valleys via large single-coordinate jumps;
//  3. keep the best (V, W) seen anywhere in the pipeline.
// Everything is deterministic; no R RNG state is touched. Pure Armadillo,
// safe inside OpenMP threads.
static arma::vec scm_multistart_core(const arma::mat& X0, const arma::vec& X1,
                                     const arma::mat& Z0_eval,
                                     const arma::vec& Z1_eval,
                                     int max_iter, double tol,
                                     arma::vec& V_out, double& best_loss) {
  const int k = X0.n_rows;
  const int n_rand = 100, n_top = 4;
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

  scm_outer_eval eval(X0, X1, Z0_eval, Z1_eval);

  arma::vec screen(n_starts);
  for (int j = 0; j < n_starts; j++) {
    screen(j) = eval(starts.col(j));
  }

  arma::uvec ord = arma::sort_index(screen);
  std::vector<arma::uword> cand;
  cand.push_back(0);  // uniform start: never-worse-than-single-start guarantee
  for (arma::uword q = 0; q < ord.n_elem; q++) {
    if ((int)cand.size() >= n_top + 1) break;
    if (ord(q) != 0 && std::isfinite(screen(ord(q)))) cand.push_back(ord(q));
  }

  best_loss = arma::datum::inf;
  arma::vec best_W, best_V;

  // Candidate pipelines are independent, so run them in parallel (inside
  // the placebo loop's existing parallel region OpenMP serialises this
  // nested level, so there is no oversubscription). Each pipeline records
  // its own best (V, W, loss) triple -- coordinate descent's reported W is
  // the warm-started QP solution recorded during its sweep, which is not
  // exactly reproducible from V alone, so re-deriving W at the end would
  // break the never-worse guarantee. The final reduction runs serially in
  // fixed candidate order with a strict '<', keeping the result
  // deterministic regardless of thread scheduling.
  const int n_cand = (int)cand.size();
  std::vector<arma::vec> slot_V(n_cand), slot_W(n_cand);
  std::vector<double>    slot_loss(n_cand, arma::datum::inf);

  // Never spawn a nested team inside the placebo loop's donor-level
  // parallel region: some runtimes create 1-thread teams (overhead) or
  // oversubscribe (38 donors x 4 candidates), which measurably slows the
  // placebo batch.
#ifdef _OPENMP
  const bool spawn_team = !omp_in_parallel();
#else
  const bool spawn_team = false;
#endif

  #pragma omp parallel for schedule(dynamic, 1) if(spawn_team)
  for (int qi = 0; qi < n_cand; qi++) {
    try {
      scm_outer_eval eval_q(X0, X1, Z0_eval, Z1_eval);
      double    loc_loss = arma::datum::inf;
      arma::vec loc_V, loc_W;
      auto keep = [&](const arma::vec& V_raw, const arma::vec& W,
                      double loss) {
        if (!std::isfinite(loss) || loss >= loc_loss) return;
        arma::vec V = arma::clamp(V_raw, 0.0, arma::datum::inf);
        double s = arma::accu(V);
        if (s < 1e-14 || !W.is_finite()) return;
        loc_loss = loss;
        loc_V    = V / s;
        loc_W    = W;
      };

      // The exact single-start path for the uniform candidate (coordinate
      // descent from V = 1/k), the NM -> coord -> NM pipeline for all.
      if (cand[qi] == 0) {
        arma::vec V = starts.col(0);
        double loss = 0.0;
        arma::vec W = scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                                             max_iter, tol, V, loss);
        keep(V, W, loss);
      }

      arma::vec W_stage;
      double f1 = 0.0;
      arma::vec V1 = nm_refine_v(eval_q, starts.col(cand[qi]), f1);
      double f1x = eval_q(V1, &W_stage);
      keep(V1, W_stage, f1x);

      arma::vec V2 = arma::clamp(V1, 0.0, arma::datum::inf);
      double s2 = arma::accu(V2);
      if (s2 >= 1e-14) {
        V2 /= s2;
        double f2 = 0.0;
        arma::vec W2 = scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                                              max_iter, tol, V2, f2);
        keep(V2, W2, f2);

        double f3 = 0.0;
        arma::vec V3 = nm_refine_v(eval_q, V2, f3);
        double f3x = eval_q(V3, &W_stage);
        keep(V3, W_stage, f3x);
      }

      slot_V[qi]    = loc_V;
      slot_W[qi]    = loc_W;
      slot_loss[qi] = loc_loss;
    } catch (...) {
      // Exceptions must not escape the OpenMP region; a failed pipeline
      // simply contributes nothing (slot_loss stays Inf).
    }
  }

  for (int qi = 0; qi < n_cand; qi++) {
    if (std::isfinite(slot_loss[qi]) && slot_loss[qi] < best_loss &&
        !slot_V[qi].is_empty() && slot_W[qi].is_finite()) {
      best_loss = slot_loss[qi];
      best_V    = slot_V[qi];
      best_W    = slot_W[qi];
    }
  }

  if (best_V.is_empty()) {
    // All pipeline stages failed (pathological inputs): fall back to the
    // plain single-start path.
    best_V = starts.col(0);
    best_W = scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                                    max_iter, tol, best_V, best_loss);
  }
  V_out = best_V;
  return best_W;
}

// Internal helper: returns only the weight vector (no Rcpp::List overhead).
// Called by scm_placebo_cpp / scm_placebo_x_cpp in inference.cpp via forward
// declaration. `multistart` selects the deterministic multi-start driver so
// placebo refits stay symmetric with a multi-start treated fit.
arma::vec scm_weights_vec_internal(const arma::mat& X0, const arma::vec& X1,
                                    const arma::mat& Z0, const arma::vec& Z1,
                                    int max_iter, double tol, bool multistart) {
  int k = X0.n_rows;
  arma::vec V_diag = arma::ones<arma::vec>(k) / k;
  double best_loss = 0.0;
  if (multistart) {
    return scm_multistart_core(X0, X1, Z0, Z1, max_iter, tol,
                               V_diag, best_loss);
  }
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
                            bool multistart = false) {
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
                          max_iter, tol, V_diag, best_loss)
    : scm_coord_descent_core(X0, X1, Z0_eval, Z1_eval,
                             max_iter, tol, V_diag, best_loss);

  // OOS mode: refit W on full pre-treatment data with selected V*
  arma::vec final_W = do_oos ? scm_inner_weights_cpp(X0, X1, V_diag) : best_W;
  // Reported loss always covers the full Z window so pre-fit quality stays
  // comparable across evaluation-window choices.
  double final_loss = arma::norm(Z1 - Z0 * final_W, 2);

  return Rcpp::List::create(
    Rcpp::Named("W")    = final_W,
    Rcpp::Named("V")    = V_diag,
    Rcpp::Named("loss") = final_loss
  );
}
