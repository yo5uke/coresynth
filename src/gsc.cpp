#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Additive fixed effects of the extended model in Xu (2017, Remark 3):
// whether the counterfactual carries a unit effect alpha_i, a time effect
// xi_t, both or neither, on top of the factors and the covariates.
struct GscFe {
  bool unit;
  bool time;
};

static GscFe gsc_parse_force(const std::string& force) {
  if (force == "none")    return {false, false};
  if (force == "unit")    return {true,  false};
  if (force == "time")    return {false, true};
  if (force == "two-way") return {true,  true};
  Rcpp::stop("force must be one of \"none\", \"unit\", \"time\", \"two-way\"");
  return {false, false};
}

// Within transformation over the control panel: remove each unit's mean over
// time (unit effect), then each period's mean across units (time effect).
// Applying both is the two-way double demeaning.
static arma::mat gsc_demean(const arma::mat& M, const GscFe& fe) {
  arma::mat out = M;
  if (fe.unit) out.each_row() -= arma::mean(out, 0);
  if (fe.time) out.each_col() -= arma::mean(out, 1);
  return out;
}

// Rank-r factors and loadings of M from its SVD; r = 0 gives empty ones.
static void gsc_factors(const arma::mat& M, int r, arma::mat& U, arma::vec& D,
                        arma::mat& V, arma::mat& F, arma::mat& L) {
  arma::svd_econ(U, D, V, M);
  if (r > 0) {
    F = U.cols(0, r - 1) * arma::diagmat(D.subvec(0, r - 1));
    L = V.cols(0, r - 1);
  } else {
    F = arma::zeros(M.n_rows, 0);
    L = arma::zeros(M.n_cols, 0);
  }
}

// Treated-unit step of Xu (2017), Step 2: each treated unit's loadings -- and
// its unit effect when the model has one -- by OLS of its pre-treatment
// outcomes, net of the covariate part, the grand mean and the time effects,
// on the factors' pre-treatment rows. Y_hat is the counterfactual over all T
// periods without the covariate contribution.
struct GscTreated {
  arma::mat  L_t;     // r x N_tr
  arma::rowvec alpha; // N_tr
  arma::mat  Y_hat;   // T x N_tr
};

static GscTreated gsc_treated_fit(const arma::mat& F, const arma::vec& beta,
                                  const arma::mat& Y_tr_pre,
                                  const arma::cube& X_tr_pre,
                                  double mu, const arma::vec& xi,
                                  bool unit_fe) {
  const arma::uword T_pre = Y_tr_pre.n_rows;
  const arma::uword n     = Y_tr_pre.n_cols;
  const bool has_xi = xi.n_elem == F.n_rows;

  arma::mat Y_dm = Y_tr_pre - mu;
  for (arma::uword j = 0; j < X_tr_pre.n_slices; j++)
    Y_dm -= X_tr_pre.slice(j) * beta(j);
  if (has_xi) Y_dm.each_col() -= xi.head(T_pre);

  arma::mat Z = F.rows(0, T_pre - 1);
  if (unit_fe) Z = arma::join_rows(arma::ones(T_pre), Z);
  arma::mat coef = Z.n_cols > 0 ? arma::mat(arma::pinv(Z) * Y_dm)
                                : arma::mat(0, n);

  GscTreated out;
  out.alpha = unit_fe ? arma::rowvec(coef.row(0)) : arma::zeros<arma::rowvec>(n);
  out.L_t   = unit_fe ? arma::mat(coef.tail_rows(F.n_cols)) : coef;
  out.Y_hat = F * out.L_t + mu;
  out.Y_hat.each_row() += out.alpha;
  if (has_xi) out.Y_hat.each_col() += xi;
  return out;
}

//' Treated-unit step of the IFE estimator for given factors (internal)
//'
//' The counterfactual of treated units whose first `nrow(Y_tr_pre)` periods
//' are untreated, from the factors, coefficients and additive effects that a
//' control-unit fit of [gsc_ife_cpp()] produced. Lets a staggered fit reuse
//' one control-unit fit for cohorts with different adoption dates.
//'
//' @return A list with `L_tr` (N_tr x r), `alpha_tr` (N_tr; zero without a
//'   unit effect) and `Y_tr_hat` (T x N_tr, without the covariate
//'   contribution).
//' @noRd
// [[Rcpp::export]]
Rcpp::List gsc_treated_cpp(const arma::mat& F, const arma::vec& beta,
                           const arma::mat& Y_tr_pre,
                           const arma::cube& X_tr_pre,
                           double mu, const arma::vec& xi, bool unit_fe) {
  if (Y_tr_pre.n_rows < F.n_cols + (unit_fe ? 1u : 0u))
    Rcpp::stop("T_pre is too short to identify the treated loadings");
  if (X_tr_pre.n_slices > beta.n_elem)
    Rcpp::stop("X_tr_pre has more slices than there are coefficients");
  GscTreated tr = gsc_treated_fit(F, beta, Y_tr_pre, X_tr_pre, mu, xi, unit_fe);
  return Rcpp::List::create(
    Rcpp::Named("L_tr")     = tr.L_t.t(),
    Rcpp::Named("alpha_tr") = tr.alpha.t(),
    Rcpp::Named("Y_tr_hat") = tr.Y_hat
  );
}

//' Fast Interactive Fixed Effects (IFE) for Generalized Synthetic Control
//'
//' Implements Xu (2017) IFE model with optional covariate adjustment and
//' optional additive fixed effects (Xu 2017, Remark 3).
//' When X_co has p > 0 slices, alternates the two steps of Bai (2009):
//'   factors: truncated SVD of Y_tilde = Y_co - X_co * beta
//'   beta:    inverse(sum_i X_i' M_F X_i) sum_i X_i' M_F Y_i, with the loadings
//'            concentrated out through the annihilator M_F of the factors
//' When X_co has 0 slices (default), falls back to the plain 3-step estimator.
//' Additive effects are swept out of the control panel by the within
//' transformation first; on a balanced panel this solves the joint least
//' squares problem exactly, and the effects are recovered from the means of
//' the residual.
//'
//' @param Y_co     Control units outcome matrix (T x N_co)
//' @param Y_tr_pre Treated units pre-treatment outcomes (T_pre x N_tr)
//' @param r        Number of latent factors (must be <= min(T, N_co)); 0 is
//'                 allowed with additive effects
//' @param X_co     Time-varying covariate cube (T x N_co x p). Pass an empty
//'                 cube (0 slices) for the covariate-free estimator.
//' @param X_tr_pre Time-varying covariate cube for treated units in the
//'                 pre-treatment window (T_pre x N_tr x p). Required for
//'                 correct Step 2 loading estimation per Xu (2017): lambda_hat
//'                 is estimated from Y_tr_pre - X_tr_pre * beta (covariate-
//'                 demeaned). Pass an empty cube (0 slices) to skip demeaning
//'                 (backward-compatible, but biased when beta != 0).
//' @param max_iter Maximum iterations of the covariate loop (default 1000)
//' @param tol      Convergence tolerance on relative beta change (default 1e-8)
//' @param force    Additive fixed effects: `"none"` (default), `"unit"`,
//'                 `"time"` or `"two-way"`.
//' @return A list with components:
//'   * `F`: estimated time factors (T x r).
//'   * `L_co`: control-unit factor loadings (N_co x r).
//'   * `L_tr`: treated-unit factor loadings (N_tr x r).
//'   * `Y_tr_hat`: estimated treated-unit counterfactual outcomes (T x N_tr),
//'     additive effects included, covariate contribution excluded.
//'   * `Y_co_hat`: fitted control outcomes on the same basis (T x N_co).
//'   * `singular_values`: singular values from the final truncated SVD.
//'   * `beta`: estimated covariate coefficients (p x 1), empty when no
//'     covariates are supplied.
//'   * `mu`, `alpha_co`, `alpha_tr`, `xi`: grand mean, unit effects of the
//'     control and treated units, and time effects (zero where `force`
//'     leaves them out).
//'   * `n_iter`: iterations of the covariate loop (0 without covariates).
//'   * `converged`: whether the covariate loop met `tol` within `max_iter`
//'     (always `TRUE` without covariates).
//' @export
// [[Rcpp::export]]
Rcpp::List gsc_ife_cpp(const arma::mat& Y_co,
                        const arma::mat& Y_tr_pre,
                        int r,
                        const arma::cube& X_co,
                        const arma::cube& X_tr_pre,
                        int max_iter = 1000,
                        double tol   = 1e-8,
                        std::string force = "none") {
  int T     = (int)Y_co.n_rows;
  int N_co  = (int)Y_co.n_cols;
  int T_pre = (int)Y_tr_pre.n_rows;
  int p     = (int)X_co.n_slices;
  const GscFe fe = gsc_parse_force(force);

  if (r < 0 || r > std::min(T, N_co))
    Rcpp::stop("r must be between 0 and min(T, N_co)");
  if (T_pre < r + (fe.unit ? 1 : 0))
    Rcpp::stop("T_pre is too short to identify the treated loadings");

  // The model is fitted on the within-transformed panel
  const arma::mat Y_d = gsc_demean(Y_co, fe);
  arma::cube X_d(X_co.n_rows, X_co.n_cols, p);
  for (int j = 0; j < p; j++) X_d.slice(j) = gsc_demean(X_co.slice(j), fe);

  arma::vec beta = arma::zeros<arma::vec>(p);
  arma::mat U, V, F, L_co;
  arma::vec D;
  int  n_iter    = 0;
  bool converged = true;

  // Factors start from the outcome panel itself (beta = 0)
  gsc_factors(Y_d, r, U, D, V, F, L_co);

  if (p > 0) {
    // Given the factors, beta solves the regression with the loadings
    // concentrated out: the annihilator M_F = I - U_r U_r' removes the factor
    // space from every covariate slice. Regressing on Y - F L' instead (plain
    // alternating least squares) reaches the same fixed point but crawls
    // towards it when the covariates load on the factors, and a capped loop
    // then stops short of it with a biased beta.
    // A tiny ridge keeps the normal matrix invertible for collinear slices.
    double ridge = 1e-8;
    converged = false;

    for (int iter = 0; iter < max_iter; iter++) {
      arma::vec beta_old = beta;

      std::vector<arma::mat> MX(p);
      for (int j = 0; j < p; j++) {
        MX[j] = X_d.slice(j);
        if (r > 0) {
          arma::mat U_r = U.cols(0, r - 1);
          MX[j] -= U_r * (U_r.t() * X_d.slice(j));
        }
      }

      arma::mat XtMX(p, p);
      arma::vec XtMY(p);
      for (int j = 0; j < p; j++) {
        for (int k = j; k < p; k++)
          XtMX(j, k) = XtMX(k, j) = arma::accu(MX[j] % X_d.slice(k));
        XtMY(j) = arma::accu(MX[j] % Y_d);
      }
      XtMX.diag() += ridge;
      arma::vec beta_new = arma::solve(XtMX, XtMY);

      // Factors given beta: truncated SVD of Y_tilde = Y_d - X_d * beta
      arma::mat Y_tilde = Y_d;
      for (int j = 0; j < p; j++)
        Y_tilde -= X_d.slice(j) * beta_new(j);
      gsc_factors(Y_tilde, r, U, D, V, F, L_co);

      double delta = arma::norm(beta_new - beta_old, 2) /
                     (1.0 + arma::norm(beta_old, 2));
      beta   = beta_new;
      n_iter = iter + 1;
      if (delta < tol) {
        converged = true;
        break;
      }
    }
  }

  // ── Additive effects from the means of the residual ───────────────────────
  // The factor part fitted to the demeaned panel has zero means along the
  // demeaned directions, so the residual's means are the effects.
  arma::mat E = Y_co - F * L_co.t();
  for (int j = 0; j < p; j++) E -= X_co.slice(j) * beta(j);
  double   mu       = (fe.unit || fe.time) ? arma::mean(arma::vectorise(E)) : 0.0;
  arma::vec alpha_co = fe.unit ? arma::vec(arma::mean(E, 0).t() - mu)
                               : arma::zeros<arma::vec>(N_co);
  arma::vec xi       = fe.time ? arma::vec(arma::mean(E, 1) - mu)
                               : arma::zeros<arma::vec>(T);

  arma::mat Y_co_hat = F * L_co.t() + mu;
  Y_co_hat.each_row() += alpha_co.t();
  Y_co_hat.each_col() += xi;

  // ── Treated loadings: Xu (2017) Step 2 ─────────────────────────────────────
  // lambda_hat_i = (F0' F0)^{-1} F0' (Y0_i - X0_i beta)
  GscTreated tr = gsc_treated_fit(F, beta, Y_tr_pre, X_tr_pre, mu, xi, fe.unit);

  return Rcpp::List::create(
    Rcpp::Named("F")               = F,
    Rcpp::Named("L_co")            = L_co,
    Rcpp::Named("L_tr")            = tr.L_t.t(),
    Rcpp::Named("Y_tr_hat")        = tr.Y_hat,
    Rcpp::Named("Y_co_hat")        = Y_co_hat,
    Rcpp::Named("singular_values") = D,
    Rcpp::Named("beta")            = beta,
    Rcpp::Named("mu")              = mu,
    Rcpp::Named("alpha_co")        = alpha_co,
    Rcpp::Named("alpha_tr")        = tr.alpha.t(),
    Rcpp::Named("xi")              = xi,
    Rcpp::Named("force")           = force,
    Rcpp::Named("n_iter")          = n_iter,
    Rcpp::Named("converged")       = converged
  );
}
