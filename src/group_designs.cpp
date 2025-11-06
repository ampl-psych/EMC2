// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

// Numerically safe dot for partial rows: sum_{k<j} L(i,k)*L(j,k)
inline double dot_partial(const arma::mat& L, arma::uword i, arma::uword j) {
  const arma::uword len = std::min(i, j);
  // Simple summation; for small p this is fine. If you expect large p and ill conditioning,
  // you can switch to Kahan summation.
  double s = 0.0;
  //const double* Li = L.memptr() + i;               // row-major access is not ideal in Armadillo,
  //const double* Lj = L.memptr() + j;               // but p is typically small here.
  //const arma::uword n_cols = L.n_cols;
  for (arma::uword k = 0; k < len; ++k) {
    // L(i,k) and L(j,k); column-major layout => access via (k*n_rows + i), so use L() API:
    s += L(i, k) * L(j, k);
  }
  return s;
}

// Lower-triangular Cholesky without LAPACK.
// On success: returns true and sets L (strictly lower + positive diag). On failure: returns false.
inline bool chol_lower_nolapack(const arma::mat& Ain,
                                arma::mat& L,
                                bool symmetrize = true,
                                double base_jitter = 0.0,   // if > 0, start by adding this jitter
                                int    max_tries  = 3)      // total attempts: base + adaptive
{
  if (Ain.n_rows != Ain.n_cols) return false;
  const arma::uword p = Ain.n_rows;

  // Start from a working copy; optionally symmetrise
  arma::mat A = Ain;
  if (symmetrize) {
    A = 0.5 * (A + A.t());
  }

  // Scale for jitter based on diagonal magnitude
  auto diag = A.diag();
  const double scale = arma::mean(arma::abs(diag)) + 1e-16;

  // Adaptive jitter schedule: 0, 1e-12, 1e-10, 1e-8 (scaled), unless base_jitter provided
  std::vector<double> jitters;
  if (base_jitter > 0.0) {
    jitters.push_back(base_jitter);
  } else {
    jitters = {0.0, 1e-12, 1e-10, 1e-8};
  }
  // keep only max_tries entries
  if ((int)jitters.size() > max_tries) jitters.resize(max_tries);

  L.zeros(p, p);
  for (size_t att = 0; att < jitters.size(); ++att) {
    const double jitter = jitters[att] * scale;

    // If jitter > 0, add to diagonal
    if (jitter > 0.0) {
      A.diag() = diag + jitter;
    } else if (att > 0) {
      // reset in case previous attempt modified diag
      A.diag() = diag;
    }

    bool ok = true;
    L.zeros();

    for (arma::uword i = 0; i < p && ok; ++i) {
      // Diagonal
      double s = 0.0;
      for (arma::uword k = 0; k < i; ++k) s += L(i, k) * L(i, k);
      double d = A(i, i) - s;
      if (!(d > 0.0) || !std::isfinite(d)) { ok = false; break; }
      double lii = std::sqrt(d);
      if (!std::isfinite(lii) || lii <= 0.0) { ok = false; break; }
      L(i, i) = lii;

      // Off-diagonals below the diagonal: j = i+1..p-1 (but we fill as "row i" in lower form),
      // equivalently iterate rows r=i+1..p-1 and set L(r,i).
      for (arma::uword r = i + 1; r < p; ++r) {
        double s2 = 0.0;
        for (arma::uword k = 0; k < i; ++k) s2 += L(r, k) * L(i, k);
        double num = A(r, i) - s2;
        double val = num / lii;
        if (!std::isfinite(val)) { ok = false; break; }
        L(r, i) = val;
      }
    }

    if (ok) return true;
  }

  // all attempts failed
  L.reset();
  return false;
}


// Multiply xb = X * beta without BLAS (manual loops).
// X: (n_subjects x p_k), beta: (p_k), xb: (n_subjects)
// Assumes column-major storage; loops over columns outermost to keep cache-friendly axpy style.
inline void matvec_noblas(const arma::mat& X, const arma::colvec& beta, arma::colvec& xb) {
  const arma::uword n = X.n_rows;
  const arma::uword m = X.n_cols;
  xb.zeros();
  for (arma::uword j = 0; j < m; ++j) {
    const double bj = beta[j];
    if (bj == 0.0) continue;
    const double* xcol = X.colptr(j);
    for (arma::uword s = 0; s < n; ++s) {
      xb[s] += bj * xcol[s];
    }
  }
}

// Build subject-specific means (p x n_subjects) without BLAS.
// group_designs: list of length p, each (n_subjects x m_k)
// params: stacked coefficients (length M = sum_k m_k)
// [[Rcpp::export]]
arma::mat calculate_subject_means(const Rcpp::List& group_designs,
                                                const arma::colvec& params) {
  const int p = group_designs.size();
  if (p == 0) Rcpp::stop("group_designs is empty.");

  arma::mat X0 = group_designs[0];
  const arma::uword n_subjects = X0.n_rows;

  arma::mat subj_mu(p, n_subjects, arma::fill::zeros);

  int par_idx = 0;
  arma::colvec xb(n_subjects, arma::fill::none);

  // k loops parameters; inside, do xb = Xk * beta_k with manual loops
  for (int k = 0; k < p; ++k) {
    arma::mat Xk = group_designs[k];                  // (n_subjects x m_k)
    const int pk = static_cast<int>(Xk.n_cols);
    arma::colvec beta_k = params.subvec(par_idx, par_idx + pk - 1);
    matvec_noblas(Xk, beta_k, xb);                    // xb = Xk * beta_k
    subj_mu.row(k) = xb.t();                          // store as row
    par_idx += pk;
  }
  return subj_mu; // (p x n_subjects)
}

// [[Rcpp::export]]
arma::cube draw_alpha_from_design(const Rcpp::List& group_designs,
                                         const arma::mat&  mu,   // (M x N), stacked [intercept + effects]
                                         const arma::cube& var)  // (p x p x N)
{
  const int p = group_designs.size();
  if (p == 0) Rcpp::stop("group_designs is empty.");

  arma::mat X0 = group_designs[0];
  const arma::uword n_subjects = X0.n_rows;

  // compute total M and sanity-shape
  int M = X0.n_cols;
  for (int k = 1; k < p; ++k) {
    arma::mat Xk = group_designs[k];
    if (Xk.n_rows != n_subjects) Rcpp::stop("All design blocks must share the same number of rows (subjects).");
    M += Xk.n_cols;
  }

  const int N = mu.n_cols;
  if ((int)mu.n_rows != M) Rcpp::stop("mu has %d rows, expected %d.", (int)mu.n_rows, M);
  if ((int)var.n_rows != p || (int)var.n_cols != p || (int)var.n_slices != N)
    Rcpp::stop("var must be (p x p x N).");

  arma::cube alpha(p, n_subjects, N, arma::fill::none);


  for (int i = 0; i < N; ++i) {
    arma::mat subj_mu = calculate_subject_means(group_designs, mu.col(i)); // (p x n_subjects)

    arma::mat Vi = var.slice(i);
    arma::mat L;

    // Try no-LAPACK Cholesky with small adaptive jitter
    if (!chol_lower_nolapack(Vi, L, /*symmetrize=*/true,
                             /*base_jitter=*/0.0,
                             /*max_tries=*/4)) {
                             Rcpp::stop("Cholesky (no-LAPACK) failed on slice %d after jitter attempts.", i + 1);
    }

    arma::mat Z = arma::randn(p, n_subjects);
    arma::mat draw = L * Z;      // OK without BLAS; Armadillo falls back to internal loops
    draw += subj_mu;

    alpha.slice(i) = std::move(draw);
  }

  return alpha; // (p x n_subjects x N)
}
