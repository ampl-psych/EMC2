// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

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
    // 1) subject means for draw i (no BLAS)
    arma::mat subj_mu = calculate_subject_means(group_designs, mu.col(i)); // (p x n_subjects)

    // 2) Cholesky once for this draw
    arma::mat L;
    arma::mat Vi = var.slice(i);
    if (!arma::chol(L, Vi, "lower")) {
      double eps = 1e-10 * arma::mean(Vi.diag());
      if (!(eps > 0)) eps = 1e-10;
      if (!arma::chol(L, Vi + eps * arma::eye(p, p), "lower")) {
        eps = std::max(eps, 1e-8);
        if (!arma::chol(L, Vi + eps * arma::eye(p, p), "lower"))
          Rcpp::stop("Cholesky failed on slice %d.", i + 1);
      }
    }

    // 3) Sample all subjects: L * Z + subj_mu (this multiply is Armadillo's internal loops)
    arma::mat Z = arma::randn(p, n_subjects);
    arma::mat draw = L * Z;      // no BLAS required; falls back to internal loops if BLAS is unavailable
    draw += subj_mu;

    alpha.slice(i) = std::move(draw);
  }

  return alpha; // (p x n_subjects x N)
}
