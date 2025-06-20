// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace arma;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

static double const log2pi = std::log(2.0 * M_PI);


/*===========================================================================
 1. calculate_subject_means  (unchanged logic, still Armadillo)
 ===========================================================================*/
// [[Rcpp::export]]
arma::mat calculate_subject_means(const Rcpp::List& group_designs,
                                  const arma::colvec& params,
                                  const int n_subjects,
                                  const int n_pars)
{
  arma::mat subj_mu(n_pars, n_subjects, arma::fill::zeros);

  for (int s = 0; s < n_subjects; ++s) {
    int par_idx = 0;
    for (int k = 0; k < n_pars; ++k) {

      arma::mat   Xk   = group_designs[k];
      arma::rowvec xsk = Xk.row(s);
      int          p_k = Xk.n_cols;

      arma::colvec beta = params.subvec(par_idx, par_idx + p_k - 1);
      subj_mu(k, s) = arma::as_scalar(xsk * beta);

      par_idx += p_k;
    }
  }
  return subj_mu;
}

/*===========================================================================
 2. standard_subj_ll  (Eigen LLT + zero-copy column maps)
 ===========================================================================*/
// [[Rcpp::export]]
arma::vec standard_subj_ll(const Rcpp::List& group_designs,
                           const arma::cube& theta_var,   // p × p × N
                           const arma::mat&  theta_mu,    // p × N
                           const arma::cube& alpha,       // p × n_subj × N
                           const int         n_subj)
{
  const int N = theta_var.n_slices;      // number of draws
  const int p = theta_mu.n_rows;         // number of parameters

  arma::vec lls(N, arma::fill::zeros);
  const MatrixXd I = MatrixXd::Identity(p, p);
  const double log2pi = std::log(2.0 * M_PI);

  // helper typedef: 1 × p const row mapped onto raw memory
  using MapRowC = Eigen::Map<
    const Eigen::Matrix<double,1,Eigen::Dynamic,Eigen::RowMajor>>;

  for (int i = 0; i < N; ++i) {

    /* 2.1  Eigen view of Σ_i and its Cholesky ---------------------- */
    const double* Si_ptr = theta_var.slice(i).memptr();
    Eigen::Map<const MatrixXd> Sigma(Si_ptr, p, p);

    Eigen::LLT<MatrixXd> llt(Sigma);
    if (llt.info() != Eigen::Success)
      Rcpp::stop("theta_var[,,%d] is not positive-definite.", i + 1);

    MatrixXd rooti = llt.matrixU().solve(I);                // U⁻¹
    const double log_const =
      rooti.diagonal().array().log().sum() - 0.5 * p * log2pi;

    /* 2.2  expected means for all subjects (p × n_subj) ------------ */
    arma::mat subj_mu = calculate_subject_means(
      group_designs, theta_mu.col(i), n_subj, p);

    /* 2.3  alpha slice for draw i  --------------------------------- */
    const double* al_ptr = alpha.slice(i).memptr();         // p × n_subj
    arma::mat alpha_i(const_cast<double*>(al_ptr), p, n_subj, false, true);

    /* 2.4  subject loop -------------------------------------------- */
    double ll_i = 0.0;
    for (int s = 0; s < n_subj; ++s) {

      MapRowC a (alpha_i.colptr(s),  p);                    // 1 × p
      MapRowC mu(subj_mu.colptr(s), p);

      RowVectorXd z = (a - mu) * rooti;
      ll_i += -0.5 * z.squaredNorm() + log_const;
    }
    lls(i) = ll_i;
  }
  return lls;
}


// -------------------------------------------------------------------
// [[Rcpp::export(name = "dmvnorm_cpp")]]
Eigen::VectorXd dmvnorm_cpp(const Eigen::MatrixXd&   X,
                            const Eigen::RowVectorXd mean,
                            const Eigen::MatrixXd&   sigma,
                            const bool               logd = false)
{
  const int n = X.rows();
  const int p = X.cols();

  // ---- 1. Cholesky Σ = Uᵀ U  --------------------------------------
  Eigen::LLT<Eigen::MatrixXd> llt(sigma);
  if (llt.info() != Eigen::Success)
    Rcpp::stop("sigma is not positive-definite.");

  // U⁻¹  (upper-triangular back-solve)
  Eigen::MatrixXd rooti =
    llt.matrixU().solve(Eigen::MatrixXd::Identity(p, p));

  // ---- 2. whiten all observations at once  ------------------------
  Eigen::MatrixXd Z = (X.rowwise() - mean) * rooti;     // n × p

  // ---- 3. quadratic forms (row-wise squared norms) ----------------
  Eigen::VectorXd quad = Z.rowwise().squaredNorm();     // length n

  // ---- 4. log-densities ------------------------------------------
  const double logdet = rooti.diagonal().array().log().sum();
  Eigen::VectorXd out =
    (-0.5 * quad.array()).matrix().array()
    + logdet
  - 0.5 * p * log2pi;

  return logd ? out : out.array().exp();
}

