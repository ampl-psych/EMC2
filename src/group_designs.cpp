// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

using namespace arma;

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
