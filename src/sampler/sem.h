#ifndef EMC2_SAMPLER_SEM_H
#define EMC2_SAMPLER_SEM_H

#include "group.h"

namespace emc2 {

class SEMSampler {
 public:
  explicit SEMSampler(const Rcpp::List& sampler);
  GroupState initialize(Random& rng, SEXP start_mu = R_NilValue,
                        SEXP start_var = R_NilValue) const;
  GroupState read_state(const Rcpp::List& last) const;
  void draw(Random& rng, const arma::mat& alpha, GroupState& state) const;
  arma::mat subject_means(const GroupState& state) const;
  double affine(GroupState& state, arma::uword parameter,
                double log_scale, double shift) const;

 private:
  arma::uword p_, n_, q_, c_;
  arma::mat x_, lambda_, b_, k_, g_, loading_precision_, loading_linear_;
  arma::mat x_cov_;
  arma::vec mean_, mean_precision_, a_e_, b_e_, a_d_, b_d_;
  arma::vec structural_precision_;
  std::vector<arma::uvec> loading_free_, structural_free_, factor_groups_;
  bool cyclic_;
  void moments(GroupState& state) const;
};

}  // namespace emc2
#endif
