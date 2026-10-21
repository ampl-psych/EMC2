#ifndef EMC2_SAMPLER_INFNT_FACTOR_H
#define EMC2_SAMPLER_INFNT_FACTOR_H

#include <RcppArmadillo.h>
#include "random.h"

namespace emc2 {

struct GroupState;

class InfiniteFactorSampler {
 public:
  explicit InfiniteFactorSampler(const Rcpp::List& sampler);
  GroupState initialize(Random& rng, SEXP start_mu = R_NilValue,
                        SEXP start_var = R_NilValue) const;
  GroupState read_state(const Rcpp::List& last) const;
  void draw(Random& rng, const arma::mat& alpha, GroupState& state) const;
  arma::mat subject_means(const GroupState& state) const;
  double affine(GroupState& state, arma::uword j,
                double log_scale, double shift) const;

 private:
  arma::uword p_, n_, factors_;
  arma::vec mean_, variance_, precision_, shape_, rate_;
  double df_, ad1_, bd1_, ad2_, bd2_;
};

}  // namespace emc2
#endif
