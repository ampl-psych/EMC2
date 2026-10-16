#ifndef EMC2_SAMPLER_GROUP_H
#define EMC2_SAMPLER_GROUP_H

#include <RcppArmadillo.h>
#include "random.h"
#include <vector>
#include <memory>

namespace emc2 {

class SEMSampler;
class InfiniteFactorSampler;

arma::mat group_inverse(const arma::mat& x);
arma::vec group_normal(Random& rng, const arma::vec& mean, const arma::mat& covariance);
arma::mat group_inverse_wishart(Random& rng, double df, const arma::mat& scale);

struct GroupState {
  arma::vec tmu, a_half, u, s, a_half_s;
  arma::mat tvar, tvinv, subj_mu;
  arma::mat lambda, eta, psi, B, K, G, delta_inv, population_var;
  arma::vec epsilon_inv, delta;
};

// Construct once per chain/stage; no R objects or callbacks enter draw().
class GroupSampler {
 public:
  explicit GroupSampler(const Rcpp::List& sampler);
  ~GroupSampler();
  GroupState initialize(Random& rng, SEXP start_mu = R_NilValue,
                        SEXP start_var = R_NilValue) const;
  GroupState read_state(const Rcpp::List& last) const;
  void draw(Random& rng, const arma::mat& alpha, GroupState& state) const;
  arma::mat subject_means(const GroupState& state) const;
  double affine(GroupState& state, arma::uword j, double log_scale, double shift) const;
  double shear(GroupState& state, arma::uword j, arma::uword k, double delta) const;
  double shear_precision(const GroupState& state, arma::uword j, arma::uword k) const;
  const std::vector<arma::uvec>& covariance_blocks() const { return blocks_; }
  bool location_move() const { return location_; }

  int parameters() const { return p_; }
  bool single() const { return single_; }

 private:
  struct RandomDesign {
    arma::uword parameter, start;
    arma::mat z, ztz;
  };
  bool single_, diagonal_gamma_ = false, location_ = false;
  std::unique_ptr<SEMSampler> sem_;
  std::unique_ptr<InfiniteFactorSampler> infnt_;
  arma::vec shape_, rate_;
  arma::uword p_, n_, random_size_ = 0;
  double v_ = 2, v_z_ = 2;
  arma::vec mean_, a_inv_sq_, az_inv_sq_, az_init_inv_sq_, block_dimensions_;
  arma::mat prior_var_, prior_precision_;
  arma::vec prior_linear_;
  std::vector<arma::mat> fixed_, crossproducts_;
  std::vector<arma::uword> starts_;
  std::vector<RandomDesign> random_;
  std::vector<arma::uvec> blocks_;
  arma::uvec unblocked_;
  std::vector<std::string> parameter_names_, mean_names_;

  arma::mat fixed_means(const arma::vec& beta) const;
  arma::mat random_means(const arma::vec& u) const;
};

}  // namespace emc2
#endif
