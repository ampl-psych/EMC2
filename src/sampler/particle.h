#ifndef EMC2_SAMPLER_PARTICLE_H
#define EMC2_SAMPLER_PARTICLE_H

#include <RcppArmadillo.h>
#include <vector>
#include "random.h"

namespace emc2 {

class Gaussian {
 public:
  Gaussian() = default;
  explicit Gaussian(const arma::mat& covariance);
  arma::mat draw(Random& rng, const arma::vec& mean, int n, double scale = 1) const;
  arma::vec log_density(const arma::mat& x, const arma::vec& mean,
                        double scale = 1) const;
 private:
  arma::mat root_, inverse_root_;
  double log_constant_ = 0;
};

class ParticleSampler {
 public:
  ParticleSampler(Random& rng, const Rcpp::List& sampler, const Rcpp::List& tune,
                  const Rcpp::List& pm_settings, const std::string& kernel);
  int block_count() const { return indices_.size(); }
  int block_component(int block) const { return shared_[block]; }
  void begin_iteration(const arma::mat& covariance);
  const arma::mat& prepare(int subject, int block, const arma::vec& current,
                          const arma::vec& group_mean, const arma::mat& prior_precision);
  arma::uword select(int subject, int block, const arma::vec& likelihood);
  const arma::mat& candidate(int subject) const { return subjects_[subject].candidates; }
  Rcpp::List settings() const;

 private:
  struct Settings {
    arma::vec mix, epsilon, proposal_counts, acc_counts;
    int n_particles, iter;
    double local_uses = 0, log_eps_sum = 0, log_ess_sum = 0,
      log_ess_n = 0, ess_target = 0, max_particles = 0;
  };
  struct Subject {
    arma::vec chain_mean, efficient_mean, likelihood_linear;
    arma::mat likelihood_precision;
    std::vector<Gaussian> chain, efficient;
    std::vector<Gaussian> conditional, local;
    std::vector<bool> has_local;
    arma::mat conditional_covariance;
    unsigned int covariance_version = 0;
    std::vector<Settings> settings;
    arma::mat candidates;
    arma::vec log_weights;
    arma::ivec counts;
  };
  Random& rng;
  std::string kernel_;
  bool exact_, frozen_;
  int proposals_, multiplier_, n0_;
  double min_local_, target_ess_, ess_scale_, max_particles_, mix_adapt_, alpha_star_;
  arma::vec target_accept_;
  std::vector<bool> local_;
  std::vector<int> shared_, labels_;
  std::vector<arma::uvec> indices_;
  std::vector<Subject> subjects_;
  Rcpp::List original_settings_;
  arma::mat covariance_;
  Gaussian group_;
  std::vector<Gaussian> group_blocks_;
  unsigned int covariance_version_ = 0;

  void update(Settings& settings, const arma::vec& weights,
              const arma::ivec& counts, int dimensions) const;
};

}  // namespace emc2
#endif
