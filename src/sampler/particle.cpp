#include "particle.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace emc2 {
namespace {

bool present(const Rcpp::List& x, const char* name) {
  return x.containsElementNamed(name) && !Rf_isNull(x[name]);
}
double number(const Rcpp::List& x, const char* name, double fallback = 0) {
  return present(x, name) ? Rcpp::as<double>(x[name]) : fallback;
}
bool flag(const Rcpp::List& x, const char* name) {
  return present(x, name) && Rcpp::as<bool>(x[name]);
}
SEXP subject_item(const Rcpp::List& x, const char* name, int subject) {
  if (!present(x, name)) return R_NilValue;
  const Rcpp::List values = x[name];
  return values[subject];
}
Rcpp::NumericVector vector(const arma::vec& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}
bool conditional_covariance(const arma::mat& likelihood, const arma::mat& prior,
                            arma::mat& covariance) {
  if (likelihood.is_empty() || prior.is_empty()) return false;
  return arma::inv_sympd(covariance, likelihood + prior) && covariance.is_finite();
}

}  // namespace

Gaussian::Gaussian(const arma::mat& covariance) {
  if (!arma::chol(root_, covariance, "lower"))
    throw std::runtime_error("singular Gaussian covariance");
  inverse_root_ = arma::inv(arma::trimatu(root_.t()));
  log_constant_ = -.5 * covariance.n_rows * std::log(2 * arma::datum::pi) -
    arma::sum(arma::log(root_.diag()));
}

arma::mat Gaussian::draw(Random& rng, const arma::vec& mean, int n, double scale) const {
  arma::mat z(n, root_.n_cols);
  for (int i = 0; i < n; ++i)
    for (arma::uword j = 0; j < z.n_cols; ++j) z(i, j) = rng.rnorm(0, 1);
  arma::mat out = scale * z * root_.t();
  out.each_row() += mean.t();
  return out;
}

arma::vec Gaussian::log_density(const arma::mat& x, const arma::vec& mean,
                               double scale) const {
  arma::mat centered = x;
  centered.each_row() -= mean.t();
  const arma::mat z = centered * inverse_root_;
  arma::vec out = log_constant_ - root_.n_cols * std::log(scale) -
    .5 * arma::sum(arma::square(z), 1) / (scale * scale);
  return out;
}

ParticleSampler::ParticleSampler(Random& rng_, const Rcpp::List& sampler, const Rcpp::List& tune,
                                const Rcpp::List& pm_settings, const std::string& kernel)
  : rng(rng_), kernel_(kernel), exact_(flag(tune, "exact")), frozen_(flag(tune, "frozen")),
    original_settings_(Rcpp::clone(pm_settings)) {
  proposals_ = kernel == "preburn" ? 2 : kernel == "sample" ? 4 : 3;
  multiplier_ = kernel == "preburn" ? 2 : 1;
  local_.assign(proposals_, false);
  local_[1] = true;
  if (kernel == "burn") local_[2] = true;
  n0_ = number(tune, "n0", 25);
  min_local_ = number(tune, "min_local");
  target_ess_ = number(tune, "target_ESS");
  ess_scale_ = number(tune, "ESS_scale", .05);
  max_particles_ = number(tune, "max_particles");
  mix_adapt_ = number(tune, "mix_adapt", .05);
  alpha_star_ = number(tune, "alphaStar", kernel == "sample" ? 2 : 3);
  target_accept_ = Rcpp::as<arma::vec>(tune["p_accept"]);
  const arma::ivec components = Rcpp::as<arma::ivec>(tune["components"]);
  const arma::ivec shared = Rcpp::as<arma::ivec>(tune["shared_ll_idx"]);
  for (const auto label : components) {
    if (std::find(labels_.begin(), labels_.end(), label) != labels_.end()) continue;
    labels_.push_back(label);
    indices_.push_back(arma::find(components == label));
    const int likelihood = shared[indices_.back()[0]];
    shared_.push_back(likelihood);
  }
  const int n = Rcpp::as<int>(sampler["n_subjects"]);
  subjects_.resize(n);
  for (int i = 0; i < n; ++i) {
    Subject& subject = subjects_[i];
    const SEXP chain_mean = subject_item(sampler, "chains_mu", i);
    const SEXP efficient_mean = subject_item(sampler, "eff_mu", i);
    if (!Rf_isNull(chain_mean)) subject.chain_mean = Rcpp::as<arma::vec>(chain_mean);
    if (!Rf_isNull(efficient_mean)) subject.efficient_mean = Rcpp::as<arma::vec>(efficient_mean);
    const SEXP chain_covariance = subject_item(sampler, "chains_var", i);
    const SEXP efficient_covariance = subject_item(sampler, "eff_var", i);
    for (const auto& idx : indices_) {
      subject.chain.emplace_back(Rf_isNull(chain_covariance) ? Gaussian() :
        Gaussian(Rcpp::as<arma::mat>(chain_covariance).submat(idx, idx)));
      subject.efficient.emplace_back(Rf_isNull(efficient_covariance) ? Gaussian() :
        Gaussian(Rcpp::as<arma::mat>(efficient_covariance).submat(idx, idx)));
    }
    SEXP likelihood = subject_item(tune, "lik_prec", i);
    if (!Rf_isNull(likelihood)) {
      Rcpp::List quadratic(likelihood);
      if (present(quadratic, "post")) quadratic = Rcpp::as<Rcpp::List>(quadratic["post"]);
      if (present(quadratic, "prec")) subject.likelihood_precision = Rcpp::as<arma::mat>(quadratic["prec"]);
      if (present(quadratic, "lin")) subject.likelihood_linear = Rcpp::as<arma::vec>(quadratic["lin"]);
    }
    const Rcpp::List settings = pm_settings[i];
    for (const auto label : labels_) {
      const Rcpp::List setting = settings[label - 1];
      Settings value;
      value.mix = Rcpp::as<arma::vec>(setting["mix"]);
      value.epsilon = Rcpp::as<arma::vec>(setting["epsilon"]);
      value.proposal_counts = Rcpp::as<arma::vec>(setting["proposal_counts"]);
      value.acc_counts = Rcpp::as<arma::vec>(setting["acc_counts"]);
      // Settings carried over from another kernel (check_sampling_settings
      // resizes them) would otherwise fail deep in the adaptation arithmetic
      const arma::uword k = proposals_;
      if (value.mix.n_elem != k || value.epsilon.n_elem != k - 1 ||
          value.proposal_counts.n_elem != k || value.acc_counts.n_elem != k)
        Rcpp::stop("particle settings of subject %d (component %d) do not match the '%s' kernel: "
                   "mix has %d entries, epsilon %d, counts %d/%d; expected %d, %d, %d/%d",
                   i + 1, label, kernel_.c_str(), (int) value.mix.n_elem, (int) value.epsilon.n_elem,
                   (int) value.proposal_counts.n_elem, (int) value.acc_counts.n_elem,
                   (int) k, (int) k - 1, (int) k, (int) k);
      value.n_particles = number(setting, "n_particles");
      value.iter = number(setting, "iter", 1);
      value.local_uses = number(setting, "local_uses");
      value.log_eps_sum = number(setting, "log_eps_sum");
      value.log_ess_sum = number(setting, "log_ess_sum");
      value.log_ess_n = number(setting, "log_ess_n");
      value.ess_target = number(setting, "ess_target");
      value.max_particles = number(setting, "max_particles");
      subject.settings.push_back(std::move(value));
    }
  }
}

void ParticleSampler::begin_iteration(const arma::mat& covariance) {
  if (covariance_.n_rows != covariance.n_rows ||
      !arma::approx_equal(covariance_, covariance, "absdiff", 0)) {
    covariance_ = covariance;
    group_ = Gaussian(covariance);
    group_blocks_.clear();
    for (const auto& block : indices_)
      group_blocks_.emplace_back(covariance.submat(block, block));
    ++covariance_version_;
  }
}

const arma::mat& ParticleSampler::prepare(int s, int b, const arma::vec& current,
                                        const arma::vec& mean,
                                        const arma::mat& prior_precision) {
  Subject& subject = subjects_[s];
  const Settings& settings = subject.settings[b];
  const auto& idx = indices_[b];
  if (exact_ && subject.covariance_version != covariance_version_) {
    subject.covariance_version = covariance_version_;
    subject.conditional.clear();
    subject.local.clear();
    subject.has_local.assign(indices_.size(), false);
    subject.conditional_covariance.reset();
    if (conditional_covariance(subject.likelihood_precision, prior_precision,
                               subject.conditional_covariance)) {
      for (const auto& block : indices_)
        subject.conditional.emplace_back(subject.conditional_covariance.submat(block, block));
    }
    for (arma::uword k = 0; k < indices_.size(); ++k) {
      const auto& block = indices_[k];
      arma::mat local_covariance;
      const bool valid = !subject.likelihood_precision.is_empty() && !prior_precision.is_empty() &&
        conditional_covariance(subject.likelihood_precision.submat(block, block),
          prior_precision.submat(block, block), local_covariance);
      subject.has_local[k] = valid;
      subject.local.emplace_back(valid ? Gaussian(local_covariance) : Gaussian());
    }
  }
  std::vector<arma::vec> centers(proposals_);
  std::vector<const Gaussian*> kernels(proposals_);
  centers[0] = mean.elem(idx);
  centers[1] = current.elem(idx);
  kernels[0] = &group_blocks_[b];
  kernels[1] = kernel_ == "preburn" || kernel_ == "burn" ?
    &group_blocks_[b] : &subject.chain[b];
  if (proposals_ >= 3) {
    centers[2] = kernel_ == "burn" ? current.elem(idx) : subject.chain_mean.elem(idx);
    kernels[2] = &subject.chain[b];
  }
  if (proposals_ == 4) {
    centers[3] = subject.efficient_mean.elem(idx);
    kernels[3] = &subject.efficient[b];
  }
  if (!subject.conditional.empty()) {
    const arma::vec conditional_mean = subject.conditional_covariance *
      (subject.likelihood_linear + prior_precision * mean);
    if (conditional_mean.is_finite()) {
      centers[2] = conditional_mean.elem(idx);
      kernels[2] = &subject.conditional[b];
    }
  }
  const double local_probability = settings.mix[1];
  const bool use_local = exact_ && rng.runif(0, 1) < local_probability;
  std::vector<int> active;
  for (int j = 0; j < proposals_; ++j)
    if (!exact_ || local_[j] == use_local) active.push_back(j);
  arma::vec mixture(active.size());
  for (arma::uword j = 0; j < active.size(); ++j) mixture[j] = settings.mix[active[j]];
  mixture /= arma::sum(mixture);
  std::vector<int> counts(active.size());
  rng.rmultinom(settings.n_particles * multiplier_, mixture.memptr(), active.size(), counts.data());
  subject.counts.zeros(proposals_);
  for (arma::uword k = 0; k < active.size(); ++k)
    subject.counts[active[k]] = exact_ ? counts[k] : std::max(1, counts[k]);
  if (use_local && subject.has_local[b]) kernels[1] = &subject.local[b];
  arma::vec scales(proposals_, arma::fill::ones);
  for (int j = 1; j < proposals_; ++j) scales[j] = settings.epsilon[j - 1];
  for (const auto j : active)
    if (exact_ && local_[j]) centers[j] = kernels[j]->draw(rng, centers[j], 1, scales[j]).row(0).t();
  const arma::uword rows = arma::sum(subject.counts) + 1;
  subject.candidates = arma::repmat(current.t(), rows, 1);
  arma::mat proposed(rows, idx.n_elem);
  proposed.row(0) = current.elem(idx).t();
  arma::uword offset = 1;
  for (const auto j : active) {
    const int count = subject.counts[j];
    if (count == 0) continue;
    proposed.rows(offset, offset + count - 1) = kernels[j]->draw(rng, centers[j], count, scales[j]);
    offset += count;
  }
  subject.candidates.cols(idx) = proposed;
  subject.log_weights = indices_.size() == 1 ?
    group_blocks_[b].log_density(proposed, mean.elem(idx)) :
    group_.log_density(subject.candidates, mean);
  if (!exact_ || !use_local) {
    arma::mat log_densities(rows, active.size());
    for (arma::uword k = 0; k < active.size(); ++k) {
      const int j = active[k];
      log_densities.col(k) = std::log(mixture[k]) + kernels[j]->log_density(proposed, centers[j], scales[j]);
    }
    arma::vec log_mixture(rows);
    for (arma::uword row = 0; row < rows; ++row) {
      const double largest = log_densities.row(row).max();
      log_mixture[row] = largest + std::log(arma::sum(arma::exp(log_densities.row(row) - largest)));
    }
    const arma::uvec finite = arma::find_finite(log_mixture);
    if (finite.n_elem)
      log_mixture.elem(arma::find_nonfinite(log_mixture)).fill(log_mixture.elem(finite).min());
    subject.log_weights -= log_mixture;
  }
  return subject.candidates;
}

arma::uword ParticleSampler::select(int s, int b, const arma::vec& likelihood) {
  Subject& subject = subjects_[s];
  arma::vec log_weights = subject.log_weights + likelihood;
  arma::vec weights = arma::exp(log_weights - log_weights.max());
  const double total = arma::sum(weights);
  if (!std::isfinite(total) || total <= 0)
    throw std::runtime_error("particle weights are not finite");
  double threshold = rng.runif(0, total);
  arma::uword selected = 0;
  while (selected + 1 < weights.n_elem && threshold > weights[selected])
    threshold -= weights[selected++];
  if (!frozen_) update(subject.settings[b], weights, subject.counts, indices_[b].n_elem);
  return selected;
}

void ParticleSampler::update(Settings& x, const arma::vec& weights,
                             const arma::ivec& counts, int d) const {
  if (++x.iter <= n0_) return;
  arma::vec rate_now(proposals_, arma::fill::value(arma::datum::nan));
  arma::uword offset = 1;
  for (int j = 0; j < proposals_; ++j) {
    x.proposal_counts[j] += counts[j];
    if (counts[j] == 0) continue;
    const double better = arma::accu(weights.subvec(offset, offset + counts[j] - 1) > weights[0]);
    x.acc_counts[j] += better;
    rate_now[j] = better / counts[j];
    offset += counts[j];
  }
  arma::vec rates(proposals_);
  for (int j = 0; j < proposals_; ++j)
    rates[j] = x.proposal_counts[j] > 0 ? x.acc_counts[j] / x.proposal_counts[j] : 0;
  for (int j = 1; j < proposals_; ++j) {
    const double target = target_accept_[j - 1];
    if (exact_) {
      if (!local_[j] || !std::isfinite(rate_now[j])) continue;
      const double signal = std::clamp((rate_now[j] - target) / target, -1.0, 3.0);
      x.epsilon[j - 1] = std::clamp(std::exp(std::log(x.epsilon[j - 1]) + .1 * signal), .01, 20.0);
      if (++x.local_uses > 20) x.log_eps_sum += std::log(x.epsilon[j - 1]);
    } else {
      const double c = (1 - 1.0 / d) * std::sqrt(2 * arma::datum::pi) *
        std::exp(alpha_star_ * alpha_star_ / 2) / (2 * alpha_star_) +
        1 / (d * target * (1 - target));
      const double signal = std::clamp((rates[j] - target) / target, -1.0, 1.0);
      x.epsilon[j - 1] = std::clamp(std::exp(std::log(x.epsilon[j - 1]) +
        c / std::max(100, x.iter) * signal), .01, 20.0);
    }
  }
  if (proposals_ > 2) {
    arma::vec performance = (rates + 1e-12) / x.mix;
    performance[0] /= arma::mean(target_accept_);
    for (int j = 1; j < proposals_; ++j) performance[j] /= target_accept_[j - 1];
    performance /= arma::sum(performance);
    x.mix = arma::clamp((1 - mix_adapt_) * x.mix + mix_adapt_ * performance, .02, arma::datum::inf);
    x.mix /= arma::sum(x.mix);
    if (min_local_ > 0 && x.mix[1] < min_local_) {
      const double previous_local = x.mix[1];
      x.mix *= (1 - min_local_) / (1 - previous_local);
      x.mix[1] = min_local_;
    }
  }
  if (proposals_ > 3) {
    const double ess = std::pow(arma::sum(weights), 2) / arma::dot(weights, weights);
    x.n_particles = std::max(25.0, std::min(max_particles_,
      std::nearbyint(x.n_particles * std::pow(target_ess_ / ess, ess_scale_))));
    x.log_ess_sum += std::log(ess / (weights.n_elem - 1));
    ++x.log_ess_n;
    x.ess_target = target_ess_;
    x.max_particles = max_particles_;
  }
}

Rcpp::List ParticleSampler::settings() const {
  Rcpp::List out = Rcpp::clone(original_settings_);
  for (arma::uword s = 0; s < subjects_.size(); ++s) {
    Rcpp::List blocks = out[s];
    for (arma::uword b = 0; b < indices_.size(); ++b) {
      const Settings& x = subjects_[s].settings[b];
      Rcpp::List value = blocks[labels_[b] - 1];
      value["mix"] = vector(x.mix);
      value["epsilon"] = vector(x.epsilon);
      value["proposal_counts"] = vector(x.proposal_counts);
      value["acc_counts"] = vector(x.acc_counts);
      value["n_particles"] = x.n_particles;
      value["iter"] = x.iter;
      if (x.local_uses > 0 || present(value, "local_uses")) {
        value["local_uses"] = x.local_uses;
        value["log_eps_sum"] = x.log_eps_sum;
      }
      if (x.log_ess_n > 0 || present(value, "log_ess_n")) {
        value["log_ess_sum"] = x.log_ess_sum;
        value["log_ess_n"] = x.log_ess_n;
        value["ess_target"] = x.ess_target;
        value["max_particles"] = x.max_particles;
      }
      blocks[labels_[b] - 1] = value;
    }
    out[s] = blocks;
  }
  return out;
}

}  // namespace emc2
