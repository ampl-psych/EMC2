#include "sem.h"
#include <algorithm>

namespace emc2 {
namespace {

arma::vec recycled(SEXP value, arma::uword size) {
  const arma::vec x = Rcpp::as<arma::vec>(value);
  arma::vec out(size);
  for (arma::uword j = 0; j < size; ++j) out[j] = x[j % x.n_elem];
  return out;
}

arma::mat loading_prior(SEXP value, arma::uword p, arma::uword q) {
  if (Rf_isMatrix(value)) return Rcpp::as<arma::mat>(value);
  const arma::vec x = recycled(value, q);
  return arma::repmat(x.t(), p, 1);
}

arma::mat covariates(SEXP value, arma::uword n) {
  if (Rf_isNull(value)) return arma::mat(n, 0);
  if (!Rf_inherits(value, "data.frame")) return Rcpp::as<arma::mat>(value);
  const Rcpp::List columns(value);
  arma::mat out(n, columns.size());
  for (int j = 0; j < columns.size(); ++j) out.col(j) = Rcpp::as<arma::vec>(columns[j]);
  return out;
}

arma::mat fixed_values(const arma::mat& x) {
  arma::mat out = x;
  out.elem(arma::find_nonfinite(out)).zeros();
  return out;
}

bool cyclic_graph(const arma::mat& b) {
  arma::umat path = b != 0;
  for (arma::uword k = 0; k < b.n_rows; ++k)
    for (arma::uword i = 0; i < b.n_rows; ++i)
      for (arma::uword j = 0; j < b.n_rows; ++j)
        path(i, j) |= path(i, k) && path(k, j);
  return arma::any(path.diag());
}

double structural_logdet(const arma::mat& b) {
  double value, sign;
  arma::log_det(value, sign, arma::eye(b.n_rows, b.n_cols) - b);
  return value;
}

}  // namespace

SEMSampler::SEMSampler(const Rcpp::List& sampler) {
  p_ = Rcpp::as<int>(sampler["n_pars"]);
  n_ = Rcpp::as<int>(sampler["n_subjects"]);
  const Rcpp::List settings = sampler["sem_settings"], prior = sampler["prior"];
  lambda_ = Rcpp::as<arma::mat>(settings["Lambda_mat"]);
  q_ = lambda_.n_cols;
  x_ = covariates(settings["covariates"], n_);
  c_ = x_.n_cols;
  b_ = q_ ? Rcpp::as<arma::mat>(settings["B_mat"]) : arma::mat(0, 0);
  k_ = Rcpp::as<arma::mat>(settings["K_mat"]);
  g_ = q_ ? Rcpp::as<arma::mat>(settings["G_mat"]) : arma::mat(0, c_);
  mean_ = recycled(prior["theta_mu_mean"], p_);
  mean_precision_ = 1 / recycled(prior["theta_mu_var"], p_);
  a_e_ = recycled(prior["a_e"], p_);
  b_e_ = recycled(prior["b_e"], p_);
  a_d_ = recycled(prior["a_d"], q_);
  b_d_ = recycled(prior["b_d"], q_);
  const arma::vec k_var = recycled(prior["K_var"], c_);
  const arma::mat lambda_var = loading_prior(prior["lambda_var"], p_, q_);
  const arma::mat lambda_mean = loading_prior(prior["lambda_mean"], p_, q_);
  loading_precision_ = 1 / arma::join_rows(arma::repmat(k_var.t(), p_, 1), lambda_var);
  loading_linear_ = loading_precision_ % arma::join_rows(arma::mat(p_, c_, arma::fill::zeros), lambda_mean);
  structural_precision_ = 1 / arma::join_cols(recycled(prior["G_var"], c_), recycled(prior["B_var"], q_));
  const arma::mat loading = arma::join_rows(k_, lambda_);
  const arma::mat structural = arma::join_rows(g_, b_);
  for (arma::uword j = 0; j < p_; ++j)
    loading_free_.push_back(arma::find_nonfinite(loading.row(j).t()));
  for (arma::uword j = 0; j < q_; ++j)
    structural_free_.push_back(arma::find_nonfinite(structural.row(j).t()));
  cyclic_ = cyclic_graph(b_);
  lambda_ = fixed_values(lambda_); b_ = fixed_values(b_);
  k_ = fixed_values(k_); g_ = fixed_values(g_);
  const Rcpp::IntegerVector labels = Rcpp::as<Rcpp::IntegerVector>(settings["factor_groups"]);
  std::vector<int> groups;
  for (arma::uword j = 0; j < q_; ++j) {
    if (std::find(groups.begin(), groups.end(), labels[j]) != groups.end()) continue;
    groups.push_back(labels[j]);
    std::vector<arma::uword> indices;
    for (arma::uword k = j; k < q_; ++k) if (labels[k] == labels[j]) indices.push_back(k);
    factor_groups_.emplace_back(indices);
    if (indices.size() > 1 && a_d_[j] <= indices.size() - 1)
      Rcpp::stop("SEM factor covariance prior degrees of freedom are too small");
  }
  x_cov_ = c_ ? arma::cov(x_) : arma::mat(0, 0);
}

arma::mat SEMSampler::subject_means(const GroupState& state) const {
  arma::mat out = state.K * x_.t() + state.lambda * state.eta.t();
  out.each_col() += state.tmu;
  return out;
}

void SEMSampler::moments(GroupState& state) const {
  state.subj_mu = subject_means(state);
  state.tvar = arma::diagmat(1 / state.epsilon_inv);
  state.tvinv = arma::diagmat(state.epsilon_inv);
  state.population_var = state.K * x_cov_ * state.K.t() + state.tvar;
  if (q_) {
    const arma::mat structural = arma::solve(arma::eye(q_, q_) - state.B, arma::eye(q_, q_));
    const arma::mat effects = state.K + state.lambda * structural * state.G;
    state.population_var = effects * x_cov_ * effects.t() + state.lambda * structural *
      group_inverse(state.delta_inv) * structural.t() * state.lambda.t() + state.tvar;
  }
}

GroupState SEMSampler::initialize(Random& rng, SEXP start_mu, SEXP start_var) const {
  GroupState state;
  state.tmu = mean_;
  if (Rf_isNull(start_mu)) {
    for (arma::uword j = 0; j < p_; ++j) state.tmu[j] += rng.rnorm(0, 1 / std::sqrt(mean_precision_[j]));
  } else state.tmu = Rcpp::as<arma::vec>(start_mu);
  state.lambda = lambda_; state.B = b_; state.K = k_; state.G = g_;
  state.eta.zeros(n_, q_); state.delta_inv.eye(q_, q_); state.epsilon_inv.ones(p_);
  moments(state);
  state.tvar = Rf_isNull(start_var) ? group_inverse_wishart(rng, 3 * p_, arma::eye(p_, p_)) :
    Rcpp::as<arma::mat>(start_var);
  state.tvinv = arma::pinv(state.tvar);
  return state;
}

GroupState SEMSampler::read_state(const Rcpp::List& last) const {
  GroupState state;
  state.tmu = Rcpp::as<arma::vec>(last["tmu"]);
  state.lambda = Rcpp::as<arma::mat>(last["lambda"]);
  state.B = Rcpp::as<arma::mat>(last["B"]);
  state.K = Rcpp::as<arma::mat>(last["K"]);
  state.G = Rcpp::as<arma::mat>(last["G"]);
  state.eta = Rcpp::as<arma::mat>(last["eta"]);
  state.delta_inv = Rcpp::as<arma::mat>(last["delta_inv"]);
  state.epsilon_inv = Rcpp::as<arma::vec>(last["epsilon_inv"]);
  moments(state);
  return state;
}

void SEMSampler::draw(Random& rng, const arma::mat& alpha, GroupState& state) const {
  GroupState next = state;
  const arma::mat y = alpha.t();
  const arma::mat residual = y - x_ * state.K.t() - state.eta * state.lambda.t();
  const arma::vec mu_var = 1 / (n_ * state.epsilon_inv + mean_precision_);
  const arma::vec mu_mean = mu_var % (state.epsilon_inv % arma::sum(residual, 0).t() + mean_precision_ % mean_);
  for (arma::uword j = 0; j < p_; ++j) next.tmu[j] = rng.rnorm(mu_mean[j], std::sqrt(mu_var[j]));
  arma::mat centered = y;
  centered.each_row() -= next.tmu.t();
  if (q_) {
    const arma::mat b0 = arma::eye(q_, q_) - state.B;
    const arma::mat eta_precision = b0.t() * state.delta_inv * b0;
    const arma::mat eta_var = group_inverse(eta_precision + state.lambda.t() *
      arma::diagmat(state.epsilon_inv) * state.lambda);
    const arma::mat eta_mean = eta_var * (state.lambda.t() * arma::diagmat(state.epsilon_inv) *
      (centered - x_ * state.K.t()).t() + b0.t() * state.delta_inv * state.G * x_.t());
    const arma::mat root = arma::chol(eta_var, "lower");
    arma::vec z(q_);
    for (arma::uword s = 0; s < n_; ++s) {
      for (double& value : z) value = rng.rnorm(0, 1);
      next.eta.row(s) = (eta_mean.col(s) + root * z).t();
    }
  }
  const arma::vec squares = arma::sum(arma::square(centered - x_ * state.K.t() - next.eta * state.lambda.t()), 0).t();
  for (arma::uword j = 0; j < p_; ++j)
    next.epsilon_inv[j] = rng.rgamma(a_e_[j] + .5 * n_, 1 / (b_e_[j] + .5 * squares[j]));
  const arma::mat design = arma::join_rows(x_, next.eta);
  arma::mat loading = arma::join_rows(state.K, state.lambda);
  for (arma::uword j = 0; j < p_; ++j) {
    const arma::uvec& free = loading_free_[j];
    if (free.empty()) continue;
    const arma::mat z = design.cols(free);
    arma::vec coefficients = loading.row(j).t();
    const arma::vec target = centered.col(j) - design * coefficients + z * coefficients.elem(free);
    arma::mat precision = next.epsilon_inv[j] * z.t() * z;
    precision.diag() += arma::vec(loading_precision_.row(j).t()).elem(free);
    const arma::mat variance = group_inverse(precision);
    const arma::vec linear = next.epsilon_inv[j] * z.t() * target + arma::vec(loading_linear_.row(j).t()).elem(free);
    coefficients.elem(free) = group_normal(rng, variance * linear, variance);
    loading.row(j) = coefficients.t();
  }
  next.K = loading.head_cols(c_); next.lambda = loading.tail_cols(q_);
  arma::mat structural = arma::join_rows(state.G, state.B);
  arma::mat eta_residual = next.eta - design * structural.t();
  for (arma::uword j = 0; j < q_; ++j) {
    const arma::uvec& free = structural_free_[j];
    if (free.empty()) continue;
    const arma::mat z = design.cols(free);
    arma::vec coefficients = structural.row(j).t();
    const double precision_j = state.delta_inv(j, j);
    const arma::vec fixed = design * coefficients - z * coefficients.elem(free);
    const arma::vec target = precision_j * (next.eta.col(j) - fixed) +
      eta_residual * state.delta_inv.col(j) - precision_j * eta_residual.col(j);
    arma::mat precision = precision_j * z.t() * z;
    precision.diag() += structural_precision_.elem(free);
    const arma::mat variance = group_inverse(precision);
    coefficients.elem(free) = group_normal(rng, variance * z.t() * target, variance);
    if (cyclic_) {
      arma::mat candidate_b = structural.tail_cols(q_);
      candidate_b.row(j) = coefficients.tail(q_).t();
      const double ratio = n_ * (structural_logdet(candidate_b) - structural_logdet(structural.tail_cols(q_)));
      if (!std::isfinite(ratio) || std::log(rng.runif(0, 1)) >= ratio) continue;
    }
    structural.row(j) = coefficients.t();
    eta_residual.col(j) = next.eta.col(j) - design * coefficients;
  }
  next.G = structural.head_cols(c_); next.B = structural.tail_cols(q_);
  for (const auto& indices : factor_groups_) {
    const arma::mat residuals = eta_residual.cols(indices);
    const arma::mat squares = residuals.t() * residuals;
    if (indices.n_elem == 1) {
      const arma::uword j = indices[0];
      next.delta_inv(j, j) = rng.rgamma(a_d_[j] + .5 * n_, 1 / (b_d_[j] + .5 * squares[0]));
    } else {
      next.delta_inv.submat(indices, indices) = group_inverse(group_inverse_wishart(rng, 
        a_d_[indices[0]] + n_, arma::diagmat(b_d_.elem(indices)) + squares));
    }
  }
  moments(next);
  state = std::move(next);
}

double SEMSampler::affine(GroupState& state, arma::uword j,
                          double log_scale, double shift) const {
  double ratio = 0;
  if (log_scale != 0) {
    const double multiplier = std::exp(-2 * log_scale);
    ratio = -2 * a_e_[j] * log_scale - b_e_[j] * state.epsilon_inv[j] * (multiplier - 1);
    state.epsilon_inv[j] *= multiplier;
    const double variance = 1 / state.epsilon_inv[j];
    state.population_var(j, j) += variance - state.tvar(j, j);
    state.tvar(j, j) = variance;
    state.tvinv(j, j) = state.epsilon_inv[j];
  }
  if (shift != 0) {
    ratio -= shift * mean_precision_[j] * (state.tmu[j] - mean_[j]) +
      .5 * shift * shift * mean_precision_[j];
    state.tmu[j] += shift;
    state.subj_mu.row(j) += shift;
  }
  return ratio;
}

}  // namespace emc2
