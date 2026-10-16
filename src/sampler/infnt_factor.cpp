#include "infnt_factor.h"
#include "group.h"

namespace emc2 {
namespace {

arma::vec recycle(SEXP x, arma::uword n) {
  const arma::vec values = Rcpp::as<arma::vec>(x);
  arma::vec out(n);
  for (arma::uword j = 0; j < n; ++j) out[j] = values[j % values.n_elem];
  return out;
}

arma::mat precision_root(const arma::mat& precision) {
  arma::mat root;
  if (!arma::chol(root, precision))
    throw std::runtime_error("singular factor precision");
  return root;
}

// Draw from N(P^-1 b, P^-1) with P = root' root: solve root' y = b, then
// root x = y + z. The non-throwing solves give the retry loop a "singular"
// error where Armadillo's would stop the stage.
arma::mat triangular_draw(const arma::mat& root, const arma::mat& linear, const arma::mat& z) {
  arma::mat mean, out;
  if (!arma::solve(mean, arma::trimatl(root.t()), linear))
    throw std::runtime_error("singular factor precision");
  const arma::mat shifted = mean + z;
  if (!arma::solve(out, arma::trimatu(root), shifted))
    throw std::runtime_error("singular factor precision");
  return out;
}

}  // namespace

InfiniteFactorSampler::InfiniteFactorSampler(const Rcpp::List& sampler)
  : p_(Rcpp::as<int>(sampler["n_pars"])),
    n_(Rcpp::as<int>(sampler["n_subjects"])),
    factors_(Rcpp::as<int>(sampler.attr("max_factors"))) {
  const Rcpp::List prior = sampler["prior"];
  mean_ = recycle(prior["theta_mu_mean"], p_);
  variance_ = recycle(prior["theta_mu_var"], p_);
  precision_ = 1 / variance_;
  shape_ = recycle(prior["as"], p_);
  rate_ = recycle(prior["bs"], p_);
  df_ = Rcpp::as<double>(prior["df"]);
  ad1_ = Rcpp::as<double>(prior["ad1"]);
  bd1_ = Rcpp::as<double>(prior["bd1"]);
  ad2_ = Rcpp::as<double>(prior["ad2"]);
  bd2_ = Rcpp::as<double>(prior["bd2"]);
}

arma::mat InfiniteFactorSampler::subject_means(const GroupState& state) const {
  return arma::repmat(state.tmu, 1, n_);
}

GroupState InfiniteFactorSampler::initialize(Random& rng, SEXP start_mu, SEXP start_var) const {
  GroupState state;
  if (Rf_isNull(start_mu)) {
    state.tmu.set_size(p_);
    for (arma::uword j = 0; j < p_; ++j)
      state.tmu[j] = rng.rnorm(mean_[j], std::sqrt(variance_[j]));
  } else state.tmu = Rcpp::as<arma::vec>(start_mu);
  state.tvar = Rf_isNull(start_var) ?
    group_inverse_wishart(rng, 3 * p_, arma::eye(p_, p_)) : Rcpp::as<arma::mat>(start_var);
  state.tvinv = arma::pinv(state.tvar,
    std::sqrt(std::numeric_limits<double>::epsilon()) * arma::norm(state.tvar, 2));
  state.lambda.zeros(p_, factors_);
  state.eta.zeros(n_, factors_);
  state.epsilon_inv.set_size(p_);
  for (arma::uword j = 0; j < p_; ++j)
    state.epsilon_inv[j] = rng.rgamma(shape_[j], 1 / rate_[j]);
  state.psi.set_size(p_, factors_);
  for (double& value : state.psi) value = rng.rgamma(df_ / 2, 2 / df_);
  state.delta.set_size(factors_);
  for (arma::uword h = 0; h < factors_; ++h)
    state.delta[h] = rng.rgamma(h == 0 ? ad1_ : ad2_, 1 / (h == 0 ? bd1_ : bd2_));
  state.subj_mu = subject_means(state);
  return state;
}

GroupState InfiniteFactorSampler::read_state(const Rcpp::List& last) const {
  GroupState state;
  state.tmu = Rcpp::as<arma::vec>(last["tmu"]);
  state.tvar = Rcpp::as<arma::mat>(last["tvar"]);
  state.tvinv = group_inverse(state.tvar);
  state.lambda = Rcpp::as<arma::mat>(last["lambda"]);
  state.eta = Rcpp::as<arma::mat>(last["eta"]);
  state.psi = Rcpp::as<arma::mat>(last["psi"]);
  state.delta = Rcpp::as<arma::vec>(last["delta"]);
  state.epsilon_inv = Rcpp::as<arma::vec>(last["epsilon_inv"]);
  state.subj_mu = subject_means(state);
  return state;
}

void InfiniteFactorSampler::draw(Random& rng, const arma::mat& alpha, GroupState& state) const {
  GroupState next = state;
  arma::vec tau = arma::cumprod(state.delta);
  arma::mat centered = alpha.each_col() - state.tmu;

  // Refresh the latent factors after the marginal participant update.
  arma::mat weighted_loadings = state.lambda.each_col() % state.epsilon_inv;
  arma::mat eta_precision = state.lambda.t() * weighted_loadings;
  eta_precision.diag() += 1;
  const arma::mat eta_root = precision_root(eta_precision);
  arma::mat noise(n_, factors_);
  for (double& value : noise) value = rng.rnorm(0, 1);
  next.eta = triangular_draw(eta_root, weighted_loadings.t() * centered, noise.t()).t();

  const arma::vec mu_variance = 1 / (n_ * state.epsilon_inv + precision_);
  const arma::vec mu_mean = mu_variance % (state.epsilon_inv %
    arma::sum(alpha - state.lambda * next.eta.t(), 1) + precision_ % mean_);
  for (arma::uword j = 0; j < p_; ++j)
    next.tmu[j] = rng.rnorm(mu_mean[j], std::sqrt(mu_variance[j]));
  centered = alpha.each_col() - next.tmu;

  const arma::mat eta_cross = next.eta.t() * next.eta;
  for (arma::uword j = 0; j < p_; ++j) {
    arma::mat precision = state.epsilon_inv[j] * eta_cross;
    precision.diag() += state.psi.row(j).t() % tau;
    const arma::vec linear = state.epsilon_inv[j] * next.eta.t() * centered.row(j).t();
    const arma::mat root = precision_root(precision);
    arma::vec z(factors_);
    for (double& value : z) value = rng.rnorm(0, 1);
    next.lambda.row(j) = triangular_draw(root, linear, z).t();
  }

  for (arma::uword h = 0; h < factors_; ++h)
    for (arma::uword j = 0; j < p_; ++j)
      next.psi(j, h) = rng.rgamma((df_ + 1) / 2,
        2 / (df_ + tau[h] * std::pow(next.lambda(j, h), 2)));

  const arma::vec loading_squares = arma::sum(next.psi % arma::square(next.lambda), 0).t();
  for (arma::uword h = 0; h < factors_; ++h) {
    const double shape = (h == 0 ? ad1_ : ad2_) + .5 * p_ * (factors_ - h);
    const double rate = (h == 0 ? bd1_ : bd2_) + .5 *
      arma::dot(tau.tail(factors_ - h), loading_squares.tail(factors_ - h)) / next.delta[h];
    const double delta = rng.rgamma(shape, 1 / rate);
    tau.tail(factors_ - h) *= delta / next.delta[h];
    next.delta[h] = delta;
  }

  const arma::vec residual_squares = arma::sum(arma::square(centered - next.lambda * next.eta.t()), 1);
  for (arma::uword j = 0; j < p_; ++j)
    next.epsilon_inv[j] = rng.rgamma(shape_[j] + .5 * n_,
                                  1 / (rate_[j] + .5 * residual_squares[j]));
  next.tvar = next.lambda * next.lambda.t();
  next.tvar.diag() += 1 / next.epsilon_inv;
  next.tvinv = group_inverse(next.tvar);
  next.subj_mu = subject_means(next);
  state = std::move(next);
}

double InfiniteFactorSampler::affine(GroupState& state, arma::uword j,
                                     double log_scale, double shift) const {
  double ratio = -precision_[j] * shift * (state.tmu[j] - mean_[j] + .5 * shift);
  if (log_scale != 0) {
    const arma::vec tau = arma::cumprod(state.delta);
    const double loading_precision = arma::dot(
      state.psi.row(j).t() % tau, arma::square(state.lambda.row(j).t()));
    ratio += (factors_ - 2 * shape_[j]) * log_scale -
      .5 * std::expm1(2 * log_scale) * loading_precision -
      rate_[j] * state.epsilon_inv[j] * std::expm1(-2 * log_scale);
    const double scale = std::exp(log_scale);
    state.lambda.row(j) *= scale;
    state.epsilon_inv[j] /= scale * scale;
    state.tvar.row(j) *= scale;
    state.tvar.col(j) *= scale;
    state.tvinv.row(j) /= scale;
    state.tvinv.col(j) /= scale;
  }
  state.tmu[j] += shift;
  state.subj_mu.row(j) += shift;
  return ratio;
}

}  // namespace emc2
