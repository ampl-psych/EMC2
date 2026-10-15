#include "group.h"
#include "sem.h"
#include "infnt_factor.h"
#include <algorithm>
#include <stdexcept>

namespace emc2 {
namespace {

bool present(const Rcpp::List& x, const char* name) {
  return x.containsElementNamed(name) && !Rf_isNull(x[name]);
}

}  // namespace

arma::mat group_inverse(const arma::mat& x) {
  arma::mat out;
  if (!arma::inv_sympd(out, x))
    throw std::runtime_error("singular group covariance or precision");
  return out;
}

arma::vec group_normal(Random& rng, const arma::vec& mean, const arma::mat& covariance) {
  arma::mat root;
  if (!arma::chol(root, covariance, "lower"))
    throw std::runtime_error("singular group covariance");
  arma::vec z(mean.n_elem);
  for (double& value : z) value = rng.rnorm(0, 1);
  return mean + root * z;
}

// Bartlett decomposition, with the same scale convention as riwish().
arma::mat group_inverse_wishart(Random& rng, double df, const arma::mat& scale) {
  const arma::uword d = scale.n_rows;
  arma::mat root;
  if (!arma::chol(root, group_inverse(scale)))
    throw std::runtime_error("singular inverse-Wishart scale");
  arma::mat z(d, d, arma::fill::zeros);
  for (arma::uword k = 0; k < d; ++k)
    z(k, k) = std::sqrt(rng.rchisq(df - k));
  for (arma::uword col = 1; col < d; ++col)
    for (arma::uword row = 0; row < col; ++row)
      z(row, col) = rng.rnorm(0, 1);
  arma::mat factor = z * root;
  return group_inverse(factor.t() * factor);
}

namespace {
arma::vec get_vector(const Rcpp::List& x, const char* name) {
  return present(x, name) ? Rcpp::as<arma::vec>(x[name]) : arma::vec();
}

std::vector<std::string> get_names(SEXP x) {
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  return Rf_isNull(names) ? std::vector<std::string>() :
    Rcpp::as<std::vector<std::string>>(names);
}

}  // namespace

GroupSampler::GroupSampler(const Rcpp::List& sampler) {
  const std::string type = Rcpp::as<std::string>(sampler["type"]);
  single_ = type == "single";
  location_ = !present(sampler, "group_designs");
  n_ = Rcpp::as<int>(sampler["n_subjects"]);
  p_ = Rcpp::as<int>(sampler["n_pars"]);
  if (type == "SEM") { sem_.reset(new SEMSampler(sampler)); return; }
  if (type == "infnt_factor") { infnt_.reset(new InfiniteFactorSampler(sampler)); return; }
  diagonal_gamma_ = type == "diagonal-gamma";
  const Rcpp::List prior = sampler["prior"];
  mean_ = Rcpp::as<arma::vec>(prior["theta_mu_mean"]);
  prior_var_ = Rcpp::as<arma::mat>(prior["theta_mu_var"]);
  mean_names_ = get_names(prior["theta_mu_mean"]);

  p_ = Rcpp::as<int>(sampler["n_pars"]);
  if (present(sampler, "par_names")) {
    Rcpp::CharacterVector names = sampler["par_names"];
    for (arma::uword k = 0; k < p_; ++k)
      parameter_names_.push_back(Rcpp::as<std::string>(names[k]));
  } else if (present(sampler, "samples")) {
    const Rcpp::List samples = sampler["samples"];
    if (present(samples, "theta_var")) {
      SEXP dims = Rf_getAttrib(samples["theta_var"], R_DimNamesSymbol);
      if (!Rf_isNull(dims)) {
        Rcpp::List dn(dims);
        if (!Rf_isNull(dn[0]))
          parameter_names_ = Rcpp::as<std::vector<std::string>>(dn[0]);
      }
    }
  }
  if (parameter_names_.empty() && mean_names_.size() == p_)
    parameter_names_ = mean_names_;

  if (single_) {
    if (mean_.n_elem != p_ || prior_var_.n_rows != p_ || prior_var_.n_cols != p_)
      Rcpp::stop("single prior dimensions do not match parameters");
    prior_precision_ = arma::pinv(prior_var_);
    return;
  }
  prior_precision_ = present(prior, "theta_mu_invar") ?
    Rcpp::as<arma::mat>(prior["theta_mu_invar"]) : arma::pinv(prior_var_);
  prior_linear_ = prior_precision_ * mean_;
  if (diagonal_gamma_) {
    shape_ = Rcpp::as<arma::vec>(prior["shape"]);
    rate_ = Rcpp::as<arma::vec>(prior["rate"]);
  } else {
    v_ = Rcpp::as<double>(prior["v"]);
    const arma::vec a = Rcpp::as<arma::vec>(prior["A"]);
    a_inv_sq_.set_size(p_);
    for (arma::uword k = 0; k < p_; ++k)
      a_inv_sq_[k] = 1 / std::pow(a[k % a.n_elem], 2);
  }

  Rcpp::List designs = present(sampler, "group_designs") ?
    Rcpp::as<Rcpp::List>(sampler["group_designs"]) : Rcpp::List();
  arma::uword offset = 0;
  for (arma::uword k = 0; k < p_; ++k) {
    SEXP design = R_NilValue;
    if (parameter_names_.size() == p_ &&
        designs.containsElementNamed(parameter_names_[k].c_str()))
      design = designs[parameter_names_[k]];
    arma::mat x(n_, 1, arma::fill::ones);
    if (!Rf_isNull(design)) {
      if (TYPEOF(design) == VECSXP && !Rf_inherits(design, "data.frame")) {
        const Rcpp::List components(design);
        x = Rcpp::as<arma::mat>(components["fixed"]);
        if (present(components, "random")) {
          const Rcpp::List parts = components["random"];
          for (R_xlen_t j = 0; j < parts.size(); ++j) {
            arma::mat z = Rcpp::as<arma::mat>(parts[j]);
            if (z.n_rows != n_ || z.n_cols == 0)
              Rcpp::stop("random group design dimensions do not match subjects");
            random_.push_back({k, random_size_, z, z.t() * z});
            random_size_ += z.n_cols;
          }
        }
      } else {
        x = Rcpp::as<arma::mat>(design);
      }
    }
    if (x.n_rows != n_ || x.n_cols == 0)
      Rcpp::stop("fixed group design dimensions do not match subjects");
    starts_.push_back(offset);
    offset += x.n_cols;
    fixed_.push_back(std::move(x));
  }
  if (offset != mean_.n_elem)
    Rcpp::stop("group prior dimensions do not match fixed group designs");
  // X_k' X_l never changes as covariance and means are sampled.
  for (arma::uword k = 0; k < p_; ++k)
    for (arma::uword l = 0; l <= k; ++l)
      crossproducts_.push_back(fixed_[k].t() * fixed_[l]);

  if (!random_.empty()) {
    if (present(prior, "v_Z")) v_z_ = Rcpp::as<double>(prior["v_Z"]);
    const arma::vec az = present(prior, "A_z") ?
      Rcpp::as<arma::vec>(prior["A_z"]) : arma::vec(1, arma::fill::value(.3));
    az_inv_sq_.set_size(random_.size());
    az_init_inv_sq_.set_size(random_.size());
    for (arma::uword k = 0; k < random_.size(); ++k) {
      az_inv_sq_[k] = 1 / std::pow(az[k < az.n_elem ? k : 0], 2);
      az_init_inv_sq_[k] = 1 / std::pow(az[k % az.n_elem], 2);
    }
  }

  Rcpp::LogicalVector blocked = present(sampler, "is_blocked") &&
    Rf_xlength(sampler["is_blocked"]) > 0 ?
    Rcpp::as<Rcpp::LogicalVector>(sampler["is_blocked"]) :
    Rcpp::LogicalVector(p_, type != "diagonal");
  Rcpp::IntegerVector groups = present(sampler, "par_group") ?
    Rcpp::as<Rcpp::IntegerVector>(sampler["par_group"]) :
    Rcpp::IntegerVector(p_, 1);
  if (blocked.size() < static_cast<R_xlen_t>(p_) || groups.size() < static_cast<R_xlen_t>(p_))
    Rcpp::stop("is_blocked (%d) and par_group (%d) must have one entry per parameter (%d)",
               (int) blocked.size(), (int) groups.size(), (int) p_);
  std::vector<int> labels;
  std::vector<std::vector<arma::uword>> indices;
  std::vector<arma::uword> unblocked;
  block_dimensions_.ones(p_);
  for (arma::uword k = 0; k < p_; ++k) {
    if (!blocked[k]) {
      unblocked.push_back(k);
      continue;
    }
    auto found = std::find(labels.begin(), labels.end(), groups[k]);
    if (found == labels.end()) {
      labels.push_back(groups[k]);
      indices.push_back({k});
    } else {
      indices[found - labels.begin()].push_back(k);
    }
  }
  for (const auto& block : indices) {
    blocks_.push_back(arma::uvec(block));
    block_dimensions_.elem(blocks_.back()).fill(block.size());
  }
  unblocked_ = arma::uvec(unblocked);
}

GroupSampler::~GroupSampler() = default;

arma::mat GroupSampler::fixed_means(const arma::vec& beta) const {
  arma::mat out(p_, n_);
  for (arma::uword k = 0; k < p_; ++k)
    out.row(k) = (fixed_[k] * beta.subvec(starts_[k],
      starts_[k] + fixed_[k].n_cols - 1)).t();
  return out;
}

arma::mat GroupSampler::random_means(const arma::vec& u) const {
  arma::mat out(p_, n_, arma::fill::zeros);
  for (const auto& part : random_)
    out.row(part.parameter) += (part.z *
      u.subvec(part.start, part.start + part.z.n_cols - 1)).t();
  return out;
}

arma::mat GroupSampler::subject_means(const GroupState& state) const {
  if (sem_) return sem_->subject_means(state);
  if (infnt_) return infnt_->subject_means(state);
  if (single_) return arma::repmat(state.tmu, 1, n_);
  arma::mat out = fixed_means(state.tmu);
  if (!random_.empty()) out += random_means(state.u);
  return out;
}

GroupState GroupSampler::initialize(Random& rng, SEXP start_mu, SEXP start_var) const {
  if (sem_) return sem_->initialize(rng, start_mu, start_var);
  if (infnt_) return infnt_->initialize(rng, start_mu, start_var);
  GroupState state;
  if (Rf_isNull(start_mu)) {
    state.tmu = mean_;
    if (!single_) {
      arma::vec values, z(mean_.n_elem);
      arma::mat vectors;
      if (!arma::eig_sym(values, vectors, prior_var_))
        throw std::runtime_error("singular initial group prior covariance");
      for (double& value : z) value = rng.rnorm(0, 1);
      state.tmu += vectors * (arma::sqrt(arma::clamp(values, 0, arma::datum::inf)) % z);
    }
  } else state.tmu = Rcpp::as<arma::vec>(start_mu);
  if (!Rf_isNull(start_var)) state.tvar = Rcpp::as<arma::mat>(start_var);
  else if (single_) state.tvar = prior_var_;
  else if (diagonal_gamma_) {
    state.tvar.zeros(p_, p_);
    for (arma::uword j = 0; j < p_; ++j) state.tvar(j, j) = 1 / rng.rgamma(10, .2);
  } else state.tvar = group_inverse_wishart(rng, 3 * p_, arma::eye(p_, p_));
  if (single_) {
    if (!arma::inv_sympd(state.tvinv, state.tvar)) state.tvinv = arma::pinv(state.tvar);
  } else {
    state.tvinv = arma::pinv(state.tvar, std::sqrt(std::numeric_limits<double>::epsilon()) *
                          arma::norm(state.tvar, 2));
  }
  if (!single_ && !diagonal_gamma_) {
    state.a_half.set_size(p_);
    for (double& value : state.a_half) value = 1 / rng.rgamma(2, 1);
    state.u.set_size(random_size_);
    state.s.set_size(random_.size());
    state.a_half_s.set_size(random_.size());
    for (arma::uword k = 0; k < random_.size(); ++k)
      state.a_half_s[k] = 1 / rng.rgamma(.5, 1 / az_init_inv_sq_[k]);
    for (arma::uword k = 0; k < random_.size(); ++k)
      state.s[k] = 1 / rng.rgamma(v_z_ / 2, state.a_half_s[k] / v_z_);
    for (arma::uword k = 0; k < random_.size(); ++k)
      for (arma::uword j = 0; j < random_[k].z.n_cols; ++j)
        state.u[random_[k].start + j] = rng.rnorm(0, std::sqrt(state.s[k]));
  }
  state.subj_mu = subject_means(state);
  return state;
}

GroupState GroupSampler::read_state(const Rcpp::List& last) const {
  if (sem_) return sem_->read_state(last);
  if (infnt_) return infnt_->read_state(last);
  GroupState state;
  state.tmu = Rcpp::as<arma::vec>(last["tmu"]);
  state.tvar = Rcpp::as<arma::mat>(last["tvar"]);
  state.tvinv = present(last, "tvinv") ? Rcpp::as<arma::mat>(last["tvinv"]) :
    (single_ ? arma::pinv(state.tvar) : group_inverse(state.tvar));
  state.a_half = get_vector(last, "a_half");
  state.u = get_vector(last, "u");
  state.s = get_vector(last, "s");
  state.a_half_s = get_vector(last, "a_half_s");
  state.subj_mu = subject_means(state);
  return state;
}

void GroupSampler::draw(Random& rng, const arma::mat& alpha, GroupState& state) const {
  if (single_) return;
  if (sem_) { sem_->draw(rng, alpha, state); return; }
  if (infnt_) { infnt_->draw(rng, alpha, state); return; }
  if (diagonal_gamma_) {
    const arma::vec precision = state.tvinv.diag();
    arma::mat mean_precision = prior_precision_;
    mean_precision.diag() += n_ * precision;
    const arma::mat variance = group_inverse(mean_precision);
    state.tmu = group_normal(rng, variance * (arma::sum(alpha, 1) % precision + prior_linear_), variance);
    state.tvar.zeros(p_, p_); state.tvinv.zeros(p_, p_);
    for (arma::uword j = 0; j < p_; ++j) {
      const arma::rowvec residual = alpha.row(j) - state.tmu[j];
      state.tvinv(j, j) = rng.rgamma(shape_[j % shape_.n_elem] + .5 * n_,
        1 / (rate_[j % rate_.n_elem] + .5 * arma::dot(residual, residual)));
      state.tvar(j, j) = 1 / state.tvinv(j, j);
    }
    state.subj_mu = subject_means(state);
    return;
  }
  // Keep the old state intact if a covariance factorization fails.
  GroupState next = state;
  arma::mat residual = alpha;
  if (!random_.empty()) residual -= random_means(state.u);
  arma::mat weighted = state.tvinv * residual;
  arma::mat precision = prior_precision_;
  arma::vec linear = prior_linear_;
  arma::uword cross_index = 0;
  for (arma::uword k = 0; k < p_; ++k) {
    const arma::uword first = starts_[k], last = first + fixed_[k].n_cols - 1;
    linear.subvec(first, last) += fixed_[k].t() * weighted.row(k).t();
    for (arma::uword l = 0; l <= k; ++l) {
      const arma::uword other = starts_[l], end = other + fixed_[l].n_cols - 1;
      const arma::mat block = state.tvinv(k, l) * crossproducts_[cross_index++];
      precision.submat(first, other, last, end) += block;
      if (l != k) precision.submat(other, first, end, last) += block.t();
    }
  }
  const arma::mat covariance = group_inverse(precision);
  next.tmu = group_normal(rng, covariance * linear, covariance);
  next.subj_mu = fixed_means(next.tmu);
  if (!random_.empty()) {
    next.subj_mu += random_means(next.u);
    residual = alpha - next.subj_mu;
    for (arma::uword j = 0; j < random_.size(); ++j) {
      const auto& part = random_[j];
      const arma::uword k = part.parameter, first = part.start;
      const arma::uword last = first + part.z.n_cols - 1;
      const double noise_precision = state.tvinv(k, k);
      arma::vec old_offset = part.z * next.u.subvec(first, last);
      arma::vec target = (state.tvinv.row(k) * residual).t() /
        noise_precision + old_offset;
      arma::mat post_precision = noise_precision * part.ztz;
      post_precision.diag() += 1 / next.s[j];
      const arma::mat post_covariance = group_inverse(post_precision);
      arma::vec u = group_normal(rng, post_covariance * (noise_precision * part.z.t() * target),
                                post_covariance);
      next.u.subvec(first, last) = u;
      arma::rowvec change = (part.z * u - old_offset).t();
      next.subj_mu.row(k) += change;
      residual.row(k) -= change;
      next.s[j] = 1 / rng.rgamma((v_z_ + u.n_elem) / 2,
        1 / (v_z_ / next.a_half_s[j] + .5 * arma::dot(u, u)));
      next.a_half_s[j] = 1 / rng.rgamma((v_z_ + 1) / 2,
        1 / (v_z_ / next.s[j] + az_inv_sq_[j]));
    }
  } else residual = alpha - next.subj_mu;
  next.tvar.zeros(p_, p_);
  for (const auto& block : blocks_) {
    arma::mat r = residual.rows(block);
    arma::mat scale = r * r.t();
    scale.diag() += 2 * v_ / state.a_half.elem(block);
    next.tvar.submat(block, block) = group_inverse_wishart(rng, v_ + block.n_elem - 1 + n_, scale);
  }
  for (const auto k : unblocked_) {
    const double rate = v_ / state.a_half[k] + .5 * arma::dot(residual.row(k), residual.row(k));
    next.tvar(k, k) = 1 / rng.rgamma((v_ + n_) / 2, 1 / rate);
  }
  next.tvinv = group_inverse(next.tvar);
  for (arma::uword k = 0; k < p_; ++k)
    next.a_half[k] = 1 / rng.rgamma((v_ + block_dimensions_[k]) / 2,
      1 / (v_ * next.tvinv(k, k) + a_inv_sq_[k]));
  state = std::move(next);
}

double GroupSampler::affine(GroupState& state, arma::uword j, double log_scale, double shift) const {
  if (sem_) return sem_->affine(state, j, log_scale, shift);
  if (infnt_) return infnt_->affine(state, j, log_scale, shift);
  if (log_scale != 0) {
    const double factor = std::exp(log_scale);
    const double ratio = diagonal_gamma_ ?
      -2 * shape_[j % shape_.n_elem] * log_scale - rate_[j % rate_.n_elem] *
        state.tvinv(j, j) * (1 / (factor * factor) - 1) :
      log_scale - (factor * factor - 1) * a_inv_sq_[j] / state.a_half[j];
    state.tvar.row(j) *= factor; state.tvar.col(j) *= factor;
    state.tvinv.row(j) /= factor; state.tvinv.col(j) /= factor;
    if (!diagonal_gamma_) state.a_half[j] /= factor * factor;
    return ratio;
  }
  const double ratio = -shift * arma::dot(prior_precision_.row(j), state.tmu - mean_) -
    .5 * shift * shift * prior_precision_(j, j);
  state.tmu[j] += shift; state.subj_mu.row(j) += shift;
  return ratio;
}

}  // namespace emc2
