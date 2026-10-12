#include "proposals.h"
#include "group.h"
#include <algorithm>
#include <numeric>

namespace emc2 {
namespace {

bool present(const Rcpp::List& x, const char* name) {
  return x.containsElementNamed(name) && !Rf_isNull(x[name]);
}

Rcpp::List shallow(SEXP x) { return Rcpp::List(Rf_shallow_duplicate(x)); }

arma::cube cube(SEXP x) {
  Rcpp::NumericVector values(x);
  Rcpp::IntegerVector d = values.attr("dim");
  return arma::cube(values.begin(), d[0], d[1], d[2]);
}

arma::uvec selected_indices(SEXP x, arma::uword n) {
  if (Rf_isNull(x)) return arma::regspace<arma::uvec>(0, n - 1);
  std::vector<arma::uword> out;
  if (TYPEOF(x) == LGLSXP) {
    Rcpp::LogicalVector keep(x);
    for (R_xlen_t k = 0; k < keep.size(); ++k) if (keep[k]) out.push_back(k);
  } else {
    Rcpp::IntegerVector keep = Rcpp::as<Rcpp::IntegerVector>(x);
    for (int k : keep) out.push_back(k - 1);
  }
  return arma::uvec(out);
}

bool positive_definite(const arma::mat& covariance) {
  arma::vec eigenvalues;
  return covariance.is_finite() && arma::eig_sym(eigenvalues, covariance) &&
    eigenvalues.min() >= 1e-8;
}

arma::mat fallback_covariance(const Rcpp::List& sampler) {
  Rcpp::List prior = sampler["prior"];
  Rcpp::NumericVector variance = prior["theta_mu_var"];
  const bool matrix = Rf_isMatrix(variance);
  arma::vec diagonal = matrix ? arma::mat(Rcpp::as<arma::mat>(variance)).diag() : Rcpp::as<arma::vec>(variance);
  Rcpp::CharacterVector names = sampler["par_names"], prior_names;
  if (matrix) {
    SEXP dn = variance.attr("dimnames");
    if (!Rf_isNull(dn)) { Rcpp::List dims(dn); if (!Rf_isNull(dims[0])) prior_names = dims[0]; }
  } else if (!Rf_isNull(variance.attr("names"))) prior_names = variance.attr("names");
  const double average = arma::mean(diagonal);
  arma::mat out(names.size(), names.size(), arma::fill::zeros);
  for (int j = 0; j < names.size(); ++j) {
    double value = average;
    if (prior_names.size()) {
      for (int k = 0; k < prior_names.size(); ++k)
        if (names[j] == prior_names[k]) { value = diagonal[k]; break; }
    } else if (diagonal.n_elem == names.size()) value = diagonal[j];
    out(j, j) = .1 * value;
  }
  return out;
}

void rescale_settings(Rcpp::List& sampler, double ratio) {
  Rcpp::List samples = shallow(sampler["samples"]);
  SEXP settings_sexp = samples.attr("pm_settings");
  if (!Rf_isNull(settings_sexp)) {
    Rcpp::List settings = Rcpp::clone(Rcpp::List(settings_sexp));
    for (R_xlen_t s = 0; s < settings.size(); ++s) {
      Rcpp::List blocks = settings[s];
      for (R_xlen_t b = 0; b < blocks.size(); ++b) {
        Rcpp::List block = blocks[b];
        if (present(block, "epsilon")) {
          Rcpp::NumericVector epsilon = block["epsilon"];
          for (double& e : epsilon) e *= ratio;
          block["epsilon"] = epsilon;
        }
        blocks[b] = block;
      }
      settings[s] = blocks;
    }
    samples.attr("pm_settings") = settings;
  }
  sampler["samples"] = samples;
}

// Cache the empirical Gaussian once; preserve condMVN's existing jitter when
// its joint covariance is nearly singular.
class Conditional {
 public:
  Conditional(const arma::mat& alpha, const arma::mat& given) : p_(alpha.n_rows), given_(given) {
    if (given_.n_rows) {
      arma::uvec varying = arma::find(arma::max(given_, 1) != arma::min(given_, 1));
      given_ = given_.rows(varying);
    }
    arma::mat draws = arma::join_cols(alpha, given_);
    mean_ = arma::mean(draws, 1);
    covariance_ = arma::cov(draws.t());
    if (!given_.n_rows) {
      variance_ = covariance_;
      return;
    }
    arma::vec eigenvalues;
    if (!arma::eig_sym(eigenvalues, covariance_)) Rcpp::stop("invalid proposal covariance");
    jitter_ = eigenvalues.min() < 1e-8;
    if (!jitter_) condition(covariance_);
  }

  Rcpp::List at(int iteration) {
    arma::vec mu = mean_.head(p_);
    if (given_.n_rows) {
      if (jitter_) {
        arma::mat covariance = covariance_;
        for (arma::uword j = 0; j < covariance.n_rows; ++j)
          covariance(j, j) += std::abs(R::rnorm(0, 1e-3));
        condition(covariance);
      }
      mu += regression_ * (given_.col(iteration) - mean_.tail(given_.n_rows));
    }
    return Rcpp::List::create(Rcpp::_ ["eff_mu"] = Rcpp::NumericVector(mu.begin(), mu.end()), Rcpp::_ ["eff_var"] = variance_);
  }

 private:
  arma::uword p_;
  bool jitter_ = false;
  arma::vec mean_;
  arma::mat given_, covariance_, regression_, variance_;

  void condition(const arma::mat& covariance) {
    const arma::uword end = covariance.n_rows - 1;
    arma::mat inverse;
    if (!arma::inv_sympd(inverse, covariance.submat(p_, p_, end, end)))
      Rcpp::stop("singular conditional proposal covariance");
    arma::mat cross = covariance.submat(0, p_, p_ - 1, end);
    regression_ = cross * inverse;
    variance_ = covariance.submat(0, 0, p_ - 1, p_ - 1) - regression_ * cross.t();
  }
};

struct ProposalBlock {
  arma::uvec global;
  std::vector<Conditional> subjects;
};

ProposalBlock make_block(const Rcpp::List& samples, const std::string& type,
                         const arma::uvec& selected, const arma::uvec& global, int only_subject = -1) {
  ProposalBlock block;
  block.global = global;
  arma::cube alpha = cube(samples["alpha"]);
  const arma::uword p = selected.n_elem, draws = alpha.n_slices;
  arma::mat unwound, means;
  arma::cube subject_means;
  if (type == "SEM" || type == "infnt_factor") {
    const arma::cube loadings = cube(samples["lambda"]), eta = cube(samples["eta"]);
    const arma::mat mu = Rcpp::as<arma::mat>(samples["theta_mu"]);
    const arma::mat epsilon = Rcpp::as<arma::mat>(samples["epsilon_inv"]);
    const arma::uword q = loadings.n_cols;
    for (arma::uword subject = only_subject < 0 ? 0 : only_subject;
         subject < (only_subject < 0 ? alpha.n_cols : only_subject + 1); ++subject) {
      arma::mat a(p, draws), given(2 * p + q + p * q, draws);
      for (arma::uword i = 0; i < draws; ++i) {
        const arma::vec values = alpha.slice(i).col(subject), errors = epsilon.col(i), means = mu.col(i);
        a.col(i) = values.elem(selected);
        arma::vec column = arma::join_cols(means.elem(selected), eta.slice(i).row(subject).t());
        column = arma::join_cols(column, type == "infnt_factor" ? arma::vec(arma::log(errors.elem(selected))) : arma::vec(errors.elem(selected)));
        given.col(i) = arma::join_cols(column, arma::vectorise(loadings.slice(i).rows(selected)));
      }
      block.subjects.emplace_back(a, given);
    }
    return block;
  }
  if (type != "single") {
    arma::cube covariance = cube(samples["theta_var"]);
    unwound.set_size(p * (p + 1) / 2, draws);
    for (arma::uword i = 0; i < draws; ++i) {
      arma::mat root;
      if (!arma::chol(root, covariance.slice(i).submat(selected, selected), "lower"))
        Rcpp::stop("singular group covariance in proposal history");
      root.diag() = arma::log(root.diag());
      arma::uword k = 0;
      for (arma::uword col = 0; col < p; ++col)
        for (arma::uword row = col; row < p; ++row) unwound(k++, i) = root(row, col);
    }
    if (present(samples, "subj_mu")) subject_means = cube(samples["subj_mu"]);
    else means = Rcpp::as<arma::mat>(samples["theta_mu"]);
  }
  for (arma::uword subject = only_subject < 0 ? 0 : only_subject;
       subject < (only_subject < 0 ? alpha.n_cols : only_subject + 1); ++subject) {
    arma::mat a(p, draws), given;
    if (type != "single") given.set_size(p + unwound.n_rows, draws);
    for (arma::uword i = 0; i < draws; ++i) {
      arma::vec value = alpha.slice(i).col(subject);
      a.col(i) = value.elem(selected);
      if (type != "single") {
        arma::vec mu = subject_means.n_elem ? subject_means.slice(i).col(subject) : means.col(i);
        given.col(i).head(p) = mu.elem(selected);
      }
    }
    if (type != "single") given.tail_rows(unwound.n_rows) = unwound;
    block.subjects.emplace_back(a, given);
  }
  return block;
}

}  // namespace

Rcpp::List chain_proposals(Rcpp::List emc, Rcpp::NumericVector alpha,
                                 Rcpp::NumericMatrix likelihood, bool do_block,
                                 std::string stage) {
  arma::cube draws = cube(alpha);
  Rcpp::List first = emc[0];
  Rcpp::List data = first["data"];
  Rcpp::IntegerVector components = data.attr("components");
  arma::mat fallback = fallback_covariance(first);
  const arma::uword p = draws.n_rows, n = draws.n_cols;
  Rcpp::List draw_names = alpha.attr("dimnames");
  SEXP par_names = draw_names.size() ? draw_names[0] : R_NilValue;
  const int count = static_cast<int>(std::nearbyint(draws.n_slices / 2.0));
  Rcpp::List out = shallow(emc);
  for (R_xlen_t ch = 0; ch < emc.size(); ++ch) {
    Rcpp::List chain = shallow(emc[ch]);
    // Match sample.int() without replacement, including R's sample.kind.
    std::vector<int> pool(draws.n_slices);
    std::iota(pool.begin(), pool.end(), 0);
    Rcpp::IntegerVector selected(count);
    for (int i = 0, remaining = pool.size(); i < count; ++i) {
      int pick = static_cast<int>(R_unif_index(remaining));
      selected[i] = pool[pick];
      pool[pick] = pool[--remaining];
    }
    Rcpp::List covariances(n), means(n);
    arma::mat average(p, p, arma::fill::zeros);
    int valid = 0;
    for (arma::uword s = 0; s < n; ++s) {
      arma::mat x(p, count);
      arma::vec weights(count);
      for (int i = 0; i < count; ++i) {
        x.col(i) = draws.slice(selected[i]).col(s);
        weights[i] = likelihood(s, selected[i]);
      }
      weights = arma::exp(weights - weights.max());
      weights /= arma::accu(weights);
      arma::vec weighted_mean = x * weights;
      Rcpp::NumericVector mu(weighted_mean.begin(), weighted_mean.end());
      if (!Rf_isNull(par_names)) mu.attr("names") = par_names;
      means[s] = mu;
      arma::mat covariance = arma::cov(x.t());
      if (do_block)
        for (arma::uword k = 0; k < p; ++k)
          for (arma::uword l = 0; l < p; ++l)
            if (components[k] != components[l]) covariance(k, l) = 0;
      if (positive_definite(covariance)) {
        Rcpp::NumericMatrix matrix = Rcpp::wrap(covariance);
        if (!Rf_isNull(par_names)) matrix.attr("dimnames") = Rcpp::List::create(par_names, par_names);
        covariances[s] = matrix;
        average += covariance;
        ++valid;
      }
    }
    if (valid) average /= valid;
    Rcpp::List previous = present(chain, "chains_var") ? Rcpp::as<Rcpp::List>(chain["chains_var"]) : Rcpp::List();
    for (arma::uword s = 0; s < n; ++s) if (Rf_isNull(covariances[s])) {
      arma::mat covariance = valid ? average :
        (s < static_cast<arma::uword>(previous.size()) && !Rf_isNull(previous[s]) ?
          Rcpp::as<arma::mat>(previous[s]) : fallback);
      covariances[s] = covariance;
      if (!valid) average += covariance / n;
    }
    const double prop_var = arma::mean(average.diag());
    SEXP old_prop_var = chain.attr("prop_var");
    const double ratio = Rf_isNull(old_prop_var) ? 2 :
      std::sqrt(Rcpp::as<double>(old_prop_var) / prop_var);
    const bool scaled = stage == "adapt" && Rcpp::as<std::string>(chain["type"]) != "single";
    if (stage != "sample" && !scaled) rescale_settings(chain, ratio);
    chain["chains_var"] = covariances;
    chain["chains_mu"] = means;
    Rf_copyMostAttrib(emc[ch], chain);
    chain.attr("prop_var") = prop_var;
    out[ch] = chain;
  }
  return out;
}

Rcpp::List efficient_proposals(Rcpp::List emc, Rcpp::List samples) {
  Rcpp::List first = emc[0], data = first["data"];
  Rcpp::IntegerVector components = data.attr("components");
  const arma::uword p = components.size();
  const int n = Rcpp::as<int>(first["n_subjects"]);
  const int iterations = Rcpp::as<int>(samples["iteration"]);
  const std::string type = Rcpp::as<std::string>(first["type"]);
  std::vector<int> labels;
  for (int component : components)
    if (std::find(labels.begin(), labels.end(), component) == labels.end()) labels.push_back(component);
  std::vector<ProposalBlock> blocks;
  for (int label : labels) {
    std::vector<arma::uword> selected;
    for (arma::uword k = 0; k < p; ++k) if (components[k] == label) selected.push_back(k);
    const arma::uvec block_indices(selected);
    blocks.push_back(make_block(samples, type, block_indices, block_indices));
  }
  Rcpp::List out = shallow(emc);
  for (R_xlen_t ch = 0; ch < emc.size(); ++ch) {
    const int iteration = static_cast<int>(std::nearbyint(iterations * (ch + 1.0) / emc.size())) - 1;
    arma::mat means(p, n, arma::fill::zeros);
    std::vector<arma::mat> covariances(n, arma::mat(p, p, arma::fill::zeros));
    for (auto& block : blocks) for (int s = 0; s < n; ++s) {
      Rcpp::List conditional = block.subjects[s].at(iteration);
      arma::vec mu = Rcpp::as<arma::vec>(conditional["eff_mu"]);
      for (arma::uword k = 0; k < block.global.n_elem; ++k) means(block.global[k], s) = mu[k];
      covariances[s].submat(block.global, block.global) = Rcpp::as<arma::mat>(conditional["eff_var"]);
    }
    Rcpp::List chain = shallow(emc[ch]), mus(n), vars(n);
    for (int s = 0; s < n; ++s) {
      mus[s] = Rcpp::NumericVector(means.colptr(s), means.colptr(s) + p);
      vars[s] = covariances[s];
    }
    chain["eff_mu"] = mus;
    chain["eff_var"] = vars;
    Rf_copyMostAttrib(emc[ch], chain);
    out[ch] = chain;
  }
  return out;
}

Rcpp::IntegerVector parameter_blocks(Rcpp::List emc, int n_blocks) {
  if (n_blocks < 1) Rcpp::stop("n_blocks must be positive");
  Rcpp::List first = emc[0], covariances = first["chains_var"], data = first["data"];
  arma::mat covariance = Rcpp::as<arma::mat>(covariances[0]);
  Rcpp::IntegerVector shared = data.attr("shared_ll_idx");
  Rcpp::IntegerVector out(shared.size());
  std::vector<int> labels;
  for (int label : shared)
    if (std::find(labels.begin(), labels.end(), label) == labels.end()) labels.push_back(label);
  int offset = 0;
  // The existing R construction uses the first chain/subject correlation.
  for (int label : labels) {
    std::vector<int> selected;
    for (int k = 0; k < shared.size(); ++k) if (shared[k] == label) selected.push_back(k);
    const int n = selected.size(), target = std::min(n_blocks, n);
    arma::mat distance(n, n);
    std::vector<std::vector<int>> members(n);
    for (int i = 0; i < n; ++i) {
      members[i].push_back(i);
      for (int j = 0; j < n; ++j)
        distance(i, j) = 1 - std::abs(covariance(selected[i], selected[j]) /
          std::sqrt(covariance(selected[i], selected[i]) * covariance(selected[j], selected[j])));
    }
    for (int remaining = n; remaining > target; --remaining) {
      double best = arma::datum::inf;
      int left = -1, right = -1;
      for (int i = 0; i < n; ++i) if (!members[i].empty())
        for (int j = i + 1; j < n; ++j) if (!members[j].empty() && distance(i, j) < best) {
          best = distance(i, j); left = i; right = j;
        }
      for (int j = 0; j < n; ++j) if (j != left && !members[j].empty())
        distance(left, j) = distance(j, left) = std::max(distance(left, j), distance(right, j));
      members[left].insert(members[left].end(), members[right].begin(), members[right].end());
      members[right].clear();
    }
    int block = 0;
    for (const auto& cluster : members) if (!cluster.empty()) {
      ++block;
      for (int index : cluster) out[selected[index]] = offset + block;
    }
    offset += block;
  }
  return out;
}

Rcpp::List conditional_proposal(Rcpp::List samples, std::string type, int subject,
                              int iteration, SEXP indices) {
  Rcpp::NumericVector alpha = samples["alpha"];
  Rcpp::IntegerVector dims = alpha.attr("dim");
  arma::uvec selected = selected_indices(indices, dims[0]);
  ProposalBlock block = make_block(samples, type, selected, selected, subject - 1);
  return block.subjects[0].at(iteration - 1);
}

Rcpp::NumericVector group_means(Rcpp::List sampler, Rcpp::NumericMatrix means,
                                     SEXP random = R_NilValue) {
  emc2::GroupSampler group(sampler);
  emc2::GroupState state;
  arma::mat beta(means.begin(), means.nrow(), means.ncol(), false);
  arma::mat u;
  if (!Rf_isNull(random)) u = Rcpp::as<arma::mat>(random);
  arma::cube out(group.parameters(), Rcpp::as<int>(sampler["n_subjects"]), means.ncol());
  for (int i = 0; i < means.ncol(); ++i) {
    state.tmu = beta.col(i);
    if (u.n_elem) state.u = u.col(i);
    out.slice(i) = group.subject_means(state);
  }
  return Rcpp::wrap(out);
}

} // namespace emc2
