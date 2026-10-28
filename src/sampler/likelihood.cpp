#include "likelihood.h"
#include "../LikelihoodContext.h"
#include "group.h"
#include <algorithm>
#include <cmath>
#include <optional>
#include <set>

namespace emc2 {
namespace {
SEXP field(const Rcpp::List& x, const char* name) {
  return x.containsElementNamed(name) ? static_cast<SEXP>(x[name]) : R_NilValue;
}
Rcpp::List optional_list(SEXP x) {
  return Rf_isNull(x) ? Rcpp::List() : Rcpp::List(x);
}
}

// =============================================================================
// Likelihood components — cached model setup and parameter-column mapping
// =============================================================================

struct Likelihoods::Part {
  int component;
  arma::uvec columns;
  bool direct;
  LikelihoodContext context;
  arma::mat coefficients;
  arma::vec values;
  Part(const Rcpp::List& spec, const Rcpp::List& model)
    : component(Rcpp::as<int>(spec["component"])),
      columns(Rcpp::as<arma::uvec>(spec["columns"]) - 1),
      direct(true),
      context(Rcpp::as<Rcpp::CharacterVector>(spec["par_names"]),
        Rcpp::as<Rcpp::DataFrame>(spec["data"]),
        Rcpp::as<Rcpp::NumericVector>(spec["constants"]),
        Rcpp::as<Rcpp::List>(spec["designs"]), model["c_name"],
        optional_list(field(model, "bound")), optional_list(field(model, "transform")),
        optional_list(field(model, "pre_transform")), std::log(1e-10),
        field(model, "trend")) {
    for (arma::uword j = 0; j < columns.n_elem; ++j) direct = direct && columns[j] == j;
  }
  Part(const Part& other)
    : component(other.component), columns(other.columns), direct(other.direct), context(other.context) {}
  void evaluate(const arma::mat& candidates, arma::vec& out) {
    // Most models use all columns in order; avoid a copy in that case.
    const double* input;
    if (direct && columns.n_elem == candidates.n_cols) input = candidates.memptr();
    else { coefficients = candidates.cols(columns); input = coefficients.memptr(); }
    values.set_size(candidates.n_rows);
    context.evaluate(input, candidates.n_rows, values.memptr());
    out += values;
  }
};

Likelihoods::Likelihoods(const Rcpp::List& inputs, int requested)
  : parts(inputs.size()), safe(inputs.size(), true), cores(std::max(1, requested)) {
  for (int s = 0; s < inputs.size(); ++s) {
    Rcpp::List specifications = inputs[s];
    for (SEXP specification : specifications) {
      Rcpp::List spec(specification);
      parts[s].emplace_back(new Part(spec, Rcpp::as<Rcpp::List>(spec["model"])));
      safe[s] = safe[s] && parts[s].back()->context.thread_safe();
    }
  }
}
Likelihoods::Likelihoods(const Likelihoods& other)
  : parts(other.parts.size()), safe(other.safe), cores(other.cores) {
  for (size_t s = 0; s < parts.size(); ++s)
    for (const auto& part : other.parts[s]) parts[s].emplace_back(new Part(*part));
}
Likelihoods::~Likelihoods() = default;
bool Likelihoods::thread_safe() const {
  return std::all_of(safe.begin(), safe.end(), [](bool value) { return value; });
}
void Likelihoods::one(int subject, const arma::mat& candidates, arma::vec& out, int component) {
  // component == 0 evaluates all tasks. A subject with no data contributes zero.
  out.zeros(candidates.n_rows);
  for (auto& part : parts[subject])
    if (!component || component == part->component) part->evaluate(candidates, out);
}

// Trial scores retain the tangent to a curved ridge even away from its maximum.
arma::mat Likelihoods::score_information(int subject, const arma::vec& centre,
                                        const arma::vec& scale) {
  const int p = centre.n_elem;
  arma::mat candidates = arma::repmat(centre.t(), 2 * p, 1);
  for (int j = 0; j < p; ++j) {
    candidates(j, j) += 1e-4 * scale[j];
    candidates(p + j, j) -= 1e-4 * scale[j];
  }
  arma::mat information(p, p, arma::fill::zeros);
  for (auto& part : parts[subject]) {
    part->coefficients = candidates.cols(part->columns);
    arma::mat values(part->context.trial_count(), 2 * p);
    part->context.evaluate(part->coefficients.memptr(), 2 * p, values.memptr(), true);
    const arma::mat score = (values.cols(0, p - 1) - values.cols(p, 2 * p - 1)) / 2e-4;
    information += score.t() * score;
  }
  return information;
}

namespace {

// =============================================================================
// Likelihood curvature — log L(x) ~= constant + linear' x - x' precision x / 2
// Used to construct proposals, never to replace the likelihood in the weights.
// =============================================================================

struct Quadratic {
  arma::mat precision;
  arma::vec linear;
};

std::optional<Quadratic> finite_difference(Likelihoods& likelihood, int subject,
                                          const arma::vec& centre, arma::vec h,
                                          double target = 1, int max_rounds = 8) {
  const int p = centre.n_elem;
  arma::vec f;
  likelihood.one(subject, centre.t(), f);
  const double f0 = f[0];
  if (!std::isfinite(f0)) return std::nullopt;
  for (double& x : h) if (!std::isfinite(x) || x <= 0) x = .1;
  arma::vec fp(p), fm(p), lo(p), hi(p);
  fp.fill(arma::datum::nan); fm.fill(arma::datum::nan);
  lo.fill(arma::datum::nan); hi.fill(arma::datum::nan);
  arma::uvec todo(p, arma::fill::ones);
  // Choose finite-difference steps giving about one log-unit of curvature.
  // These are new likelihood evaluations, centred on the retained draw mean.
  for (int round = 0; round < max_rounds; ++round) {
    arma::uvec k = arma::find(todo);
    const int n = k.n_elem;
    arma::mat candidates = arma::repmat(centre.t(), 2 * n, 1);
    for (int i = 0; i < n; ++i) {
      candidates(i, k[i]) += h[k[i]];
      candidates(i + n, k[i]) -= h[k[i]];
    }
    likelihood.one(subject, candidates, f);
    for (int i = 0; i < n; ++i) {
      const int j = k[i];
      fp[j] = f[i]; fm[j] = f[i + n];
      const double drop = f0 - (fp[j] + fm[j]) / 2;
      if (std::isfinite(drop) && drop > target / 3 && drop < target * 3) todo[j] = 0;
      else if (round + 1 < max_rounds) {
        const bool small = std::isfinite(drop) && drop <= target / 3;
        if (small) lo[j] = std::isfinite(lo[j]) ? std::max(lo[j], h[j]) : h[j];
        else hi[j] = std::isfinite(hi[j]) ? std::min(hi[j], h[j]) : h[j];
        const double multiplier = !std::isfinite(drop) ? .25 : drop <= 0 ? 4 :
          std::min(10., std::max(.1, std::sqrt(target / drop)));
        h[j] = std::isfinite(lo[j]) && std::isfinite(hi[j]) ? std::sqrt(lo[j] * hi[j]) : h[j] * multiplier;
      }
    }
    if (arma::accu(todo) == 0) break;
  }
  arma::mat H = arma::diagmat((2 * f0 - fp - fm) / arma::square(h));
  // Four corners for each mixed derivative, evaluated in one batch.
  const int pairs = p * (p - 1) / 2;
  if (pairs) {
    arma::mat candidates = arma::repmat(centre.t(), 4 * pairs, 1);
    int i = 0;
    for (int j = 0; j < p; ++j) for (int k = j + 1; k < p; ++k, ++i) {
      candidates(i, j) += h[j]; candidates(i, k) += h[k];
      candidates(i + pairs, j) += h[j]; candidates(i + pairs, k) -= h[k];
      candidates(i + 2 * pairs, j) -= h[j]; candidates(i + 2 * pairs, k) += h[k];
      candidates(i + 3 * pairs, j) -= h[j]; candidates(i + 3 * pairs, k) -= h[k];
    }
    likelihood.one(subject, candidates, f);
    i = 0;
    for (int j = 0; j < p; ++j) for (int k = j + 1; k < p; ++k, ++i)
      H(j, k) = H(k, j) = -(f[i] - f[i + pairs] - f[i + 2 * pairs] +
        f[i + 3 * pairs]) / (4 * h[j] * h[k]);
  }
  arma::vec gradient = (fp - fm) / (2 * h);
  for (double& x : H) if (!std::isfinite(x)) x = 0;
  for (int j = 0; j < p; ++j) {
    if (todo[j] || !std::isfinite(fp[j]) || !std::isfinite(fm[j])) {
      H.row(j).zeros(); H.col(j).zeros(); gradient[j] = 0;
    } else if (!std::isfinite(gradient[j])) gradient[j] = 0;
  }
  arma::vec values; arma::mat vectors;
  if (!arma::eig_sym(values, vectors, H)) return std::nullopt;
  // Flat/negative curvature contributes no precision; the group prior supplies it.
  values = arma::clamp(values, 0, arma::datum::inf);
  H = vectors * arma::diagmat(values) * vectors.t();
  H = (H + H.t()) / 2;
  return Quadratic{H, H * centre + gradient};
}

// coda::effectiveSize: Yule-Walker AR spectrum, with AIC-selected order.
double ar_effective_size(const arma::vec& x) {
  const int n = x.n_elem;
  if (n < 2 || !x.is_finite()) return 1;
  arma::vec y = x - arma::mean(x);
  arma::vec time = arma::regspace<arma::vec>(1, n) - (n + 1.) / 2;
  arma::vec residual = y - time * (arma::dot(time, y) / arma::dot(time, time));
  if (arma::stddev(residual) <= std::sqrt(std::numeric_limits<double>::epsilon())) return 0;
  const int max_order = std::min(n - 1, static_cast<int>(std::floor(10 * std::log10(n))));
  arma::vec covariance(max_order + 1);
  for (int lag = 0; lag <= max_order; ++lag)
    covariance[lag] = arma::dot(y.head(n - lag), y.tail(n - lag)) / n;
  if (covariance[0] <= 0) return 1;
  arma::vec coefficients(max_order, arma::fill::zeros);
  double innovation = covariance[0], best_variance = innovation, best_sum = 0;
  double best_aic = n * std::log(innovation);
  int best_order = 0;
  for (int m = 1; m <= max_order; ++m) {
    double numerator = covariance[m];
    for (int j = 1; j < m; ++j) numerator -= coefficients[j - 1] * covariance[m - j];
    const double reflection = numerator / innovation;
    arma::vec previous = coefficients;
    for (int j = 1; j < m; ++j)
      coefficients[j - 1] = previous[j - 1] - reflection * previous[m - j - 1];
    coefficients[m - 1] = reflection;
    innovation *= 1 - reflection * reflection;
    if (!(innovation > 0)) return 1;
    const double aic = n * std::log(innovation) + 2 * m;
    if (aic < best_aic) {
      best_aic = aic; best_variance = innovation; best_order = m;
      best_sum = arma::accu(coefficients.head(m));
    }
  }
  const double spectrum = best_variance * n / (n - best_order - 1) / std::pow(1 - best_sum, 2);
  return spectrum == 0 ? 0 : n * arma::var(x) / spectrum;
}

// ---------------------------------------------------------------------------
// Draw-based refinement — subtract the group precision from posterior precision
// Work in prior-standardised coordinates to compare directions on the same scale.
// ---------------------------------------------------------------------------
std::optional<Quadratic> draw_precision(const arma::mat& draws, const arma::mat& P,
                                        const arma::vec& mu, const Quadratic& fd,
                                        int n_chains, double tau = .5) {
  const int p = draws.n_rows, N = draws.n_cols;
  std::set<double> unique;
  for (int i = 0; i < N; ++i) unique.insert(draws(0, i));
  if (unique.size() < static_cast<size_t>(std::max(50, 5 * p))) return std::nullopt;
  arma::vec effective(p, arma::fill::zeros);
  const int block = N / n_chains;
  for (int j = 0; j < p; ++j) for (int chain = 0; chain < n_chains; ++chain)
    effective[j] += ar_effective_size(draws.row(j).cols(chain * block, (chain + 1) * block - 1).t());
  tau = std::max(tau, std::pow(1 + std::sqrt(p / std::max(arma::median(effective), static_cast<double>(p))), 2) - 1);
  const arma::vec mean = arma::mean(draws, 1);
  arma::mat inverse_covariance;
  if (!arma::inv_sympd(inverse_covariance, arma::cov(draws.t()))) return std::nullopt;
  arma::vec values; arma::mat vectors;
  if (!arma::eig_sym(values, vectors, P) || arma::any(values <= 0)) return std::nullopt;
  const arma::mat Ph = vectors * arma::diagmat(arma::sqrt(values)) * vectors.t();
  const arma::mat Pih = vectors * arma::diagmat(1 / arma::sqrt(values)) * vectors.t();
  arma::mat M = Pih * inverse_covariance * Pih;
  if (!arma::eig_sym(values, vectors, (M + M.t()) / 2)) return std::nullopt;
  const arma::vec lambda = values - 1;
  // Use empirical curvature only where it clearly exceeds the group precision.
  const arma::uvec trusted = arma::find(lambda > tau), other = arma::find(lambda <= tau);
  const arma::mat U = vectors.cols(trusted);
  const arma::vec mt = Ph * mean, mut = Ph * mu;
  arma::mat Lw = U * arma::diagmat(lambda.elem(trusted)) * U.t();
  arma::vec lw = U * (values.elem(trusted) % (U.t() * mt) - U.t() * mut);
  // In the remaining directions use capped finite-difference curvature.
  if (other.n_elem) {
    const arma::mat C = vectors.cols(other), Fw = Pih * fd.precision * Pih;
    arma::mat Fc = C.t() * Fw * C;
    if (!arma::eig_sym(values, vectors, (Fc + Fc.t()) / 2)) return std::nullopt;
    const arma::mat Vc = C * vectors;
    arma::vec phi = arma::clamp(values, 0, arma::datum::inf);
    const arma::uvec capped = arma::find(phi > tau);
    phi = arma::clamp(phi, 0, tau);
    Lw += Vc * arma::diagmat(phi) * Vc.t();
    arma::vec lfd = Vc.t() * (Pih * fd.linear - Fw * (U * (U.t() * mt)));
    const arma::vec empirical = (phi + 1) % (Vc.t() * mt) - Vc.t() * mut;
    lfd.elem(capped) = empirical.elem(capped);
    lw += Vc * lfd;
  }
  arma::mat L = Ph * Lw * Ph;
  L = (L + L.t()) / 2;
  arma::vec linear = Ph * lw;
  if (!L.is_finite() || !linear.is_finite()) return std::nullopt;
  return Quadratic{L, linear};
}

Rcpp::List wrap_quadratic(const Quadratic& x, const Rcpp::CharacterVector& names) {
  Rcpp::NumericMatrix precision = Rcpp::wrap(x.precision);
  precision.attr("dimnames") = Rcpp::List::create(names, names);
  Rcpp::NumericVector linear(x.linear.begin(), x.linear.end());
  linear.attr("names") = names;
  return Rcpp::List::create(Rcpp::Named("prec") = precision, Rcpp::Named("lin") = linear);
}
} // namespace

// ---------------------------------------------------------------------------
// Pool the same history window for participant and group draws
// ---------------------------------------------------------------------------
Rcpp::List likelihood_precision(Rcpp::List emc, Rcpp::IntegerVector history_idx,
                           Rcpp::List ll_inputs, int n_cores) {
  const Rcpp::List first = emc[0];
  const int p = Rcpp::as<int>(first["n_pars"]), n = Rcpp::as<int>(first["n_subjects"]);
  const int window = history_idx.size(), chains = emc.size(), N = window * chains;
  std::vector<arma::mat> draws(n, arma::mat(p, N));
  std::vector<std::vector<int>> history(chains);
  for (int chain = 0; chain < chains; ++chain) {
    Rcpp::List sampler = emc[chain], samples = sampler["samples"];
    Rcpp::NumericVector alpha = samples["alpha"];
    const double* values = alpha.begin();
    Rcpp::CharacterVector stages = samples["stage"];
    std::vector<int> selected;
    for (int i = 0; i < stages.size(); ++i)
      if (stages[i] == "preburn" || stages[i] == "burn" || stages[i] == "adapt" || stages[i] == "sample") selected.push_back(i);
    // history_idx indexes the filtered history, not the full sample arrays.
    for (int i : history_idx) history[chain].push_back(selected[i - 1]);
    for (int i = 0; i < window; ++i) for (int s = 0; s < n; ++s)
      std::copy_n(values + static_cast<size_t>(history[chain][i]) * p * n + s * p,
                  p, draws[s].colptr(chain * window + i));
  }
  arma::mat group_precision, group_mean, group_covariance;
  bool group = Rcpp::as<std::string>(first["type"]) == "standard" && N >= 100;
  if (group) {
    try {
      arma::mat covariance(p, p, arma::fill::zeros);
      arma::vec precision_mean(p, arma::fill::zeros), precision_m2(p, arma::fill::zeros);
      int count = 0;
      group_mean.zeros(p, n);
      GroupSampler population(first);
      for (int chain = 0; chain < chains; ++chain) {
        Rcpp::List sampler = emc[chain], samples = sampler["samples"];
        GroupState state;
        Rcpp::NumericMatrix means = samples["theta_mu"];
        const double* mean_values = means.begin();
        state.tmu.zeros(means.nrow());
        const SEXP random_sexp = field(samples, "theta_u");
        Rcpp::NumericMatrix random = Rf_isNull(random_sexp) ? Rcpp::NumericMatrix(0, 0) : Rcpp::NumericMatrix(random_sexp);
        state.u.zeros(random.nrow());
        Rcpp::NumericVector variance = samples["theta_var"];
        const double* variance_values = variance.begin();
        for (int i : history[chain]) {
          const arma::mat cov(const_cast<double*>(variance_values + static_cast<size_t>(i) * p * p), p, p, false, true);
          covariance += cov;
          const arma::vec precision = 1 / cov.diag(), delta = precision - precision_mean;
          precision_mean += delta / ++count;
          precision_m2 += delta % (precision - precision_mean);
          state.tmu += arma::vec(const_cast<double*>(mean_values + static_cast<size_t>(i) * means.nrow()), means.nrow(), false, true);
          if (random.nrow()) state.u += arma::vec(random.begin() + static_cast<size_t>(i) * random.nrow(), random.nrow(), false, true);
        }
        state.tmu /= window; state.u /= window;
        group_mean += population.subject_means(state);
      }
      // Group precision is the inverse of the average covariance.
      covariance /= N; group_mean /= chains;
      group_covariance = covariance;
      // Use draw-based curvature only when all marginal precisions have CV < .5.
      const arma::vec cv = arma::sqrt(precision_m2 / (N - 1)) / precision_mean;
      group = cv.is_finite() && cv.max() < .5 && arma::inv_sympd(group_precision, covariance) &&
        group_mean.is_finite() && group_mean.n_rows == group_precision.n_rows;
    } catch (const std::exception&) { group = false; }
  }
  // Finite differences remain available when group history is not stable enough.
  Likelihoods likelihood(ll_inputs, n_cores);
  // Learn affine velocities b + G z from this chain's history. They stay fixed
  // during sampling, along with the other proposal geometry.
  const arma::ivec groups = Rf_isNull(field(first, "par_group")) ? arma::ivec() :
    Rcpp::as<arma::ivec>(first["par_group"]);
  arma::mat centres(p, n);
  for (int s = 0; s < n; ++s) centres.col(s) = arma::mean(draws[s], 1);
  const arma::vec centre = arma::mean(centres, 1), scale = arma::stddev(centres, 0, 1);
  const bool affine = Rcpp::as<std::string>(first["type"]) == "standard" &&
    Rf_isNull(field(first, "group_designs")) && groups.n_elem &&
    arma::all(groups == groups[0]) && scale.is_finite() && arma::all(scale > 0);
  std::vector<arma::mat> scores(n);
  std::vector<std::optional<Quadratic>> finite(n), posterior(n);
  likelihood.for_each([&](int s) {
    try {
      finite[s] = finite_difference(likelihood, s, centres.col(s), arma::stddev(draws[s], 0, 1));
      if (group && finite[s]) posterior[s] = draw_precision(draws[s], group_precision, group_mean.col(s), *finite[s], chains);
      if (affine) scores[s] = likelihood.score_information(s, centres.col(s), scale);
    } catch (const std::exception&) {}
  });
  Rcpp::CharacterVector names = first["par_names"];
  Rcpp::List result(n);
  for (int s = 0; s < n; ++s) if (finite[s]) {
    Rcpp::List precision = wrap_quadratic(*finite[s], names);
    if (posterior[s]) precision["post"] = wrap_quadratic(*posterior[s], names);
    result[s] = precision;
  }
  if (affine) {
    arma::mat information(p * (p + 1), p * (p + 1), arma::fill::zeros);
    arma::mat positions(p + 1, p + 1, arma::fill::zeros);
    for (int s = 0; s < n; ++s) if (!scores[s].is_empty()) {
      arma::vec position(p + 1); position[0] = 1;
      position.tail(p) = (centres.col(s) - centre) / scale;
      information += arma::kron(position * position.t(), scores[s]);
      positions += position * position.t();
    }
    auto generator = [&](const arma::vec& field) {
      const arma::mat velocity = arma::reshape(field, p, p + 1);
      arma::mat out(p + 1, p + 1, arma::fill::zeros);
      out.submat(0, 0, p - 1, p - 1) = arma::diagmat(scale) *
        velocity.cols(1, p) * arma::diagmat(1 / scale);
      out.submat(0, p, p - 1, p) = scale % velocity.col(0) -
        out.submat(0, 0, p - 1, p - 1) * centre;
      return arma::mat(out / arma::abs(out).max());
    };
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    if (information.is_finite() && arma::eig_sym(eigenvalues, eigenvectors, information)) {
      Rcpp::List fields = Rcpp::List::create(generator(eigenvectors.col(0)));
      // Retain fields whose likelihood curvature is weaker than the group
      // precision. Whiten the field norm so nearly constant predictors do not
      // acquire an artificially small information eigenvalue.
      arma::mat precision, root;
      if (!group_covariance.is_empty() && arma::inv_sympd(precision, group_covariance) &&
          arma::chol(root, arma::kron(positions, arma::diagmat(scale) * precision * arma::diagmat(scale)))) {
        const arma::mat inverse_root = arma::inv(arma::trimatu(root));
        const arma::mat whitened = inverse_root.t() * information * inverse_root;
        if (arma::eig_sym(eigenvalues, eigenvectors, arma::symmatu(whitened))) {
          const arma::uvec weak = arma::find(eigenvalues < 1);
          if (weak.n_elem > 1) {
            fields = Rcpp::List(weak.n_elem);
            for (arma::uword k = 0; k < weak.n_elem; ++k)
              fields[k] = Rcpp::wrap(generator(inverse_root * eigenvectors.col(weak[k])));
          }
        }
      }
      result.attr("affine_generators") = fields;
    }
  }
  return result;
}

} // namespace emc2
