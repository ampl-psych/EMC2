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
  out.zeros(candidates.n_rows);
  for (auto& part : parts[subject])
    if (!component || component == part->component) part->evaluate(candidates, out);
}

namespace {
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
  const arma::uvec trusted = arma::find(lambda > tau), other = arma::find(lambda <= tau);
  const arma::mat U = vectors.cols(trusted);
  const arma::vec mt = Ph * mean, mut = Ph * mu;
  arma::mat Lw = U * arma::diagmat(lambda.elem(trusted)) * U.t();
  arma::vec lw = U * (values.elem(trusted) % (U.t() * mt) - U.t() * mut);
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
    for (int i : history_idx) history[chain].push_back(selected[i - 1]);
    for (int i = 0; i < window; ++i) for (int s = 0; s < n; ++s)
      std::copy_n(values + static_cast<size_t>(history[chain][i]) * p * n + s * p,
                  p, draws[s].colptr(chain * window + i));
  }
  arma::mat group_precision, group_mean;
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
      covariance /= N; group_mean /= chains;
      // Use draw-based curvature only when all marginal precisions have CV < .5.
      const arma::vec cv = arma::sqrt(precision_m2 / (N - 1)) / precision_mean;
      group = cv.is_finite() && cv.max() < .5 && arma::inv_sympd(group_precision, covariance) &&
        group_mean.is_finite() && group_mean.n_rows == group_precision.n_rows;
    } catch (const std::exception&) { group = false; }
  }
  Likelihoods likelihood(ll_inputs, n_cores);
  std::vector<std::optional<Quadratic>> finite(n), posterior(n);
  likelihood.for_each([&](int s) {
    try {
      finite[s] = finite_difference(likelihood, s, arma::mean(draws[s], 1), arma::stddev(draws[s], 0, 1));
      if (group && finite[s]) posterior[s] = draw_precision(draws[s], group_precision, group_mean.col(s), *finite[s], chains);
    } catch (const std::exception&) {}
  });
  Rcpp::CharacterVector names = first["par_names"];
  Rcpp::List result(n);
  for (int s = 0; s < n; ++s) if (finite[s]) {
    Rcpp::List precision = wrap_quadratic(*finite[s], names);
    if (posterior[s]) precision["post"] = wrap_quadratic(*posterior[s], names);
    result[s] = precision;
  }
  return result;
}

} // namespace emc2
