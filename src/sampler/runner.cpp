#include "runner.h"
#include "likelihood.h"
#include "group.h"
#include "particle.h"
#include <algorithm>
#include <exception>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace emc2 {
namespace {
bool present(const Rcpp::List& x, const char* key) {
  return x.containsElementNamed(key) && !Rf_isNull(x[key]);
}

// =============================================================================
// SampleStore — R arrays allocated at the boundary, filled without R calls
// =============================================================================

class SampleStore {
 public:
  Rcpp::List data;
  explicit SampleStore(SEXP x) : data(Rf_isNull(x) ? Rcpp::List() :
                                     Rcpp::clone(Rcpp::as<Rcpp::List>(x))) {
    for (int k = 0; k < 17; ++k) if (present(data, keys[k])) {
      Rcpp::NumericVector values = data[keys[k]];
      fields[k] = values.begin();
    }
  }
  int index() const { return present(data, "idx") ? Rcpp::as<int>(data["idx"]) - 1 : 0; }
  arma::mat alpha(int p, int n) const {
    Rcpp::NumericVector values = data["alpha"];
    return arma::mat(values.begin() + static_cast<size_t>(index()) * p * n, p, n);
  }
  arma::vec likelihood(int n) const {
    Rcpp::NumericVector values = data["subj_ll"];
    return arma::vec(values.begin() + static_cast<size_t>(index()) * n, n);
  }
  Rcpp::List last() const {
    // Stored names -> GroupState names. The last dimension is always iteration.
    Rcpp::List out;
    const int j = index();
    const char* names[] = {"tmu", "tvar", "a_half", "u", "s", "a_half_s", "alpha", "subj_ll",
                           "lambda", "eta", "epsilon_inv", "psi", "delta", "B", "K", "G", "delta_inv"};
    for (int k = 0; k < 17; ++k) if (k != 6 && k != 7 && present(data, keys[k])) {
      Rcpp::NumericVector values = data[keys[k]];
      Rcpp::IntegerVector dims = values.attr("dim");
      int size = 1;
      for (R_xlen_t d = 0; d + 1 < dims.size(); ++d) size *= dims[d];
      Rcpp::NumericVector slice(values.begin() + static_cast<size_t>(j) * size,
                                values.begin() + static_cast<size_t>(j + 1) * size);
      if (dims.size() > 2) slice.attr("dim") = Rcpp::IntegerVector(dims.begin(), dims.end() - 1);
      out[names[k]] = slice;
    }
    if (present(data, "last_theta_var_inv")) out["tvinv"] = data["last_theta_var_inv"];
    return out;
  }
  void group(const GroupState& state, int j) {
    put(0, state.tmu, j); put(1, state.population_var.n_elem ? state.population_var : state.tvar, j);
    put(2, state.a_half, j); put(3, state.u, j);
    put(4, state.s, j); put(5, state.a_half_s, j);
    put(8, state.lambda, j); put(9, state.eta, j); put(10, state.epsilon_inv, j);
    put(11, state.psi, j); put(12, state.delta, j); put(13, state.B, j);
    put(14, state.K, j); put(15, state.G, j); put(16, state.delta_inv, j);
  }
  void particles(const arma::mat& alpha, const arma::vec& ll, int j) {
    put(6, alpha, j); put(7, ll, j);
  }
  void finish(const GroupState& state, int j) {
    data["idx"] = j + 1;
    if (present(data, "theta_var")) {
      Rcpp::NumericMatrix inv = Rcpp::wrap(state.tvinv);
      Rcpp::NumericVector var = data["theta_var"];
      SEXP names = var.attr("dimnames");
      if (!Rf_isNull(names)) {
        Rcpp::List dn(names);
        inv.attr("dimnames") = Rcpp::List::create(dn[0], dn[1]);
      }
      data["last_theta_var_inv"] = inv;
    }
  }
 private:
  const char* keys[17] = {"theta_mu", "theta_var", "a_half", "theta_u", "theta_s", "a_half_s", "alpha", "subj_ll",
                          "lambda", "eta", "epsilon_inv", "psi", "delta", "B", "K", "G", "delta_inv"};
  double* fields[17] = {};
  template<class T> void put(int field, const T& value, int j) {
    if (fields[field])
      std::copy(value.begin(), value.end(), fields[field] + static_cast<size_t>(j) * value.n_elem);
  }
};

struct Population {
  GroupSampler main;
  GroupState state;
  explicit Population(const Rcpp::List& sampler) : main(sampler) {}
};

// =============================================================================
// Joint moves — scale subject deviations, or translate subjects and means together
// =============================================================================

class JointMoves {
  Rcpp::List settings;
  arma::vec step, loc, sum_step, sum_loc, mean, curvature, shift;
  arma::mat directions, prior_directions;
  std::vector<arma::mat> proposals;
  std::vector<arma::vec> values;
  int uses;
  bool do_location;
 public:
  JointMoves(const GroupSampler& group, const Rcpp::List& init,
             const Rcpp::List& sampler, const Rcpp::List& tune, int n)
    : settings(Rcpp::clone(init)), step(Rcpp::as<arma::vec>(init["step"])),
      loc(Rcpp::as<arma::vec>(init["loc"])),
      sum_step(Rcpp::as<arma::vec>(init["log_step_sum"])),
      sum_loc(Rcpp::as<arma::vec>(init["log_loc_sum"])), proposals(n), values(n),
      uses(Rcpp::as<int>(init["uses"])) {
    do_location = group.location_move();
    if (!do_location) return;
    const int p = loc.n_elem;
    Rcpp::List prior = sampler["prior"];
    Rcpp::NumericVector variance = prior["theta_mu_var"];
    Rcpp::NumericVector centre = prior["theta_mu_mean"];
    mean.set_size(p); shift.set_size(p);
    for (int j = 0; j < p; ++j) mean[j] = centre[j % centre.size()];
    arma::mat covariance(p, p, arma::fill::zeros);
    if (variance.hasAttribute("dim")) covariance = Rcpp::as<arma::mat>(variance);
    else for (int j = 0; j < p; ++j) covariance(j, j) = variance[j % variance.size()];
    // Use the saved directions throughout a sampling stage. Otherwise estimate
    // curvature in prior-standardised coordinates and shrink the stiff directions.
    if (present(init, "directions")) directions = Rcpp::as<arma::mat>(init["directions"]);
    else {
      const arma::mat root = arma::chol(covariance, "lower");
      arma::mat precision(p, p, arma::fill::zeros);
      if (present(tune, "lik_prec")) {
        Rcpp::List likelihood = tune["lik_prec"];
        for (SEXP x : likelihood) if (!Rf_isNull(x))
          precision += Rcpp::as<arma::mat>(Rcpp::List(x)["prec"]);
      }
      precision = root.t() * precision * root;
      precision.diag() += 1;
      arma::vec eigenvalues;
      arma::mat eigenvectors;
      arma::eig_sym(eigenvalues, eigenvectors, arma::symmatu(precision));
      directions = root * eigenvectors *
        arma::diagmat(1 / arma::sqrt(arma::clamp(eigenvalues, 1, arma::datum::inf)));
      loc.ones();
    }
    arma::mat precision;
    if (present(prior, "theta_mu_invar")) {
      Rcpp::NumericVector inverse = prior["theta_mu_invar"];
      if (inverse.hasAttribute("dim")) precision = Rcpp::as<arma::mat>(inverse);
      else {
        precision.zeros(p, p);
        for (int j = 0; j < p; ++j) precision(j, j) = inverse[j % inverse.size()];
      }
    } else precision = group_inverse(covariance);
    // Cache the linear/quadratic terms of the normal mean-prior ratio.
    prior_directions = precision * directions;
    curvature = arma::sum(directions % prior_directions, 0).t();
  }
  void draw(Random& rng, Population& pop, arma::mat& alpha, arma::vec& ll, Likelihoods& likelihood, bool frozen) {
    const int p = pop.state.tvar.n_rows, n = alpha.n_cols;
    arma::vec accepted(p, arma::fill::zeros), accepted_loc(p, arma::fill::zeros);
    for (int k = 0; k < (do_location ? 2 * p : p); ++k) {
      const bool scale = k < p;
      const int j = k % p;
      const double delta = rng.rnorm(0, scale ? step[j] : loc[j]);
      const double factor = scale ? std::exp(delta) : 1;
      GroupState proposed;
      double lp = 0;
      if (scale) {
        proposed = pop.state;
        lp = pop.main.affine(proposed, j, delta, 0);
      } else {
        // Translating alpha and mu equally leaves the conditional group density
        // unchanged. Only the mean prior and subject likelihoods enter the ratio.
        shift = delta * directions.col(j);
        lp = -delta * arma::dot(prior_directions.col(j), pop.state.tmu - mean) -
          .5 * delta * delta * curvature[j];
      }
      for (int s = 0; s < n; ++s) {
        proposals[s] = alpha.col(s).t();
        if (scale) proposals[s](0, j) = pop.state.subj_mu(j, s) +
          factor * (alpha(j, s) - pop.state.subj_mu(j, s));
        else proposals[s].row(0) += shift.t();
      }
      // Exact likelihood check for all subjects; curvature only sets the proposal.
      likelihood.evaluate([&](int s) -> const arma::mat& { return proposals[s]; }, values);
      double lr = lp;
      for (int s = 0; s < n; ++s) lr += values[s][0] - ll[s];
      if (std::isfinite(lr) && std::log(rng.runif(0, 1)) < lr) {
        for (int s = 0; s < n; ++s) ll[s] = values[s][0];
        if (scale) {
          for (int s = 0; s < n; ++s) alpha(j, s) = proposals[s](0, j);
          accepted[j] = 1;
          pop.state = std::move(proposed);
        } else {
          alpha.each_col() += shift;
          pop.state.tmu += shift;
          pop.state.subj_mu.each_col() += shift;
          accepted_loc[j] = 1;
        }
      }
    }
    // Warmup tuning; R freezes the geometric average of the adaptation tail.
    if (!frozen) {
      step %= arma::exp(.1 * (accepted - .3));
      if (do_location) loc %= arma::exp(.1 * (accepted_loc - .3));
      if (++uses > 20) { sum_step += arma::log(step); sum_loc += arma::log(loc); }
    }
  }
  Rcpp::List finish() {
    auto save = [&](const char* key, const arma::vec& x) {
      Rcpp::NumericVector target = settings[key];
      std::copy(x.begin(), x.end(), target.begin());
    };
    save("step", step); save("loc", loc); save("log_step_sum", sum_step); save("log_loc_sum", sum_loc);
    settings["uses"] = uses;
    if (do_location) settings["directions"] = Rcpp::wrap(directions);
    return settings;
  }
};

void attach_samples(Rcpp::List& sampler, SampleStore& store, Population& pop, int j) {
  store.finish(pop.state, j); sampler["samples"] = store.data;
}
} // namespace

// =============================================================================
// Initial particles — draw from the group distribution and weight by likelihood
// =============================================================================

Rcpp::List initialize_sampler(Rcpp::List sampler, SEXP start_mu, SEXP start_var,
                      int particles, Likelihoods& likelihood) {
  sampler = Rcpp::List(Rf_shallow_duplicate(sampler));
  Population pop(sampler);
  SampleStore store(static_cast<SEXP>(sampler["samples"]));
  Random rng;
  pop.state = pop.main.initialize(rng, start_mu, start_var);

  const int n = Rcpp::as<int>(sampler["n_subjects"]), p = Rcpp::as<int>(sampler["n_pars"]);
  std::vector<arma::mat> candidates(n);
  std::vector<arma::vec> values(n);
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues, eigenvectors, pop.state.tvar))
    Rcpp::stop("invalid initial covariance");
  const arma::mat root = eigenvectors * arma::diagmat(arma::sqrt(arma::clamp(eigenvalues, 0, arma::datum::inf)));
  for (int s = 0; s < n; ++s) {
    arma::mat z(particles, p);
    for (int i = 0; i < particles; ++i)
      for (int j = 0; j < p; ++j) z(i, j) = rng.rnorm(0, 1);
    candidates[s] = z * root.t();
    candidates[s].each_row() += pop.state.subj_mu.col(s).t();
  }
  likelihood.evaluate([&](int s) -> const arma::mat& { return candidates[s]; }, values);
  arma::mat alpha(p, n); arma::vec ll(n);
  for (int s = 0; s < n; ++s) {
    arma::vec weights = arma::exp(values[s] - values[s].max());
    const double total = arma::sum(weights);
    if (!std::isfinite(total) || total <= 0) Rcpp::stop("initial particle weights are not finite");
    double threshold = rng.runif(0, total);
    arma::uword selected = 0;
    while (selected + 1 < weights.n_elem && threshold >= weights[selected])
      threshold -= weights[selected++];
    while (weights[selected] == 0) --selected;
    alpha.col(s) = candidates[s].row(selected).t(); ll[s] = values[s][selected];
  }
  store.group(pop.state, 0); store.particles(alpha, ll, 0);
  attach_samples(sampler, store, pop, 0);
  store.data.attr("rng") = rng.state();
  sampler["init"] = true;
  return sampler;
}

namespace {

// =============================================================================
// Chain — one batch of group, joint and particle updates
// Construction/finish handle R storage; advance() performs the sampling updates.
// =============================================================================

class Chain {
  Rcpp::List sampler;
  SampleStore store;
  Random rng;
  Population pop;
  Likelihoods& likelihood;
  ParticleSampler particle;
  std::unique_ptr<JointMoves> joint;
  int start, retries;
  bool frozen;
  arma::mat alpha, component_ll;
  arma::vec ll;
  std::vector<arma::vec> values;
 public:
  int iteration = 0;
  Chain(const Rcpp::List& args, Likelihoods& context)
    : sampler(Rf_shallow_duplicate(args["sampler"])), store(static_cast<SEXP>(sampler["samples"])),
      rng(static_cast<SEXP>(store.data.attr("rng"))), pop(sampler), likelihood(context),
      particle(rng, sampler, args["tune"], args["pm_settings"],
        Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(args["tune"])["kernel"])) {
    pop.state = pop.main.single() ? pop.main.initialize(rng) : pop.main.read_state(store.last());
    Rcpp::List tune = args["tune"], on_singular = args["on_singular"];
    if (Rcpp::as<bool>(tune["do_scale"]))
      joint.reset(new JointMoves(pop.main, args["scale_settings"], sampler, tune,
                                pop.state.subj_mu.n_cols));
    start = store.index();
    retries = Rcpp::as<int>(on_singular["max_retries"]);
    frozen = Rcpp::as<bool>(tune["frozen"]);
    const int n = pop.state.subj_mu.n_cols, p = pop.state.tvar.n_rows;
    alpha = store.alpha(p, n);
    ll = store.likelihood(n);
    values.resize(n);
    component_ll.zeros(n, particle.block_count());
  }
  void advance(int end) {
    while (iteration < end) {
      ++iteration;
      for (int attempt = 0; ; ++attempt) {
        try { pop.main.draw(rng, alpha, pop.state); break; }
        catch (const std::runtime_error& e) {
          if (attempt >= retries || std::string(e.what()).find("singular") == std::string::npos) throw;
        }
      }
      if (joint) joint->draw(rng, pop, alpha, ll, likelihood, frozen);

      particle.begin_iteration(pop.state.tvar);
      component_ll.zeros();
      // Blocks are sequential: prepare() sees alpha after all earlier updates.
      for (int b = 0; b < particle.block_count(); ++b) {
        for (arma::uword s = 0; s < alpha.n_cols; ++s)
          particle.prepare(s, b, alpha.col(s), pop.state.subj_mu.col(s), pop.state.tvinv);
        const int component = particle.block_component(b);
        likelihood.evaluate([&](int s) -> const arma::mat& { return particle.candidate(s); }, values, component);
        for (arma::uword s = 0; s < alpha.n_cols; ++s) {
          const int selected = particle.select(s, b, values[s]);
          alpha.col(s) = particle.candidate(s).row(selected).t();
          // Shared likelihoods are overwritten by their last block, not added twice.
          component_ll(s, component - 1) = values[s][selected];
        }
      }
      ll = arma::sum(component_ll, 1);
      const int j = start + iteration;
      store.group(pop.state, j); store.particles(alpha, ll, j);
    }
  }
  Rcpp::List finish() {
    store.data.attr("pm_settings") = particle.settings();
    if (joint) store.data.attr("scale_move") = joint->finish();
    store.data.attr("rng") = rng.state();
    attach_samples(sampler, store, pop, start + iteration);
    return sampler;
  }
};
} // namespace

// =============================================================================
// Chain scheduling — shared fixed likelihood inputs, separate mutable workspaces
// =============================================================================

SamplerContexts::SamplerContexts(const Rcpp::List& inputs, int n_chains, int cores)
  : cores_per_chain(std::max(1, cores)) {
  chains.emplace_back(new Likelihoods(inputs, cores_per_chain));
  for (int c = 1; c < n_chains; ++c) chains.emplace_back(new Likelihoods(*chains[0]));
}
SamplerContexts::~SamplerContexts() = default;

Rcpp::List sample_chains(Rcpp::List args, SamplerContexts& contexts,
                        std::string stage, int iter, int cores_for_chains, bool verbose) {
  const int n = args.size(), concurrent = std::min(n, std::max(1, cores_for_chains));
#ifdef _OPENMP
  const int threads = concurrent * contexts.cores_per_chain;
  const bool parallel = concurrent > 1 && contexts.chains[0]->thread_safe();
#else
  if (!contexts.warned_openmp && (concurrent > 1 || contexts.cores_per_chain > 1)) {
    Rcpp::warning("OpenMP is not available in this EMC2 build. Sampling will run serially, ignoring the requested parallel core settings.");
    contexts.warned_openmp = true;
  }
#endif
  // Allocate R storage and restore RNGs before entering the parallel region.
  std::vector<std::unique_ptr<Chain>> chains;
  for (int c = 0; c < n; ++c) chains.emplace_back(new Chain(args[c], *contexts.chains[c]));
  std::vector<std::exception_ptr> errors(n);

  // Return to the R thread every ten iterations for interrupts and progress.
  for (int end = 0; end < iter;) {
    end = std::min(iter, end + 10);
    Rcpp::checkUserInterrupt();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) if(parallel)
#pragma omp master
#endif
    {
      // Chain and subject tasks share one OpenMP team, avoiding nested teams.
      for (int first = 0; first < n; first += concurrent) {
#ifdef _OPENMP
#pragma omp taskgroup
#endif
        for (int c = first; c < std::min(n, first + concurrent); ++c) {
#ifdef _OPENMP
#pragma omp task firstprivate(c, end) if(parallel)
#endif
          try { chains[c]->advance(end); }
          catch (...) { errors[c] = std::current_exception(); }
        }
      }
    }
    for (int c = 0; c < n; ++c) if (errors[c]) {
      try { std::rethrow_exception(errors[c]); }
      catch (const std::exception& error) {
        Rcpp::stop("%s\n  [context: chain %d, '%s' stage, iteration %d of %d]",
                   error.what(), c + 1, stage, chains[c]->iteration, iter);
      }
    }
    if (verbose) Rcpp::Rcout << "\r" << stage << ": " << end << "/" << iter << std::flush;
  }
  if (verbose) Rcpp::Rcout << "\n";
  Rcpp::List out(n);
  out.attr("names") = args.attr("names");
  for (int c = 0; c < n; ++c) out[c] = chains[c]->finish();
  return out;
}

// Group-only sampling conditional on supplied subject coefficients (run_hyper).
Rcpp::List sample_group(Rcpp::List sampler, Rcpp::NumericVector alpha, int iter) {
  sampler = Rcpp::List(Rf_shallow_duplicate(sampler));
  GroupSampler group(sampler);
  Random rng;
  GroupState state = group.initialize(rng);
  SampleStore store(static_cast<SEXP>(sampler["samples"]));
  Rcpp::IntegerVector dims = alpha.attr("dim");
  const int p = dims[0], n = dims[1], draws = dims[2];
  store.group(state, 0);
  for (int i = 1; i < iter; ++i) {
    Rcpp::checkUserInterrupt();
    arma::mat current(alpha.begin() + static_cast<size_t>(draws == 1 ? 0 : i) * p * n, p, n, false);
    group.draw(rng, current, state);
    store.group(state, i);
  }
  store.finish(state, iter - 1);
  sampler["samples"] = store.data;
  return sampler;
}

} // namespace emc2
