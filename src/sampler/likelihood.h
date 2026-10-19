#ifndef EMC2_SAMPLER_LIKELIHOOD_H
#define EMC2_SAMPLER_LIKELIHOOD_H

#include <RcppArmadillo.h>
#include <exception>
#include <memory>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace emc2 {

// One set of likelihood components per subject. Copies share fixed model inputs;
// coefficient, trend and likelihood buffers belong to the individual chain.
class Likelihoods {
  struct Part;
  std::vector<std::vector<std::unique_ptr<Part>>> parts;
  std::vector<bool> safe;
  int cores;
public:
  Likelihoods(const Rcpp::List& inputs, int requested);
  Likelihoods(const Likelihoods& other);
  ~Likelihoods();
  bool thread_safe() const;
  void one(int subject, const arma::mat& candidates, arma::vec& out, int component = 0);

  template<class Work> void for_each(Work work) {
    std::exception_ptr error;
    const int n = parts.size();
    auto one = [&](int s) {
      if (!safe[s]) return;
      try { work(s); }
      catch (...) {
#ifdef _OPENMP
#pragma omp critical(emc_likelihood_error)
#endif
        { if (!error) error = std::current_exception(); }
      }
    };
#ifdef _OPENMP
    if (omp_in_parallel()) {
      // Already in the chain team; schedule subject tasks on that same team.
#pragma omp taskgroup
      for (int k = 0; k < std::min(cores, n); ++k) {
#pragma omp task firstprivate(k) shared(one)
        for (int s = k; s < n; s += cores) one(s);
      }
    } else {
#pragma omp parallel for num_threads(cores) schedule(static) if(cores > 1)
      for (int s = 0; s < n; ++s) one(s);
    }
#else
    for (int s = 0; s < n; ++s) one(s);
#endif
    // Task errors go back to the caller; APIs requiring R run on the R thread.
    if (error) std::rethrow_exception(error);
    for (int s = 0; s < n; ++s) if (!safe[s]) work(s);
  }

  template<class Get> void evaluate(Get candidates, std::vector<arma::vec>& out, int component = 0) {
    for_each([&](int s) { one(s, candidates(s), out[s], component); });
  }
};

Rcpp::List likelihood_precision(Rcpp::List emc, Rcpp::IntegerVector history_idx,
                                Rcpp::List ll_inputs, int n_cores);

} // namespace emc2
#endif
