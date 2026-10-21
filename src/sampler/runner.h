#ifndef EMC2_SAMPLER_RUNNER_H
#define EMC2_SAMPLER_RUNNER_H
#include <RcppArmadillo.h>
#include <memory>
#include <vector>
namespace emc2 {
class Likelihoods;

// Owned by an R external pointer; reused across fit() batches.
struct SamplerContexts {
  std::vector<std::unique_ptr<Likelihoods>> chains;
  int cores_per_chain;
#ifndef _OPENMP
  bool warned_openmp = false;
#endif
  SamplerContexts(const Rcpp::List& inputs, int n_chains, int cores);
  ~SamplerContexts();
};
Rcpp::List sample_chains(Rcpp::List args, SamplerContexts& contexts,
                        std::string stage, int iter, int cores_for_chains, bool verbose);
Rcpp::List initialize_sampler(Rcpp::List sampler, SEXP start_mu, SEXP start_var,
                             int particles, Likelihoods& likelihood);
Rcpp::List sample_group(Rcpp::List sampler, Rcpp::NumericVector alpha, int iter);
}
#endif
