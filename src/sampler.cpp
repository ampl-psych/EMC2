#include "sampler/runner.h"
#include "sampler/proposals.h"
#include "sampler/likelihood.h"
#include "openmp.h"

// =============================================================================
// Rcpp entry points — R objects and OpenMP pool lifetime are handled here
// Sampling, proposal construction and likelihood caching live in sampler/.
// =============================================================================

// [[Rcpp::export]]
Rcpp::List native_run_hyper(Rcpp::List sampler, Rcpp::NumericVector alpha, int iter) {
  return emc2::sample_group(sampler, alpha, iter);
}

// [[Rcpp::export]]
Rcpp::List native_chain_proposals(Rcpp::List emc, Rcpp::NumericVector alpha,
    Rcpp::NumericMatrix likelihood, bool do_block, std::string stage) {
  return emc2::chain_proposals(emc, alpha, likelihood, do_block, stage);
}

// [[Rcpp::export]]
Rcpp::List native_eff_proposals(Rcpp::List emc, Rcpp::List samples) {
  return emc2::efficient_proposals(emc, samples);
}

// [[Rcpp::export]]
Rcpp::IntegerVector native_sub_blocking(Rcpp::List emc, int n_blocks) {
  return emc2::parameter_blocks(emc, n_blocks);
}

// [[Rcpp::export]]
Rcpp::List native_conditionals(Rcpp::List samples, std::string type, int subject,
    int iteration, SEXP indices) {
  return emc2::conditional_proposal(samples, type, subject, iteration, indices);
}

// [[Rcpp::export]]
Rcpp::NumericVector native_group_means(Rcpp::List sampler, Rcpp::NumericMatrix means,
    SEXP random = R_NilValue) {
  return emc2::group_means(sampler, means, random);
}

// [[Rcpp::export]]
Rcpp::List native_lik_prec(Rcpp::List emc, Rcpp::IntegerVector history_idx,
    Rcpp::List ll_inputs, int n_cores) {
  emc2::OpenMPPoolRelease release;
  return emc2::likelihood_precision(emc, history_idx, ll_inputs, n_cores);
}

// [[Rcpp::export]]
SEXP native_sampler_contexts(Rcpp::List inputs, int n_chains, int n_cores) {
  return Rcpp::XPtr<emc2::SamplerContexts>(new emc2::SamplerContexts(inputs, n_chains, n_cores), true);
}

// [[Rcpp::export]]
Rcpp::List native_run_chains(Rcpp::List args, SEXP cache, std::string stage,
                           int iter, int cores_for_chains, bool verboseProgress) {
  emc2::OpenMPPoolRelease release;
  Rcpp::XPtr<emc2::SamplerContexts> contexts(cache);
  return emc2::sample_chains(args, *contexts, stage, iter, cores_for_chains, verboseProgress);
}

// [[Rcpp::export]]
Rcpp::List native_init_chains(Rcpp::List emc, SEXP cache, SEXP start_mu,
                            SEXP start_var, int particles, bool force) {
  emc2::OpenMPPoolRelease release;
  Rcpp::XPtr<emc2::SamplerContexts> contexts(cache);
  Rcpp::List out(emc.size());
  out.attr("names") = emc.attr("names");
  for (int c = 0; c < emc.size(); ++c)
    out[c] = force || !Rcpp::as<bool>(Rcpp::as<Rcpp::List>(emc[c])["init"]) ?
      emc2::initialize_sampler(emc[c], start_mu, start_var, particles, *contexts->chains[c]) : emc[c];
  return out;
}
