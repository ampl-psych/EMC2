#ifndef EMC2_SAMPLER_PROPOSALS_H
#define EMC2_SAMPLER_PROPOSALS_H

#include <RcppArmadillo.h>

namespace emc2 {

Rcpp::List chain_proposals(Rcpp::List emc, Rcpp::NumericVector alpha,
                                 Rcpp::NumericMatrix likelihood, bool do_block,
                                 std::string stage);
Rcpp::List efficient_proposals(Rcpp::List emc, Rcpp::List samples);
Rcpp::IntegerVector parameter_blocks(Rcpp::List emc, int n_blocks);
Rcpp::List conditional_proposal(Rcpp::List samples, std::string type, int subject,
                              int iteration, SEXP indices);
Rcpp::NumericVector group_means(Rcpp::List sampler, Rcpp::NumericMatrix means,
                                     SEXP random);

} // namespace emc2
#endif
