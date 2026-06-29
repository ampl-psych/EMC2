// Rcpp.h must come before headers that pull in <R.h> (model_RDM.h chain does),
// to keep the Intel/NCI build happy (see RcppCore/Rcpp#1410).
#include <Rcpp.h>
#include "race_gng.h"
using namespace Rcpp;

// Test entry point: P(no-go accumulator finishes first in [lower, upper]) for a
// single RDM trial. v/B/A/t0/s are per-accumulator (length n_acc); nogo is the
// 1-based index of the no-go accumulator.
// [[Rcpp::export]]
double rdm_gng_withheld_R(NumericVector v, NumericVector B, NumericVector A,
                          NumericVector t0, NumericVector s, int nogo,
                          double lower, double upper) {
  const int n_acc = v.size();
  return rdm_gng_withheld_prob(v.begin(), B.begin(), A.begin(), t0.begin(),
                               s.begin(), n_acc, nogo - 1, lower, upper);
}
