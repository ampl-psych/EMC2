#pragma once
#include <Rcpp.h>

// Signature for a custom trend kernel mapping function provided by the user.
// The function receives:
//  - trend_pars: per-trial parameter matrix for the kernel (base parameters are not included)
//  - input:      per-trial input matrix (covariates and/or par_input columns)
// and returns a NumericVector of length nrow(input) with the kernel output.
using userfun_t = Rcpp::NumericVector (*)(Rcpp::NumericMatrix /*trend_pars*/,
                                          Rcpp::NumericMatrix /*input*/);

// Helper macro to create and return an external pointer to the given function.
// Ensure you also add `// [[Rcpp::export]]` before the declaration line created by this macro
// in your user file so that R can call it

#define EMC2_MAKE_PTR(FUNCNAME)                                       \
SEXP EMC2_make_##FUNCNAME##_ptr() {                                   \
  auto p = new userfun_t(&FUNCNAME);                                  \
  Rcpp::XPtr<userfun_t> xp(p, true);                                  \
  return xp;                                                          \
}
