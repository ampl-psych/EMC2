#include <Rcpp.h>
#include "EMC2/userfun.hpp"

// [[Rcpp::export]]
Rcpp::NumericVector EMC2_call_custom_trend(Rcpp::NumericMatrix trend_pars,
                                           Rcpp::NumericMatrix input,
                                           SEXP funptrSEXP) {
  Rcpp::XPtr<userfun_t> funptr(funptrSEXP);
  if (funptr.get() == nullptr) Rcpp::stop("Null function pointer.");
  userfun_t f = *funptr;
  if (!f) Rcpp::stop("Invalid function pointer.");
  Rcpp::NumericVector res = f(trend_pars, input);
  if (res.size() != input.nrow()) {
    Rcpp::stop("Custom trend function must return a vector of length nrow(input).");
  }
  return res;
}
