#include <Rcpp.h>
#include "EMC2/userfun.hpp"   // our typedef

using namespace Rcpp;

// [[Rcpp::export]]
double apply_and_square(double a, double b, SEXP funptrSEXP) {
  // Recover function pointer from external pointer
  XPtr<userfun_t> funptr(funptrSEXP);
  userfun_t f = *funptr;
  if (!f) stop("Null function pointer.");

  // Call user-supplied C++ function
  double result = f(a, b);

  // Square the result
  return result * result;
}
