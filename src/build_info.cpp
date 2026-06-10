// src/build_info.cpp
#include <Rcpp.h>
#include "build_info.h"

// [[Rcpp::export(name = ".emc2_build_info")]]
void emc2_build_info() {
  Rcpp::Rcout
  << "EMC2 build configuration\n"
  << "  PNORM_MODE  : " << EMC2_PNORM_MODE     << "\n"
  << "  fast-math   : " << EMC2_FAST_MATH      << "\n"
  << "  native      : " << EMC2_ENABLE_NATIVE  << "\n"
  << "  CXXFLAGS    : " << EMC2_EXTRA_CXXFLAGS << "\n"
  << "  Accelerate  : " << EMC2_ACCELERATE     << "\n";
}
