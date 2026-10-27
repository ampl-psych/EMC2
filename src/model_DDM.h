// This DDM implementation is taken from the wienR package for R (https://github.com/RaphaelHartmann/WienR), which is distributed under a GNU GPL-2 license
//
// Original authors of the source code:
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann
//
// The code has been modified to work with EMC2s internal structures

#ifndef MODEL_DDM_H
#define MODEL_DDM_H

#include <Rcpp.h>
using namespace Rcpp;
#include "utility_functions.h"
#include "pnorm_utils.h"
#include "hcubature.h"
#include "RaceSpec.h"

// Struct for cubature
struct my_params {
  double t;
  int low_or_up; // density boundary, -1 or +1
  double a;
  double v;
  double t0;
  double w;
  double sw;
  double sv;
  double st;
  };

// Filler function for the DDM
void fill_ddm(const double* rts,
              const int* R,
              const ParamTable& pt,
              const RaceSpec& spec,
              const std::vector<int>& idx_all,
              double* __restrict__ ll_row,
              const std::vector<int>& is_ok);

double ddm_log_response_interval(const ParamTable& pt,const RaceSpec& spec,int row,
                                 int response,double lower,double upper);

NumericVector d_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars);

#endif // MODEL_DDM_H
