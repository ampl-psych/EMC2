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
//#include "pdf_fncs.h"
//#include "fncs_seven.h"
//#include "tools.h"
#include "pnorm_utils.h"
#include "hcubature.h"
#include "RaceSpec.h"

// Struct for cubature
struct my_params {
  double t;
  int low_or_up;
  double a;
  double v;
  double t0;
  double w;
  double sw;
  double sv;
  double st;
  };

// Filler function for the DDM
void fill_ddm(const NumericVector& rts,
              const IntegerVector& R,
              const ParamTable& pt,
              const RaceSpec& spec,
              const std::vector<int>& idx_all,
              double* __restrict__ ll_row);

// Survival function for the DDM
void ddm_survivor(const std::vector<int>& idx,
                  const std::vector<double>& bound,
                  const ParamTable& pt,
                  const RaceSpec& spec,
                  double* __restrict__ out,
                  RaceScratch& scratch);


NumericVector d_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars, std::vector<int> is_ok);

#endif // MODEL_DDM_H
