#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "RaceSpec.h"
#include "math_utils.h"  // must be before Rcpp
#include "wald_functions.h"
#include "ParamTable.h"
#include "hcubature.h"

using namespace Rcpp;

// This new filling function checks whether A==0, if so --> runs digt0 and pigt0
void drdm_prdm_fast(const double*           rt,
                    const ParamTable&       pt,
                    const RaceSpec&         spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__    ll_row,
                    RaceScratch&            scratch,
                    const std::vector<int>& is_ok, bool log_p = false);

void rdm_survivor(const std::vector<int>& idx,
                  const std::vector<double>& bound,
                  const ParamTable& pt,
                  const RaceSpec& spec,
                  double* __restrict__ out,
                  RaceScratch& scratch,
                  const std::vector<int>& is_ok, bool log_p = false);

#endif // rdm_h
