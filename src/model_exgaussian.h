#ifndef exgaussian_h
#define exgaussian_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"  // must be before Rcpp
#include "exgaussian_functions.h"
#include "ParamTable.h"
#include "hcubature.h"

using namespace Rcpp;

void dexg_pexg_fast(const double*           rt,
                    const ParamTable&       pt,
                    const RaceSpec&         spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__    ll_row,
                    RaceScratch&            scratch);

void exg_survivor(const std::vector<int>&    idx,
                  const std::vector<double>& bound,
                  const ParamTable&          pt,
                  const RaceSpec&            spec,
                  double* __restrict__       out,
                  RaceScratch&               scratch);


void exg_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out);

#endif
