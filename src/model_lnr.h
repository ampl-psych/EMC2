#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"  // must be before Rcpp
#include "pnorm_utils.h"
#include "ParamTable.h"
#include "hcubature.h"

using namespace Rcpp;

// Hot path: gather → compute → scatter
// Uses scratch.v for m, scratch.s for s — B, A, sv left unused.
void dlnr_plnr_fast(const double*           rt,
                    const ParamTable&       pt,
                    const RaceSpec&         spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__    ll_row,
                    RaceScratch&            scratch,
                    const std::vector<int>& is_ok, bool log_p = false);

void lnr_survivor(const std::vector<int>&     idx,
                  const std::vector<double>&  bound,
                  const ParamTable&           pt,
                  const RaceSpec&             spec,
                  double* __restrict__        out,
                  RaceScratch&                scratch,
                  const std::vector<int>& is_ok, bool log_p = false);

#endif
