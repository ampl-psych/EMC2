#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"  // must be before Rcpp
#include "pnorm_utils.h"
#include "ParamTable.h"
// #include "CensorSpec.h"

using namespace Rcpp;

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok);

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok);

// void dlnr_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const LogicalVector& winner,
//                double* ll_row);

// void plnr_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const LogicalVector& winner,
//                double* ll_row);


// Hot path: gather → compute → scatter
void dlnr_plnr_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch);

void plnr_fast(const NumericVector&    rts,
               const ParamTable&       pt,
               const RaceSpec&         spec,
               const std::vector<int>& idx,
               double* __restrict__    ll_row,
               RaceScratch&            scratch);

void lnr_survivor(const std::vector<int>&     idx,
                  const std::vector<double>&  bound,
                  const ParamTable&           pt,
                  const RaceSpec&             spec,
                  double* __restrict__        out,
                  RaceScratch&                scratch);


#endif
