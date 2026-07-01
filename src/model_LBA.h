#ifndef model_LBA_h
#define model_LBA_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"
#include "pnorm_utils.h"
#include "ParamTable.h"
#include "hcubature.h"

constexpr double LBA_A_ASYMPTOTIC = 1e-10;

// R exports — declared here, defined in model_LBA.cpp
Rcpp::NumericVector dlba(Rcpp::NumericVector t,
                         Rcpp::NumericVector A, Rcpp::NumericVector b,
                         Rcpp::NumericVector v, Rcpp::NumericVector sv,
                         bool posdrift);

Rcpp::NumericVector plba(Rcpp::NumericVector t,
                         Rcpp::NumericVector A, Rcpp::NumericVector b,
                         Rcpp::NumericVector v, Rcpp::NumericVector sv,
                         bool posdrift);

// ---------------------------------------------------------------------------
// Split compute functions — small, always inlined, no ODR issue in header
// ---------------------------------------------------------------------------

[[gnu::always_inline]] inline double dlba_core(double t, double A, double b,
                                               double v, double sv, double denom)
{
  const double zs     = t * sv;
  const double cmz    = b - t * v;
  const double cz     = cmz / zs;
  const double cz_max = (cmz - A) / zs;
  return (v * (PNORM_STD(cz,     true, false) - PNORM_STD(cz_max, true, false))
            + sv * (DNORM_STD(cz_max) - DNORM_STD(cz))) / (A * denom);
}

[[gnu::always_inline]] inline double dlba_noA(double t, double b,
                                              double v, double sv, double denom)
{
  return DNORM(b / t, v, sv) * b / (t * t * denom);
}

[[gnu::always_inline]] inline double plba_core(double t, double A, double b,
                                               double v, double sv, double denom)
{
  const double zs     = t * sv;
  const double cmz    = b - t * v;
  const double xx     = cmz - A;
  const double cz     = cmz / zs;
  const double cz_max = xx / zs;
  return (1.0 + (zs * (DNORM_STD(cz_max) - DNORM_STD(cz))
                   + xx  * PNORM_STD(cz_max, true, false)
                   - cmz * PNORM_STD(cz,     true, false)) / A) / denom;
}

[[gnu::always_inline]] inline double plba_noA(double t, double b,
                                              double v, double sv, double denom)
{
  return PNORM_STD((b / t - v) / sv, false, false) / denom;
}

// ---------------------------------------------------------------------------
// Function declarations — defined in model_LBA.cpp
// ---------------------------------------------------------------------------

void dlba_plba_fast(const Rcpp::NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch);

void plba_fast(const Rcpp::NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const std::vector<int>& idx,
               double* __restrict__ ll_row,
               RaceScratch& scratch);

void lba_survivor(const std::vector<int>& idx,
                  const std::vector<double>& bound,
                  const ParamTable& pt,
                  const RaceSpec& spec,
                  double* __restrict__ out,
                  RaceScratch& scratch);

void lba_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out);

#endif // model_LBA_h
