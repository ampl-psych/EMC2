#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "TruncSpec.h"
#include "math_utils.h"

using namespace Rcpp;

// -----------------------------------------------------------------------------
// calculate_normalization_constant()
//
// 1. Calls setup->fill_truncate for LT and UT bounds into p_lower / p_upper
//    (raw survivors per accumulator row, length n_rows)
// 2. Reduces across n_acc accumulators per trial:
//      prod_LT = prod_k p_lower[base + k]
//      prod_UT = prod_k p_upper[base + k]
// 3. Returns log(prod_LT - prod_UT) per trial (length n_trials)
// 4. SM TO DO - *EXCEPT* for censored trials! The normalization constant there should be 0 (so log(1)), because that case is handled by CensorSpec.
//    but we *do* need S_RACE(UT) and S_RACE(LT) there.. although we could just have censor overwrite ll_trial there instead?
// -----------------------------------------------------------------------------
std::vector<double> TruncSpec::calculate_normalization_constant() const
{
  // Fill p_lower: S_k(LT - t0) for active LT rows, 1.0 elsewhere (default)
  std::fill(p_lower.begin(), p_lower.end(), 1.0);
  std::fill(S_lower.begin(), S_lower.end(), 1.0);
  setup->fill_survivor(idx_LT, LT, *pt, setup->spec, p_lower.data(), *scratch);

  // Fill p_upper: S_k(UT - t0) for active UT rows, 0.0 elsewhere (default)
  // Default 0.0 because S(Inf) = 0 — contributes nothing to the product difference
  std::fill(p_upper.begin(), p_upper.end(), 0.0);
  std::fill(S_upper.begin(), S_upper.end(), 0.0);
  setup->fill_survivor(idx_UT, UT, *pt, setup->spec, p_upper.data(), *scratch);

  std::vector<double> log_Z(n_trials);

  for (int t = 0; t < n_trials; ++t) {
    const int base = t * n_acc;
    double prod_LT = 1.0, prod_UT = 1.0;
    for (int k = 0; k < n_acc; ++k) {
      prod_LT *= p_lower[base + k];
      prod_UT *= p_upper[base + k];
    }
    // store trialwise RACE
    S_lower[t] = prod_LT;
    S_upper[t] = prod_UT;

    const double Z = prod_LT - prod_UT;
    // how to deal with Z ~ 0? if left to 0, Z = log(0) = -Inf, which is *subtracted* from the trialwise LL, which then become huge.
    // Some clamping here - but at what value?
    // that is - what probability of the race ending between UT and LT do we minimally require? 1e-2?
    // log_Z[t] = (Z > 1e-2) ? std::log(Z) : 0.0;
    log_Z[t] = (Z > 1e-2) ? Z : 1.0;  // NB: log(1.0) = 0 --> no truncation when LT-UT ~= 0.
  }
  // fast log
  vec_log(log_Z.data(), n_trials);

  return log_Z;
}

// -----------------------------------------------------------------------------
// make_trunc_spec()
// -----------------------------------------------------------------------------
TruncSpec make_trunc_spec(const DataFrame& data,
                          int n_trials,
                          int n_acc,
                          const RaceModelSetup& setup,
                          const ParamTable& pt,
                          RaceScratch& scratch)
{
  TruncSpec trunc;
  trunc.n_trials = n_trials;
  trunc.n_acc    = n_acc;
  trunc.n_rows   = n_trials * n_acc;
  trunc.setup    = &setup;
  trunc.pt       = &pt;
  trunc.scratch  = &scratch;

  std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(data.names());
  const bool has_LT = std::find(names.begin(), names.end(), "LT") != names.end();
  const bool has_UT = std::find(names.begin(), names.end(), "UT") != names.end();

  if (!has_LT && !has_UT) {
    // Fill buffers with default values for censoring
    trunc.S_lower.resize(trunc.n_rows);
    trunc.S_upper.resize(trunc.n_rows);
    std::fill(trunc.S_lower.begin(), trunc.S_lower.end(), 1.0); // at t=0, probability of surviving equals 1
    std::fill(trunc.S_upper.begin(), trunc.S_upper.end(), 0.0); // at t=inf, probability of surviving equals 0

    return trunc;
  }

  NumericVector lt_tmp, ut_tmp;
  const double* LT_ptr = nullptr;
  const double* UT_ptr = nullptr;

  if (has_LT) { lt_tmp = data["LT"]; LT_ptr = REAL(lt_tmp); }
  if (has_UT) { ut_tmp = data["UT"]; UT_ptr = REAL(ut_tmp); }

  trunc.idx_LT.reserve(trunc.n_rows);
  trunc.idx_UT.reserve(trunc.n_rows);

  // Bounds are stored at row (accumulator) level — each trial's n_acc rows
  // share the same LT/UT value, tiled from the trial-level input
  trunc.LT.resize(trunc.n_rows, 0.0);
  trunc.UT.resize(trunc.n_rows, R_PosInf);

  for (int t = 0; t < trunc.n_trials; ++t) {
    const double lt = (LT_ptr != nullptr) ? LT_ptr[t] : 0.0;
    const double ut = (UT_ptr != nullptr) ? UT_ptr[t] : R_PosInf;

    const bool active_LT = (lt > 0.0) && std::isfinite(lt);
    const bool active_UT = std::isfinite(ut);

    for (int k = 0; k < n_acc; ++k) {
      const int row = t * n_acc + k;
      if (active_LT) {
        trunc.LT[row] = lt;
        trunc.idx_LT.push_back(row);
      }
      if (active_UT) {
        trunc.UT[row] = ut;
        trunc.idx_UT.push_back(row);
      }
    }
  }

  // Allocate working buffers
  trunc.p_lower.resize(trunc.n_rows);
  trunc.p_upper.resize(trunc.n_rows);
  trunc.S_lower.resize(trunc.n_rows);
  trunc.S_upper.resize(trunc.n_rows);
  // fill defaults
  std::fill(trunc.S_lower.begin(), trunc.S_lower.end(), 1.0); // at t=0, probability of surviving equals 1
  std::fill(trunc.S_upper.begin(), trunc.S_upper.end(), 0.0); // at t=inf, probability of surviving equals 0

  return trunc;
}
