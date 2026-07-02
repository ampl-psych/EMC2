#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "TruncSpec.h"
#include "math_utils.h"

using namespace Rcpp;

// -----------------------------------------------------------------------------
// calculate_normalization_constant()
//
// 1. Calls setup->fill_survivor for LT and UT bounds into p_lower / p_upper
//    (raw survivors per accumulator row, length n_rows)
// 2. Reduces across n_acc accumulators per trial:
//      prod_LT = prod_k p_lower[base + k]
//      prod_UT = prod_k p_upper[base + k]
// 3. Stores S_RACE(LT) and S_RACE(UT) per trial in S_lower / S_upper
//    (used by CensorSpec::fill_censored_rows)
// 4. Writes log(prod_LT - prod_UT) per trial into pre-allocated log_Z member
// -----------------------------------------------------------------------------
void TruncSpec::calculate_normalization_constant() const
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

  for (int t = 0; t < n_trials; ++t) {
    const int base = t * n_acc;
    double prod_LT = 1.0, prod_UT = 1.0;
    for (int k = 0; k < n_acc; ++k) {
      prod_LT *= p_lower[base + k];
      prod_UT *= p_upper[base + k];
    }
    S_lower[t] = prod_LT;
    S_upper[t] = prod_UT;

    const double Z = prod_LT - prod_UT;
    // Clamp: require at least some probability mass in the truncation window.
    // Z ~ 0 means the window [LT, UT] captures almost no mass under current
    // parameters — log(1) = 0 means no truncation correction, which is the
    // least harmful fallback.
    log_Z[t] = (Z > 1e-2) ? Z : 1.0;
  }

  // Fast vectorised log in-place
  vec_log(log_Z.data(), n_trials);
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

  // Always allocate all buffers — CensorSpec may read S_lower/S_upper even
  // when there is no truncation (defaults: S_lower=1, S_upper=0)
  trunc.p_lower.resize(trunc.n_rows);
  trunc.p_upper.resize(trunc.n_rows);
  trunc.S_lower.resize(trunc.n_trials, 1.0);
  trunc.S_upper.resize(trunc.n_trials, 0.0);
  trunc.log_Z.resize(trunc.n_trials);

  if (!has_LT && !has_UT) return trunc;

  NumericVector lt_tmp, ut_tmp;
  const double* LT_ptr = nullptr;
  const double* UT_ptr = nullptr;

  if (has_LT) { lt_tmp = data["LT"]; LT_ptr = REAL(lt_tmp); }
  if (has_UT) { ut_tmp = data["UT"]; UT_ptr = REAL(ut_tmp); }

  trunc.idx_LT.reserve(trunc.n_rows);
  trunc.idx_UT.reserve(trunc.n_rows);

  // Bounds stored at row (accumulator) level — each trial's n_acc rows
  // share the same LT/UT value, tiled from the trial-level input
  trunc.LT.resize(trunc.n_rows, 0.0);
  trunc.UT.resize(trunc.n_rows, R_PosInf);

  for (int t = 0; t < n_trials; ++t) {
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

  return trunc;
}
