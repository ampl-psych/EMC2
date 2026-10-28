#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "TruncSpec.h"
#include "composite_functions.h"

using namespace Rcpp;

// Retain interval mass on the log scale through likelihood normalization.
void TruncSpec::calculate_normalization_constant(const std::vector<int>& is_ok) const
{
  if (setup->log_response_interval) {
    for (int trial=0;trial<n_trials;++trial) {
      const int base=trial*n_acc;
      log_Z[trial]=0.0;
      if (is_ok[base] && participating[base])
        log_Z[trial]=setup->log_response_interval(*pt,setup->spec,base,0,LT[base],UT[base]);
    }
    return;
  }
  std::fill(log_p_lower.begin(), log_p_lower.end(), 0.0);
  std::fill(log_p_upper.begin(), log_p_upper.end(), neg_inf());
  std::fill(log_S_lower.begin(), log_S_lower.end(), 0.0);
  std::fill(log_S_upper.begin(), log_S_upper.end(), neg_inf());
  setup->fill_survivor(idx_LT, LT, *pt, setup->spec, log_p_lower.data(), *scratch, is_ok, true);
  setup->fill_survivor(idx_UT, UT, *pt, setup->spec, log_p_upper.data(), *scratch, is_ok, true);
  for (int t = 0; t < n_trials; ++t) {
    const int base=t*n_acc;
    log_Z[t]=0.0;
    if (!is_ok[base]) continue;
    double lower=0.0, upper=0.0;
    bool any_participating=false;
    for (int k=0; k<n_acc; ++k) {
      const int row=base+k;
      if (!participating[row]) continue;
      any_participating=true;
      lower+=log_p_lower[row];
      upper+=log_p_upper[row];
    }
    if (!any_participating) continue;
    log_S_lower[t]=lower;
    log_S_upper[t]=upper;
    log_Z[t]=log_diff_exp(lower,upper);
  }
}

// -----------------------------------------------------------------------------
// make_trunc_spec()
// -----------------------------------------------------------------------------
TruncSpec make_trunc_spec(const DataFrame& data,
                          int n_trials,
                          int n_acc,
                          const std::vector<bool>& participating,
                          const RaceModelSetup* setup,
                          const ParamTable& pt,
                          RaceScratch& scratch)
{
  TruncSpec trunc;
  trunc.n_trials      = n_trials;
  trunc.n_acc         = n_acc;
  trunc.n_rows        = n_trials * n_acc;
  trunc.setup         = setup;
  trunc.participating = participating;
  trunc.pt       = &pt;
  trunc.scratch  = &scratch;

  std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(data.names());
  const bool has_LT = std::find(names.begin(), names.end(), "LT") != names.end();
  const bool has_UT = std::find(names.begin(), names.end(), "UT") != names.end();

  // Always allocate all buffers — CensorSpec may read log_S_lower/log_S_upper even
  // when there is no truncation (defaults: log_S_lower=0, log_S_upper=-Inf)
  trunc.log_p_lower.resize(trunc.n_rows);
  trunc.log_p_upper.resize(trunc.n_rows);
  trunc.log_S_lower.resize(trunc.n_trials, 0.0);
  trunc.log_S_upper.resize(trunc.n_trials, neg_inf());
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
    const double lt = (LT_ptr != nullptr) ? LT_ptr[t*n_acc] : 0.0;
    const double ut = (UT_ptr != nullptr) ? UT_ptr[t*n_acc] : R_PosInf;

    const bool active_LT = (lt > 0.0) && std::isfinite(lt);
    const bool active_UT = std::isfinite(ut);

    for (int k = 0; k < n_acc; ++k) {
      const int row = t * n_acc + k;
      if (active_LT) {
        trunc.LT[row] = lt;
        if (participating[row]) trunc.idx_LT.push_back(row);
      }
      if (active_UT) {
        trunc.UT[row] = ut;
        if (participating[row]) trunc.idx_UT.push_back(row);
      }
    }
  }

  return trunc;
}
