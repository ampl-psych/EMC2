#include <Rcpp.h>
#include <algorithm>
#include "CensorSpec.h"

using namespace Rcpp;


// -----------------------------------------------------------------------------
// fill_censored_rows()
//
// 1. Calls setup->fill_survivor for LC and UC bounds into p_lower / p_upper
//    (raw survivors per accumulator row, length n_rows)
// 2. Reduces across n_acc accumulators per trial:
//      prod_LC = prod_k p_lower[base + k]
//      prod_UC = prod_k p_upper[base + k]
// 3. Uses truncation survivors to calculate S_lower = S_RACE(LT) - S_RACE(LC)
//    and S_upper = S_RACE(UC) - S_RACE(UT)
// 4. Fills S_upper in ll_trial
// -----------------------------------------------------------------------------
void CensorSpec::fill_censored_rows(const TruncSpec& trunc,
                                    //std::vector<double>& S_race_UT,
                                    //std::vector<double>& S_race_LT,
                                    NumericVector& ll_trial,
                                    const double min_ll) const
{
  // Fill p_lower: S_k(LC - t0) for active LC rows, 1.0 elsewhere (default -- at t=0, S=1)
  std::fill(p_lower.begin(), p_lower.end(), 1.0);
  setup->fill_survivor(idx_L, LC, *pt, setup->spec, p_lower.data(), *scratch);

  // Fill p_upper: S_k(UC - t0) for active UC rows, 0.0 elsewhere (default -- at t=inf, S=0)
  std::fill(p_upper.begin(), p_upper.end(), 0.0);
  setup->fill_survivor(idx_U, UC, *pt, setup->spec, p_upper.data(), *scratch);

  // ToDo: For cases with known responses, fill p_upper and p_lower with integration methods here

  // Minimal LL
  auto clamp = [min_ll](double v) {
    return (v > min_ll) ? v : min_ll;
  };

  // Lower censoring
  for (int row : idx_L_tr) {
    const int base = row * n_acc;
    double prod_LC = 1.0;
    for(int k = 0; k < n_acc; ++k) {
      prod_LC *= p_lower[base + k];
    }
    const double p = trunc.S_lower[row] - prod_LC; // P(LT <= T <= LC) = S_RACE(LT) - S_RACE(LC)
    ll_trial[row] = std::log(clamp(p));
    }

  // Upper censoring
  for(int row : idx_U_tr) {
    const int base = row * n_acc;
    double prod_UC = 1.0;
    for(int k = 0; k < n_acc; ++k) {
      prod_UC *= p_upper[base + k];
    }
    const double p = prod_UC - trunc.S_upper[row]; // P(UC <= T <= UT) = S_RACE(UC) - S_RACE(UT)
    ll_trial[row] = std::log(clamp(p));
  }

  // Both censoring -- doesn't need truncation
  for(int row : idx_B_tr) {
    const int base = row * n_acc;
    double prod_LC = 1.0, prod_UC = 1.0;
    for(int k = 0; k < n_acc; ++k) {
      prod_UC *= p_upper[base + k];
      prod_LC *= p_lower[base + k];
    }
    const double p = (1-prod_LC) + prod_UC;        // P(T <= LC || P >= UC) = (1-S_RACE(LC)) + S_RACE(UC)
    ll_trial[row] = std::log(clamp(p));

  }

  // -------------------------------------------------------------------------
  // Known-response cases — numerical integration
  // Integrates f_winner(t) * prod_{j != winner} S_j(t) over the censoring window.
  // Bounds are per-trial vectors built from trunc.LT/UT and LC/UC.
  // -------------------------------------------------------------------------

  // Lower-censored, known response: integrate over [LT, LC]
  if (!idx_L_known.empty()) {
    const int n = (int)idx_L_known.size();
    std::vector<double> lo(n), hi(n), out(n);
    for (int j = 0; j < n; ++j) {
      const int base = idx_L_known[j];
      lo[j] = trunc.LT.empty() ? 0.0 : trunc.LT[base];
      hi[j] = LC[base];
    }
    setup->fill_survivor_with_response(idx_L_known, winner_L_known, lo, hi,
                                       n_acc, *pt, setup->spec, out.data());
    for (int j = 0; j < n; ++j)
      ll_trial[idx_L_known[j] / n_acc] = std::log(clamp(out[j]));
  }

  // Upper-censored, known response: integrate over [UC, UT]
  if (!idx_U_known.empty()) {
    const int n = (int)idx_U_known.size();
    std::vector<double> lo(n), hi(n), out(n);

    for (int j = 0; j < n; ++j) {
      const int base = idx_U_known[j];
      lo[j] = UC[base];
      // Don't integrate to inf but to CENS_UPPER_CAP
      hi[j] = std::isfinite(trunc.UT[base]) ? trunc.UT[base] : CENS_UPPER_CAP;
    }
    setup->fill_survivor_with_response(idx_U_known, winner_U_known, lo, hi,
                                       n_acc, *pt, setup->spec, out.data());
    for (int j = 0; j < n; ++j) {
      ll_trial[idx_U_known[j] / n_acc] = std::log(clamp(out[j]));
    }
  }

  // Both-censored, known response: sum integrals over [LT, LC] and [UC, UT]
  if (!idx_B_known.empty()) {
    const int n = (int)idx_B_known.size();
    std::vector<double> lo1(n), hi1(n), lo2(n), hi2(n), out1(n), out2(n);
    for (int j = 0; j < n; ++j) {
      const int base = idx_B_known[j];
      lo1[j] = trunc.LT.empty() ? 0.0           : trunc.LT[base];
      hi1[j] = LC[base];
      lo2[j] = UC[base];
      hi2[j] = std::isfinite(trunc.UT[base]) ? trunc.UT[base] : CENS_UPPER_CAP;
    }
    setup->fill_survivor_with_response(idx_B_known, winner_B_known, lo1, hi1,
                                       n_acc, *pt, setup->spec, out1.data());
    setup->fill_survivor_with_response(idx_B_known, winner_B_known, lo2, hi2,
                                       n_acc, *pt, setup->spec, out2.data());
    for (int j = 0; j < n; ++j)
      ll_trial[idx_B_known[j] / n_acc] = std::log(clamp(out1[j] + out2[j]));
  }
}

// -----------------------------------------------------------------------------
// make_censor_spec()
//
// missingness codes:
//   1 = lower-censored
//   2 = upper-censored
//   3 = both-censored
//
// For each censored trial, data$R determines the path:
//   NA  → analytical (no known response)
//   1,2,... → integration (known response; winner = R - 1, 0-based)
// -----------------------------------------------------------------------------
CensorSpec make_censor_spec(const DataFrame& data,
                            int n_trials,
                            int n_acc,
                            const RaceModelSetup& setup,
                            const ParamTable& pt,
                            RaceScratch& scratch)
{
  CensorSpec censor;
  censor.n_trials = n_trials;
  censor.n_acc    = n_acc;
  censor.n_rows   = n_trials * n_acc;
  censor.setup    = &setup;
  censor.pt       = &pt;
  censor.scratch  = &scratch;

  std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(data.names());
  const bool has_LC_col      = std::find(names.begin(), names.end(), "LC")          != names.end();
  const bool has_UC_col      = std::find(names.begin(), names.end(), "UC")          != names.end();
  const bool has_missingness = std::find(names.begin(), names.end(), "missingness") != names.end();

  if (!has_missingness) return censor;

  IntegerVector missingness_tmp = data["missingness"];
  const int*    missingness_ptr = INTEGER(missingness_tmp);

  IntegerVector R_tmp = data["R"];
  const int*    R_ptr = INTEGER(R_tmp);   // NA_integer_ sentinel for unknown response

  NumericVector lc_tmp, uc_tmp;
  if (has_LC_col) {
    lc_tmp = data["LC"];
    censor.LC.assign(REAL(lc_tmp), REAL(lc_tmp) + censor.n_rows);
    censor.idx_L.reserve(censor.n_rows);
    censor.idx_L_tr.reserve(censor.n_trials);
    censor.idx_L_known.reserve(censor.n_trials);
    censor.winner_L_known.reserve(censor.n_trials);
  }
  if (has_UC_col) {
    uc_tmp = data["UC"];
    censor.UC.assign(REAL(uc_tmp), REAL(uc_tmp) + censor.n_rows);
    censor.idx_U.reserve(censor.n_rows);
    censor.idx_U_tr.reserve(censor.n_trials);
    censor.idx_B_tr.reserve(censor.n_trials);
    censor.idx_U_known.reserve(censor.n_trials);
    censor.idx_B_known.reserve(censor.n_trials);
    censor.winner_U_known.reserve(censor.n_trials);
    censor.winner_B_known.reserve(censor.n_trials);
  }

  // Row-level loop: populate per-accumulator index lists (analytical cases only)
  for (int i = 0; i < censor.n_rows; ++i) {
    switch (missingness_ptr[i]) {
    case 1: censor.idx_L.push_back(i); break;
    case 2: censor.idx_U.push_back(i); break;
    case 3: censor.idx_L.push_back(i); censor.idx_U.push_back(i); break;
    default: break;
    }
  }

  // Trial-level loop (step by n_acc): populate trial index lists
  for (int base = 0; base < censor.n_rows; base += n_acc) {
    const int miss = missingness_ptr[base];
    if (miss == 0) continue;

    const int R_val  = R_ptr[base];
    const bool known = (R_val != NA_INTEGER);
    const int winner = known ? (R_val - 1) : -1;  // 0-based; -1 unused for analytical

    switch (miss) {
    case 1:
      if (known) { censor.idx_L_known.push_back(base); censor.winner_L_known.push_back(winner); }
      else         censor.idx_L_tr.push_back(base / n_acc);
      break;
    case 2:
      if (known) { censor.idx_U_known.push_back(base); censor.winner_U_known.push_back(winner); }
      else         censor.idx_U_tr.push_back(base / n_acc);
      break;
    case 3:
      if (known) { censor.idx_B_known.push_back(base); censor.winner_B_known.push_back(winner); }
      else         censor.idx_B_tr.push_back(base / n_acc);
      break;
    default: break;
    }
  }

  // Validate consistency
  if (!censor.idx_L.empty() && censor.LC.empty())
    Rcpp::stop("missingness contains lower-censored rows (1 or 3) but no LC column found in data");
  if (!censor.idx_U.empty() && censor.UC.empty())
    Rcpp::stop("missingness contains upper-censored rows (2 or 3) but no UC column found in data");

  // Allocate working buffers for row-wise survivors
  censor.p_lower.resize(censor.n_rows);
  censor.p_upper.resize(censor.n_rows);

  return censor;
}


