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
                                    std::vector<double>& ll_trial,
                                    const double min_ll) const
{
  // Fill p_lower: S_k(LC - t0) for active LC rows, 1.0 elsewhere (default -- at t=0, S=1)
  std::fill(p_lower.begin(), p_lower.end(), 1.0);
  setup->fill_survivor(idx_L, LC, *pt, setup->spec, p_lower.data(), *scratch);

  // Survivor implementations leave t <= t0 untouched, where S(t - t0) = 1.
  std::fill(p_upper.begin(), p_upper.end(), 1.0);
  setup->fill_survivor(idx_U, UC, *pt, setup->spec, p_upper.data(), *scratch);

  // also fill UC for mixed rows (with silent accumulator)
  setup->fill_survivor(idx_mixed_rows, UC, *pt, setup->spec, p_upper.data(), *scratch);

  // Minimal LL
  auto clamp = [min_ll](double v) {
    return (v > min_ll) ? v : min_ll;
  };

  // Lower censoring
  for (int row : idx_L_tr) {
    const int base = row * n_acc;
    double prod_LC = 1.0;
    for(int k = 0; k < n_acc; ++k) {
      if (!participating[base + k]) continue;  // skip non-participating accumulator
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
      if (!participating[base + k]) continue;  // skip non-participating accumulator
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
      if (!participating[base + k]) continue;  // skip non-participating accumulator
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
    for (int j = 0; j < n; ++j) {
      const int base = idx_L_known[j];
      lo_L[j] = trunc.LT.empty() ? 0.0 : trunc.LT[base];
      hi_L[j] = LC[base];
    }
    setup->fill_survivor_with_response(idx_L_known, winner_L_known, lo_L, hi_L,
                                       n_acc, *pt, setup->spec, out_L.data());
    for (int j = 0; j < n; ++j)
      ll_trial[idx_L_known[j] / n_acc] = std::log(clamp(out_L[j]));
  }

  // Upper-censored, known response: integrate over [UC, UT]
  if (!idx_U_known.empty()) {
    const int n = (int)idx_U_known.size();
    for (int j = 0; j < n; ++j) {
      const int base = idx_U_known[j];
      // Don't integrate to inf but to CENS_UPPER_CAP
      lo_U[j] = UC[base];
      hi_U[j] = (!trunc.UT.empty() && std::isfinite(trunc.UT[base])) ? trunc.UT[base] : CENS_UPPER_CAP;
    }
    setup->fill_survivor_with_response(idx_U_known, winner_U_known, lo_U, hi_U,
                                       n_acc, *pt, setup->spec, out_U.data());
    for (int j = 0; j < n; ++j) {
      ll_trial[idx_U_known[j] / n_acc] = std::log(clamp(out_U[j]));
    }
  }

  // Both-censored, known response: sum integrals over [LT, LC] and [UC, UT]
  if (!idx_B_known.empty()) {
    const int n = (int)idx_B_known.size();
    for (int j = 0; j < n; ++j) {
      const int base = idx_B_known[j];
      lo1_B[j] = trunc.LT.empty() ? 0.0 : trunc.LT[base];
      hi1_B[j] = LC[base];
      lo2_B[j] = UC[base];
      hi2_B[j] = (!trunc.UT.empty() && std::isfinite(trunc.UT[base])) ? trunc.UT[base] : CENS_UPPER_CAP;
    }
    setup->fill_survivor_with_response(idx_B_known, winner_B_known, lo1_B, hi1_B,
                                       n_acc, *pt, setup->spec, out1_B.data());
    setup->fill_survivor_with_response(idx_B_known, winner_B_known, lo2_B, hi2_B,
                                       n_acc, *pt, setup->spec, out2_B.data());
    for (int j = 0; j < n; ++j)
      ll_trial[idx_B_known[j] / n_acc] = std::log(clamp(out1_B[j] + out2_B[j]));
  }

  // Silent accumulators
  if (!idx_mixed_trials.empty()) {
    const int n = (int)silent_acc.size();

    // Build integration bounds — [0, UC] per silent accumulator
    for (int j = 0; j < n; ++j) {
      const int base = silent_trial[j] * n_acc;
      lo_mixed[j] = 0.0;
      hi_mixed[j] = UC[base + silent_acc[j]];
    }

    // One batched integration call over all silent accumulators
    setup->fill_survivor_with_response(
        silent_trial, silent_acc,
        lo_mixed, hi_mixed,
        n_acc, *pt, setup->spec, out_mixed.data()
    );

    // Zero ll_trial for mixed trials before accumulation
    for (int trial : idx_mixed_trials)
      ll_trial[trial] = 0.0;

    // Sum integration results per trial
    for (int j = 0; j < n; ++j)
      ll_trial[silent_trial[j]] += out_mixed[j];

    // Add S_race(UC) analytically and apply log
    for (int trial : idx_mixed_trials) {
      const int base = trial * n_acc;
      double prod_UC = 1.0;
      for (int k = 0; k < n_acc; ++k) {
        if (!participating[base + k]) continue;
        prod_UC *= p_upper[base + k];
      }
      ll_trial[trial] = std::log(clamp(ll_trial[trial] + prod_UC));
    }
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
                            const std::vector<bool>& participating,
                            const RaceModelSetup* setup,
                            const ParamTable& pt,
                            RaceScratch& scratch)
{
  CensorSpec censor;
  censor.n_trials = n_trials;
  censor.n_acc    = n_acc;
  censor.n_rows   = n_trials * n_acc;
  censor.setup    = setup;
  censor.participating = participating;
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
  const int*    R_ptr = INTEGER(R_tmp);

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
    // Aggregate missingness across all accumulator rows for this trial
    bool has_L = false, has_U = false;
    int  R_val = NA_INTEGER;
    for (int k = 0; k < n_acc; ++k) {
      const int m = missingness_ptr[base + k];
      if (m == 1 || m == 3) has_L = true;
      if (m == 2 || m == 3) has_U = true;
      // Take first non-NA R value as the winner
      if (R_val == NA_INTEGER && R_ptr[base + k] != NA_INTEGER)
        R_val = R_ptr[base + k];
    }

    const int miss = (has_L ? 1 : 0) + (has_U ? 2 : 0);
    if (miss == 0) continue;

    const bool known = (R_val != NA_INTEGER);
    const int winner = known ? (R_val - 1) : -1;

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

  // Mixed trials -- go-nogo, sst have silent accumulators
  censor.idx_mixed_trials.reserve(n_trials);
  censor.idx_mixed_rows.reserve(n_trials);
  censor.silent_trial.reserve(n_trials);
  censor.silent_acc.reserve(n_trials);

  for (int base = 0; base < censor.n_rows; base += n_acc) {
    bool has_silent   = false;
    bool has_observed = false;
    for (int k = 0; k < n_acc; ++k) {
      const int m = missingness_ptr[base + k];
      if (m == 4)          has_silent   = true;
      if (m != 4)          has_observed = true;  // NA, 0, 1, 2, 3 all count
      }
    if (!has_silent) continue;

    if (!has_observed)
      Rcpp::stop("Trial %d has all accumulators marked silent (missingness=4) "
                   "with no other accumulators present", base / n_acc);

    const int trial = base / n_acc;
    censor.idx_mixed_trials.push_back(trial);
    censor.idx_mixed_rows.push_back(base);

    for (int k = 0; k < n_acc; ++k) {
      if (missingness_ptr[base + k] == 4) {
        censor.silent_trial.push_back(trial);
        censor.silent_acc.push_back(k);
      }
    }
  }

  // Size integration buffers
  const int n_mixed = (int)censor.silent_acc.size();
  censor.lo_mixed.resize(n_mixed);
  censor.hi_mixed.resize(n_mixed);
  censor.out_mixed.resize(n_mixed);

  // Validate consistency
  if (!censor.idx_L.empty() && censor.LC.empty())
    Rcpp::stop("missingness contains lower-censored rows (1 or 3) but no LC column found in data");
  if (!censor.idx_U.empty() && censor.UC.empty())
    Rcpp::stop("missingness contains upper-censored rows (2 or 3) but no UC column found in data");
  if (!censor.silent_acc.empty() && censor.UC.empty())
    Rcpp::stop("missingness contains silent accumulators (4) but no UC column found in data");

  // Working buffers for analytical survivor products
  censor.p_lower.resize(censor.n_rows);
  censor.p_upper.resize(censor.n_rows);

  // Pre-allocate integration buffers — sized to each known-response index list.
  // These are reused every particle iteration with no further allocation.
  const int n_L = (int)censor.idx_L_known.size();
  const int n_U = (int)censor.idx_U_known.size();
  const int n_B = (int)censor.idx_B_known.size();

  censor.lo_L.resize(n_L);  censor.hi_L.resize(n_L);  censor.out_L.resize(n_L);
  censor.lo_U.resize(n_U);  censor.hi_U.resize(n_U);  censor.out_U.resize(n_U);
  censor.lo1_B.resize(n_B); censor.hi1_B.resize(n_B); censor.out1_B.resize(n_B);
  censor.lo2_B.resize(n_B); censor.hi2_B.resize(n_B); censor.out2_B.resize(n_B);

  return censor;
}



