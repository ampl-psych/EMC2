#include <Rcpp.h>
#include <algorithm>
#include "CensorSpec.h"

using namespace Rcpp;

// CensorSpec make_censor_spec(const DataFrame& data, int n_trials)
// {
//   CensorSpec censor;
//
//   std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(data.names());
//   const bool has_LC_col      = std::find(names.begin(), names.end(), "LC")          != names.end();
//   const bool has_UC_col      = std::find(names.begin(), names.end(), "UC")          != names.end();
//   const bool has_missingness = std::find(names.begin(), names.end(), "missingness") != names.end();
//
//   if (!has_missingness) return censor;
//
//   // Raw pointers — kept alive by the Rcpp wrappers below
//   IntegerVector missingness_tmp = data["missingness"];
//   const int* missingness_ptr    = INTEGER(missingness_tmp);
//
//   NumericVector lc_tmp, uc_tmp;
//   if (has_LC_col) {
//     lc_tmp = data["LC"];
//     censor.LC.assign(REAL(lc_tmp), REAL(lc_tmp) + n_trials);
//     censor.idx_L.reserve(n_trials);
//   }
//   if (has_UC_col) {
//     uc_tmp = data["UC"];
//     censor.UC.assign(REAL(uc_tmp), REAL(uc_tmp) + n_trials);
//     censor.idx_U.reserve(n_trials);
//     censor.idx_B.reserve(n_trials);
//   }
//
//   for (int i = 0; i < n_trials; ++i) {
//     switch (missingness_ptr[i]) {
//     case 1: censor.idx_L.push_back(i); break;
//     case 2: censor.idx_U.push_back(i); break;
//     case 3: censor.idx_B.push_back(i); break;
//     default: break;
//     }
//   }
//
//   // Validate consistency
//   if (!censor.idx_L.empty() && censor.LC.empty())
//     Rcpp::stop("missingness contains lower-censored rows (1 or 3) but no LC column found in data");
//   if (!censor.idx_U.empty() && censor.UC.empty())
//     Rcpp::stop("missingness contains upper-censored rows (2 or 3) but no UC column found in data");
//
//   return censor;
// }

// -----------------------------------------------------------------------------
// fill_censored_rows()
//
// 1. Calls setup->fill_truncate for LC and UC bounds into p_lower / p_upper
//    (raw survivors per accumulator row, length n_rows)
// 2. Reduces across n_acc accumulators per trial:
//      prod_LC = prod_k p_lower[base + k]
//      prod_UC = prod_k p_upper[base + k]
// 3. Uses truncation survivors to calculate S_lower = S_RACE(LT) - S_RACE(LC)
//    and S_upper = S_RACE(UC) - S_RACE(UT)
// 4. Fills S_upper in ll_trial
// -----------------------------------------------------------------------------
void CensorSpec::fill_censored_rows(std::vector<double>& S_race_UT,
                                    std::vector<double>& S_race_LT,
                                    NumericVector& ll_trial,
                                    const double min_ll) const
{
  // Fill p_lower: S_k(LC - t0) for active LC rows, 1.0 elsewhere (default -- at t=0, S=1)
  std::fill(p_lower.begin(), p_lower.end(), 1.0);
  setup->fill_truncate(idx_L, LC, *pt, setup->spec, p_lower.data(), *scratch);

  // Fill p_upper: S_k(UC - t0) for active UC rows, 0.0 elsewhere (default -- at t=inf, S=0)
  std::fill(p_upper.begin(), p_upper.end(), 0.0);
  setup->fill_truncate(idx_U, UC, *pt, setup->spec, p_upper.data(), *scratch);

  // ToDo: For cases with known responses, fill p_upper and p_lower with integration methods here


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
    const double p = S_race_LT[row] - prod_LC; // P(LT <= T <= LC) = S_RACE(LT) - S_RACE(LC)
    ll_trial[row] = std::log(clamp(p));
    }

  // Upper censoring
  for(int row : idx_U_tr) {
    const int base = row * n_acc;
    double prod_UC = 1.0;
    for(int k = 0; k < n_acc; ++k) {
      prod_UC *= p_upper[base + k];
    }
    const double p = prod_UC - S_race_UT[row]; // P(UC <= T <= UT) = S_RACE(UC) - S_RACE(UT)
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
    //
    // Rcpp::Rcout << "trial " << row << ", "
    //             << " LC = " << LC[row] << ", "
    //             << " UC = " << UC[row] << ", "
    //             << " prod(LC) = " << prod_LC << ", "
    //             << " prod(UC) = " << prod_UC << ", "
    //             << " p = " << (1-prod_LC) + prod_UC << '\n';

  }

  // additional cases here...
}

// -----------------------------------------------------------------------------
// make_censor_spec()
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

  // Raw pointers — kept alive by the Rcpp wrappers below
  IntegerVector missingness_tmp = data["missingness"];
  const int* missingness_ptr    = INTEGER(missingness_tmp);

  NumericVector lc_tmp, uc_tmp;
  if (has_LC_col) {
    lc_tmp = data["LC"];
    censor.LC.assign(REAL(lc_tmp), REAL(lc_tmp) + censor.n_rows);
    censor.idx_L.reserve(censor.n_rows);
    censor.idx_L_tr.reserve(censor.n_trials);
  }
  if (has_UC_col) {
    uc_tmp = data["UC"];
    censor.UC.assign(REAL(uc_tmp), REAL(uc_tmp) + censor.n_rows);
    censor.idx_U.reserve(censor.n_rows);
    // censor.idx_B.reserve(censor.n_rows);
    censor.idx_U_tr.reserve(censor.n_trials);
    censor.idx_B_tr.reserve(censor.n_trials);
  }

  // fill row indices used to calculate per-accumulator survival S_k
  for (int i = 0; i < censor.n_rows; ++i) {
    switch (missingness_ptr[i]) {
    case 1: censor.idx_L.push_back(i); break;
    case 2: censor.idx_U.push_back(i); break;
    case 3:
      censor.idx_L.push_back(i);
      censor.idx_U.push_back(i); break;
    default: break;
    }
  }

  // Fill trialwise indices used to calculate per-trial survival S_RACE
  for (int i = 0; i < censor.n_rows; i += n_acc) {
    switch (missingness_ptr[i]) {
    case 1: censor.idx_L_tr.push_back(i / n_acc); break;
    case 2: censor.idx_U_tr.push_back(i / n_acc); break;
    case 3: censor.idx_B_tr.push_back(i / n_acc); break;
    default: break;
    }
  }

  // Validate consistency
  if (!censor.idx_L.empty() && censor.LC.empty())
    Rcpp::stop("missingness contains lower-censored rows (1 or 3) but no LC column found in data");
  if (!censor.idx_U.empty() && censor.UC.empty())
    Rcpp::stop("missingness contains upper-censored rows (2 or 3) but no UC column found in data");


  // Allocate working buffers for row-wise survivors (lower, upper) and race-wise (trialwise) surivors
  censor.p_lower.resize(censor.n_rows);
  censor.p_upper.resize(censor.n_rows);
  // censor.S_lower.resize(censor.n_rows);
  // censor.S_upper.resize(censor.n_rows);

  return censor;
}
