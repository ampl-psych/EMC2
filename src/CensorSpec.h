#ifndef CENSORSPEC_H
#define CENSORSPEC_H

#include <vector>
#include "RaceSetup.h"
#include "ParamTable.h"


struct CensorSpec {
  std::vector<int>    idx_L;   // lower-censored rows
  std::vector<int>    idx_U;   // upper-censored rows
  std::vector<double> LC;      // lower bounds (length n_trials, or empty)
  std::vector<double> UC;      // upper bounds (length n_trials, or empty)

  std::vector<int>    idx_L_tr;   // lower-censored trials
  std::vector<int>    idx_U_tr;   // upper-censored trials
  std::vector<int>    idx_B_tr;   // both-censored trials

  // Go/no-go withheld trials (missingness == 4): the likelihood is the integral
  // P(no-go accumulator finishes first in the window), not a survivor difference.
  std::vector<int>    idx_G_tr;          // go/nogo withheld trials
  int                 nogo_idx = -1;     // 0-based no-go accumulator within a trial
  std::vector<double> gng_lower;         // integration window lower bound (length n_trials)
  std::vector<double> gng_upper;         // integration window upper bound (length n_trials)

  int n_rows    = 0;
  int n_trials  = 0;
  int n_acc     = 0;

  // Working buffers — length n_rows, filled by fill_truncate calls
  mutable std::vector<double> p_lower;   // S_k(LT - t0) per accumulator row
  mutable std::vector<double> p_upper;   // S_k(UT - t0) per accumulator row

  // References to model objects — set at construction time
  const RaceModelSetup* setup  = nullptr;
  const ParamTable*     pt     = nullptr;
  RaceScratch*          scratch = nullptr;

  bool any() const { return !idx_L.empty() || !idx_U.empty() || !idx_G_tr.empty(); }
  void fill_censored_rows(std::vector<double>& S_race_UT,
                          std::vector<double>& S_race_LT,
                          NumericVector& ll_trial,
                          const double min_ll) const;
};

// CensorSpec make_censor_spec(const Rcpp::DataFrame& data,
//                             int n_trials);

CensorSpec make_censor_spec(const Rcpp::DataFrame& data,
                            int n_trials,
                            int n_acc,
                            const RaceModelSetup& setup,
                            const ParamTable& pt,
                            RaceScratch& scratch);


#endif //CENSORSPEC_H
