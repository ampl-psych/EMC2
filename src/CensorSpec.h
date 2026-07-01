#ifndef CENSORSPEC_H
#define CENSORSPEC_H

#include <vector>
#include "RaceSetup.h"
#include "ParamTable.h"
#include "TruncSpec.h"


struct CensorSpec {
  std::vector<int>    idx_L;   // lower-censored rows, no response
  std::vector<int>    idx_U;   // upper-censored rows, no response
  std::vector<int>    idx_B;   // both-censored rows,  no response
  std::vector<double> LC;      // lower bounds (length n_trials, or empty)
  std::vector<double> UC;      // upper bounds (length n_trials, or empty)

  std::vector<int>    idx_L_tr;   // lower-censored trials
  std::vector<int>    idx_U_tr;   // upper-censored trials
  std::vector<int>    idx_B_tr;   // both-censored trials

  // Index lists — known-response cases, trial base rows (t * n_acc)
  // winner_*_known are parallel vectors holding the within-trial winner index [0, n_acc-1]
  std::vector<int> idx_L_known;      std::vector<int> winner_L_known;
  std::vector<int> idx_U_known;      std::vector<int> winner_U_known;
  std::vector<int> idx_B_known;      std::vector<int> winner_B_known;

  // Upper limit of integration for censoring
  double CENS_UPPER_CAP = 30.0;

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

  bool any() const {
    return !idx_L.empty() || !idx_U.empty() || !idx_B.empty() ||
      !idx_L_known.empty() || !idx_U_known.empty() || !idx_B_known.empty();
  }
  void fill_censored_rows(const TruncSpec& trunc,
                          //std::vector<double>& S_race_UT,
                          //std::vector<double>& S_race_LT,
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
