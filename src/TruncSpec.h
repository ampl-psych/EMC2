#ifndef TRUNCSPEC_H
#define TRUNCSPEC_H

#include <vector>
#include <Rcpp.h>
#include "ParamTable.h"
#include "RaceSpec.h"
#include "RaceSetup.h"

struct TruncSpec {
  // Bound vectors — length n_rows (one entry per accumulator row)
  std::vector<double> LT;
  std::vector<double> UT;

  // Index lists — which rows have an active LT / UT bound
  std::vector<int> idx_LT;
  std::vector<int> idx_UT;

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

  bool any() const { return !idx_LT.empty() || !idx_UT.empty(); }

  // Fills p_lower / p_upper, reduces across accumulators, returns log Z per trial
  // (length n_trials, ready to subtract from ll_trial)
  std::vector<double> calculate_normalization_constant() const;
};

TruncSpec make_trunc_spec(const Rcpp::DataFrame& data,
                          int n_trials,
                          int n_acc,
                          const RaceModelSetup& setup,
                          const ParamTable& pt,
                          RaceScratch& scratch);

#endif // TRUNCSPEC_H
