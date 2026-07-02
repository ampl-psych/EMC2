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

  mutable std::vector<double> S_lower;   // S_RACE(LT - t0) per trial
  mutable std::vector<double> S_upper;   // S_RACE(UT - t0) per trial


  // Pre-allocated output buffer — filled by calculate_normalization_constant(),
  // length n_trials. Avoids heap allocation in the particle loop.
  mutable std::vector<double> log_Z;

  // References to model objects — set at construction time, valid for
  // the lifetime of the particle loop (pt is mutated each iteration but
  // the pointer itself remains stable)
  const RaceModelSetup* setup   = nullptr;
  const ParamTable*     pt      = nullptr;
  RaceScratch*          scratch = nullptr;

  bool any() const { return !idx_LT.empty() || !idx_UT.empty(); }

  // Fills p_lower / p_upper, reduces across accumulators, returns log Z per trial
  // (length n_trials, ready to subtract from ll_trial)
  void calculate_normalization_constant() const;
};

TruncSpec make_trunc_spec(const Rcpp::DataFrame& data,
                          int n_trials,
                          int n_acc,
                          const RaceModelSetup& setup,
                          const ParamTable& pt,
                          RaceScratch& scratch);

#endif // TRUNCSPEC_H
