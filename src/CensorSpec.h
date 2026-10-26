#ifndef CENSORSPEC_H
#define CENSORSPEC_H

#include <vector>
#include "RaceSetup.h"
#include "ParamTable.h"
#include "TruncSpec.h"


struct CensorSpec {
  std::vector<int>    idx_L;   // lower-censored rows, no response
  std::vector<int>    idx_U;   // upper-censored rows, no response
  std::vector<double> LC;      // lower bounds (length n_rows, or empty)
  std::vector<double> UC;      // upper bounds (length n_rows, or empty)

  std::vector<int>    idx_L_tr;   // lower-censored trials
  std::vector<int>    idx_U_tr;   // upper-censored trials
  std::vector<int>    idx_B_tr;   // both-censored trials

  // Index lists — known-response cases, trial base rows (t * n_acc)
  // winner_*_known are parallel vectors holding the within-trial winner index [0, n_acc-1]
  std::vector<int> idx_L_known;      std::vector<int> winner_L_known;
  std::vector<int> idx_U_known;      std::vector<int> winner_U_known;
  std::vector<int> idx_B_known;      std::vector<int> winner_B_known;
  std::vector<bool> participating;

  std::vector<int>  idx_silent;       // trial base rows for silent accumulators
  std::vector<int>  winner_silent;    // parallel: accumulator index [0, n_acc-1]
  int n_rows    = 0;
  int n_trials  = 0;
  int n_acc     = 0;

  // Log survival buffers, also reused by the known-response integrator
  mutable std::vector<double> log_p_lower;   // log S_k(LC - t0) per accumulator row
  mutable std::vector<double> log_p_upper;   // log S_k(UC - t0) per accumulator row

  // Shared integration work, allocated once for each likelihood worker.
  mutable std::vector<double> bound_work, breaks;
  mutable std::vector<int> active_rows, winner_row;

  double log_interval(int base, int winner, double lower, double upper,
                      const std::vector<int>& is_ok) const;

  // References to model objects — set at construction time, valid for
  // the lifetime of the particle loop (pt is mutated each iteration but
  // the pointer itself remains stable)
  const RaceModelSetup* setup   = nullptr;
  const ParamTable*     pt      = nullptr;
  RaceScratch*          scratch = nullptr;

  bool any() const {
    return !idx_L.empty() || !idx_U.empty() ||
      !idx_L_known.empty() || !idx_U_known.empty() || !idx_B_known.empty() ||
      !idx_silent.empty();
  }

  void fill_censored_rows(const TruncSpec& trunc,
                          const std::vector<int>& is_ok,
                          std::vector<double>& ll_trial) const;

  // Update paramtable for multithreading
  void rebind(ParamTable& pt, RaceScratch& scratch) {
    this->pt      = &pt;
    this->scratch = &scratch;
  }
};


CensorSpec make_censor_spec(const Rcpp::DataFrame& data,
                            int n_trials,
                            int n_acc,
                            const std::vector<bool>& participating,
                            const RaceModelSetup* setup,
                            const ParamTable& pt,
                            RaceScratch& scratch);


#endif //CENSORSPEC_H
