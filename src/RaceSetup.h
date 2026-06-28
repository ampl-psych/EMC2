#ifndef RACE_SETUP_H
#define RACE_SETUP_H

#include "ParamTable.h"
#include "RaceSpec.h"    // ← must come before model headers
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
using namespace Rcpp;

// ---------------------------------------------------------------------------
// RaceSpec — parameter column indices in ParamTable
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Dispatcher types
// ---------------------------------------------------------------------------

// Combined pdf+cdf fill — hot path, takes pre-allocated scratch
using race_combined_fn = void(*)(const NumericVector&,
                              const ParamTable&,
                              const RaceSpec&,
                              const std::vector<int>& idx_win,
                              const std::vector<int>& idx_los,
                              double* __restrict__ raw,
                              RaceScratch& scratch);

// using censor_fn = void(*)(const CensorSpec& censor,
//                           const ParamTable& pt,
//                           const RaceSpec& spec,
//                           double* __restrict__ ll_row,
//                           RaceScratch& scratch);

using survivor_fn = void(*)(const std::vector<int>& idx,
                            const std::vector<double>& bound,
                            const ParamTable& pt,
                            const RaceSpec& spec,
                            double* __restrict__ out,
                            RaceScratch& scratch);


// ---------------------------------------------------------------------------
// RaceModelSetup
// ---------------------------------------------------------------------------

struct RaceModelSetup {
  RaceSpec          spec;
  race_combined_fn  fill_both;   // single-pass pdf+cdf — use this in hot path
  survivor_fn       fill_survivor;
};

// ---------------------------------------------------------------------------
// make_race_setup
// ---------------------------------------------------------------------------

inline RaceModelSetup make_race_setup(const String& type, const ParamTable& pt)
{
  RaceModelSetup s;

  if (type == "RDM" || type == "RDM-A0") {
    s.spec.col_v        = pt.base_index_for("v");
    s.spec.col_B        = pt.base_index_for("B");
    s.spec.col_A        = pt.base_index_for("A");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.spec.col_s        = pt.base_index_for("s");
    s.fill_both         = drdm_prdm_fast;
    s.fill_survivor     = rdm_survivor;
  } else if (type == "LBA") {
    s.spec.col_v        = pt.base_index_for("v");
    s.spec.col_sv       = pt.base_index_for("sv");
    s.spec.col_B        = pt.base_index_for("B");
    s.spec.col_A        = pt.base_index_for("A");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.fill_both         = dlba_plba_fast;
    s.fill_survivor     = lba_survivor;
  } else { // LNR
    s.spec.col_m        = pt.base_index_for("m");
    s.spec.col_s        = pt.base_index_for("s");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.fill_both         = dlnr_plnr_fast;
    s.fill_survivor     = lnr_survivor;
  }

  return s;
}

#endif // RACE_SETUP_H
