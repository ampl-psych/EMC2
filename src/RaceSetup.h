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

// Legacy per-accumulator fill (kept for drdm_fast / prdm_fast etc.)
typedef void (*race_fast_fn)(const NumericVector& rts,
              const ParamTable& pt,
              const RaceSpec& spec,
              const LogicalVector& winner,
              double* raw);

// Combined pdf+cdf fill — hot path, takes pre-allocated scratch
using race_combined_fn = void(*)(const NumericVector&,
                              const ParamTable&,
                              const RaceSpec&,
                              const std::vector<int>& idx_win,
                              const std::vector<int>& idx_los,
                              double* __restrict__ raw,
                              RaceScratch& scratch);

// ---------------------------------------------------------------------------
// RaceModelSetup
// ---------------------------------------------------------------------------

struct RaceModelSetup {
  RaceSpec         spec;
  race_fast_fn     fill_pdf;
  race_fast_fn     fill_cdf;
  race_combined_fn fill_both;   // single-pass pdf+cdf — use this in hot path
  int              col_na_marker;
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
    s.col_na_marker     = s.spec.col_v;
    s.fill_pdf          = drdm_fast;
    s.fill_cdf          = prdm_fast;
    s.fill_both         = drdm_prdm_fast;
  } else if (type == "LBA") {
    s.spec.col_v        = pt.base_index_for("v");
    s.spec.col_sv       = pt.base_index_for("sv");
    s.spec.col_B        = pt.base_index_for("B");
    s.spec.col_A        = pt.base_index_for("A");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.col_na_marker     = s.spec.col_v;
    s.fill_pdf          = dlba_fast;
    s.fill_cdf          = plba_fast;
    s.fill_both         = dlba_plba_fast;
  } else { // LNR
    s.spec.col_m        = pt.base_index_for("m");
    s.spec.col_s        = pt.base_index_for("s");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.col_na_marker     = s.spec.col_m;
    s.fill_pdf          = dlnr_fast;
    s.fill_cdf          = plnr_fast;
    s.fill_both         = dlnr_plnr_fast;
  }

  return s;
}

#endif // RACE_SETUP_H
