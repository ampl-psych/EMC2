#ifndef RACE_SETUP_H
#define RACE_SETUP_H

#include "ParamTable.h"
#include "RaceSpec.h"    // ← must come before model headers
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_exgaussian.h"
using namespace Rcpp;

// Combined density and survival fill, using pre-allocated scratch.
using race_combined_fn = void(*)(const double*           rts_ptr,
                                 const ParamTable&       pt,
                                 const RaceSpec&         spec,
                                 const std::vector<int>& idx_win,
                                 const std::vector<int>& idx_los,
                                 double* __restrict__    raw,
                                 RaceScratch&            scratch,
                                 const std::vector<int>& is_ok, bool log_p);

// Survivor function, needed for censoring & truncation
using survivor_fn = void(*)(const std::vector<int>& idx,
                            const std::vector<double>& bound,
                            const ParamTable& pt,
                            const RaceSpec& spec,
                            double* __restrict__ out,
                            RaceScratch& scratch,
                            const std::vector<int>& is_ok, bool log_p);

// Diffusion interval probability: response 0 includes both boundaries.
using response_interval_fn = double(*)(const ParamTable&, const RaceSpec&, int, int,
                                       double, double);

// ---------------------------------------------------------------------------
// RaceModelSetup
// ---------------------------------------------------------------------------

struct RaceModelSetup {
  RaceSpec                 spec;
  race_combined_fn         fill_both;
  survivor_fn              fill_survivor;
  response_interval_fn log_response_interval = nullptr;
};

inline std::pair<int,int> race_scratch_slots(const String& type) {
  if (type == "RDM" || type == "RDM-A0") return {RDM::N_DBL,  RDM::N_INT};
  if (type == "LBA")                      return {LBA::N_DBL,  LBA::N_INT};
  if (type == "LNR")                      return {LNR::N_DBL,  LNR::N_INT};
  if (type == "REXG")                     return {REXG::N_DBL, REXG::N_INT};
  return {1, 1};  // DDM: scratch unused, minimal allocation
}

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
    s.fill_both                   = drdm_prdm_fast;
    s.fill_survivor               = rdm_survivor;
  } else if (type == "LBA") {
    s.spec.col_v        = pt.base_index_for("v");
    s.spec.col_sv       = pt.base_index_for("sv");
    s.spec.col_B        = pt.base_index_for("B");
    s.spec.col_A        = pt.base_index_for("A");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.fill_both                   = dlba_plba_fast;
    s.fill_survivor               = lba_survivor;
  } else if (type == "REXG") {
    s.spec.col_mu       = pt.base_index_for("mu");
    s.spec.col_sigma    = pt.base_index_for("sigma");
    s.spec.col_tau      = pt.base_index_for("tau");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.fill_both                   = dexg_pexg_fast;
    s.fill_survivor               = exg_survivor;
  } else if(type == "DDM") {
    s.spec.col_s        = pt.base_index_for("s");
    s.spec.col_a        = pt.base_index_for("a");
    s.spec.col_v        = pt.base_index_for("v");
    s.spec.col_sv       = pt.base_index_for("sv");
    s.spec.col_Z        = pt.base_index_for("Z");
    s.spec.col_SZ       = pt.base_index_for("SZ");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.spec.col_st0      = pt.base_index_for("st0");
    // no fill function, directly called from the branch in calc_ll instead
    // signature mismatch -- no idx_win / idx_los in the DDM
    s.fill_both                   = nullptr;
    s.fill_survivor               = nullptr;
    s.log_response_interval       = ddm_log_response_interval;
  } else { // LNR
    s.spec.col_m        = pt.base_index_for("m");
    s.spec.col_s        = pt.base_index_for("s");
    s.spec.col_t0       = pt.base_index_for("t0");
    s.fill_both                   = dlnr_plnr_fast;
    s.fill_survivor               = lnr_survivor;
  }

  return s;
}

inline MRISpec make_mri_spec(const ParamTable& pt,
                      const Rcpp::CharacterVector& keep_names,
                      bool is_ar1)
{
  MRISpec spec;
  const int m = keep_names.size();

  // Last column = sigma, second-to-last = rho (AR1 only)
  spec.col_sigma = pt.base_index_for(Rcpp::as<std::string>(keep_names[m - 1]));
  spec.col_rho   = is_ar1 ? pt.base_index_for(Rcpp::as<std::string>(keep_names[m - 2])) : -1;

  const int n_mean = is_ar1 ? m - 2 : m - 1;
  spec.col_means.resize(n_mean);
  spec.n_mean_cols = n_mean;
  for (int j = 0; j < n_mean; ++j)
    spec.col_means[j] = pt.base_index_for(Rcpp::as<std::string>(keep_names[j]));

  return spec;
}

inline ChoiceOnlySpec make_choice_only_spec(const ParamTable&      pt,
                                            const std::string&     type)
{
  ChoiceOnlySpec spec;

  if (type == "MULTINOMIAL_LOGIT" || type == "MULTINOMIAL_PROBIT") {
    spec.col_utility = pt.base_index_for("utility");
  } else {
    spec.col_location = pt.base_index_for("location");
    spec.col_scale    = pt.base_index_for("scale");
    spec.col_cut      = pt.base_index_for("cut");
  }
  return spec;
}


#endif // RACE_SETUP_H
