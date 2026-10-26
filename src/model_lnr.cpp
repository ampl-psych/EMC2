#include "model_lnr.h"

using namespace Rcpp;

// Hot path: gather → compute → scatter
// Uses scratch.v for m, scratch.s for s — B, A, sv left unused.
void dlnr_plnr_fast(const double* rt,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch,
                    const std::vector<int>& is_ok, bool log_p)
{
  const double* __restrict__ m  = &pt.base(0, spec.col_m);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  double* __restrict__ sc_t_eff = scratch.dslot(LNR::t_eff);
  double* __restrict__ sc_m = scratch.dslot(LNR::m);
  double* __restrict__ sc_s = scratch.dslot(LNR::s);
  double* __restrict__ sc_out = scratch.dslot(LNR::out);

  // --- Winners: gather (m → scratch.v, s → scratch.s) ---
  for (int j = 0; j < n_win; ++j) {
    const int i      = idx_win[j];
    sc_t_eff[j] = rt[i] - t0[i];
    sc_m[j]     = m[i];
    sc_s[j]     = s[i];
  }

  // --- Winners: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_win; ++j) {
    if (!is_ok[idx_win[j]] || std::isnan(sc_m[j]) || sc_t_eff[j] <= 0.0) {
      sc_out[j] = log_p ? neg_inf() : 0.0;
      continue;
    }
    const double lt=std::log(sc_t_eff[j]);
    double val = log_p ? LOG_DNORM(lt,sc_m[j],sc_s[j])-lt
                      : DNORM_STD((lt-sc_m[j])/sc_s[j])/(sc_t_eff[j]*sc_s[j]);
    sc_out[j] = log_p ? val : ((std::isfinite(val) && val >= 0.0) ? val : 0.0);
  }

  // --- Winners: scatter ---
  for (int j = 0; j < n_win; ++j) ll_row[idx_win[j]] = sc_out[j];

  // --- Losers: gather ---
  for (int j = 0; j < n_los; ++j) {
    const int i      = idx_los[j];
    sc_t_eff[j] = rt[i] - t0[i];
    sc_m[j]     = m[i];
    sc_s[j]     = s[i];
  }

  // --- Losers: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_los; ++j) {
    if (!is_ok[idx_los[j]] || std::isnan(sc_m[j]) || sc_t_eff[j] <= 0.0) {
      sc_out[j] = log_p ? 0.0 : 1.0;
      continue;
    }
    double val = PNORM_STD((std::log(sc_t_eff[j])-sc_m[j])/sc_s[j],false,log_p);
    if      (!log_p && (!std::isfinite(val) || val < 0.0)) val = 0.0;
    else if (!log_p && val > 1.0)                         val = 1.0;
    sc_out[j] = val;
  }

  // --- Losers: scatter ---
  // Store the directly evaluated survivor.
  for (int j = 0; j < n_los; ++j) ll_row[idx_los[j]] = sc_out[j];
}

// =============================================================================
// LNR — gather-scatter survivor filler
// Fills raw survivors S_k(bound - t0) into out[i] for each i in idx.
// Rows with teff <= 0 get S=1 directly (no CDF call needed).
// =============================================================================
void lnr_survivor(const std::vector<int>&     idx,
                  const std::vector<double>&  bound,
                  const ParamTable&           pt,
                  const RaceSpec&             spec,
                  double* __restrict__        out,
                  RaceScratch&                scratch,
                  const std::vector<int>& is_ok, bool log_p)
{

  const double* __restrict__ bound_ = bound.data();
  const double* __restrict__ m_     = &pt.base(0, spec.col_m);
  const double* __restrict__ t0_    = &pt.base(0, spec.col_t0);
  const double* __restrict__ s_     = &pt.base(0, spec.col_s);

  double* __restrict__ sc_teff = scratch.dslot(LNR::t_eff);
  double* __restrict__ sc_m = scratch.dslot(LNR::m);
  double* __restrict__ sc_s = scratch.dslot(LNR::s);
  double* __restrict__ sc_out = scratch.dslot(LNR::out);
  int*    __restrict__ sc_idx  = scratch.islot(ScratchInt::idx_win);

  // --- Gather: filter teff <= 0, pack valid rows contiguously ---
  int n = 0;
  for (int j = 0; j < (int)idx.size(); ++j) {
    const int i       = idx[j];
    if (!is_ok[i]) continue;
    const double teff = bound_[i] - t0_[i];
    if (teff <= 0.0) { out[i] = log_p ? 0.0 : 1.0; continue; }  // Before the accumulator starts, survival is one.
    sc_teff[n] = teff;
    sc_m   [n] = m_[i];
    sc_s   [n] = s_[i];
    sc_idx [n] = i;
    n++;
  }

  // --- Compute: PLNORM over contiguous scratch — no branching, vectorisable ---
#pragma omp simd
  for (int j = 0; j < n; ++j) {
    sc_out[j] = PNORM_STD((std::log(sc_teff[j])-sc_m[j])/sc_s[j],false,log_p);
  }

  for (int j = 0; j < n; ++j) out[sc_idx[j]] = sc_out[j];
}
