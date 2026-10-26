#include "model_exgaussian.h"

// Race fill functions for density and survivor
// ---------------------------------------------------------------------------
// dexg_pexg_fast
//
// Fills ll_row for winners (PDF) and losers (survival = 1 - CDF).
// Single flat gather → compute → scatter pass; no A branching.
//
// Parameters read from ParamTable via spec:
//   col_mu   → mu  (Gaussian mean component)
//   col_sigma → sigma (Gaussian SD component)
//   col_tau  → tau  (exponential rate component)
//   col_t0   → t0  (non-decision time)
//
// Convention (mirrors drdm_prdm_fast):
//   winners → ll_row[i] = f(t - t0)          (density,  >= 0)
//   losers  → ll_row[i] = 1 - F(t - t0)      (survival, in [0,1])
// ---------------------------------------------------------------------------
void dexg_pexg_fast(const double*           rt,
                    const ParamTable&       pt,
                    const RaceSpec&         spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__    ll_row,
                    RaceScratch&            scratch,
                    const std::vector<int>& is_ok, bool log_p)
{
  const double* __restrict__ mu    = &pt.base(0, spec.col_mu);
  const double* __restrict__ sigma = &pt.base(0, spec.col_sigma);
  const double* __restrict__ tau   = &pt.base(0, spec.col_tau);
  const double* __restrict__ t0    = &pt.base(0, spec.col_t0);

  // Scratch aliases — single buffer set (no A branching)
  double* __restrict__ sc_teff  = scratch.dslot(REXG::t_eff);
  double* __restrict__ sc_mu    = scratch.dslot(REXG::mu);
  double* __restrict__ sc_sigma = scratch.dslot(REXG::sigma);
  double* __restrict__ sc_tau   = scratch.dslot(REXG::tau);
  double* __restrict__ sc_out   = scratch.dslot(REXG::out);
  int*    __restrict__ sc_idx   = scratch.islot(ScratchInt::idx_win);

  // =========================================================================
  // WINNERS — density
  // =========================================================================

  // --- Gather ---
  int n_win = 0;
  for (int j = 0; j < (int)idx_win.size(); ++j) {
    const int i        = idx_win[j];
    if (!is_ok[i]) continue;
    const double teff  = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = log_p ? neg_inf() : 0.0; continue; }

    sc_teff [n_win] = teff;
    sc_mu   [n_win] = mu[i];
    sc_sigma[n_win] = sigma[i];
    sc_tau  [n_win] = tau[i];
    sc_idx  [n_win] = i;
    n_win++;
  }

  // --- Compute ---
#pragma omp simd
  for (int j = 0; j < n_win; ++j) {
    sc_out[j] = dexg(sc_teff[j], sc_mu[j], sc_sigma[j], sc_tau[j], /*log_d=*/log_p);
  }

  // --- Scatter ---
  for (int j = 0; j < n_win; ++j) {
    const double val = sc_out[j];
    ll_row[sc_idx[j]] = log_p ? val : ((std::isfinite(val) && val >= 0.0) ? val : 0.0);
  }

  // =========================================================================
  // LOSERS — survival (1 - CDF)
  // =========================================================================

  // Reuse scratch idx slot for losers
  int* __restrict__ sc_idx_los = scratch.islot(ScratchInt::idx_los);

  // --- Gather ---
  int n_los = 0;
  for (int j = 0; j < (int)idx_los.size(); ++j) {
    const int i        = idx_los[j];
    if (!is_ok[i]) continue;
    const double teff  = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = log_p ? 0.0 : 1.0; continue; }   // S(t <= 0) = 1

    sc_teff [n_los] = teff;
    sc_mu   [n_los] = mu[i];
    sc_sigma[n_los] = sigma[i];
    sc_tau  [n_los] = tau[i];
    sc_idx_los[n_los] = i;
    n_los++;
  }

  // --- Compute ---
#pragma omp simd
  for (int j = 0; j < n_los; ++j) {
    // upper tail: 1 - F(t) = pexg(..., lower_tail=false)
    sc_out[j] = pexg(sc_teff[j], sc_mu[j], sc_sigma[j], sc_tau[j],
                     /*lower_tail=*/false, /*log_p=*/log_p);
  }

  // --- Scatter ---
  for (int j = 0; j < n_los; ++j) {
    double val = sc_out[j];
    if (!log_p && (!std::isfinite(val) || val < 0.0)) val = 0.0;
    else if (!log_p && val > 1.0)                   val = 1.0;
    ll_row[sc_idx_los[j]] = val;
  }
}

// ---------------------------------------------------------------------------
// exg_survivor
//
// Fills out[i] = S(bound[i] - t0[i]) = 1 - F(bound[i] - t0[i]) for each
// i in idx, optionally on the log scale. Active pre-start bounds have S=1.
//
// Call twice at the call site:
//   exg_survivor(idx_LT + idx_both, LT,  ll_lower, pt, spec, scratch)
//   exg_survivor(idx_UT + idx_both, UT,  ll_upper, pt, spec, scratch)
//
// Convention (mirrors rdm_survivor):
//   out[i] = 1 - F(bound[i] - t0[i])    (survival, in [0,1])
// ---------------------------------------------------------------------------
void exg_survivor(const std::vector<int>&    idx,
                  const std::vector<double>& bound,
                  const ParamTable&          pt,
                  const RaceSpec&            spec,
                  double* __restrict__       out,
                  RaceScratch&               scratch,
                  const std::vector<int>& is_ok, bool log_p)
{
  const double* __restrict__ mu    = &pt.base(0, spec.col_mu);
  const double* __restrict__ sigma = &pt.base(0, spec.col_sigma);
  const double* __restrict__ tau   = &pt.base(0, spec.col_tau);
  const double* __restrict__ t0    = &pt.base(0, spec.col_t0);

  // Scratch aliases — single buffer set
  double* __restrict__ sc_teff  = scratch.dslot(REXG::t_eff);
  double* __restrict__ sc_mu    = scratch.dslot(REXG::mu);
  double* __restrict__ sc_sigma = scratch.dslot(REXG::sigma);
  double* __restrict__ sc_tau   = scratch.dslot(REXG::tau);
  double* __restrict__ sc_out   = scratch.dslot(REXG::out);
  int*    __restrict__ sc_idx   = scratch.islot(ScratchInt::idx_win);

  // --- Gather ---
  int n = 0;
  for (int j = 0; j < (int)idx.size(); ++j) {
    const int i       = idx[j];
    if (!is_ok[i]) continue;
    const double teff = bound[i] - t0[i];
    if (teff <= 0.0) { out[i] = log_p ? 0.0 : 1.0; continue; }  // S(t <= 0) = 1.

    sc_teff [n] = teff;
    sc_mu   [n] = mu[i];
    sc_sigma[n] = sigma[i];
    sc_tau  [n] = tau[i];
    sc_idx  [n] = i;
    n++;
  }

  // --- Compute: upper tail directly to avoid 1 - F cancellation ---
#pragma omp simd
  for (int j = 0; j < n; ++j) {
    sc_out[j] = pexg(sc_teff[j], sc_mu[j], sc_sigma[j], sc_tau[j],
                     /*lower_tail=*/false, /*log_p=*/log_p);
  }

  for (int j = 0; j < n; ++j) out[sc_idx[j]] = sc_out[j];
}
