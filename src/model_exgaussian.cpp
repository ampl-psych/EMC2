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
                    RaceScratch&            scratch)
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
    const double teff  = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = 0.0; continue; }

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
    sc_out[j] = dexg(sc_teff[j], sc_mu[j], sc_sigma[j], sc_tau[j], /*log_d=*/false);
  }

  // --- Scatter ---
  for (int j = 0; j < n_win; ++j) {
    const double val = sc_out[j];
    ll_row[sc_idx[j]] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
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
    const double teff  = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = 1.0; continue; }   // S(t <= 0) = 1

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
                     /*lower_tail=*/false, /*log_p=*/false);
  }

  // --- Scatter ---
  for (int j = 0; j < n_los; ++j) {
    double val = sc_out[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                   val = 1.0;
    ll_row[sc_idx_los[j]] = val;
  }
}


// ---------------------------------------------------------------------------
// exg_survivor
//
// Fills out[i] = S(bound[i] - t0[i]) = 1 - F(bound[i] - t0[i]) for each
// i in idx. Unvisited rows (teff <= 0) are left untouched so the caller's
// default fill (1.0 for lower bound, 0.0 for upper bound) is preserved.
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
                  RaceScratch&               scratch)
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
    const double teff = bound[i] - t0[i];
    if (teff <= 0.0) continue;  // S(t <= 0) = 1; leave out[i] at caller default

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
                     /*lower_tail=*/false, /*log_p=*/false);
  }

  // --- Scatter ---
  for (int j = 0; j < n; ++j) {
    double val = sc_out[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                   val = 1.0;
    out[sc_idx[j]] = val;
  }
}


// ---------------------------------------------------------------------------
// EXGIntegrand — per-trial parameter block for the hcubature callback.
// All arrays are length n_acc, caller-owned (no allocation inside callback).
// ---------------------------------------------------------------------------
struct EXGIntegrand {
  int           n_acc;
  int           winner;   // 0-based within-trial index
  const double* mu;
  const double* sigma;
  const double* tau;
  const double* t0;
};

// hcubature callback: x[0] = absolute time t
// retval[0] = f_winner(t - t0_winner) * prod_{j != winner} S_j(t - t0_j)
inline int int_exg_survivor(unsigned /*dim*/, const double* x, void* p,
                            unsigned /*fdim*/, double* retval)
{
  const EXGIntegrand* P = static_cast<const EXGIntegrand*>(p);
  const double t       = x[0];
  const int    w       = P->winner;
  const double t_eff_w = t - P->t0[w];

  // Density of winner — early exit if non-positive effective time
  if (t_eff_w <= 0.0) { retval[0] = 0.0; return 0; }
  const double dens = dexg(t_eff_w, P->mu[w], P->sigma[w], P->tau[w], /*log_d=*/false);
  if (!(dens > 0.0)) { retval[0] = 0.0; return 0; }

  // Survival of each loser: S_j(t - t0_j) = pexg(..., lower_tail=false)
  // t_eff_j <= 0 means accumulator j hasn't started yet: S_j = 1, skip.
  double out = dens;
  for (int j = 0; j < P->n_acc; ++j) {
    if (j == w) continue;
    const double t_eff_j = t - P->t0[j];
    if (t_eff_j <= 0.0) continue;
    out *= pexg(t_eff_j, P->mu[j], P->sigma[j], P->tau[j],
                /*lower_tail=*/false, /*log_p=*/false);
  }
  retval[0] = out;
  return 0;
}

// Single-trial integration: f_winner * prod S_losers over [lower, upper].
inline double exg_survivor_scalar(const double* mu,  const double* sigma,
                                  const double* tau, const double* t0,
                                  int n_acc, int winner,
                                  double lower, double upper,
                                  double abstol = 1e-8, double reltol = 1e-6,
                                  int maxeval = 6000)
{
  EXGIntegrand ig{ n_acc, winner, mu, sigma, tau, t0 };
  double val = 0.0, err = 0.0;
  hcubature(int_exg_survivor, static_cast<void*>(&ig), 1, &lower, &upper,
            maxeval, abstol, reltol, &val, &err);
  if (!std::isfinite(val) || val < 0.0) return 0.0;
  if (val > 1.0)                         return 1.0;
  return val;
}

// ---------------------------------------------------------------------------
// exg_survivor_with_response
// ---------------------------------------------------------------------------
void exg_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out)
{
  std::vector<double> mu_buf(n_acc), sigma_buf(n_acc),
  tau_buf(n_acc), t0_buf(n_acc);

  const int n = (int)idx.size();
  for (int j = 0; j < n; ++j) {
    const int base = idx[j];  // = trial * n_acc

    for (int acc = 0; acc < n_acc; ++acc) {
      const int row  = base + acc;
      mu_buf   [acc] = pt.base(row, spec.col_mu);
      sigma_buf[acc] = pt.base(row, spec.col_sigma);
      tau_buf  [acc] = pt.base(row, spec.col_tau);
      t0_buf   [acc] = pt.base(row, spec.col_t0);
    }

    out[j] = exg_survivor_scalar(mu_buf.data(), sigma_buf.data(),
                                 tau_buf.data(), t0_buf.data(),
                                 n_acc, winner[j], lower[j], upper[j]);
  }
}

