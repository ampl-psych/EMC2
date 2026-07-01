#include "model_LBA.h"

using namespace Rcpp;


// R exports

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)
{
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    double denom = 1.0;
    if (posdrift) {
      denom = PNORM_STD(v[i] / sv[i], /*lower=*/true, /*logp=*/false);
      if (denom < 1e-10) denom = 1e-10;
    }
    double val = (A[i] > LBA_A_ASYMPTOTIC)
      ? dlba_core(t[i], A[i], b[i], v[i], sv[i], denom)
        : dlba_noA (t[i],       b[i], v[i], sv[i], denom);
    if (val < 0.0) val = 0.0;
    pdf[i] = val;
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)
{
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    double denom = 1.0;
    if (posdrift) {
      denom = PNORM_STD(v[i] / sv[i], /*lower=*/true, /*logp=*/false);
      if (denom < 1e-10) denom = 1e-10;
    }
    double val = (A[i] > LBA_A_ASYMPTOTIC)
      ? plba_core(t[i], A[i], b[i], v[i], sv[i], denom)
        : plba_noA (t[i],       b[i], v[i], sv[i], denom);
    if      (val < 0.0) val = 0.0;
    else if (val > 1.0) val = 1.0;
    cdf[i] = val;
  }
  return cdf;
}


// ---------------------------------------------------------------------------
// Hot path: one-pass gather split on A, core dominates -> primary scratch
// ---------------------------------------------------------------------------
void dlba_plba_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch)
{
  // Note that for the losers, the *SURVIVAL* probability is filled, *NOT* the CDF
  const double* __restrict__ rt  = rts.begin();
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // primary scratch — core entries (A > LBA_A_ASYMPTOTIC)
  double* __restrict__ sc_teff  = scratch.t_eff.data();
  double* __restrict__ sc_v     = scratch.v.data();
  double* __restrict__ sc_sv    = scratch.sv.data();
  double* __restrict__ sc_B     = scratch.B.data();
  double* __restrict__ sc_A     = scratch.A.data();
  double* __restrict__ sc_out   = scratch.out.data();
  int*    __restrict__ sc_idx   = scratch.idx_win0.data();

  // secondary scratch — noA entries (A <= LBA_A_ASYMPTOTIC)
  // note: sc_sv_c reuses s_c since noA path needs sv but not s
  double* __restrict__ sc_teff_c = scratch.t_eff_c.data();
  double* __restrict__ sc_v_c    = scratch.v_c.data();
  double* __restrict__ sc_sv_c   = scratch.s_c.data();
  double* __restrict__ sc_B_c    = scratch.B_c.data();
  double* __restrict__ sc_out_c  = scratch.out_c.data();
  int*    __restrict__ sc_idx_c  = scratch.idx_win_c.data();

  // =========================================================================
  // WINNERS
  // =========================================================================

  // --- One-pass gather: core -> primary, noA -> secondary ---
  int n_win_core = 0, n_win_noA = 0;
  for (int j = 0; j < n_win; ++j) {
    const int i       = idx_win[j];
    const double teff = rt[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { ll_row[i] = 0.0; continue; }
    if (A[i] > LBA_A_ASYMPTOTIC) {
      sc_teff[n_win_core] = teff;
      sc_v   [n_win_core] = v[i];
      sc_sv  [n_win_core] = sv[i];
      sc_B   [n_win_core] = B[i];
      sc_A   [n_win_core] = A[i];
      sc_idx [n_win_core] = i;
      n_win_core++;
    } else {
      sc_teff_c[n_win_noA] = teff;
      sc_v_c   [n_win_noA] = v[i];
      sc_sv_c  [n_win_noA] = sv[i];
      sc_B_c   [n_win_noA] = B[i];
      sc_idx_c [n_win_noA] = i;
      n_win_noA++;
    }
  }

  // --- Compute: dlba_core winners ---
  for (int j = 0; j < n_win_core; ++j) {
    double denom = PNORM_STD(sc_v[j] / sc_sv[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = dlba_core(sc_teff[j], sc_A[j], sc_B[j] + sc_A[j],
                           sc_v[j], sc_sv[j], denom);
    sc_out[j] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
  }

  // --- Compute: dlba_noA winners ---
  for (int j = 0; j < n_win_noA; ++j) {
    double denom = PNORM_STD(sc_v_c[j] / sc_sv_c[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = dlba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], denom);
    sc_out_c[j] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
  }

  // --- Scatter ---
  for (int j = 0; j < n_win_core; ++j) ll_row[sc_idx  [j]] = sc_out  [j];
  for (int j = 0; j < n_win_noA;  ++j) ll_row[sc_idx_c[j]] = sc_out_c[j];

  // =========================================================================
  // LOSERS
  // =========================================================================

  sc_idx   = scratch.idx_los0.data();
  sc_idx_c = scratch.idx_los_c.data();

  // --- One-pass gather: core -> primary, noA -> secondary ---
  int n_los_core = 0, n_los_noA = 0;
  for (int j = 0; j < n_los; ++j) {
    const int i       = idx_los[j];
    const double teff = rt[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { ll_row[i] = 1.0; continue; }
    if (A[i] > LBA_A_ASYMPTOTIC) {
      sc_teff[n_los_core] = teff;
      sc_v   [n_los_core] = v[i];
      sc_sv  [n_los_core] = sv[i];
      sc_B   [n_los_core] = B[i];
      sc_A   [n_los_core] = A[i];
      sc_idx [n_los_core] = i;
      n_los_core++;
    } else {
      sc_teff_c[n_los_noA] = teff;
      sc_v_c   [n_los_noA] = v[i];
      sc_sv_c  [n_los_noA] = sv[i];
      sc_B_c   [n_los_noA] = B[i];
      sc_idx_c [n_los_noA] = i;
      n_los_noA++;
    }
  }

  // --- Compute: plba_core losers ---
  for (int j = 0; j < n_los_core; ++j) {
    double denom = PNORM_STD(sc_v[j] / sc_sv[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = plba_core(sc_teff[j], sc_A[j], sc_B[j] + sc_A[j],
                           sc_v[j], sc_sv[j], denom);
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    sc_out[j] = 1.0 - val;
  }

  // --- Compute: plba_noA losers ---
  for (int j = 0; j < n_los_noA; ++j) {
    double denom = PNORM_STD(sc_v_c[j] / sc_sv_c[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = plba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], denom);
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    sc_out_c[j] = 1.0 - val;
  }

  // --- Scatter ---
  for (int j = 0; j < n_los_core; ++j) ll_row[sc_idx  [j]] = sc_out  [j];
  for (int j = 0; j < n_los_noA;  ++j) ll_row[sc_idx_c[j]] = sc_out_c[j];
}




// Censoring, truncation
// =============================================================================
// plba_fast — gather-scatter, core/noA split
// denom = PNORM_STD(v/sv) computed in gather, stored in scratch.
// Compute loop is branch-free and vectorisable.
// =============================================================================
void plba_fast(const NumericVector&    rts,
               const ParamTable&       pt,
               const RaceSpec&         spec,
               const std::vector<int>& idx,
               double* __restrict__    ll_row,
               RaceScratch&            scratch)
{
  const double* __restrict__ rt  = rts.begin();
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);

  const int n = (int)idx.size();

  double* __restrict__ sc_teff  = scratch.t_eff.data();
  double* __restrict__ sc_v     = scratch.v.data();
  double* __restrict__ sc_sv    = scratch.sv.data();
  double* __restrict__ sc_B     = scratch.B.data();
  double* __restrict__ sc_A     = scratch.A.data();
  double* __restrict__ sc_denom = scratch.denom.data();   // new
  double* __restrict__ sc_out   = scratch.out.data();
  int*    __restrict__ sc_idx   = scratch.idx_win0.data();

  double* __restrict__ sc_teff_c  = scratch.t_eff_c.data();
  double* __restrict__ sc_v_c     = scratch.v_c.data();
  double* __restrict__ sc_sv_c    = scratch.s_c.data();
  double* __restrict__ sc_B_c     = scratch.B_c.data();
  double* __restrict__ sc_denom_c = scratch.denom_c.data(); // new
  double* __restrict__ sc_out_c   = scratch.out_c.data();
  int*    __restrict__ sc_idx_c   = scratch.idx_win_c.data();

  // --- Gather: compute teff, denom; split core / noA ---
  int n_core = 0, n_noA = 0;
  for (int j = 0; j < n; ++j) {
    const int i       = idx[j];
    const double teff = rt[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { ll_row[i] = 1.0; continue; }
    double denom = PNORM_STD(v[i] / sv[i], true, false);
    if (denom < 1e-10) denom = 1e-10;
    if (A[i] > LBA_A_ASYMPTOTIC) {
      sc_teff [n_core] = teff;  sc_v[n_core]     = v[i];  sc_sv[n_core]    = sv[i];
      sc_B    [n_core] = B[i];  sc_A[n_core]     = A[i];  sc_denom[n_core] = denom;
      sc_idx  [n_core] = i;
      n_core++;
    } else {
      sc_teff_c [n_noA] = teff; sc_v_c[n_noA]    = v[i];  sc_sv_c[n_noA]   = sv[i];
      sc_B_c    [n_noA] = B[i]; sc_denom_c[n_noA] = denom; sc_idx_c[n_noA] = i;
      n_noA++;
    }
  }

  // --- Compute: branch-free, vectorisable ---
  for (int j = 0; j < n_core; ++j)
    sc_out[j] = plba_core(sc_teff[j], sc_A[j], sc_B[j] + sc_A[j], sc_v[j], sc_sv[j], sc_denom[j]);
  for (int j = 0; j < n_noA; ++j)
    sc_out_c[j] = plba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], sc_denom_c[j]);

  // --- Scatter: clamp, write survivor ---
  for (int j = 0; j < n_core; ++j) {
    double val = sc_out[j];
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    ll_row[sc_idx[j]] = 1.0 - val;
  }
  for (int j = 0; j < n_noA; ++j) {
    double val = sc_out_c[j];
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    ll_row[sc_idx_c[j]] = 1.0 - val;
  }
}


// =============================================================================
// lba_fill_truncate — gather-scatter, core/noA split
// Fills raw survivors S_k(bound - t0) into out[] for each i in idx.
// =============================================================================
void lba_survivor(const std::vector<int>&    idx,
                  const std::vector<double>&  bound,
                  const ParamTable&           pt,
                  const RaceSpec&             spec,
                  double* __restrict__        out,
                  RaceScratch&                scratch)
{
  const double* __restrict__ bound_ = bound.data();
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);

  const int n = (int)idx.size();

  double* __restrict__ sc_teff   = scratch.t_eff.data();
  double* __restrict__ sc_v      = scratch.v.data();
  double* __restrict__ sc_sv     = scratch.sv.data();
  double* __restrict__ sc_B      = scratch.B.data();
  double* __restrict__ sc_A      = scratch.A.data();
  double* __restrict__ sc_denom  = scratch.denom.data();
  double* __restrict__ sc_out    = scratch.out.data();
  int*    __restrict__ sc_idx    = scratch.idx_win0.data();

  double* __restrict__ sc_teff_c  = scratch.t_eff_c.data();
  double* __restrict__ sc_v_c     = scratch.v_c.data();
  double* __restrict__ sc_sv_c    = scratch.s_c.data();
  double* __restrict__ sc_B_c     = scratch.B_c.data();
  double* __restrict__ sc_denom_c = scratch.denom_c.data();
  double* __restrict__ sc_out_c   = scratch.out_c.data();
  int*    __restrict__ sc_idx_c   = scratch.idx_win_c.data();

  // --- Gather ---
  int n_core = 0, n_noA = 0;
  for (int j = 0; j < n; ++j) {
    const int i       = idx[j];
    const double teff = bound_[i] - t0[i];
    // if (std::isnan(v[i]) || teff <= 0.0) { out[i] = 1.0; continue; }
    if (std::isnan(v[i]) || teff <= 0.0) { continue; }  // skip and leave out untouched - reverts to default value
    double denom = PNORM_STD(v[i] / sv[i], true, false);
    if (denom < 1e-10) denom = 1e-10;
    if (A[i] > LBA_A_ASYMPTOTIC) {
      sc_teff[n_core]  = teff; sc_v[n_core]      = v[i];  sc_sv[n_core]    = sv[i];
      sc_B[n_core]     = B[i]; sc_A[n_core]      = A[i];  sc_denom[n_core] = denom;
      sc_idx[n_core]   = i;    n_core++;
    } else {
      sc_teff_c[n_noA]  = teff; sc_v_c[n_noA]    = v[i];  sc_sv_c[n_noA]   = sv[i];
      sc_B_c[n_noA]     = B[i]; sc_denom_c[n_noA] = denom; sc_idx_c[n_noA] = i;
      n_noA++;
    }
  }

  // --- Compute: branch-free, vectorisable ---
  for (int j = 0; j < n_core; ++j)
    sc_out[j] = plba_core(sc_teff[j], sc_A[j], sc_B[j] + sc_A[j], sc_v[j], sc_sv[j], sc_denom[j]);
  for (int j = 0; j < n_noA; ++j)
    sc_out_c[j] = plba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], sc_denom_c[j]);

  // --- Scatter: clamp, write survivor ---
  for (int j = 0; j < n_core; ++j) {
    double val = sc_out[j];
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    out[sc_idx[j]] = 1.0 - val;
  }
  for (int j = 0; j < n_noA; ++j) {
    double val = sc_out_c[j];
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    out[sc_idx_c[j]] = 1.0 - val;
  }
}

// -----------------------------------------------------------------------------
// LBA survivor with known response — numerical integration
//
// Computes P(winner finishes first in [lower, upper]) for one trial by
// integrating f_winner(t) * prod_{j != winner} S_j(t) over [lower, upper].
// Parameters: v, sv, B, A, t0 — core/noA split on A > LBA_A_ASYMPTOTIC.
// -----------------------------------------------------------------------------

struct LBAIntegrand {
  int           n_acc;
  int           winner;   // 0-based within-trial index
  const double* v;
  const double* sv;
  const double* B;        // threshold = B + A (stored as raw B from ParamTable)
  const double* A;
  const double* t0;
};

// hcubature callback: x[0] = absolute time t
// retval[0] = f_winner(t - t0_winner) * prod_{j != winner} S_j(t - t0_j)
inline int int_lba_survivor(unsigned /*dim*/, const double* x, void* p,
                            unsigned /*fdim*/, double* retval)
{
  const LBAIntegrand* P = static_cast<const LBAIntegrand*>(p);
  const double t       = x[0];
  const int    w       = P->winner;
  const double t_eff_w = t - P->t0[w];

  // Density of winner — early exit if t <= t0
  if (t_eff_w <= 0.0) { retval[0] = 0.0; return 0; }
  double denom_w = PNORM_STD(P->v[w] / P->sv[w], true, false);
  if (denom_w < 1e-10) denom_w = 1e-10;
  const double dens = (P->A[w] > LBA_A_ASYMPTOTIC)
    ? dlba_core(t_eff_w, P->A[w], P->B[w] + P->A[w], P->v[w], P->sv[w], denom_w)
      : dlba_noA (t_eff_w,          P->B[w],             P->v[w], P->sv[w], denom_w);
  if (!(dens > 0.0)) { retval[0] = 0.0; return 0; }

  // Survivor of each loser: 1 - plba(t - t0_j, ...)
  double out = dens;
  for (int j = 0; j < P->n_acc; ++j) {
    if (j == w) continue;
    const double t_eff_j = t - P->t0[j];
    if (t_eff_j <= 0.0) continue;  // S(t <= t0) = 1, product unaffected
    double denom_j = PNORM_STD(P->v[j] / P->sv[j], true, false);
    if (denom_j < 1e-10) denom_j = 1e-10;
    double cdf = (P->A[j] > LBA_A_ASYMPTOTIC)
      ? plba_core(t_eff_j, P->A[j], P->B[j] + P->A[j], P->v[j], P->sv[j], denom_j)
        : plba_noA (t_eff_j,          P->B[j],             P->v[j], P->sv[j], denom_j);
    if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
    else if (cdf > 1.0)                         cdf = 1.0;
    out *= (1.0 - cdf);
  }
  retval[0] = out;
  return 0;
}

// Integrate f_winner * prod S_losers over [lower, upper] for one trial.
// All parameter arrays are length n_acc, caller-owned.
inline double lba_survivor_scalar(const double* v, const double* sv,
                                  const double* B, const double* A,
                                  const double* t0,
                                  int n_acc, int winner,
                                  double lower, double upper,
                                  double abstol = 1e-8, double reltol = 1e-6,
                                  int maxeval = 6000)
{
  LBAIntegrand ig{ n_acc, winner, v, sv, B, A, t0 };
  double val = 0.0, err = 0.0;
  hcubature(int_lba_survivor, static_cast<void*>(&ig), 1, &lower, &upper,
            maxeval, abstol, reltol, &val, &err);
  if (!std::isfinite(val) || val < 0.0) return 0.0;
  if (val > 1.0)                         return 1.0;
  return val;
}

// -----------------------------------------------------------------------------
// lba_survivor_with_response()
//
// For each entry j in idx:
//   - idx[j]    = base row (t * n_acc) into ParamTable
//   - winner[j] = within-trial winner index [0, n_acc-1]
//   - lower[j]  = lower integration bound for this trial
//   - upper[j]  = upper integration bound for this trial
//   - out[j]    = integral result, one value per trial
//
// Scratch buffers allocated once (length n_acc), reused across trials.
// denom is computed inside the integrand callback to avoid storing it.
// -----------------------------------------------------------------------------
void lba_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out)
{
  std::vector<double> v_buf(n_acc), sv_buf(n_acc), B_buf(n_acc),
  A_buf(n_acc), t0_buf(n_acc);

  const int n = (int)idx.size();
  for (int j = 0; j < n; ++j) {
    const int base = idx[j];  // = t * n_acc

    for (int acc = 0; acc < n_acc; ++acc) {
      const int row = base + acc;
      v_buf[acc]    = pt.base(row, spec.col_v);
      sv_buf[acc]   = pt.base(row, spec.col_sv);
      B_buf[acc]    = pt.base(row, spec.col_B);
      A_buf[acc]    = pt.base(row, spec.col_A);
      t0_buf[acc]   = pt.base(row, spec.col_t0);
    }

    out[j] = lba_survivor_scalar(v_buf.data(), sv_buf.data(),
                                 B_buf.data(), A_buf.data(), t0_buf.data(),
                                 n_acc, winner[j], lower[j], upper[j]);
  }
}


