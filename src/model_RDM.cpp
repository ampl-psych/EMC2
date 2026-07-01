#include "model_RDM.h"

using namespace Rcpp;


// ---------------------------------------------------------------------------
// R-exported scalar functions
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0)
{
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    double t_eff = t[i] - t0[i];
    pdf[i] = (t_eff <= 0.0) ? 0.0 : digt(t_eff, B[i] + 0.5 * A[i], v[i], 0.5 * A[i]);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0)
{
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    double t_eff = t[i] - t0[i];
    cdf[i] = (t_eff <= 0.0) ? 0.0 : pigt(t_eff, B[i] + 0.5 * A[i], v[i], 0.5 * A[i]);
  }
  return cdf;
}


// ---------------------------------------------------------------------------
// Fast ParamTable-based functions
// ---------------------------------------------------------------------------
// void drdm_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const std::vector<int>& idx,
//                double* __restrict__ ll_row)
// {
//   const double* __restrict__ rt = rts.begin();
//   const double* __restrict__ v  = &pt.base(0, spec.col_v);
//   const double* __restrict__ B  = &pt.base(0, spec.col_B);
//   const double* __restrict__ A  = &pt.base(0, spec.col_A);
//   const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
//   const double* __restrict__ s  = &pt.base(0, spec.col_s);
//
//   for (int j = 0; j < (int)idx.size(); ++j) {
//     const int i = idx[j];
//
//     const double t_eff = rt[i] - t0[i];
//     if (t_eff <= 0.0) { ll_row[i] = 0.0; continue; }
//
//     const double inv_s = 1.0 / s[i];
//     double pdf;
//     if (A[i] < A_EPS) {
//       double l = v[i] * inv_s;
//       double k = B[i] * inv_s;
//       clamp_l(l);
//       pdf = digt0(t_eff, k, l);
//     } else {
//       double a = 0.5 * A[i] * inv_s;
//       double l = v[i] * inv_s;
//       double k = B[i] * inv_s + a;
//       clamp_a_l(a, l);
//       pdf = digt_core(t_eff, k, l, a);
//     }
//     ll_row[i] = (std::isfinite(pdf) && pdf >= 0.0) ? pdf : 0.0;
//   }
// }


// void prdm_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const std::vector<int>& idx,
//                double* __restrict__ ll_row)
// {
//   const double* __restrict__ rt = rts.begin();
//   const double* __restrict__ v  = &pt.base(0, spec.col_v);
//   const double* __restrict__ B  = &pt.base(0, spec.col_B);
//   const double* __restrict__ A  = &pt.base(0, spec.col_A);
//   const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
//   const double* __restrict__ s  = &pt.base(0, spec.col_s);
//
//   for (int j = 0; j < (int)idx.size(); ++j) {
//     const int i = idx[j];
//
//     const double t_eff = rt[i] - t0[i];
//     if (t_eff <= 0.0) { ll_row[i] = 1.0; continue; }  // survival = 1
//
//     const double inv_s = 1.0 / s[i];
//     double cdf;
//     if (A[i] < A_EPS) {
//       double l = v[i] * inv_s;
//       double k = B[i] * inv_s;
//       clamp_l(l);
//       cdf = pigt0(t_eff, k, l);
//     } else {
//       double a = 0.5 * A[i] * inv_s;
//       double l = v[i] * inv_s;
//       double k = B[i] * inv_s + a;
//       clamp_a_l(a, l);
//       cdf = pigt_core(t_eff, k, l, a);
//     }
//     if (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
//     else if (cdf > 1.0)                    cdf = 1.0;
//     ll_row[i] = 1.0 - cdf;  // survival
//   }
// }


// This new filling function checks whether A==0, if so --> runs digt0 and pigt0
void drdm_prdm_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch)
{
  // Note that for the losers, the *SURVIVAL* probability is filled, *NOT* the CDF
  const double* __restrict__ rt = rts.begin();
  const double* __restrict__ v  = &pt.base(0, spec.col_v);
  const double* __restrict__ B  = &pt.base(0, spec.col_B);
  const double* __restrict__ A  = &pt.base(0, spec.col_A);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // noA scratch (primary)
  double* __restrict__ sc_teff = scratch.t_eff.data();
  double* __restrict__ sc_v    = scratch.v.data();
  double* __restrict__ sc_B    = scratch.B.data();
  double* __restrict__ sc_s    = scratch.s.data();
  double* __restrict__ sc_out  = scratch.out.data();
  int*    __restrict__ sc_idx  = scratch.idx_win0.data();

  // core scratch (secondary)
  double* __restrict__ sc_teff_c = scratch.t_eff_c.data();
  double* __restrict__ sc_v_c    = scratch.v_c.data();
  double* __restrict__ sc_B_c    = scratch.B_c.data();
  double* __restrict__ sc_A_c    = scratch.A_c.data();
  double* __restrict__ sc_s_c    = scratch.s_c.data();
  double* __restrict__ sc_out_c  = scratch.out_c.data();
  int*    __restrict__ sc_idx_c  = scratch.idx_win_c.data();

  // =========================================================================
  // WINNERS
  // =========================================================================

  // --- One-pass gather: split into noA (primary, most likely) and core (secondary, most likely) ---
  int n_win_noA = 0, n_win_core = 0;
  for (int j = 0; j < n_win; ++j) {
    const int i       = idx_win[j];

    // first guard against teff < 0 (non-decision time can't be > rt)
    const double teff = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = 0.0; continue; }

    // Check whether A equals 0
    if (A[i] < A_EPS) {
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { ll_row[i] = 0.0; continue; }  // degenerate
      sc_teff[n_win_noA] = teff;
      sc_v   [n_win_noA] = v[i];
      sc_B   [n_win_noA] = B[i];
      sc_s   [n_win_noA] = s[i];
      sc_idx [n_win_noA] = i;
      n_win_noA++;
    } else {
      sc_teff_c[n_win_core] = teff;
      sc_v_c   [n_win_core] = v[i];
      sc_B_c   [n_win_core] = B[i];
      sc_A_c   [n_win_core] = A[i];
      sc_s_c   [n_win_core] = s[i];
      sc_idx_c [n_win_core] = i;
      n_win_core++;
    }
  }

  // --- Compute: digt0 winners ---
#pragma omp simd
  for (int j = 0; j < n_win_noA; ++j) {
    const double inv_s = 1.0 / sc_s[j];
    double l = sc_v[j] * inv_s;
    double k = sc_B[j] * inv_s;
    clamp_l(l);
    sc_out[j] = digt0(sc_teff[j], k, l);
  }

  // --- Compute: digt_core winners ---
#pragma omp simd
  for (int j = 0; j < n_win_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l); // clamp a and l
    sc_out_c[j] = digt_core(sc_teff_c[j], k, l, a);
  }

  // --- Scatter: noA winners ---
  for (int j = 0; j < n_win_noA; ++j) {
    const double val = sc_out[j];
    ll_row[sc_idx[j]] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;    // guard against pdf < 0
  }

  // --- Scatter: core winners ---
  for (int j = 0; j < n_win_core; ++j) {
    const double val = sc_out_c[j];
    ll_row[sc_idx_c[j]] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;  // guard against pdf < 0
  }

  // =========================================================================
  // LOSERS
  // =========================================================================

  sc_idx  = scratch.idx_los0.data();
  sc_idx_c = scratch.idx_los_c.data();

  // --- One-pass gather: split into noA (primary, most likely) and core (secondary, most likely) ---
  int n_los_noA = 0, n_los_core = 0;
  for (int j = 0; j < n_los; ++j) {
    const int i       = idx_los[j];

    // first guard against teff < 0 (non-decision time can't be > rt)
    const double teff = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = 1.0; continue; }

    // fill scratch depending on whether we have A tiny or not
    if (A[i] < A_EPS) {
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { ll_row[i] = 1.0; continue; }  // survival = 1
      sc_teff[n_los_noA] = teff;
      sc_v   [n_los_noA] = v[i];
      sc_B   [n_los_noA] = B[i];
      sc_s   [n_los_noA] = s[i];
      sc_idx [n_los_noA] = i;
      n_los_noA++;
    } else {
      sc_teff_c[n_los_core] = teff;
      sc_v_c   [n_los_core] = v[i];
      sc_B_c   [n_los_core] = B[i];
      sc_A_c   [n_los_core] = A[i];
      sc_s_c   [n_los_core] = s[i];
      sc_idx_c [n_los_core] = i;
      n_los_core++;
    }
  }

  // --- Compute: pigt0 losers ---
#pragma omp simd
  for (int j = 0; j < n_los_noA; ++j) {
    const double inv_s = 1.0 / sc_s[j];
    double l = sc_v[j] * inv_s;
    double k = sc_B[j] * inv_s;
    clamp_l(l);
    sc_out[j] = pigt0(sc_teff[j], k, l);
  }

  // --- Compute: pigt_core losers ---
#pragma omp simd
  for (int j = 0; j < n_los_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l);
    sc_out_c[j] = pigt_core(sc_teff_c[j], k, l, a);
  }

  // --- Scatter: noA losers ---
  for (int j = 0; j < n_los_noA; ++j) {
    double val = sc_out[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;    // guard against cdf < 0
    else if (val > 1.0)                    val = 1.0;   // guard against cdf > 1
    ll_row[sc_idx[j]] = 1.0 - val;                      // fill in survival (not CDF)
  }

  // --- Scatter: core losers ---
  for (int j = 0; j < n_los_core; ++j) {
    double val = sc_out_c[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;    // guard against cdf < 0
    else if (val > 1.0)                    val = 1.0;   // guard against cdf > 1
    ll_row[sc_idx_c[j]] = 1.0 - val;                    // fill in survival (not CDF)
  }
}

// Survivor function. Gather-scatter will probably pay off here since requires computation on all trials
// rdm_fill_truncate: fills raw survivors S_k(bound - t0) for the given idx set
// Call twice at the call site — once with (idx_LT+idx_both, LT, ll_lower)
//                             — once with (idx_UT+idx_both, UT, ll_upper)
// Unvisited rows retain their std::fill default (1.0 for lower, 0.0 for upper)
// =============================================================================
void rdm_survivor(const std::vector<int>& idx,
                  const std::vector<double>& bound,
                  const ParamTable& pt,
                  const RaceSpec& spec,
                  double* __restrict__ out,
                  RaceScratch& scratch)
{
  const int n = (int)idx.size();

  const double* __restrict__ v  = &pt.base(0, spec.col_v);
  const double* __restrict__ B  = &pt.base(0, spec.col_B);
  const double* __restrict__ A  = &pt.base(0, spec.col_A);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);

  // noA scratch (primary)
  double* __restrict__ sc_teff  = scratch.t_eff.data();
  double* __restrict__ sc_v     = scratch.v.data();
  double* __restrict__ sc_B     = scratch.B.data();
  double* __restrict__ sc_s     = scratch.s.data();
  double* __restrict__ sc_out   = scratch.out.data();
  int*    __restrict__ sc_idx   = scratch.idx_win0.data();

  // core scratch (secondary)
  double* __restrict__ sc_teff_c = scratch.t_eff_c.data();
  double* __restrict__ sc_v_c    = scratch.v_c.data();
  double* __restrict__ sc_B_c    = scratch.B_c.data();
  double* __restrict__ sc_A_c    = scratch.A_c.data();
  double* __restrict__ sc_s_c    = scratch.s_c.data();
  double* __restrict__ sc_out_c  = scratch.out_c.data();
  int*    __restrict__ sc_idx_c  = scratch.idx_win_c.data();

  // --- One-pass gather: split noA (primary) / core (secondary) ---
  int n_noA = 0, n_core = 0;
  for (int j = 0; j < n; ++j) {
    const int i        = idx[j];
    const double teff  = bound[i] - t0[i];
    //if (teff <= 0.0) { out[i] = 1.0; continue; }  // S(t<=0) = 1, skip
    if (teff <= 0.0) { continue; }  // S(t<=0) = 1, skip but leave out untouched - reverts to default value

    if (A[i] < A_EPS) {
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { out[i] = 1.0; continue; }  // degenerate
      sc_teff[n_noA] = teff;
      sc_v   [n_noA] = v[i];
      sc_B   [n_noA] = B[i];
      sc_s   [n_noA] = s[i];
      sc_idx [n_noA] = i;
      n_noA++;
    } else {
      sc_teff_c[n_core] = teff;
      sc_v_c   [n_core] = v[i];
      sc_B_c   [n_core] = B[i];
      sc_A_c   [n_core] = A[i];
      sc_s_c   [n_core] = s[i];
      sc_idx_c [n_core] = i;
      n_core++;
    }
  }

  // --- Compute: pigt0 noA ---
#pragma omp simd
  for (int j = 0; j < n_noA; ++j) {
    const double inv_s = 1.0 / sc_s[j];
    double l = sc_v[j] * inv_s;
    double k = sc_B[j] * inv_s;
    clamp_l(l);
    sc_out[j] = pigt0(sc_teff[j], k, l);
  }

  // --- Compute: pigt_core core ---
#pragma omp simd
  for (int j = 0; j < n_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l);
    sc_out_c[j] = pigt_core(sc_teff_c[j], k, l, a);
  }

  // --- Scatter: noA — store survival (1 - CDF) ---
  for (int j = 0; j < n_noA; ++j) {
    double val = sc_out[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                   val = 1.0;
    out[sc_idx[j]] = 1.0 - val;
  }

  // --- Scatter: core — store survival (1 - CDF) ---
  for (int j = 0; j < n_core; ++j) {
    double val = sc_out_c[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                    val = 1.0;
    out[sc_idx_c[j]] = 1.0 - val;
  }
}

// -----------------------------------------------------------------------------
// RDM survivor with known response — numerical integration
//
// Computes P(winner finishes first in [lower, upper]) for one trial by
// integrating f_winner(t) * prod_{j != winner} S_j(t) over [lower, upper].
// -----------------------------------------------------------------------------

// Per-trial parameter block, pointing into caller-owned contiguous arrays.
// All arrays are length n_acc.
struct RDMIntegrand {
  int           n_acc;
  int           winner;   // 0-based within-trial index
  const double* k;        // (B + 0.5*A) / s per accumulator
  const double* l;        // v / s
  const double* a;        // 0.5 * A / s  (0 when A < A_EPS)
  const double* t0;
};

// hcubature callback: x[0] = absolute time t
// retval[0] = f_winner(t - t0_winner) * prod_{j != winner} S_j(t - t0_j)
inline int int_rdm_survivor(unsigned /*dim*/, const double* x, void* p,
                            unsigned /*fdim*/, double* retval)
{
  const RDMIntegrand* P = static_cast<const RDMIntegrand*>(p);
  const double t       = x[0];
  const int    w       = P->winner;
  const double t_eff_w = t - P->t0[w];

  // Density of winner — early exit if zero (t <= t0 or degenerate)
  const double dens = (P->a[w] < A_EPS)
    ? digt0(t_eff_w, P->k[w], P->l[w])
      : digt_core(t_eff_w, P->k[w], P->l[w], P->a[w]);
  if (!(dens > 0.0)) { retval[0] = 0.0; return 0; }

  // Survivor of each loser
  double out = dens;
  for (int j = 0; j < P->n_acc; ++j) {
    if (j == w) continue;
    const double t_eff_j = t - P->t0[j];
    const double cdf = (P->a[j] < A_EPS)
      ? pigt0(t_eff_j, P->k[j], P->l[j])
        : pigt_core(t_eff_j, P->k[j], P->l[j], P->a[j]);
    out *= (1.0 - cdf);
  }
  retval[0] = out;
  return 0;
}

// Integrate f_winner * prod S_losers over [lower, upper] for one trial.
// k/l/a/t0 are pre-scaled, length n_acc, caller-owned — no allocation here.
inline double rdm_survivor_scalar(const double* k, const double* l,
                                  const double* a, const double* t0,
                                  int n_acc, int winner,
                                  double lower, double upper,
                                  double abstol = 1e-8, double reltol = 1e-6,
                                  int maxeval = 6000)
{
  RDMIntegrand ig{ n_acc, winner, k, l, a, t0 };
  double val = 0.0, err = 0.0;
  hcubature(int_rdm_survivor, static_cast<void*>(&ig), 1, &lower, &upper,
            maxeval, abstol, reltol, &val, &err);
  if (!std::isfinite(val) || val < 0.0) return 0.0;
  if (val > 1.0)                         return 1.0;
  return val;
}

// -----------------------------------------------------------------------------
// rdm_survivor_with_response()
//
// For each entry j in idx:
//   - idx[j]    = base row (t * n_acc) into ParamTable
//   - winner[j] = within-trial winner index [0, n_acc-1]
//   - lower[j]  = lower integration bound for this trial
//   - upper[j]  = upper integration bound for this trial
//   - out[j]    = integral result, one value per trial
//
// Scratch buffers are allocated once (length n_acc) and reused across trials.
// -----------------------------------------------------------------------------
void rdm_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out)
{
  std::vector<double> k_buf(n_acc), l_buf(n_acc), a_buf(n_acc), t0_buf(n_acc);

  const int n = (int)idx.size();
  for (int j = 0; j < n; ++j) {
    const int base = idx[j];  // = t * n_acc

    // Gather and pre-scale parameters for all accumulators in this trial
    for (int acc = 0; acc < n_acc; ++acc) {
      const int    row   = base + acc;
      const double inv_s = 1.0 / pt.base(row, spec.col_s);
      const double A_acc = pt.base(row, spec.col_A);
      l_buf[acc]  = pt.base(row, spec.col_v) * inv_s;
      a_buf[acc]  = (A_acc < A_EPS) ? 0.0 : 0.5 * A_acc * inv_s;
      k_buf[acc]  = pt.base(row, spec.col_B) * inv_s + a_buf[acc];
      t0_buf[acc] = pt.base(row, spec.col_t0);
    }

    out[j] = rdm_survivor_scalar(k_buf.data(), l_buf.data(),
                                 a_buf.data(), t0_buf.data(),
                                 n_acc, winner[j], lower[j], upper[j]);
  }
}
