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


// Censoring
// =============================================================================
// rdm_censor — gather-scatter, noA/core split
// Three passes: idx_L, idx_U, idx_B.
//
//   idx_L: ll_row[i] = CDF(LC - t0)
//   idx_U: ll_row[i] = 1 - CDF(UC - t0)
//   idx_B: ll_row[i] = CDF(LC - t0) + 1 - CDF(UC - t0)
//
// Scaled params (k, l, a) computed in gather — compute loop is branch-free.
// noA path uses primary scratch; core path uses secondary scratch.
// =============================================================================
void rdm_censor(const CensorSpec&    censor,
                const ParamTable&    pt,
                const RaceSpec&      spec,
                double* __restrict__ ll_row,
                RaceScratch&         scratch)
{
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);
  const double* __restrict__ s   = &pt.base(0, spec.col_s);
  const double* __restrict__ LC  = censor.LC.data();
  const double* __restrict__ UC  = censor.UC.data();

  // noA scratch (primary) — scaled params: k = B/s, l = v/s
  double* __restrict__ sc_teff  = scratch.t_eff.data();
  double* __restrict__ sc_l     = scratch.v.data();    // l = v/s
  double* __restrict__ sc_k     = scratch.B.data();    // k = B/s
  double* __restrict__ sc_out   = scratch.out.data();
  int*    __restrict__ sc_idx   = scratch.idx_win0.data();

  // core scratch (secondary) — scaled params: k = B/s + a, l = v/s, a = 0.5*A/s
  double* __restrict__ sc_teff_c = scratch.t_eff_c.data();
  double* __restrict__ sc_l_c    = scratch.v_c.data();   // l = v/s
  double* __restrict__ sc_k_c    = scratch.B_c.data();   // k = B/s + a
  double* __restrict__ sc_a_c    = scratch.A_c.data();   // a = 0.5*A/s
  double* __restrict__ sc_out_c  = scratch.out_c.data();
  int*    __restrict__ sc_idx_c  = scratch.idx_win_c.data();

  // ---- Pass 1: lower-censored — ll_row[i] = CDF(LC - t0) ----
  {
    int n_noA = 0, n_core = 0;
    for (int j = 0; j < (int)censor.idx_L.size(); ++j) {
      const int i       = censor.idx_L[j];
      const double teff = LC[i] - t0[i];
      if (teff <= 0.0) { ll_row[i] = 0.0; continue; }
      const double inv_s = 1.0 / s[i];
      if (A[i] < A_EPS) {
        double l = v[i] * inv_s;  double k = B[i] * inv_s;
        if (k > K_MAX || k < 0.0) { ll_row[i] = 0.0; continue; }
        clamp_l(l);
        sc_teff[n_noA] = teff;  sc_l[n_noA] = l;  sc_k[n_noA] = k;
        sc_idx [n_noA] = i;     n_noA++;
      } else {
        double l = v[i] * inv_s;  double a = 0.5 * A[i] * inv_s;  double k = B[i] * inv_s + a;
        clamp_a_l(a, l);
        sc_teff_c[n_core] = teff;  sc_l_c[n_core] = l;  sc_k_c[n_core] = k;  sc_a_c[n_core] = a;
        sc_idx_c [n_core] = i;     n_core++;
      }
    }

    // compute
#pragma omp simd
    for (int j = 0; j < n_noA;  ++j) sc_out  [j] = pigt0     (sc_teff  [j], sc_k  [j], sc_l  [j]);
#pragma omp simd
    for (int j = 0; j < n_core; ++j) sc_out_c[j] = pigt_core (sc_teff_c[j], sc_k_c[j], sc_l_c[j], sc_a_c[j]);

    // scatter
    for (int j = 0; j < n_noA; ++j) {
      double val = sc_out[j];
      if (!std::isfinite(val) || val < 0.0) val = 0.0; else if (val > 1.0) val = 1.0;
      ll_row[sc_idx[j]] = val;
    }
    for (int j = 0; j < n_core; ++j) {
      double val = sc_out_c[j];
      if (!std::isfinite(val) || val < 0.0) val = 0.0; else if (val > 1.0) val = 1.0;
      ll_row[sc_idx_c[j]] = val;
    }
  }

  // ---- Pass 2: upper-censored — ll_row[i] = 1 - CDF(UC - t0) ----
  {
    int n_noA = 0, n_core = 0;
    for (int j = 0; j < (int)censor.idx_U.size(); ++j) {
      const int i       = censor.idx_U[j];
      const double teff = UC[i] - t0[i];
      if (teff <= 0.0) { ll_row[i] = 1.0; continue; }
      const double inv_s = 1.0 / s[i];
      if (A[i] < A_EPS) {
        double l = v[i] * inv_s;  double k = B[i] * inv_s;
        if (k > K_MAX || k < 0.0) { ll_row[i] = 1.0; continue; }
        clamp_l(l);
        sc_teff[n_noA] = teff;  sc_l[n_noA] = l;  sc_k[n_noA] = k;
        sc_idx [n_noA] = i;     n_noA++;
      } else {
        double l = v[i] * inv_s;  double a = 0.5 * A[i] * inv_s;  double k = B[i] * inv_s + a;
        clamp_a_l(a, l);
        sc_teff_c[n_core] = teff;  sc_l_c[n_core] = l;  sc_k_c[n_core] = k;  sc_a_c[n_core] = a;
        sc_idx_c [n_core] = i;     n_core++;
      }
    }

    // compute
#pragma omp simd
    for (int j = 0; j < n_noA;  ++j) sc_out  [j] = pigt0     (sc_teff  [j], sc_k  [j], sc_l  [j]);
#pragma omp simd
    for (int j = 0; j < n_core; ++j) sc_out_c[j] = pigt_core (sc_teff_c[j], sc_k_c[j], sc_l_c[j], sc_a_c[j]);

    // scatter
    for (int j = 0; j < n_noA; ++j) {
      double val = sc_out[j];
      if (!std::isfinite(val) || val < 0.0) val = 0.0; else if (val > 1.0) val = 1.0;
      ll_row[sc_idx[j]] = 1.0 - val;
    }
    for (int j = 0; j < n_core; ++j) {
      double val = sc_out_c[j];
      if (!std::isfinite(val) || val < 0.0) val = 0.0; else if (val > 1.0) val = 1.0;
      ll_row[sc_idx_c[j]] = 1.0 - val;
    }
  }

  // ---- Pass 3: both-censored — ll_row[i] = CDF(LC - t0) + 1 - CDF(UC - t0) ----
  // Single gather packs both LC and UC teff; scaled params computed once per row.
  {
    double* __restrict__ sc_teff_UC  = scratch.k.data();   // UC offsets, noA
    double* __restrict__ sc_teff_UC_c = scratch.l.data();  // UC offsets, core

    int n_noA = 0, n_core = 0;
    for (int j = 0; j < (int)censor.idx_B.size(); ++j) {
      const int i        = censor.idx_B[j];
      const double inv_s = 1.0 / s[i];
      if (A[i] < A_EPS) {
        double l = v[i] * inv_s;  double k = B[i] * inv_s;
        if (k > K_MAX || k < 0.0) { ll_row[i] = 0.0; continue; }
        clamp_l(l);
        sc_teff   [n_noA] = LC[i] - t0[i];  sc_teff_UC[n_noA] = UC[i] - t0[i];
        sc_l[n_noA] = l;  sc_k[n_noA] = k;  sc_idx[n_noA] = i;
        n_noA++;
      } else {
        double l = v[i] * inv_s;  double a = 0.5 * A[i] * inv_s;  double k = B[i] * inv_s + a;
        clamp_a_l(a, l);
        sc_teff_c   [n_core] = LC[i] - t0[i];  sc_teff_UC_c[n_core] = UC[i] - t0[i];
        sc_l_c[n_core] = l;  sc_k_c[n_core] = k;  sc_a_c[n_core] = a;  sc_idx_c[n_core] = i;
        n_core++;
      }
    }
    for (int j = 0; j < n_noA; ++j) {
      double cdf_LC = (sc_teff   [j] > 0.0) ? pigt0(sc_teff   [j], sc_k[j], sc_l[j]) : 0.0;
      double cdf_UC = (sc_teff_UC[j] > 0.0) ? pigt0(sc_teff_UC[j], sc_k[j], sc_l[j]) : 0.0;
      if (!std::isfinite(cdf_LC) || cdf_LC < 0.0) cdf_LC = 0.0; else if (cdf_LC > 1.0) cdf_LC = 1.0;
      if (!std::isfinite(cdf_UC) || cdf_UC < 0.0) cdf_UC = 0.0; else if (cdf_UC > 1.0) cdf_UC = 1.0;
      ll_row[sc_idx[j]] = cdf_LC + 1.0 - cdf_UC;
    }
    for (int j = 0; j < n_core; ++j) {
      double cdf_LC = (sc_teff_c   [j] > 0.0) ? pigt_core(sc_teff_c   [j], sc_k_c[j], sc_l_c[j], sc_a_c[j]) : 0.0;
      double cdf_UC = (sc_teff_UC_c[j] > 0.0) ? pigt_core(sc_teff_UC_c[j], sc_k_c[j], sc_l_c[j], sc_a_c[j]) : 0.0;
      if (!std::isfinite(cdf_LC) || cdf_LC < 0.0) cdf_LC = 0.0; else if (cdf_LC > 1.0) cdf_LC = 1.0;
      if (!std::isfinite(cdf_UC) || cdf_UC < 0.0) cdf_UC = 0.0; else if (cdf_UC > 1.0) cdf_UC = 1.0;
      ll_row[sc_idx_c[j]] = cdf_LC + 1.0 - cdf_UC;
    }
  }
}



// Truncation. Gather-scatter will probably pay off here since requires computation on all trials
// rdm_fill_truncate: fills raw survivors S_k(bound - t0) for the given idx set
// Call twice at the call site — once with (idx_LT+idx_both, LT, ll_lower)
//                             — once with (idx_UT+idx_both, UT, ll_upper)
// Unvisited rows retain their std::fill default (1.0 for lower, 0.0 for upper)
// =============================================================================
void rdm_truncate(const std::vector<int>& idx,
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
