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

// This new filling function checks whether A==0, if so --> runs digt0 and pigt0
void drdm_prdm_fast(const double* rt,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch,
                    const std::vector<int>& is_ok, bool log_p)
{
  // Note that for the losers, the *SURVIVAL* probability is filled, *NOT* the CDF
  const double* __restrict__ v  = &pt.base(0, spec.col_v);
  const double* __restrict__ B  = &pt.base(0, spec.col_B);
  const double* __restrict__ A  = &pt.base(0, spec.col_A);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // noA scratch (primary)
  double* __restrict__ sc_teff = scratch.dslot(RDM::t_eff);
  double* __restrict__ sc_v    = scratch.dslot(RDM::v);
  double* __restrict__ sc_B    = scratch.dslot(RDM::B);
  double* __restrict__ sc_s    = scratch.dslot(RDM::s);
  double* __restrict__ sc_out  = scratch.dslot(RDM::out);
  int*    __restrict__ sc_idx  = scratch.islot(ScratchInt::idx_win);

  // core scratch (secondary)
  double* __restrict__ sc_teff_c = scratch.dslot(RDM::t_eff_c);
  double* __restrict__ sc_v_c    = scratch.dslot(RDM::v_c);
  double* __restrict__ sc_B_c    = scratch.dslot(RDM::B_c);
  double* __restrict__ sc_A_c    = scratch.dslot(RDM::A_c);
  double* __restrict__ sc_s_c    = scratch.dslot(RDM::s_c);
  double* __restrict__ sc_out_c  = scratch.dslot(RDM::out_c);
  int*    __restrict__ sc_idx_c  = scratch.islot(ScratchInt::idx_win_c);

  // =========================================================================
  // WINNERS
  // =========================================================================

  // --- One-pass gather: split into noA (primary, most likely) and core (secondary, most likely) ---
  int n_win_noA = 0, n_win_core = 0;
  for (int j = 0; j < n_win; ++j) {
    const int i       = idx_win[j];
    if (!is_ok[i]) continue;

    // first guard against teff < 0 (non-decision time can't be > rt)
    const double teff = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = log_p ? neg_inf() : 0.0; continue; }

    // Check whether A equals 0
    if (A[i] < A_EPS) {
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { ll_row[i] = log_p ? neg_inf() : 0.0; continue; }  // degenerate
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
    sc_out[j] = digt0(sc_teff[j], k, l, log_p);
  }

  // --- Compute: digt_core winners ---
#pragma omp simd
  for (int j = 0; j < n_win_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l); // clamp a and l
    sc_out_c[j] = digt_core(sc_teff_c[j], k, l, a, log_p);
  }

  // --- Scatter: noA winners ---
  for (int j = 0; j < n_win_noA; ++j) {
    const double val = sc_out[j];
    ll_row[sc_idx[j]] = log_p ? val : ((std::isfinite(val) && val >= 0.0) ? val : 0.0);    // guard against pdf < 0
  }

  // --- Scatter: core winners ---
  for (int j = 0; j < n_win_core; ++j) {
    const double val = sc_out_c[j];
    ll_row[sc_idx_c[j]] = log_p ? val : ((std::isfinite(val) && val >= 0.0) ? val : 0.0);  // guard against pdf < 0
  }

  // =========================================================================
  // LOSERS
  // =========================================================================
  sc_idx  = scratch.islot(ScratchInt::idx_los);
  sc_idx_c = scratch.islot(ScratchInt::idx_los_c);

  // --- One-pass gather: split into noA (primary, most likely) and core (secondary, most likely) ---
  int n_los_noA = 0, n_los_core = 0;
  for (int j = 0; j < n_los; ++j) {
    const int i       = idx_los[j];
    if (!is_ok[i]) continue;

    // first guard against teff < 0 (non-decision time can't be > rt)
    const double teff = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = log_p ? 0.0 : 1.0; continue; }

    // fill scratch depending on whether we have A tiny or not
    if (A[i] < A_EPS) {
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { ll_row[i] = log_p ? 0.0 : 1.0; continue; }  // survival = 1
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
    sc_out[j] = pigt0(sc_teff[j], k, l, false, log_p);
  }

  // --- Compute: pigt_core losers ---
#pragma omp simd
  for (int j = 0; j < n_los_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l);
    sc_out_c[j] = pigt_core(sc_teff_c[j], k, l, a, false, log_p);
  }

  // --- Scatter: noA losers ---
  for (int j = 0; j < n_los_noA; ++j) {
    double val = sc_out[j];
    if (!log_p && (!std::isfinite(val) || val < 0.0)) val = 0.0;    // guard against survival < 0
    else if (!log_p && val > 1.0)                    val = 1.0;   // guard against survival > 1
    ll_row[sc_idx[j]] = val;                      // direct survival
  }

  // --- Scatter: core losers ---
  for (int j = 0; j < n_los_core; ++j) {
    double val = sc_out_c[j];
    if (!log_p && (!std::isfinite(val) || val < 0.0)) val = 0.0;    // guard against survival < 0
    else if (!log_p && val > 1.0)                    val = 1.0;   // guard against survival > 1
    ll_row[sc_idx_c[j]] = val;                    // direct survival
  }
}

// Survivor function. Gather-scatter will probably pay off here since requires computation on all trials
// rdm_fill_truncate: fills raw survivors S_k(bound - t0) for the given idx set
// Call twice at the call site — once with (idx_LT+idx_both, LT, ll_lower)
//                             — once with (idx_UT+idx_both, UT, ll_upper)
// Unvisited rows retain the caller's default; active pre-start bounds have S=1.
// =============================================================================
void rdm_survivor(const std::vector<int>& idx,
                  const std::vector<double>& bound,
                  const ParamTable& pt,
                  const RaceSpec& spec,
                  double* __restrict__ out,
                  RaceScratch& scratch,
                  const std::vector<int>& is_ok, bool log_p)
{
  const int n = (int)idx.size();

  const double* __restrict__ v  = &pt.base(0, spec.col_v);
  const double* __restrict__ B  = &pt.base(0, spec.col_B);
  const double* __restrict__ A  = &pt.base(0, spec.col_A);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);

  // noA scratch (primary)
  double* __restrict__ sc_teff = scratch.dslot(RDM::t_eff);
  double* __restrict__ sc_v    = scratch.dslot(RDM::v);
  double* __restrict__ sc_B    = scratch.dslot(RDM::B);
  double* __restrict__ sc_s    = scratch.dslot(RDM::s);
  double* __restrict__ sc_out  = scratch.dslot(RDM::out);
  int*    __restrict__ sc_idx  = scratch.islot(ScratchInt::idx_win);

  // core scratch (secondary)
  double* __restrict__ sc_teff_c = scratch.dslot(RDM::t_eff_c);
  double* __restrict__ sc_v_c    = scratch.dslot(RDM::v_c);
  double* __restrict__ sc_B_c    = scratch.dslot(RDM::B_c);
  double* __restrict__ sc_A_c    = scratch.dslot(RDM::A_c);
  double* __restrict__ sc_s_c    = scratch.dslot(RDM::s_c);
  double* __restrict__ sc_out_c  = scratch.dslot(RDM::out_c);
  int*    __restrict__ sc_idx_c  = scratch.islot(ScratchInt::idx_win_c);

  // --- One-pass gather: split noA (primary) / core (secondary) ---
  int n_noA = 0, n_core = 0;
  for (int j = 0; j < n; ++j) {
    const int i        = idx[j];
    if (!is_ok[i]) continue;
    const double teff  = bound[i] - t0[i];
    if (teff <= 0.0) { out[i] = log_p ? 0.0 : 1.0; continue; }  // S(t<=0) = 1.

    if (A[i] < A_EPS) {
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { out[i] = log_p ? 0.0 : 1.0; continue; }  // degenerate
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
    sc_out[j] = pigt0(sc_teff[j], k, l, false, log_p);
  }

  // --- Compute: pigt_core core ---
#pragma omp simd
  for (int j = 0; j < n_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l);
    sc_out_c[j] = pigt_core(sc_teff_c[j], k, l, a, false, log_p);
  }

  for (int j = 0; j < n_noA; ++j) out[sc_idx[j]] = sc_out[j];
  for (int j = 0; j < n_core; ++j) out[sc_idx_c[j]] = sc_out_c[j];
}
