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
    if (t[i]<=0.0 || is_pos_inf(t[i])) { pdf[i]=0.0; continue; }
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
    if (t[i]<=0.0) { cdf[i]=0.0; continue; }
    double denom = 1.0;
    if (posdrift) {
      denom = PNORM_STD(v[i] / sv[i], /*lower=*/true, /*logp=*/false);
      if (denom < 1e-10) denom = 1e-10;
    }
    double val = (A[i] > LBA_A_ASYMPTOTIC)
      ? plba_core(t[i], A[i], b[i], v[i], sv[i], denom, true, false, posdrift)
        : plba_noA (t[i],       b[i], v[i], sv[i], denom, true, false, posdrift);
    if      (val < 0.0) val = 0.0;
    else if (val > 1.0) val = 1.0;
    cdf[i] = val;
  }
  return cdf;
}

// ---------------------------------------------------------------------------
// Hot path: one-pass gather split on A, core dominates -> primary scratch
// ---------------------------------------------------------------------------
void dlba_plba_fast(const double* rt,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch,
                    const std::vector<int>& is_ok, bool log_p)
{
  // Note that for the losers, the *SURVIVAL* probability is filled, *NOT* the CDF
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // primary scratch — core entries (A > LBA_A_ASYMPTOTIC)
  double* __restrict__ sc_teff  = scratch.dslot(LBA::t_eff);
  double* __restrict__ sc_v     = scratch.dslot(LBA::v);
  double* __restrict__ sc_sv    = scratch.dslot(LBA::sv);
  double* __restrict__ sc_B     = scratch.dslot(LBA::B);
  double* __restrict__ sc_A     = scratch.dslot(LBA::A);
  double* __restrict__ sc_out   = scratch.dslot(LBA::out);
  int*    __restrict__ sc_idx   = scratch.islot(ScratchInt::idx_win);

  // secondary scratch — noA entries (A <= LBA_A_ASYMPTOTIC)
  double* __restrict__ sc_teff_c = scratch.dslot(LBA::t_eff_c);
  double* __restrict__ sc_v_c    = scratch.dslot(LBA::v_c);
  double* __restrict__ sc_sv_c   = scratch.dslot(LBA::sv_c);
  double* __restrict__ sc_B_c    = scratch.dslot(LBA::B_c);
  double* __restrict__ sc_out_c  = scratch.dslot(LBA::out_c);
  int*    __restrict__ sc_idx_c  = scratch.islot(ScratchInt::idx_win_c);

  // =========================================================================
  // WINNERS
  // =========================================================================

  // --- One-pass gather: core -> primary, noA -> secondary ---
  int n_win_core = 0, n_win_noA = 0;
  for (int j = 0; j < n_win; ++j) {
    const int i       = idx_win[j];
    if (!is_ok[i]) continue;
    const double teff = rt[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { ll_row[i] = log_p ? neg_inf() : 0.0; continue; }
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
                           sc_v[j], sc_sv[j], denom, log_p);
    sc_out[j] = log_p ? val : ((std::isfinite(val) && val >= 0.0) ? val : 0.0);
  }

  // --- Compute: dlba_noA winners ---
  for (int j = 0; j < n_win_noA; ++j) {
    double denom = PNORM_STD(sc_v_c[j] / sc_sv_c[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = dlba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], denom, log_p);
    sc_out_c[j] = log_p ? val : ((std::isfinite(val) && val >= 0.0) ? val : 0.0);
  }

  // --- Scatter ---
  for (int j = 0; j < n_win_core; ++j) ll_row[sc_idx  [j]] = sc_out  [j];
  for (int j = 0; j < n_win_noA;  ++j) ll_row[sc_idx_c[j]] = sc_out_c[j];

  // =========================================================================
  // LOSERS
  // =========================================================================

  sc_idx   = scratch.islot(ScratchInt::idx_los);
  sc_idx_c = scratch.islot(ScratchInt::idx_los_c);

  // --- One-pass gather: core -> primary, noA -> secondary ---
  int n_los_core = 0, n_los_noA = 0;
  for (int j = 0; j < n_los; ++j) {
    const int i       = idx_los[j];
    if (!is_ok[i]) continue;
    const double teff = rt[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { ll_row[i] = log_p ? 0.0 : 1.0; continue; }
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
                           sc_v[j], sc_sv[j], denom, false, log_p);
    if      (!log_p && (!std::isfinite(val) || val < 0.0)) val = 0.0;
    else if (!log_p && val > 1.0)                         val = 1.0;
    sc_out[j] = val;
  }

  // --- Compute: plba_noA losers ---
  for (int j = 0; j < n_los_noA; ++j) {
    double denom = PNORM_STD(sc_v_c[j] / sc_sv_c[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = plba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], denom, false, log_p);
    if      (!log_p && (!std::isfinite(val) || val < 0.0)) val = 0.0;
    else if (!log_p && val > 1.0)                         val = 1.0;
    sc_out_c[j] = val;
  }

  // --- Scatter ---
  for (int j = 0; j < n_los_core; ++j) ll_row[sc_idx  [j]] = sc_out  [j];
  for (int j = 0; j < n_los_noA;  ++j) ll_row[sc_idx_c[j]] = sc_out_c[j];
}

// Censoring, truncation
void lba_survivor(const std::vector<int>&    idx,
                  const std::vector<double>&  bound,
                  const ParamTable&           pt,
                  const RaceSpec&             spec,
                  double* __restrict__        out,
                  RaceScratch&                scratch,
                  const std::vector<int>& is_ok, bool log_p)
{
  const double* __restrict__ bound_ = bound.data();
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);

  const int n = (int)idx.size();

  double* __restrict__ sc_teff  = scratch.dslot(LBA::t_eff);
  double* __restrict__ sc_v     = scratch.dslot(LBA::v);
  double* __restrict__ sc_sv    = scratch.dslot(LBA::sv);
  double* __restrict__ sc_B     = scratch.dslot(LBA::B);
  double* __restrict__ sc_A     = scratch.dslot(LBA::A);
  double* __restrict__ sc_denom = scratch.dslot(LBA::denom);   // new
  double* __restrict__ sc_out   = scratch.dslot(LBA::out);
  int*    __restrict__ sc_idx   = scratch.islot(ScratchInt::idx_win);

  double* __restrict__ sc_teff_c  = scratch.dslot(LBA::t_eff_c);
  double* __restrict__ sc_v_c     = scratch.dslot(LBA::v_c);
  double* __restrict__ sc_sv_c    = scratch.dslot(LBA::sv_c);
  double* __restrict__ sc_B_c     = scratch.dslot(LBA::B_c);
  double* __restrict__ sc_denom_c = scratch.dslot(LBA::denom_c);   // new
  double* __restrict__ sc_out_c   = scratch.dslot(LBA::out_c);
  int*    __restrict__ sc_idx_c   = scratch.islot(ScratchInt::idx_win_c);

  // --- Gather ---
  int n_core = 0, n_noA = 0;
  for (int j = 0; j < n; ++j) {
    const int i       = idx[j];
    if (!is_ok[i]) continue;
    const double teff = bound_[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { out[i] = log_p ? 0.0 : 1.0; continue; }  // Before the accumulator starts, survival is one.
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
    sc_out[j] = plba_core(sc_teff[j], sc_A[j], sc_B[j] + sc_A[j], sc_v[j], sc_sv[j], sc_denom[j], false, log_p);
  for (int j = 0; j < n_noA; ++j)
    sc_out_c[j] = plba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], sc_denom_c[j], false, log_p);

  for (int j = 0; j < n_core; ++j) out[sc_idx[j]] = sc_out[j];
  for (int j = 0; j < n_noA; ++j) out[sc_idx_c[j]] = sc_out_c[j];
}
