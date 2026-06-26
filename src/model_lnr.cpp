#include "model_lnr.h"

using namespace Rcpp;

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  // 0 = m, 1 = s, 2 = t0
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = PLNORM(rts[i] - pars(i,2), pars(i,0), pars(i,1));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return out;
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  // 0 = m, 1 = s, 2 = t0
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = DLNORM(rts[i] - pars(i,2), pars(i,0), pars(i,1));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return out;
}

// void dlnr_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const LogicalVector& winner,
//                double* ll_row)
// {
//   const int N = rts.size();
//
//   const double* rt = rts.begin();
//   const double* m  = &pt.base(0, spec.col_m);
//   const double* s  = &pt.base(0, spec.col_s);
//   const double* t0 = &pt.base(0, spec.col_t0);
//
//   int* win_ptr = LOGICAL(winner);
//
// #pragma omp simd
//   for (int i = 0; i < N; ++i) {
//     if (!win_ptr[i]) continue;
//
//     if (std::isnan(m[i])) { ll_row[i] = 0.0; continue; }
//
//     const double t_eff = rt[i] - t0[i];
//     if (t_eff <= 0.0)    { ll_row[i] = 0.0; continue; }
//
//     double pdf = DLNORM(t_eff, m[i], s[i]);
//     ll_row[i] = (std::isfinite(pdf) && pdf >= 0.0) ? pdf : 0.0;
//   }
// }
//
// void plnr_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const LogicalVector& winner,
//                double* ll_row)
// {
//   const int N = rts.size();
//
//   const double* rt = rts.begin();
//   const double* m  = &pt.base(0, spec.col_m);
//   const double* s  = &pt.base(0, spec.col_s);
//   const double* t0 = &pt.base(0, spec.col_t0);
//
//   int* win_ptr = LOGICAL(winner);
//
// #pragma omp simd
//   for (int i = 0; i < N; ++i) {
//     if (win_ptr[i]) continue;
//
//     if (std::isnan(m[i])) { ll_row[i] = 0.0; continue; }
//
//     const double t_eff = rt[i] - t0[i];
//     if (t_eff <= 0.0)    { ll_row[i] = 0.0; continue; }
//
//     double cdf = PLNORM(t_eff, m[i], s[i]);
//     if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
//     else if (cdf > 1.0)                         cdf = 1.0;
//     ll_row[i] = cdf;
//   }
// }


// Hot path: gather → compute → scatter
// Uses scratch.v for m, scratch.s for s — B, A, sv left unused.
void dlnr_plnr_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch)
{
  const double* __restrict__ rt = rts.begin();
  const double* __restrict__ m  = &pt.base(0, spec.col_m);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // --- Winners: gather (m → scratch.v, s → scratch.s) ---
  for (int j = 0; j < n_win; ++j) {
    const int i      = idx_win[j];
    scratch.t_eff[j] = rt[i] - t0[i];
    scratch.v[j]     = m[i];
    scratch.s[j]     = s[i];
  }

  // --- Winners: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_win; ++j) {
    if (std::isnan(scratch.v[j]) || scratch.t_eff[j] <= 0.0) {
      scratch.out[j] = 0.0;
      continue;
    }
    double val = DLNORM(scratch.t_eff[j], scratch.v[j], scratch.s[j]);
    scratch.out[j] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
  }

  // --- Winners: scatter ---
  for (int j = 0; j < n_win; ++j) ll_row[idx_win[j]] = scratch.out[j];

  // --- Losers: gather ---
  for (int j = 0; j < n_los; ++j) {
    const int i      = idx_los[j];
    scratch.t_eff[j] = rt[i] - t0[i];
    scratch.v[j]     = m[i];
    scratch.s[j]     = s[i];
  }

  // --- Losers: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_los; ++j) {
    if (std::isnan(scratch.v[j]) || scratch.t_eff[j] <= 0.0) {
      scratch.out[j] = 0.0;
      continue;
    }
    double val = PLNORM(scratch.t_eff[j], scratch.v[j], scratch.s[j]);
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    scratch.out[j] = val;
  }

  // --- Losers: scatter ---
  // SURVIVOR! 1-CDF, not CDF
  for (int j = 0; j < n_los; ++j) ll_row[idx_los[j]] = 1-scratch.out[j];
}

// =============================================================================
// plnr_fast — gather-scatter
// Fills survivor S_k(rt - t0) into ll_row for each i in idx.
// Used for loser accumulators in the main likelihood.
// =============================================================================
void plnr_fast(const NumericVector&    rts,
               const ParamTable&       pt,
               const RaceSpec&         spec,
               const std::vector<int>& idx,
               double* __restrict__    ll_row,
               RaceScratch&            scratch)
{
  const double* __restrict__ rt = rts.begin();
  const double* __restrict__ m_ = &pt.base(0, spec.col_m);
  const double* __restrict__ s_ = &pt.base(0, spec.col_s);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);

  double* __restrict__ sc_teff = scratch.t_eff.data();
  double* __restrict__ sc_m    = scratch.v.data();
  double* __restrict__ sc_s    = scratch.s.data();
  double* __restrict__ sc_out  = scratch.out.data();
  int*    __restrict__ sc_idx  = scratch.idx_win0.data();

  // --- Gather ---
  int n = 0;
  for (int j = 0; j < (int)idx.size(); ++j) {
    const int i       = idx[j];
    const double teff = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = 1.0; continue; }
    sc_teff[n] = teff;
    sc_m   [n] = m_[i];
    sc_s   [n] = s_[i];
    sc_idx [n] = i;
    n++;
  }

  // --- Compute ---
#pragma omp simd
  for (int j = 0; j < n; ++j) {
    sc_out[j] = PLNORM(sc_teff[j], sc_m[j], sc_s[j]);
  }

  // --- Scatter: 1 - cdf ---
  for (int j = 0; j < n; ++j) {
    double cdf = sc_out[j];
    if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
    else if (cdf > 1.0)                         cdf = 1.0;
    ll_row[sc_idx[j]] = 1.0 - cdf;
  }
}

// // =============================================================================
// // lnr_censor — gather-scatter
// // Three passes: lower-censored (idx_L), upper-censored (idx_U), both (idx_B).
// //
// //   idx_L: ll_row[i] = CDF(LC - t0)
// //   idx_U: ll_row[i] = 1 - CDF(UC - t0)
// //   idx_B: ll_row[i] = CDF(LC - t0) + 1 - CDF(UC - t0)
// //
// // For idx_B, LC and UC are gathered into separate scratch buffers (t_eff / k),
// // two PLNORM passes are run, then combined in the scatter.
// // =============================================================================
// void lnr_censor(const CensorSpec&    censor,
//                 const ParamTable&    pt,
//                 const RaceSpec&      spec,
//                 double* __restrict__ ll_row,
//                 RaceScratch&         scratch)
// {
//   const double* __restrict__ m_  = &pt.base(0, spec.col_m);
//   const double* __restrict__ s_  = &pt.base(0, spec.col_s);
//   const double* __restrict__ t0_ = &pt.base(0, spec.col_t0);
//   const double* __restrict__ LC_ = censor.LC.data();
//   const double* __restrict__ UC_ = censor.UC.data();
//
//   double* __restrict__ sc_teff = scratch.t_eff.data();
//   double* __restrict__ sc_m    = scratch.v.data();
//   double* __restrict__ sc_s    = scratch.s.data();
//   double* __restrict__ sc_out  = scratch.out.data();
//   int*    __restrict__ sc_idx  = scratch.idx_win0.data();
//
//   // ---- Pass 1: lower-censored — ll_row[i] = CDF(LC - t0) ----
//   {
//     int n = 0;
//     for (int j = 0; j < (int)censor.idx_L.size(); ++j) {
//       const int i       = censor.idx_L[j];
//       const double teff = LC_[i] - t0_[i];
//       if (teff <= 0.0) { ll_row[i] = 0.0; continue; }
//       sc_teff[n] = teff;
//       sc_m   [n] = m_[i];
//       sc_s   [n] = s_[i];
//       sc_idx [n] = i;
//       n++;
//     }
// #pragma omp simd
//     for (int j = 0; j < n; ++j)
//       sc_out[j] = PLNORM(sc_teff[j], sc_m[j], sc_s[j]);
//
//     for (int j = 0; j < n; ++j) {
//       double cdf = sc_out[j];
//       if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
//       else if (cdf > 1.0)                        cdf = 1.0;
//       ll_row[sc_idx[j]] = cdf;
//     }
//   }
//
//   // ---- Pass 2: upper-censored — ll_row[i] = 1 - CDF(UC - t0) ----
//   {
//     int n = 0;
//     for (int j = 0; j < (int)censor.idx_U.size(); ++j) {
//       const int i       = censor.idx_U[j];
//       const double teff = UC_[i] - t0_[i];
//       if (teff <= 0.0) { ll_row[i] = 1.0; continue; }
//       sc_teff[n] = teff;
//       sc_m   [n] = m_[i];
//       sc_s   [n] = s_[i];
//       sc_idx [n] = i;
//       n++;
//     }
// #pragma omp simd
//     for (int j = 0; j < n; ++j)
//       sc_out[j] = PLNORM(sc_teff[j], sc_m[j], sc_s[j]);
//
//     for (int j = 0; j < n; ++j) {
//       double cdf = sc_out[j];
//       if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
//       else if (cdf > 1.0)                         cdf = 1.0;
//       ll_row[sc_idx[j]] = 1.0 - cdf;
//     }
//   }
//
//   // ---- Pass 3: both-censored — ll_row[i] = CDF(LC - t0) + 1 - CDF(UC - t0) ----
//   // LC and UC gathered into separate time buffers (t_eff and k),
//   // two compute passes, combined in scatter.
//   {
//     double* __restrict__ sc_teff_LC = scratch.t_eff.data();  // LC offsets
//     double* __restrict__ sc_teff_UC = scratch.k.data();      // UC offsets (spare buffer)
//     double* __restrict__ sc_out_LC  = scratch.out.data();
//     double* __restrict__ sc_out_UC  = scratch.out_c.data();  // spare output buffer
//
//     int n = 0;
//     for (int j = 0; j < (int)censor.idx_B.size(); ++j) {
//       const int i = censor.idx_B[j];
//       sc_teff_LC[n] = LC_[i] - t0_[i];
//       sc_teff_UC[n] = UC_[i] - t0_[i];
//       sc_m      [n] = m_[i];
//       sc_s      [n] = s_[i];
//       sc_idx    [n] = i;
//       n++;
//     }
//
// #pragma omp simd
//     for (int j = 0; j < n; ++j)
//       sc_out_LC[j] = (sc_teff_LC[j] > 0.0) ? PLNORM(sc_teff_LC[j], sc_m[j], sc_s[j]) : 0.0;
//
// #pragma omp simd
//     for (int j = 0; j < n; ++j)
//       sc_out_UC[j] = (sc_teff_UC[j] > 0.0) ? PLNORM(sc_teff_UC[j], sc_m[j], sc_s[j]) : 0.0;
//
//     for (int j = 0; j < n; ++j) {
//       double cdf_LC = sc_out_LC[j];
//       double cdf_UC = sc_out_UC[j];
//       if      (!std::isfinite(cdf_LC) || cdf_LC < 0.0) cdf_LC = 0.0;
//       else if (cdf_LC > 1.0)                            cdf_LC = 1.0;
//       if      (!std::isfinite(cdf_UC) || cdf_UC < 0.0) cdf_UC = 0.0;
//       else if (cdf_UC > 1.0)                            cdf_UC = 1.0;
//       ll_row[sc_idx[j]] = cdf_LC + 1.0 - cdf_UC;
//     }
//   }
// }


// =============================================================================
// LNR — gather-scatter truncation filler
// Fills raw survivors S_k(bound - t0) into out[i] for each i in idx.
// Rows with teff <= 0 get S=1 directly (no CDF call needed).
// =============================================================================
void lnr_truncate(const std::vector<int>&     idx,
                  const std::vector<double>&  bound,
                  const ParamTable&           pt,
                  const RaceSpec&             spec,
                  double* __restrict__        out,
                  RaceScratch&                scratch)
{

  const double* __restrict__ bound_ = bound.data();
  const double* __restrict__ m_     = &pt.base(0, spec.col_m);
  const double* __restrict__ t0_    = &pt.base(0, spec.col_t0);
  const double* __restrict__ s_     = &pt.base(0, spec.col_s);

  double* __restrict__ sc_teff = scratch.t_eff.data();
  double* __restrict__ sc_m    = scratch.v.data();    // m reuses v buffer
  double* __restrict__ sc_s    = scratch.s.data();
  double* __restrict__ sc_out  = scratch.out.data();
  int*    __restrict__ sc_idx  = scratch.idx_win0.data();

  // --- Gather: filter teff <= 0, pack valid rows contiguously ---
  int n = 0;
  for (int j = 0; j < (int)idx.size(); ++j) {
    const int i       = idx[j];
    const double teff = bound_[i] - t0_[i];
    // if (teff <= 0.0) { out[i] = 1.0; continue; }  // S(t<=0) = 1, no CDF call needed
    if (teff <= 0.0) { continue; }  // Impossible parameter set - leave untouched, fall back to default value in out
    sc_teff[n] = teff;
    sc_m   [n] = m_[i];
    sc_s   [n] = s_[i];
    sc_idx [n] = i;
    n++;
  }

  // --- Compute: PLNORM over contiguous scratch — no branching, vectorisable ---
#pragma omp simd
  for (int j = 0; j < n; ++j) {
    sc_out[j] = PLNORM(sc_teff[j], sc_m[j], sc_s[j]);
  }

  // --- Scatter: clamp cdf to [0,1], write survivor 1 - cdf ---
  for (int j = 0; j < n; ++j) {
    double cdf = sc_out[j];
    if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
    else if (cdf > 1.0)                        cdf = 1.0;
    out[sc_idx[j]] = 1.0 - cdf;
  }
}
