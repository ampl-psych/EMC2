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

// -----------------------------------------------------------------------------
// LNR survivor with known response — numerical integration
//
// Computes P(winner finishes first in [lower, upper]) for one trial by
// integrating f_winner(t) * prod_{j != winner} S_j(t) over [lower, upper].
// -----------------------------------------------------------------------------

struct LNRIntegrand {
  int           n_acc;
  int           winner;   // 0-based within-trial index
  const double* m;        // lognormal mean per accumulator
  const double* s;        // lognormal sd per accumulator
  const double* t0;
};

// hcubature callback: x[0] = absolute time t
// retval[0] = f_winner(t - t0_winner) * prod_{j != winner} S_j(t - t0_j)
inline int int_lnr_survivor(unsigned /*dim*/, const double* x, void* p,
                            unsigned /*fdim*/, double* retval)
{
  const LNRIntegrand* P = static_cast<const LNRIntegrand*>(p);
  const double t       = x[0];
  const int    w       = P->winner;
  const double t_eff_w = t - P->t0[w];

  // Density of winner — early exit if t <= t0
  if (t_eff_w <= 0.0) { retval[0] = 0.0; return 0; }
  const double dens = DLNORM(t_eff_w, P->m[w], P->s[w]);
  if (!(dens > 0.0))  { retval[0] = 0.0; return 0; }

  // Survivor of each loser: 1 - plnorm(t - t0_j, m_j, s_j)
  double out = dens;
  for (int j = 0; j < P->n_acc; ++j) {
    if (j == w) continue;
    const double t_eff_j = t - P->t0[j];
    const double cdf = (t_eff_j <= 0.0) ? 0.0 : PLNORM(t_eff_j, P->m[j], P->s[j]);
    out *= (1.0 - cdf);
  }
  retval[0] = out;
  return 0;
}

// Integrate f_winner * prod S_losers over [lower, upper] for one trial.
// m/s/t0 are raw (no pre-scaling needed for LNR), length n_acc, caller-owned.
inline double lnr_survivor_scalar(const double* m, const double* s,
                                  const double* t0,
                                  int n_acc, int winner,
                                  double lower, double upper,
                                  double abstol = 1e-8, double reltol = 1e-6,
                                  int maxeval = 6000)
{
  LNRIntegrand ig{ n_acc, winner, m, s, t0 };
  double val = 0.0, err = 0.0;
  hcubature(int_lnr_survivor, static_cast<void*>(&ig), 1, &lower, &upper,
            maxeval, abstol, reltol, &val, &err);
  if (!std::isfinite(val) || val < 0.0) return 0.0;
  if (val > 1.0)                         return 1.0;
  return val;
}

// -----------------------------------------------------------------------------
// lnr_survivor_with_response()
//
// For each entry j in idx:
//   - idx[j]    = base row (t * n_acc) into ParamTable
//   - winner[j] = within-trial winner index [0, n_acc-1]
//   - lower[j]  = lower integration bound for this trial
//   - upper[j]  = upper integration bound for this trial
//   - out[j]    = integral result, one value per trial
//
// Scratch buffers allocated once (length n_acc), reused across trials.
// No pre-scaling needed — LNR parameters enter DLNORM/PLNORM directly.
// -----------------------------------------------------------------------------
void lnr_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out)
{
  std::vector<double> m_buf(n_acc), s_buf(n_acc), t0_buf(n_acc);

  const int n = (int)idx.size();
  for (int j = 0; j < n; ++j) {
    const int base = idx[j];  // = t * n_acc

    for (int acc = 0; acc < n_acc; ++acc) {
      const int row  = base + acc;
      m_buf[acc]     = pt.base(row, spec.col_m);
      s_buf[acc]     = pt.base(row, spec.col_s);
      t0_buf[acc]    = pt.base(row, spec.col_t0);
    }

    out[j] = lnr_survivor_scalar(m_buf.data(), s_buf.data(), t0_buf.data(),
                                 n_acc, winner[j], lower[j], upper[j]);
  }
}
