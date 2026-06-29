#ifndef RACE_GNG_H
#define RACE_GNG_H

// Race-model go/nogo: a withheld (inhibited) trial contributes the probability
// that the NO-GO accumulator finishes first within the response window,
//   P = INT_{lower}^{upper} f_nogo(t) * prod_{go} S_go(t) dt,
// with integrand f(t) = density_nogo(t) * prod_{losers}(1 - cdf_go(t)). There is
// no closed form, so we integrate numerically with hcubature (the same 1D
// integrator ddiff/pdiff already use). The integrator here is generic; only the
// per-accumulator scalar pdf/cdf is model-specific (RDM first).

#include <vector>
#include <cmath>
#include "model_RDM.h"   // scalar digt() / pigt() + A_EPS
#include "model_LBA.h"   // scalar dlba_core/noA, plba_core/noA + LBA_A_ASYMPTOTIC
#include "gauss.h"        // hcubature
#include "ParamTable.h"
#include "RaceSpec.h"

// Finite upper bound used when the integration window is unbounded above
// (no upper truncation). The withheld integrand f_nogo(t) * prod S_go(t) decays
// to 0 well before this for realistic RT-model parameters, so capping keeps the
// adaptive integrator on a finite domain with negligible truncation error.
static const double GNG_UPPER_CAP = 30.0;

// Per-trial RDM accumulator parameters, pre-scaled to the digt/pigt convention:
//   k = (B + 0.5*A)/s, l = v/s, a = 0.5*A/s, with decision time t - t0.
struct RDMGngIntegrand {
  int n_acc;
  int nogo;            // 0-based index of the no-go accumulator
  const double* k;     // length n_acc
  const double* l;
  const double* a;
  const double* t0;
};

// hcubature integrand: x[0] = absolute time t; retval[0] = f_nogo(t) * prod S_go(t).
inline int int_rdm_gng_withheld(unsigned /*dim*/, const double* x, void* p,
                                unsigned /*fdim*/, double* retval) {
  const RDMGngIntegrand* P = static_cast<const RDMGngIntegrand*>(p);
  const double t  = x[0];
  const int    ng = P->nogo;
  const double dens = digt(t - P->t0[ng], P->k[ng], P->l[ng], P->a[ng]);  // 0 if t<=t0
  if (!(dens > 0.0)) { retval[0] = 0.0; return 0; }
  double out = dens;
  for (int j = 0; j < P->n_acc; ++j) {
    if (j == ng) continue;
    const double Fj = pigt(t - P->t0[j], P->k[j], P->l[j], P->a[j]);       // 0 if t<=t0
    out *= (1.0 - Fj);
  }
  retval[0] = out;
  return 0;
}

// P(no-go accumulator finishes first in [lower, upper]) for one RDM trial, from
// raw per-accumulator parameters (length n_acc each). nogo is 0-based.
inline double rdm_gng_withheld_prob(const double* v, const double* B,
                                    const double* A, const double* t0,
                                    const double* s, int n_acc, int nogo,
                                    double lower, double upper,
                                    double abstol = 1e-8, double reltol = 1e-6,
                                    int maxeval = 6000) {
  std::vector<double> k(n_acc), l(n_acc), a(n_acc), t0v(n_acc);
  for (int j = 0; j < n_acc; ++j) {
    const double inv_s = 1.0 / s[j];
    a[j]   = 0.5 * A[j] * inv_s;
    l[j]   = v[j] * inv_s;
    k[j]   = (B[j] + 0.5 * A[j]) * inv_s;
    t0v[j] = t0[j];
  }
  RDMGngIntegrand ig{n_acc, nogo, k.data(), l.data(), a.data(), t0v.data()};
  double val = 0.0, err = 0.0, lo = lower, hi = upper;
  hcubature(int_rdm_gng_withheld, static_cast<void*>(&ig), 1, &lo, &hi,
            maxeval, abstol, reltol, &val, &err);
  if (!std::isfinite(val) || val < 0.0) val = 0.0;
  else if (val > 1.0) val = 1.0;
  return val;
}

// ParamTable-based backend (the generic gng_withheld_fn slot in RaceModelSetup).
// Gathers the trial's per-accumulator RDM params (rows base_row .. base_row+n_acc-1)
// from the natural-scale ParamTable and returns P(no-go finishes first in
// [lower, upper]). nogo is 0-based.
inline double rdm_gng_withheld(const ParamTable& pt, const RaceSpec& spec,
                               int base_row, int n_acc, int nogo,
                               double lower, double upper) {
  std::vector<double> v(n_acc), B(n_acc), A(n_acc), t0(n_acc), s(n_acc);
  for (int k = 0; k < n_acc; ++k) {
    const int r = base_row + k;
    v[k]  = pt.base(r, spec.col_v);
    B[k]  = pt.base(r, spec.col_B);
    A[k]  = pt.base(r, spec.col_A);
    t0[k] = pt.base(r, spec.col_t0);
    s[k]  = pt.base(r, spec.col_s);
  }
  return rdm_gng_withheld_prob(v.data(), B.data(), A.data(), t0.data(), s.data(),
                               n_acc, nogo, lower, upper);
}

// ===========================================================================
// LBA go/nogo withheld integrand (posdrift = TRUE, matching the LBA dfun/pfun)
// ===========================================================================
// Single-accumulator LBA pdf/cdf at decision time td, mirroring dlba()/plba():
// b = B + A, denom = pnorm(v/sv), A>LBA_A_ASYMPTOTIC selects the core/noA branch.
inline double lba_d_one(double td, double A, double b, double v, double sv, double denom) {
  if (td <= 0.0) return 0.0;
  double val = (A > LBA_A_ASYMPTOTIC) ? dlba_core(td, A, b, v, sv, denom)
                                      : dlba_noA (td,    b, v, sv, denom);
  return (val > 0.0) ? val : 0.0;
}
inline double lba_p_one(double td, double A, double b, double v, double sv, double denom) {
  if (td <= 0.0) return 0.0;
  double val = (A > LBA_A_ASYMPTOTIC) ? plba_core(td, A, b, v, sv, denom)
                                      : plba_noA (td,    b, v, sv, denom);
  if (val < 0.0) return 0.0;
  return (val > 1.0) ? 1.0 : val;
}

struct LBAGngIntegrand {
  int n_acc, nogo;
  const double* A; const double* b; const double* v; const double* sv;
  const double* t0; const double* denom;
};
inline int int_lba_gng_withheld(unsigned, const double* x, void* p, unsigned, double* retval) {
  const LBAGngIntegrand* P = static_cast<const LBAGngIntegrand*>(p);
  const double t = x[0]; const int ng = P->nogo;
  const double dens = lba_d_one(t - P->t0[ng], P->A[ng], P->b[ng], P->v[ng], P->sv[ng], P->denom[ng]);
  if (!(dens > 0.0)) { retval[0] = 0.0; return 0; }
  double out = dens;
  for (int j = 0; j < P->n_acc; ++j) {
    if (j == ng) continue;
    out *= (1.0 - lba_p_one(t - P->t0[j], P->A[j], P->b[j], P->v[j], P->sv[j], P->denom[j]));
  }
  retval[0] = out;
  return 0;
}
inline double lba_gng_withheld(const ParamTable& pt, const RaceSpec& spec,
                               int base_row, int n_acc, int nogo,
                               double lower, double upper) {
  std::vector<double> A(n_acc), b(n_acc), v(n_acc), sv(n_acc), t0(n_acc), denom(n_acc);
  for (int k = 0; k < n_acc; ++k) {
    const int r = base_row + k;
    A[k]  = pt.base(r, spec.col_A);
    v[k]  = pt.base(r, spec.col_v);
    sv[k] = pt.base(r, spec.col_sv);
    t0[k] = pt.base(r, spec.col_t0);
    b[k]  = pt.base(r, spec.col_B) + A[k];                 // b = B + A
    double d = R::pnorm(v[k] / sv[k], 0.0, 1.0, 1, 0);     // posdrift normaliser
    denom[k] = (d < 1e-10) ? 1e-10 : d;
  }
  LBAGngIntegrand ig{n_acc, nogo, A.data(), b.data(), v.data(), sv.data(), t0.data(), denom.data()};
  double val = 0.0, err = 0.0, lo = lower, hi = upper;
  hcubature(int_lba_gng_withheld, static_cast<void*>(&ig), 1, &lo, &hi, 6000, 1e-8, 1e-6, &val, &err);
  if (!std::isfinite(val) || val < 0.0) val = 0.0; else if (val > 1.0) val = 1.0;
  return val;
}

// ===========================================================================
// LNR go/nogo withheld integrand (lognormal race)
// ===========================================================================
struct LNRGngIntegrand {
  int n_acc, nogo;
  const double* m; const double* s; const double* t0;
};
inline int int_lnr_gng_withheld(unsigned, const double* x, void* p, unsigned, double* retval) {
  const LNRGngIntegrand* P = static_cast<const LNRGngIntegrand*>(p);
  const double t = x[0]; const int ng = P->nogo;
  const double td = t - P->t0[ng];
  if (td <= 0.0) { retval[0] = 0.0; return 0; }
  double out = R::dlnorm(td, P->m[ng], P->s[ng], 0);
  if (!(out > 0.0)) { retval[0] = 0.0; return 0; }
  for (int j = 0; j < P->n_acc; ++j) {
    if (j == ng) continue;
    const double tdj = t - P->t0[j];
    const double Fj = (tdj <= 0.0) ? 0.0 : R::plnorm(tdj, P->m[j], P->s[j], 1, 0);
    out *= (1.0 - Fj);
  }
  retval[0] = out;
  return 0;
}
inline double lnr_gng_withheld(const ParamTable& pt, const RaceSpec& spec,
                               int base_row, int n_acc, int nogo,
                               double lower, double upper) {
  std::vector<double> m(n_acc), s(n_acc), t0(n_acc);
  for (int k = 0; k < n_acc; ++k) {
    const int r = base_row + k;
    m[k]  = pt.base(r, spec.col_m);
    s[k]  = pt.base(r, spec.col_s);
    t0[k] = pt.base(r, spec.col_t0);
  }
  LNRGngIntegrand ig{n_acc, nogo, m.data(), s.data(), t0.data()};
  double val = 0.0, err = 0.0, lo = lower, hi = upper;
  hcubature(int_lnr_gng_withheld, static_cast<void*>(&ig), 1, &lo, &hi, 6000, 1e-8, 1e-6, &val, &err);
  if (!std::isfinite(val) || val < 0.0) val = 0.0; else if (val > 1.0) val = 1.0;
  return val;
}

#endif // RACE_GNG_H
