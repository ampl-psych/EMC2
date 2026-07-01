#ifndef ss_exg_analytic_h
#define ss_exg_analytic_h

// ---------------------------------------------------------------------------
// Exact closed form for the SS-EXG stop-success integral with ONE go runner,
// plus a self-contained bivariate normal CDF (port of Alan Genz's BVND).
//
// Derivation and pure-R validation (relerr ~1e-16 against 1e-12-tolerance
// adaptive integration): WorkingTests/stop_success_methods.R (Section 1). Both the
// ex-Gaussian density and survivor are (exponential x normal-CDF) objects, so
// f_S(x) * S_G(x + SSD) expands over the go survivor's two branches into terms
// C_A exp(gamma x) Phi(a_0 + b_0 x) Phi(a_1 + b_1 x) with gamma < 0; one
// integration by parts per term gives, for L = lbS (or -Inf), U = Inf:
//
//   J_A(L) = exp(gamma L)/|gamma| * Phi(a_0 + b_0 L) Phi(a_1 + b_1 L)
//          + (1/|gamma|) sum_k sign(b_k) exp(delta_k) Q_k(L)
//
//   xbar_k = -a_k/b_k,  mustar_k = xbar_k + gamma/b_k^2,
//   delta_k = gamma xbar_k + gamma^2/(2 b_k^2)
//   Q_k(L)  = INT_{vL}^{Inf} phi(v) Phi(alpha_j + beta_j v) dv,  j != k,
//             vL = |b_k| (L - mustar_k), beta_j = b_j/|b_k|,
//             alpha_j = a_j + b_j mustar_k
//
// L = -Inf (full line):  Q_k = Phi(alpha/sqrt(1+beta^2))      -> 4 pnorm total
// L finite (truncated):  Q_k = Phi2(-vL, alpha/s; rho=beta/s) -> 4 Phi2 +
//                        2 boundary terms,  s = sqrt(1+beta^2)
//
// Guards (deterministic, computed before any CDF call; on failure the caller
// falls back to GL — the small discontinuity at the switch is harmless for
// particle Metropolis, which uses no gradients):
//   * sigS > ANALYTIC_MAX_SIGS_OVER_TAUS * tauS: two reasons converge here.
//     (1) dexg()'s Mills-ratio tail branch (z < -8) starts to cover
//         non-negligible integrand mass, so the quadrature routes would
//         disagree with this exact form — but dexg now carries the
//         -1/z^2 + 2.5/z^4 corrections (exgaussian_functions.h), shrinking that
//         disagreement to ~1e-6, which let the threshold relax from 4 to 10.
//     (2) As tauS -> 0, logC0 (= -log tauS + muS/tauS + sigS^2/2tauS^2) and the
//         per-term exponent gain large cancelling magnitudes; the net stays
//         moderate but precision erodes, so we keep a finite ceiling and hand
//         the (now sharp, near-Gaussian) peak to the density-bumped GL route
//         (gl_auto_nodes), which resolves it well.
//   * any term exponent logC_A + delta_k (or logC_A + gamma*L) > 600, or any
//     non-finite intermediate: exp() overflow regime.
//   * kinked domain (lbS + SSD < lbG, or lbS = -Inf with finite lbG): the
//     truncated go survivor is piecewise on the integration domain.
//
// cens port: this header must be included AFTER exgaussian_functions.h (uses
// pexg() for the truncation normaliser so it matches the quadrature integrand
// exactly) and after the SS compat shim (emc2_isfinite). We include
// exgaussian_functions.h directly so the header is self-sufficient.
// ---------------------------------------------------------------------------

#include <cmath>
#include <algorithm>
#include "exgaussian_functions.h"

// Largest sigS/tauS for which the exact n_go==1 closed form is trusted; above
// it "auto" hands off to the density-bumped GL route. Relaxed from 4 to 10 once
// dexg()'s tail branch gained the higher-order Mills corrections (see the guard
// note above). Mirrored by stop_success_texg_analytic1_R in R/stop_success_gl.R.
constexpr double ANALYTIC_MAX_SIGS_OVER_TAUS = 10.0;

// ---------------------------------------------------------------------------
// Bivariate normal CDF: port of Alan Genz's BVND (tvpack.f, Fortran),
// self-contained except for the univariate normal CDF from Rmath.
// Maximum absolute error ~5e-16. Mirrored in R by pbvn_R() (stop_success_gl.R)
// so the R and C++ analytic routes agree to machine precision.
// ---------------------------------------------------------------------------

inline double bvn_phid(double z) { return R::pnorm(z, 0.0, 1.0, 1, 0); }

// P(X > dh, Y > dk) for standard bivariate normal with correlation r.
inline double bvn_genz_bvnu(double dh, double dk, double r) {
  // Gauss-Legendre half-rules: 6-point (|r|<.3), 12-point (<.75), 20-point.
  static const double XW6[3][2] = {
    {-0.9324695142031522, 0.1713244923791705},
    {-0.6612093864662647, 0.3607615730481384},
    {-0.2386191860831970, 0.4679139345726904}};
  static const double XW12[6][2] = {
    {-0.9815606342467191, 0.04717533638651177},
    {-0.9041172563704750, 0.1069393259953183},
    {-0.7699026741943050, 0.1600783285433464},
    {-0.5873179542866171, 0.2031674267230659},
    {-0.3678314989981802, 0.2334925365383547},
    {-0.1252334085114692, 0.2491470458134029}};
  static const double XW20[10][2] = {
    {-0.9931285991850949, 0.01761400713915212},
    {-0.9639719272779138, 0.04060142980038694},
    {-0.9122344282513259, 0.06267204833410906},
    {-0.8391169718222188, 0.08327674157670475},
    {-0.7463319064601508, 0.1019301198172404},
    {-0.6360536807265150, 0.1181945319615184},
    {-0.5108670019508271, 0.1316886384491766},
    {-0.3737060887154196, 0.1420961093183821},
    {-0.2277858511416451, 0.1491729864726037},
    {-0.07652652113349733, 0.1527533871307259}};
  const double TWOPI = 6.283185307179586;

  if (std::isnan(dh) || std::isnan(dk) || std::isnan(r)) return NA_REAL;
  if (dh == R_PosInf || dk == R_PosInf) return 0.0;
  if (dh == R_NegInf) return (dk == R_NegInf) ? 1.0 : bvn_phid(-dk);
  if (dk == R_NegInf) return bvn_phid(-dh);

  const double (*xw)[2];
  int lg;
  if (std::fabs(r) < 0.3)       { xw = XW6;  lg = 3; }
  else if (std::fabs(r) < 0.75) { xw = XW12; lg = 6; }
  else                          { xw = XW20; lg = 10; }

  double h = dh, k = dk, hk = h * k, bvn = 0.0;
  if (std::fabs(r) < 0.925) {
    if (std::fabs(r) > 0.0) {
      const double hs = (h * h + k * k) / 2.0;
      const double asr = std::asin(r);
      for (int i = 0; i < lg; ++i) {
        for (int is = -1; is <= 1; is += 2) {
          const double sn = std::sin(asr * (is * xw[i][0] + 1.0) / 2.0);
          bvn += xw[i][1] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
        }
      }
      bvn = bvn * asr / (2.0 * TWOPI);
    }
    bvn += bvn_phid(-h) * bvn_phid(-k);
  } else {
    if (r < 0.0) { k = -k; hk = -hk; }
    if (std::fabs(r) < 1.0) {
      const double as = (1.0 - r) * (1.0 + r);
      double a = std::sqrt(as);
      const double bs = (h - k) * (h - k);
      const double c = (4.0 - hk) / 8.0;
      const double d = (12.0 - hk) / 16.0;
      double asr = -(bs / as + hk) / 2.0;
      if (asr > -100.0) {
        bvn = a * std::exp(asr) *
          (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as * as / 5.0);
      }
      if (-hk < 100.0) {
        const double b = std::sqrt(bs);
        bvn -= std::exp(-hk / 2.0) * std::sqrt(TWOPI) * bvn_phid(-b / a) * b *
          (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
      }
      a /= 2.0;
      for (int i = 0; i < lg; ++i) {
        for (int is = -1; is <= 1; is += 2) {
          double xs = a * (is * xw[i][0] + 1.0);
          xs = xs * xs;
          const double rs = std::sqrt(1.0 - xs);
          asr = -(bs / xs + hk) / 2.0;
          if (asr > -100.0) {
            bvn += a * xw[i][1] * std::exp(asr) *
              (std::exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs -
               (1.0 + c * xs * (1.0 + d * xs)));
          }
        }
      }
      bvn = -bvn / TWOPI;
    }
    if (r > 0.0) {
      bvn += bvn_phid(-std::max(h, k));
    } else {
      bvn = -bvn;
      if (k > h) {
        // Phi(k) - Phi(h); use upper tails when both positive for accuracy
        if (h >= 0.0) bvn += bvn_phid(-h) - bvn_phid(-k);
        else          bvn += bvn_phid(k) - bvn_phid(h);
      }
    }
  }
  if (bvn < 0.0) bvn = 0.0;
  if (bvn > 1.0) bvn = 1.0;
  return bvn;
}

// P(X <= h, Y <= k), correlation r.
inline double bvn_cdf(double h, double k, double r) {
  return bvn_genz_bvnu(-h, -k, r);
}

// ---------------------------------------------------------------------------
// n_go == 1 stop-success probability (NOT log), matching the semantics of
// ss_texg_stop_success_lpdf with upper = Inf. Sets ok = false (and returns 0)
// when a guard trips and the caller should use GL instead.
// ---------------------------------------------------------------------------
inline double ss_texg_stop_success_analytic1(
    double SSD,
    double muG, double sigG, double tauG, double lbG,
    double muS, double sigS, double tauS, double lbS,
    bool& ok
) {
  ok = false;
  if (!(sigS > 0.0) || !(tauS > 0.0) || !(sigG > 0.0) || !(tauG > 0.0)) return 0.0;
  if (sigS > ANALYTIC_MAX_SIGS_OVER_TAUS * tauS) return 0.0;  // see guard note

  const bool full_line = (lbS == R_NegInf) && (lbG == R_NegInf);
  if (!full_line) {
    if (!emc2_isfinite(lbS)) return 0.0;        // lbS = -Inf, finite lbG: kink
    if (!(lbS + SSD >= lbG)) return 0.0;        // kink inside the domain
  }

  const double m = muG - SSD;
  const double a0 = -muS / sigS - sigS / tauS;
  const double b0 = 1.0 / sigS;
  const double logC0 = -std::log(tauS) + muS / tauS + sigS * sigS / (2.0 * tauS * tauS);

  double total = 0.0;
  for (int A = 0; A < 2; ++A) {
    double gam, logC, a1, b1;
    if (A == 0) {           // go survivor's Gaussian branch
      gam = -1.0 / tauS;
      logC = logC0;
      a1 = m / sigG;
      b1 = -1.0 / sigG;
    } else {                // go survivor's exponential branch
      gam = -1.0 / tauS - 1.0 / tauG;
      logC = logC0 + m / tauG + sigG * sigG / (2.0 * tauG * tauG);
      a1 = -m / sigG - sigG / tauG;
      b1 = 1.0 / sigG;
    }
    const double inv_abs_gam = -1.0 / gam;      // gam < 0
    const double a[2] = {a0, a1};
    const double b[2] = {b0, b1};

    if (!full_line) {       // boundary term at L = lbS
      const double e = logC + gam * lbS;
      if (!(e <= 600.0)) return 0.0;            // overflow guard (catches NaN)
      total += std::exp(e) * inv_abs_gam *
        bvn_phid(a0 + b0 * lbS) * bvn_phid(a1 + b1 * lbS);
    }

    for (int kk = 0; kk < 2; ++kk) {
      const double ak = a[kk], bk = b[kk];
      const double xbar = -ak / bk;
      const double mustar = xbar + gam / (bk * bk);
      const double e = logC + gam * xbar + gam * gam / (2.0 * bk * bk);
      if (!(e <= 600.0)) return 0.0;            // overflow guard (catches NaN)
      const int j = 1 - kk;
      const double absbk = std::fabs(bk);
      const double beta = b[j] / absbk;
      const double alpha = a[j] + b[j] * mustar;
      const double s = std::sqrt(1.0 + beta * beta);
      double Q;
      if (full_line) {
        Q = bvn_phid(alpha / s);
      } else {
        const double vL = absbk * (lbS - mustar);
        Q = bvn_cdf(-vL, alpha / s, beta / s);
      }
      const double sgn = (bk > 0.0) ? 1.0 : -1.0;
      total += sgn * std::exp(e) * inv_abs_gam * Q;
    }
  }

  if (!emc2_isfinite(total) || total <= 0.0) return 0.0;
  double p = total;
  if (!full_line) {
    // truncation normaliser: matches the dtexg/ptexg integrand of the
    // quadrature routes (uses the same pexg)
    const double NS = pexg(lbS, muS, sigS, tauS, false, false);
    const double NG = emc2_isfinite(lbG)
      ? pexg(lbG, muG, sigG, tauG, false, false) : 1.0;
    const double N = NS * NG;
    if (!emc2_isfinite(N) || !(N > 0.0)) return 0.0;
    p = total / N;
  }
  if (!emc2_isfinite(p) || p <= 0.0 || p > 1.0 + 1e-10) return 0.0;
  if (p > 1.0) p = 1.0;
  ok = true;
  return p;
}

#endif
