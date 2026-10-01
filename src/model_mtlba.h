#ifndef MODEL_MTLBA_H
#define MODEL_MTLBA_H

// Multiple-threshold LBA (MTLBA). C++ twin of R/model_MTLBA.R; see there for
// the model and the numerical method.
//
// Joint per-trial kernel: each trial has two contiguous accumulator rows. A
// response (R, RR) is looked up as (stopping accumulator s, its threshold j,
// the other accumulator's band lo):
//   ll = log f_LBA(dt; A_s, T_{s,j}, v_s, sv_s) + log P(T_{o,lo} < e_o < T_{o,lo+1} | v_o > 0),
// dt = rt - t0_s (t0 from the STOPPING accumulator's row). BE: s = R, j = K,
// lo = K - RR. TC: folded position p = 1 .. 2K, p odd = accumulator 1 stops at
// K - (p-1)/2, p even = accumulator 2 stops at p/2, lo = K - j. A don't-know
// trial (RR = 1) is the log-sum of the (R = 1, RR = 1) and (R = 2, RR = 1)
// kernels. Thresholds T_j = A + (b - A) d_j, b = A + B, d_{K-k} =
// exp(-(c1 + ... + ck)), built here from c1 .. c{K-1} (Ttransform is not
// applied in C++). The stopping rule and the don't-know flag are read from the
// model's c_options (attr(dadm, "c_options"), set by design_model), K from the
// number of criteria.
//
// The band is formed on the log scale: sentinels H(lowest) = Phi(-v/sv),
// S(lowest) = Phi(v/sv); band * Phi(v/sv) = S(lo) - S(hi) = H(hi) - H(lo),
// taking the route with the smaller leading term; log G by the Mills ratio
// (z < 0) and Laplace's continued fraction (z < -5, 40 terms). Log tails use
// R::pnorm(log.p) (Rmath), not the package's Hart pnorm, which cancels there;
// Rmath's pnorm/dnorm are pure arithmetic for a standard normal, so safe in
// the OpenMP loop.
//
// Censoring and truncation (EMC2's make_missing coding: missingness 1/2/3,
// LC/UC/LT/UT columns, read from each trial's first row): a censored trial
// (response not retained, R = NA) is the log probability that the race stops
// in its censored window(s); a truncated trial is divided by P(LT < rt < UT).
// Both come from the stopping-time survivor (mtlba_logstop), which needs t0
// equal for the two accumulators of such a trial.
//
// Correlated drift rates (rho != 0, read from accumulator 1's row; the two
// drifts bivariate normal truncated to both positive): the kernel is a
// one-dimensional integral over the stopping accumulator's drift, done by the
// R default rule mtlba_quad_rule(): 16 Gauss-Legendre nodes per panel, the
// window of mtlba_quad_window (half-width 8 sds), one panel per 4 widths of
// the integrand's narrowest feature, 1 to 64 panels (mtlba_logkernel_rho),
// divided by P(both > 0) from Genz's bivariate normal (mtlba_bvu). The
// censoring/truncation terms are one integral over accumulator 1's drift per
// distinct parameter row and bound (mtlba_stop_rho1, the R reference rule:
// half-width 12 sds, one panel per feature width, at most 1024 panels per
// piece), cached per likelihood call. Trials with rho = 0 take the closed forms.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "ParamTable.h"
#include "model_mtlnr.h"   // mtlnr_log_diff_exp, mtlnr_log_pnorm_interval

struct MTLBASpec {
  int  K       = 1;
  bool TC      = false;   // threshold count (else balance of evidence)
  bool DK      = false;   // rating 1 is one shared don't-know response
  int  col_v   = -1;
  int  col_sv  = -1;
  int  col_B   = -1;
  int  col_A   = -1;
  int  col_t0  = -1;
  int  col_rho = -1;
  std::vector<int>    col_c;   // c1 .. c{K-1}
  std::vector<double> rt;      // per trial
  std::vector<int>    RR;      // per trial, 1 .. K
  std::vector<int>    R;       // per trial, choice 1 / 2 (accumulator level)
  std::vector<int>    row1;    // per trial, row of accumulator 1 (first lR level)
  // censoring / truncation, per trial
  std::vector<int>    miss;    // 0 observed, 1 lower, 2 upper, 3 both censored
  std::vector<char>   trunc;   // LT > 0 or UT finite
  std::vector<double> LT, UT, LC, UC;
  bool any_missing = false;    // any censored or truncated trial
};

// Upper limit on the number of ratings (stack buffers in mtlba_logstop).
constexpr int MTLBA_MAX_K = 32;

// Data, model flags and column lookups, once per call (main thread). Refuses
// what the likelihood does not (yet) define: other than two accumulators,
// censoring with the response retained, missingness code 4 and missing
// responses that are not censored (contamination).
inline MTLBASpec make_mtlba_spec(const Rcpp::DataFrame& data, const ParamTable& pt, int n_lR)
{
  if (n_lR != 2) Rcpp::stop("MTLBA: needs exactly two accumulators per trial (found %d)", n_lR);
  MTLBASpec s;
  s.col_v   = pt.base_index_for("v");
  s.col_sv  = pt.base_index_for("sv");
  s.col_B   = pt.base_index_for("B");
  s.col_A   = pt.base_index_for("A");
  s.col_t0  = pt.base_index_for("t0");
  s.col_rho = pt.base_index_for("rho");
  for (int k = 1; ; ++k) {
    auto it = pt.name_to_base_idx.find("c" + std::to_string(k));
    if (it == pt.name_to_base_idx.end()) break;
    s.col_c.push_back(it->second);
  }
  s.K = (int)s.col_c.size() + 1;

  // rule and don't know: the model's c_options, stored on the data by design_model()
  SEXP opt = data.attr("c_options");
  if (Rf_isNull(opt))
    Rcpp::stop("MTLBA: the data have no c_options attribute (remake them with design_model() or make_emc())");
  Rcpp::NumericVector o(opt);
  Rcpp::CharacterVector on = o.names();
  auto opt_val = [&](const char* x) {
    for (int i = 0; i < on.size(); ++i) if (on[i] == x) return (double)o[i];
    Rcpp::stop("MTLBA: c_options has no '%s'", x);
    return 0.0;
  };
  s.TC = opt_val("TC") != 0;
  s.DK = opt_val("dont_know") != 0;
  if ((int)opt_val("n_ratings") != s.K)
    Rcpp::stop("MTLBA: the model has %d ratings but the parameters give %d",
               (int)opt_val("n_ratings"), s.K);

  const int n = data.nrow();
  if (n % 2 != 0) Rcpp::stop("MTLBA: odd number of data rows");
  Rcpp::CharacterVector nm = data.names();
  auto has = [&](const char* x) { for (int i = 0; i < nm.size(); ++i) if (nm[i] == x) return true; return false; };

  if (!has("RR") && s.K > 1) Rcpp::stop("MTLBA: the data have no rating column RR");
  if (s.K > MTLBA_MAX_K) Rcpp::stop("MTLBA: at most %d ratings", MTLBA_MAX_K);

  auto num_col = [&](const char* x, double def) {
    return has(x) ? Rcpp::as<Rcpp::NumericVector>(data[x]) : Rcpp::NumericVector(n, def);
  };
  Rcpp::NumericVector rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  Rcpp::IntegerVector lR = data["lR"];
  Rcpp::IntegerVector Rv = data["R"];
  Rcpp::LogicalVector win = data["winner"];
  Rcpp::NumericVector rr = num_col("RR", 1.0);
  Rcpp::IntegerVector ms = has("missingness") ? Rcpp::as<Rcpp::IntegerVector>(data["missingness"])
                                              : Rcpp::IntegerVector(n, NA_INTEGER);
  Rcpp::NumericVector lt = num_col("LT", 0.0), ut = num_col("UT", R_PosInf);
  Rcpp::NumericVector lc = num_col("LC", 0.0), uc = num_col("UC", R_PosInf);
  const bool has_LC = has("LC"), has_UC = has("UC");
  const int nt = n / 2;
  s.rt.resize(nt); s.RR.resize(nt); s.R.resize(nt); s.row1.resize(nt);
  s.miss.resize(nt); s.trunc.resize(nt);
  s.LT.resize(nt); s.UT.resize(nt); s.LC.resize(nt); s.UC.resize(nt);
  for (int t = 0; t < nt; ++t) {
    const int a = 2 * t, b = a + 1;
    if (!((lR[a] == 1 && lR[b] == 2) || (lR[a] == 2 && lR[b] == 1)))
      Rcpp::stop("MTLBA: trial %d does not have one row per accumulator", t + 1);
    s.row1[t] = (lR[a] == 1) ? a : b;
    s.LT[t] = lt[a]; s.UT[t] = ut[a]; s.LC[t] = lc[a]; s.UC[t] = uc[a];
    s.trunc[t] = (lt[a] > 0 && std::isfinite(lt[a])) || std::isfinite(ut[a]);
    const int m = ms[a];
    if (m != NA_INTEGER) {
      if (m < 1 || m > 3)
        Rcpp::stop("MTLBA: missingness code %d is not supported (row %d)", m, a + 1);
      if (Rv[a] != NA_INTEGER || Rv[b] != NA_INTEGER)
        Rcpp::stop("MTLBA: censoring with the response retained (LCresponse/UCresponse = TRUE) is not supported; censored trials need R = NA");
      if ((m == 1 || m == 3) && !has_LC)
        Rcpp::stop("missingness contains lower-censored rows (1 or 3) but no LC column found in data");
      if ((m == 2 || m == 3) && !has_UC)
        Rcpp::stop("missingness contains upper-censored rows (2 or 3) but no UC column found in data");
      s.miss[t] = m;
      s.R[t] = 0; s.RR[t] = 0; s.rt[t] = NA_REAL;
      continue;
    }
    s.miss[t] = 0;
    if (Rcpp::LogicalVector::is_na(win[a]) || Rcpp::LogicalVector::is_na(win[b]) ||
        (bool)win[a] == (bool)win[b] || Rv[a] == NA_INTEGER || Rv[b] == NA_INTEGER ||
        !std::isfinite(rt[a]))
      Rcpp::stop("MTLBA: trial %d has a missing response or rt that is not censored (contamination is not supported)", t + 1);
    const double v = rr[a];
    if (!(v >= 1 && v <= s.K && v == std::floor(v)))
      Rcpp::stop("MTLBA: RR must be an integer in 1 .. %d (row %d)", s.K, a + 1);
    s.R[t]    = win[a] ? lR[a] : lR[b];
    s.RR[t]   = (int)v;
    s.rt[t]   = rt[a];
  }
  for (int t = 0; t < nt; ++t) if (s.miss[t] || s.trunc[t]) { s.any_missing = true; break; }
  return s;
}

// R's pmax for two scalars: NaN if either is NaN.
inline double mtlba_pmax(double a, double b)
{
  if (std::isnan(a) || std::isnan(b)) return R_NaN;
  return a > b ? a : b;
}

// log(exp(la) - exp(lb)): NaN if either is NaN, -Inf unless la > lb.
inline double mtlba_ldiff(double la, double lb)
{
  if (std::isnan(la) || std::isnan(lb)) return R_NaN;
  if (!(la > lb)) return R_NegInf;
  return mtlnr_log_diff_exp(la, lb);
}

// log(exp(a) + exp(b))
inline double mtlba_lse(double a, double b)
{
  if (std::isnan(a) || std::isnan(b)) return R_NaN;
  const double m = a > b ? a : b;
  if (m == R_NegInf) return R_NegInf;
  return m + std::log1p(std::exp(-std::fabs(a - b)));
}

// log G(z), G(z) = z Phi(z) + phi(z), accurate in the lower tail.
inline double mtlba_logG(double z)
{
  if (std::isnan(z)) return R_NaN;
  if (z >= 0) return std::log(z * R::pnorm(z, 0.0, 1.0, 1, 0) + R::dnorm(z, 0.0, 1.0, 0));
  const double x = -z;
  const double ld = R::dnorm(x, 0.0, 1.0, 1);
  if (x < 5) {
    const double xR = x * std::exp(R::pnorm(-x, 0.0, 1.0, 1, 1) - ld);
    return ld + std::log1p(-xR);
  }
  double tt = x;
  for (int k = 40; k >= 2; --k) tt = x + k / tt;   // as the R twin (40 = 120 to the bit)
  const double q = 1.0 / tt;
  return ld + std::log(q) - std::log(x + q);
}

// log of the mean of Phi (lower) or 1 - Phi over the window [z1 - al, z1] by
// 4-node Gauss-Legendre (the R twin's mtlba_gl(4)): positive terms.
inline double mtlba_logmean4(double z1, double al, bool lower)
{
  static const double U4[4]  = {0.069431844202973714, 0.33000947820757187,
                                0.66999052179242813, 0.93056815579702623};
  static const double LW4[4] = {-1.7491171784568635, -1.1206352922365659,
                                -1.1206352922365659, -1.7491171784568635};   // log(w / 2)
  const double zb = z1 - al;
  double L[4], M = R_NegInf;
  for (int i = 0; i < 4; ++i) {
    L[i] = R::pnorm(zb + al * U4[i], 0.0, 1.0, lower, 1) + LW4[i];
    if (L[i] > M) M = L[i];
  }
  if (M == R_NegInf) return R_NegInf;
  double sum = 0.0;
  for (int i = 0; i < 4; ++i) sum += std::exp(L[i] - M);
  return M + std::log(sum);
}

// log H(y) and log S(y) for evidence U(0, A) + v t, v ~ N(v, sv) untruncated,
// y >= A (separately: the band needs only one route). A narrow window,
// al max(1, |z1|) < 0.1, takes the 4-node mean (the difference of G loses
// digits as al -> 0).
inline double mtlba_logH(double y, double t, double A, double v, double sv)
{
  const double z1 = (y - t * v) / (t * sv);
  if (A <= 1e-10) return R::pnorm(z1, 0.0, 1.0, 1, 1);
  const double al = A / (t * sv);
  if (al * mtlba_pmax(1.0, std::fabs(z1)) < 0.1) return mtlba_logmean4(z1, al, true);
  return mtlba_ldiff(mtlba_logG(z1), mtlba_logG(z1 - al)) - std::log(al);
}
inline double mtlba_logS(double y, double t, double A, double v, double sv)
{
  const double z1 = (y - t * v) / (t * sv);
  if (A <= 1e-10) return R::pnorm(z1, 0.0, 1.0, 0, 1);
  const double al = A / (t * sv);
  if (al * mtlba_pmax(1.0, std::fabs(z1)) < 0.1) return mtlba_logmean4(z1, al, false);
  return mtlba_ldiff(mtlba_logG(-(z1 - al)), mtlba_logG(-z1)) - std::log(al);
}

// log P(lo_y < e < hi_y, v > 0) at decision time t; lowest: no threshold below.
// The route (S or H) is chosen by the smaller leading term, S(lo) vs H(hi).
inline double mtlba_logband_un(double t, double A, double v, double sv,
                               double lo_y, double hi_y, bool lowest)
{
  const double zv  = v / sv;
  const double hH  = mtlba_logH(hi_y, t, A, v, sv);
  const double lS  = lowest ? R::pnorm(zv, 0.0, 1.0, 1, 1) : mtlba_logS(lo_y, t, A, v, sv);
  if (lS < hH) return mtlba_ldiff(lS, mtlba_logS(hi_y, t, A, v, sv));
  const double lH = lowest ? R::pnorm(-zv, 0.0, 1.0, 1, 1) : mtlba_logH(lo_y, t, A, v, sv);
  return mtlba_ldiff(hH, lH);
}

// log P(lo_y < e < hi_y | v > 0): the band over Phi(v/sv), floored at 1e-10.
inline double mtlba_logband(double t, double A, double v, double sv,
                            double lo_y, double hi_y, bool lowest)
{
  return mtlba_logband_un(t, A, v, sv, lo_y, hi_y, lowest) -
         mtlba_pmax(R::pnorm(v / sv, 0.0, 1.0, 1, 1), std::log(1e-10));
}

// log positive-drift LBA density of reaching threshold bs at t.
inline double mtlba_logf(double t, double A, double bs, double v, double sv)
{
  const double lden = mtlba_pmax(R::pnorm(v / sv, 0.0, 1.0, 1, 1), std::log(1e-10));
  if (A <= 1e-10)
    return R::dnorm(bs / t, v, sv, 1) + std::log(bs) - 2.0 * std::log(t) - lden;
  const double cz  = (bs - t * v) / (t * sv);
  const double czm = cz - A / (t * sv);
  const double lp1 = R::dnorm(cz, 0.0, 1.0, 1);
  const double lp0 = R::dnorm(czm, 0.0, 1.0, 1);
  const double lpi = mtlnr_log_pnorm_interval(czm, cz);
  const double M   = mtlba_pmax(mtlba_pmax(lp1, lp0), lpi);   // no term overflows
  const double dP  = std::exp(lpi - M);
  const double inner = v * dP + sv * (std::exp(lp0 - M) - std::exp(lp1 - M));
  const double lo = (inner > 0) ? std::log(inner) + M : R_NegInf;  // NaN -> -Inf
  return lo - std::log(A) - lden;
}

// ---- correlated drift rates (rho != 0) ----

// R's pmin for two scalars: NaN if either is NaN.
inline double mtlba_pmin(double a, double b)
{
  if (std::isnan(a) || std::isnan(b)) return R_NaN;
  return a < b ? a : b;
}

// 16-node Gauss-Legendre rule on (-1, 1): the R twin's mtlba_gl(16), with
// u = (x + 1) / 2 and log(w / 2) for a panel of unit width.
constexpr int MTLBA_GL_N = 16;
struct MTLBAGL {
  double u[MTLBA_GL_N], w2[MTLBA_GL_N], lw2[MTLBA_GL_N];
  MTLBAGL() {
    static const double x[MTLBA_GL_N] = {
      -0.98940093499164994, -0.9445750230732326, -0.86563120238783176, -0.755404408355003,
      -0.61787624440264377, -0.45801677765722737, -0.28160355077925892, -0.095012509837637441,
       0.095012509837637441, 0.28160355077925892, 0.45801677765722737, 0.61787624440264377,
       0.755404408355003, 0.86563120238783176, 0.9445750230732326, 0.98940093499164994};
    static const double w[MTLBA_GL_N] = {
      0.027152459411754058, 0.062253523938647776, 0.095158511682492897, 0.12462897125553395,
      0.14959598881657685, 0.16915651939500256, 0.18260341504492361, 0.1894506104550685,
      0.1894506104550685, 0.18260341504492361, 0.16915651939500256, 0.14959598881657685,
      0.12462897125553395, 0.095158511682492897, 0.062253523938647776, 0.027152459411754058};
    for (int i = 0; i < MTLBA_GL_N; ++i) {
      u[i] = (x[i] + 1.0) / 2.0;
      w2[i] = w[i] / 2.0;
      lw2[i] = std::log(w2[i]);
    }
  }
};
inline const MTLBAGL& mtlba_gl16()
{
  static const MTLBAGL g;
  return g;
}

// The fit rule of the correlated kernel (R: mtlba_quad_rule()) and the rule of
// the stopping-time terms (R: mtlba_quad_reference()).
constexpr double MTLBA_QUAD_HW   = 8.0;    // window half-width, sds of the stopping drift
constexpr double MTLBA_QUAD_PW   = 4.0;    // feature widths per panel
constexpr int    MTLBA_QUAD_MAXP = 64;
constexpr double MTLBA_STOP_HW   = 12.0;
constexpr double MTLBA_STOP_PW   = 1.0;
constexpr int    MTLBA_STOP_MAXP = 1024;

// Number of panels for a window of the given width: one per pw widths of the
// narrowest feature sig, 1 .. maxp; 0 if undefined.
inline int mtlba_quad_panels(double width, double sig, double pw, int maxp)
{
  const double q = std::ceil(width / (pw * sig));
  if (std::isnan(q)) return 0;
  if (q >= maxp) return maxp;
  return q < 1 ? 1 : (int)q;
}

// log of the correlated kernel times P(both drifts > 0): the stopping
// accumulator (A_s, threshold bs, v_s, sv_s) reaches bs at dt while the other
// (A_o, v_o, sv_o) lies in its band (lo_y, hi_y, lowest). As the R twin
// mtlba_logkernel_rho with the default rule: log-scale integrand over the
// stopping drift x (the band at the conditional mean and sd, unnormalised),
// window as mtlba_quad_window, summed by log-sum-exp. A_s = 0: no integral.
inline double mtlba_logkernel_rho(double dt, double As, double bs, double vs, double svs,
                                  double Ao, double vo, double svo,
                                  double lo_y, double hi_y, bool lowest, double rho)
{
  const double sc    = svo * std::sqrt(1.0 - rho * rho);
  const double slope = rho * svo / svs;
  auto lint = [&](double x) {
    const double mc = vo + slope * (x - vs);
    const double l = std::log(x) + R::dnorm(x, vs, svs, 1) +
                     mtlba_logband_un(dt, Ao, mc, sc, lo_y, hi_y, lowest);
    return std::isnan(l) ? R_NegInf : l;
  };
  if (As <= 1e-10) {
    const double l = lint(bs / dt) - std::log(dt);
    return std::isnan(l) ? R_NegInf : l;
  }
  // window: the range of x clipped to within r of the interval (v_s, m)
  const double lo_x = (bs - As) / dt;
  const double hi_x = bs / dt;
  const double yl = lowest ? 0.0 : mtlba_pmax(0.0, (lo_y - Ao) / dt);
  const double m  = vs + rho * svs / svo * (mtlba_pmin(mtlba_pmax(vo, yl), hi_y / dt) - vo);
  const double c_lo = mtlba_pmin(vs, m), c_hi = mtlba_pmax(vs, m);
  const double d = mtlba_pmax(mtlba_pmax(lo_x - c_hi, c_lo - hi_x), 0.0);
  const double r = std::sqrt(d * d + (MTLBA_QUAD_HW * svs) * (MTLBA_QUAD_HW * svs));
  const double a = mtlba_pmax(lo_x, c_lo - r);
  const double width = mtlba_pmin(hi_x, c_hi + r) - a;
  const double sig = svs * mtlba_pmin(1.0, std::sqrt(1.0 - rho * rho) / std::fabs(rho));
  const int P = mtlba_quad_panels(width, sig, MTLBA_QUAD_PW, MTLBA_QUAD_MAXP);
  if (P == 0 || !(width > 0)) return R_NegInf;
  const MTLBAGL& g = mtlba_gl16();
  const double h = width / P;
  double L[MTLBA_QUAD_MAXP * MTLBA_GL_N];
  double M = R_NegInf;
  int n = 0;
  for (int k = 0; k < P; ++k)
    for (int i = 0; i < MTLBA_GL_N; ++i) {
      const double l = lint(a + h * (g.u[i] + k)) + g.lw2[i];
      L[n++] = l;
      if (l > M) M = l;
    }
  if (M == R_NegInf) return R_NegInf;
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += std::exp(L[i] - M);
  const double out = M + std::log(sum) + std::log(h) - std::log(As);
  return std::isnan(out) ? R_NegInf : out;
}

// P(Z1 > h, Z2 > k) for a standard bivariate normal with correlation r,
// |r| < 1: Genz (2004) BVND, as the R twin mtlba_bvu.
inline double mtlba_bvu(double h, double k, double r)
{
  static const double W3[3]  = {0.1713244923791705, 0.3607615730481384, 0.4679139345726904};
  static const double X3[3]  = {-0.9324695142031522, -0.6612093864662647, -0.2386191860831970};
  static const double W6[6]  = {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                                0.2031674267230659, 0.2334925365383547, 0.2491470458134029};
  static const double X6[6]  = {-0.9815606342467191, -0.9041172563704750, -0.7699026741943050,
                                -0.5873179542866171, -0.3678314989981802, -0.1252334085114692};
  static const double W10[10] = {0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
                                 0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
                                 0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
                                 0.1527533871307259};
  static const double X10[10] = {-0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
                                 -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
                                 -0.5108670019508271, -0.3737060887154196, -0.2277858511416451,
                                 -0.07652652113349733};
  if (!std::isfinite(h) || !std::isfinite(k) || !std::isfinite(r) || !(std::fabs(r) < 1)) return R_NaN;
  const double ar = std::fabs(r);
  const double* W = ar < 0.3 ? W3 : (ar < 0.75 ? W6 : W10);
  const double* X = ar < 0.3 ? X3 : (ar < 0.75 ? X6 : X10);
  const int ng = ar < 0.3 ? 3 : (ar < 0.75 ? 6 : 10);
  const double TWO_PI = 2.0 * M_PI;
  if (ar < 0.925) {
    const double hk = h * k, hs = (h * h + k * k) / 2.0, asr = std::asin(r);
    double sum = 0.0;
    for (int sgn = 0; sgn < 2; ++sgn)
      for (int i = 0; i < ng; ++i) {
        const double x  = sgn ? -X[i] : X[i];
        const double sn = std::sin(asr * ((x + 1.0) / 2.0));
        sum += std::exp((sn * hk - hs) / (1.0 - sn * sn)) * W[i];
      }
    return sum * asr / (2.0 * TWO_PI) + R::pnorm(-h, 0.0, 1.0, 1, 0) * R::pnorm(-k, 0.0, 1.0, 1, 0);
  }
  const bool neg = r < 0;
  const double kk = neg ? -k : k;
  const double hk = h * kk;
  const double as = (1.0 - r) * (1.0 + r);
  const double a  = std::sqrt(as);
  const double bs = (h - kk) * (h - kk);
  const double cc = (4.0 - hk) / 8.0;
  const double d  = (12.0 - hk) / 16.0;
  const double asr = -(bs / as + hk) / 2.0;
  double bvn = (asr > -100) ?
    a * std::exp(asr) * (1.0 - cc * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 + cc * d * as * as / 5.0) : 0.0;
  const double b = std::sqrt(bs);
  if (-hk < 100)
    bvn -= std::exp(-hk / 2.0) * std::sqrt(TWO_PI) * R::pnorm(-b / a, 0.0, 1.0, 1, 0) * b *
           (1.0 - cc * bs * (1.0 - d * bs / 5.0) / 3.0);
  double sum = 0.0;
  for (int sgn = 0; sgn < 2; ++sgn)
    for (int i = 0; i < ng; ++i) {
      const double x  = sgn ? -X[i] : X[i];
      const double q  = a / 2.0 * (1.0 + x);
      const double xs = q * q;
      const double rs = std::sqrt(1.0 - xs);
      const double as2 = -(bs / xs + hk) / 2.0;
      if (as2 > -100)
        sum += W[i] * std::exp(as2) *
               (std::exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs - (1.0 + cc * xs * (1.0 + d * xs)));
    }
  bvn = -(bvn + (a / 2.0) * sum) / TWO_PI;
  if (neg)
    return -bvn + (kk > h ? R::pnorm(kk, 0.0, 1.0, 1, 0) - R::pnorm(h, 0.0, 1.0, 1, 0) : 0.0);
  return bvn + R::pnorm(-(h > kk ? h : kk), 0.0, 1.0, 1, 0);
}

// log P(both drifts > 0)
inline double mtlba_logppos(double v1, double sv1, double v2, double sv2, double rho)
{
  return std::log(mtlba_pmax(mtlba_bvu(-v1 / sv1, -v2 / sv2, rho), 0.0));
}

// Correlated twin of mtlba_logstop for one parameter row and a finite bound
// y > 0, as the R twin mtlba_stop_rho1: one integral over accumulator 1's
// drift x > 0, composite Gauss-Legendre between the kinks of accumulator 1's
// crossing counts, all terms positive, over P(both > 0). T1, T2: thresholds
// 1 .. K of the two accumulators. x, w: work space.
inline void mtlba_stop_rho1(double y, double A1, double v1, double s1, const double* T1,
                            double A2, double v2, double s2, const double* T2,
                            double rho, int K, bool TC, double& lP, double& lS,
                            std::vector<double>& x, std::vector<double>& w)
{
  lP = lS = R_NaN;
  bool fin = std::isfinite(y) && std::isfinite(A1) && std::isfinite(v1) && std::isfinite(s1) &&
             std::isfinite(A2) && std::isfinite(v2) && std::isfinite(s2) && std::isfinite(rho);
  for (int n = 0; fin && n < K; ++n) fin = std::isfinite(T1[n]) && std::isfinite(T2[n]);
  if (!fin || !(std::fabs(rho) < 1) || !(s1 > 0) || !(s2 > 0)) return;
  // window: v2 > 0 drags x towards m when v2's mean is negative
  const double m  = v1 + rho * s1 / s2 * ((v2 > 0 ? v2 : 0.0) - v2);
  const double hi = v1 > m ? v1 : m, lw = v1 < m ? v1 : m;
  const double d  = -hi > 0 ? -hi : 0.0;
  const double r  = std::sqrt(d * d + (MTLBA_STOP_HW * s1) * (MTLBA_STOP_HW * s1));
  const double a  = (lw - r > 0) ? lw - r : 0.0;
  const double b  = hi + r;
  double brk[2 * MTLBA_MAX_K + 2];
  int nb = 0;
  brk[nb++] = a; brk[nb++] = b;
  for (int n = 0; n < K; ++n)
    for (double kink : {(T1[n] - A1) / y, T1[n] / y})
      if (kink > a && kink < b) brk[nb++] = kink;
  std::sort(brk, brk + nb);
  nb = (int)(std::unique(brk, brk + nb) - brk);
  const double sig = s1 * mtlba_pmin(1.0, std::sqrt(1.0 - rho * rho) / std::fabs(rho));
  const MTLBAGL& g = mtlba_gl16();
  x.clear(); w.clear();
  for (int i = 0; i + 1 < nb; ++i) {
    const double width = brk[i + 1] - brk[i];
    const int p = mtlba_quad_panels(width, sig, MTLBA_STOP_PW, MTLBA_STOP_MAXP);
    const double h = width / p;
    for (int k = 0; k < p; ++k)
      for (int q = 0; q < MTLBA_GL_N; ++q) {
        x.push_back(brk[i] + h * (g.u[q] + k));
        w.push_back(h * g.w2[q]);
      }
  }
  const double sc = s2 * std::sqrt(1.0 - rho * rho);
  const double slope = rho * s2 / s1;
  double sumP = 0.0, sumS = 0.0;
  double Fle[MTLBA_MAX_K + 1];
  Fle[0] = 0.0;
  const size_t nx = x.size();
  for (size_t i = 0; i < nx; ++i) {
    const double xi = x[i];
    const double wi = w[i] * R::dnorm(xi, v1, s1, 0);
    const double mc = v2 + slope * (xi - v1);
    // P(e_1 < T_{1,n} | x), n = 1 .. K
    for (int n = 1; n <= K; ++n) {
      const double T = T1[n - 1];
      if (A1 <= 1e-10) Fle[n] = (T > xi * y) ? 1.0 : 0.0;
      else {
        const double f = (T - xi * y) / A1;
        Fle[n] = f < 0 ? 0.0 : (f > 1 ? 1.0 : f);
      }
    }
    double Sint = 0.0, Pint = (1.0 - Fle[K]) * R::pnorm(mc / sc, 0.0, 1.0, 1, 0);
    for (int n1 = TC ? 0 : K - 1; n1 < K; ++n1) {
      // accumulator 1 at n1 crossings (BE: below b), accumulator 2 against T_{2,mm}
      const double wn = TC ? Fle[n1 + 1] - Fle[n1] : Fle[K];
      if (wn == 0) continue;
      const int mm = TC ? K - n1 : K;
      Sint += wn * std::exp(mtlba_logband_un(y, A2, mc, sc, 0.0, T2[mm - 1], true));
      Pint += wn * std::exp(mtlba_logS(T2[mm - 1], y, A2, mc, sc));
    }
    sumP += wi * Pint;
    sumS += wi * Sint;
  }
  const double lpp = mtlba_logppos(v1, s1, v2, s2, rho);
  lP = std::log(sumP) - lpp;
  lS = std::log(sumS) - lpp;
}

// Per-trial column pointers
struct MTLBACols {
  const double *v, *sv, *B, *A, *t0, *rho;
  std::vector<const double*> c;
};

// Threshold j (1 .. K; K = b) of the accumulator in row r.
inline double mtlba_threshold(const MTLBACols& p, int K, int r, int j)
{
  const double A = p.A[r];
  const double b = p.B[r] + A;
  if (j >= K) return b;
  double C = 0.0;
  for (int k = 0; k < K - j; ++k) C += p.c[k][r];
  return A + (b - A) * std::exp(-C);
}

// log P(both > 0) is the same for all trials of a design cell: a small cache,
// local to one likelihood call, keyed on the exact (v_s, sv_s, v_o, sv_o, rho).
struct MTLBAPposCache {
  static constexpr int MAX_ENTRIES = 16;
  int n = 0, next = 0;
  double keys[MAX_ENTRIES][5], vals[MAX_ENTRIES];
  double get(double vs, double svs, double vo, double svo, double rho) {
    for (int e = 0; e < n; ++e) {
      const double* k = keys[e];
      if (k[0] == vs && k[1] == svs && k[2] == vo && k[3] == svo && k[4] == rho) return vals[e];
    }
    const double val = mtlba_pmax(mtlba_logppos(vs, svs, vo, svo, rho), std::log(1e-10));
    const int e = (n < MAX_ENTRIES) ? n++ : next++ % MAX_ENTRIES;
    keys[e][0] = vs; keys[e][1] = svs; keys[e][2] = vo; keys[e][3] = svo; keys[e][4] = rho;
    vals[e] = val;
    return val;
  }
};

// log kernel: accumulator in row rs stops at its threshold j while the one in
// row ro lies in band lo; rho the trial's drift-rate correlation (0: closed
// form). -Inf for undefined values (non-finite parameters, dt <= 0, rho not
// inside (-1, 1)).
inline double mtlba_logkernel(const MTLBACols& p, int K, double rt, int rs, int ro, int j, int lo,
                              double rho, MTLBAPposCache& ppos)
{
  const double dt = rt - p.t0[rs];
  if (!(std::isfinite(dt) && dt > 0)) return R_NegInf;
  for (int r : {rs, ro})
    if (!std::isfinite(p.v[r]) || !std::isfinite(p.sv[r]) || !std::isfinite(p.A[r]) ||
        !std::isfinite(p.B[r] + p.A[r])) return R_NegInf;
  const double Ts = mtlba_threshold(p, K, rs, j);
  const double hi = mtlba_threshold(p, K, ro, lo + 1);
  const double ly = (lo == 0) ? 0.0 : mtlba_threshold(p, K, ro, lo);
  if (!std::isfinite(Ts) || !std::isfinite(hi) || !std::isfinite(ly)) return R_NegInf;
  if (!(std::fabs(rho) < 1)) return R_NegInf;   // NaN too
  if (rho == 0)
    return mtlba_logf(dt, p.A[rs], Ts, p.v[rs], p.sv[rs]) +
           mtlba_logband(dt, p.A[ro], p.v[ro], p.sv[ro], ly, hi, lo == 0);
  return mtlba_logkernel_rho(dt, p.A[rs], Ts, p.v[rs], p.sv[rs], p.A[ro], p.v[ro], p.sv[ro],
                             ly, hi, lo == 0, rho) -
         ppos.get(p.v[rs], p.sv[rs], p.v[ro], p.sv[ro], rho);
}

// Response lookup -> (stopping accumulator 1/2, threshold j, other's band lo).
inline void mtlba_lookup(int R, int RR, int K, bool TC, int& s, int& j, int& lo)
{
  if (!TC) {
    s = R; j = K; lo = K - RR;
    return;
  }
  const int p = (R == 1) ? K - RR + 1 : K + RR;
  if (p % 2 == 1) { s = 1; j = K - (p - 1) / 2; }
  else            { s = 2; j = p / 2; }
  lo = K - j;
}

// log P(dt <= y) (lP) and log P(dt > y) (lS) for the race's stopping time at
// decision-time bound y; rows r1 / r2 = accumulators 1 / 2. As the R twin
// mtlba_logstop: with n_a the thresholds crossed by accumulator a at y,
// F_a(n) = P(n_a <= n), Fc_a(n) = 1 - F_a(n), p_a(n) = P(n_a = n) (v_a > 0),
//   BE: S = F_1(K-1) F_2(K-1),   P = Fc_1(K-1) + F_1(K-1) Fc_2(K-1)
//   TC: S = sum_{n1<K} p_1(n1) F_2(K-1-n1),
//       P = Fc_1(K-1) + sum_{n1<K} p_1(n1) Fc_2(K-1-n1),
// all terms positive (log-sum-exp). NaN for undefined parameters.
inline void mtlba_logstop(const MTLBACols& p, int K, bool TC, double y, int r1, int r2,
                          double& lP, double& lS)
{
  if (std::isnan(y)) { lP = lS = R_NaN; return; }
  if (y == R_PosInf) { lP = 0.0; lS = R_NegInf; return; }
  if (!(y > 0)) { lP = R_NegInf; lS = 0.0; return; }
  double F[2][MTLBA_MAX_K], Fc[2][MTLBA_MAX_K], band[2][MTLBA_MAX_K];
  const int rows[2] = {r1, r2};
  for (int a = 0; a < 2; ++a) {
    const int r = rows[a];
    const double A = p.A[r], v = p.v[r], sv = p.sv[r];
    if (!std::isfinite(v) || !std::isfinite(sv) || !std::isfinite(A) ||
        !std::isfinite(p.B[r] + A)) { lP = lS = R_NaN; return; }
    const double lden = mtlba_pmax(R::pnorm(v / sv, 0.0, 1.0, 1, 1), std::log(1e-10));
    double lo_y = 0.0;
    for (int nn = 0; nn < K; ++nn) {
      const double hi = mtlba_threshold(p, K, r, nn + 1);
      if (!std::isfinite(hi)) { lP = lS = R_NaN; return; }
      F[a][nn]  = mtlba_logband(y, A, v, sv, 0.0, hi, true);
      Fc[a][nn] = mtlba_logS(hi, y, A, v, sv) - lden;
      band[a][nn] = (nn == 0) ? F[a][0] : mtlba_logband(y, A, v, sv, lo_y, hi, false);
      lo_y = hi;
    }
  }
  if (!TC) {
    lS = F[0][K - 1] + F[1][K - 1];
    lP = mtlba_lse(Fc[0][K - 1], F[0][K - 1] + Fc[1][K - 1]);
  } else {
    lS = band[0][0] + F[1][K - 1];
    lP = mtlba_lse(Fc[0][K - 1], band[0][0] + Fc[1][K - 1]);
    for (int n1 = 1; n1 < K; ++n1) {
      lS = mtlba_lse(lS, band[0][n1] + F[1][K - 1 - n1]);
      lP = mtlba_lse(lP, band[0][n1] + Fc[1][K - 1 - n1]);
    }
  }
}

// log P(a < dt <= b) from the (lP, lS) pairs at a and b: the route with the
// smaller leading term.
inline double mtlba_logwindow(double Pa, double Sa, double Pb, double Sb)
{
  if (Sa < Pb) return mtlba_ldiff(Sa, Sb);   // false for NaN
  return mtlba_ldiff(Pb, Pa);
}

// The stopping-time terms at LT and UT depend only on the trial's parameter
// rows and the two bounds, not on the rt, so trials of the same design cell
// share them. A small cache, local to one likelihood call, keyed on the exact
// values (v, sv, A, B, t0, criteria of both rows, LT, UT); NaN never matches.
struct MTLBAStopCache {
  static constexpr int MAX_ENTRIES = 64;
  int n_key = 0, n = 0, next = 0;
  std::vector<double> keys, vals;   // MAX_ENTRIES x n_key, MAX_ENTRIES x 4
  std::vector<double> key;
  void init(int K) {
    n_key = 2 * (5 + K - 1) + 2;
    keys.resize((size_t)MAX_ENTRIES * n_key);
    vals.resize((size_t)MAX_ENTRIES * 4);
    key.resize(n_key);
  }
  void get(const MTLBACols& p, int K, bool TC, int r1, int r2, double LT, double UT,
           double& PLT, double& SLT, double& PUT, double& SUT) {
    const double t0 = p.t0[r1];
    int k = 0;
    for (int r : {r1, r2}) {
      key[k++] = p.v[r]; key[k++] = p.sv[r]; key[k++] = p.A[r]; key[k++] = p.B[r];
      key[k++] = p.t0[r];
      for (int c = 0; c < K - 1; ++c) key[k++] = p.c[c][r];
    }
    key[k++] = LT; key[k++] = UT;
    for (int e = 0; e < n; ++e) {
      const double* ke = keys.data() + (size_t)e * n_key;
      int i = 0;
      while (i < n_key && ke[i] == key[i]) ++i;
      if (i == n_key) {
        const double* v = vals.data() + (size_t)e * 4;
        PLT = v[0]; SLT = v[1]; PUT = v[2]; SUT = v[3];
        return;
      }
    }
    mtlba_logstop(p, K, TC, LT - t0, r1, r2, PLT, SLT);
    mtlba_logstop(p, K, TC, UT - t0, r1, r2, PUT, SUT);
    const int e = (n < MAX_ENTRIES) ? n++ : next++ % MAX_ENTRIES;   // round-robin when full
    std::copy(key.begin(), key.end(), keys.begin() + (size_t)e * n_key);
    double* v = vals.data() + (size_t)e * 4;
    v[0] = PLT; v[1] = SLT; v[2] = PUT; v[3] = SUT;
  }
};

// Correlated stopping-time terms (rho != 0): log P(dt <= y) and log P(dt > y)
// at one bound y, as the R twin mtlba_logstop_any. The integral depends only on
// the parameter rows and y, so it is cached per likelihood call, keyed on the
// exact values (y, rho, then v, sv, A, B and the criteria of both rows).
struct MTLBAStopRhoCache {
  static constexpr int MAX_ENTRIES = 64;
  int n_key = 0, n = 0, next = 0;
  std::vector<double> keys, vals;   // MAX_ENTRIES x n_key, MAX_ENTRIES x 2
  std::vector<double> key, x, w;
  void init(int K) {
    n_key = 2 + 2 * (4 + K - 1);
    keys.resize((size_t)MAX_ENTRIES * n_key);
    vals.resize((size_t)MAX_ENTRIES * 2);
    key.resize(n_key);
  }
  void get(const MTLBACols& p, int K, bool TC, double y, int r1, int r2, double rho,
           double& lP, double& lS) {
    if (std::isnan(y)) { lP = lS = R_NaN; return; }
    if (y == R_PosInf) { lP = 0.0; lS = R_NegInf; return; }
    if (!(y > 0)) { lP = R_NegInf; lS = 0.0; return; }
    if (n_key == 0) init(K);
    int k = 0;
    key[k++] = y; key[k++] = rho;
    for (int r : {r1, r2}) {
      key[k++] = p.v[r]; key[k++] = p.sv[r]; key[k++] = p.A[r]; key[k++] = p.B[r];
      for (int c = 0; c < K - 1; ++c) key[k++] = p.c[c][r];
    }
    for (int e = 0; e < n; ++e) {
      const double* ke = keys.data() + (size_t)e * n_key;
      int i = 0;
      while (i < n_key && ke[i] == key[i]) ++i;
      if (i == n_key) { lP = vals[2 * e]; lS = vals[2 * e + 1]; return; }
    }
    double T1[MTLBA_MAX_K], T2[MTLBA_MAX_K];
    for (int j = 1; j <= K; ++j) {
      T1[j - 1] = mtlba_threshold(p, K, r1, j);
      T2[j - 1] = mtlba_threshold(p, K, r2, j);
    }
    mtlba_stop_rho1(y, p.A[r1], p.v[r1], p.sv[r1], T1, p.A[r2], p.v[r2], p.sv[r2], T2,
                    rho, K, TC, lP, lS, x, w);
    const int e = (n < MAX_ENTRIES) ? n++ : next++ % MAX_ENTRIES;   // round-robin when full
    std::copy(key.begin(), key.end(), keys.begin() + (size_t)e * n_key);
    vals[2 * e] = lP; vals[2 * e + 1] = lS;
  }
};

// Return codes of mtlba_trial_ll (the caller stops, outside any parallel region).
enum { MTLBA_OK = 0, MTLBA_T0 = 2 };
inline const char* mtlba_error_message(int code)
{
  (void)code;
  return "MTLBA: censoring and truncation need t0 equal for both accumulators (t0 must not depend on lR or lM)";
}

// Trial log-likelihoods into ll_buf (indexed by trial), floored at min_ll.
// Observed trials: the response kernel; censored trials: the probability of
// stopping in the censored window(s); truncated trials: minus log P(LT < rt < UT).
// rho is read from accumulator 1's row (a trial-level parameter). Returns
// MTLBA_T0 if a censored or truncated trial has unequal t0 (the caller stops,
// outside any parallel region).
inline int mtlba_trial_ll(const MTLBASpec& s, const ParamTable& pt, double min_ll,
                          double* __restrict__ ll_buf)
{
  MTLBACols p;
  p.v   = pt.base.colptr(s.col_v);
  p.sv  = pt.base.colptr(s.col_sv);
  p.B   = pt.base.colptr(s.col_B);
  p.A   = pt.base.colptr(s.col_A);
  p.t0  = pt.base.colptr(s.col_t0);
  p.rho = pt.base.colptr(s.col_rho);
  for (int col : s.col_c) p.c.push_back(pt.base.colptr(col));
  const int K = s.K;
  int code = MTLBA_OK;
  MTLBAStopCache cache;
  MTLBAStopRhoCache rcache;
  MTLBAPposCache ppos;
  if (s.any_missing) cache.init(K);

  const int nt = (int)s.rt.size();
  for (int t = 0; t < nt; ++t) {
    const int r1 = s.row1[t];
    const int r2 = (r1 % 2 == 0) ? r1 + 1 : r1 - 1;
    const double rho = p.rho[r1];
    const int miss = s.miss[t];
    double ll = R_NegInf;
    if (miss == 0) {
      const int RR = s.RR[t];
      int a, j, lo;
      if (s.DK && RR == 1) {
        mtlba_lookup(1, 1, K, s.TC, a, j, lo);
        const double l1 = mtlba_logkernel(p, K, s.rt[t], a == 1 ? r1 : r2, a == 1 ? r2 : r1, j, lo, rho, ppos);
        mtlba_lookup(2, 1, K, s.TC, a, j, lo);
        const double l2 = mtlba_logkernel(p, K, s.rt[t], a == 1 ? r1 : r2, a == 1 ? r2 : r1, j, lo, rho, ppos);
        ll = mtlba_lse(l1, l2);
      } else {
        mtlba_lookup(s.R[t], RR, K, s.TC, a, j, lo);
        ll = mtlba_logkernel(p, K, s.rt[t], a == 1 ? r1 : r2, a == 1 ? r2 : r1, j, lo, rho, ppos);
      }
    }
    if (miss != 0 || s.trunc[t]) {
      const double t0 = p.t0[r1];
      if (!std::isnan(t0) && !std::isnan(p.t0[r2]) && t0 != p.t0[r2]) {
        if (code == MTLBA_OK) code = MTLBA_T0;
        ll_buf[t] = min_ll;
        continue;
      }
      const bool corr = !(rho == 0);   // NaN: the correlated route, which returns NaN
      auto stop_at = [&](double y, double& P, double& S) {
        if (corr) rcache.get(p, K, s.TC, y, r1, r2, rho, P, S);
        else mtlba_logstop(p, K, s.TC, y, r1, r2, P, S);
      };
      double PLT, SLT, PUT, SUT;
      if (corr) {
        stop_at(s.LT[t] - t0, PLT, SLT);
        stop_at(s.UT[t] - t0, PUT, SUT);
      } else {
        cache.get(p, K, s.TC, r1, r2, s.LT[t], s.UT[t], PLT, SLT, PUT, SUT);
      }
      if (miss != 0) {
        double P, S;
        if (miss == 1 || miss == 3) {
          stop_at(s.LC[t] - t0, P, S);
          ll = mtlba_logwindow(PLT, SLT, P, S);
        }
        if (miss == 2 || miss == 3) {
          stop_at(s.UC[t] - t0, P, S);
          ll = mtlba_lse(ll, mtlba_logwindow(P, S, PUT, SUT));
        }
      }
      if (s.trunc[t]) {
        const double lZ = mtlba_logwindow(PLT, SLT, PUT, SUT);
        ll = (lZ == R_NegInf) ? R_NegInf : ll - lZ;
      }
    }
    ll_buf[t] = (ll > min_ll) ? ll : min_ll;   // NaN -> min_ll
  }
  return code;
}

#endif // MODEL_MTLBA_H
