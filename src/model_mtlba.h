#ifndef MODEL_MTLBA_H
#define MODEL_MTLBA_H

// Multiple-threshold LBA (MTLBA), uncorrelated (rho = 0). C++ twin of
// R/model_MTLBA.R; see there for the model and the numerical method.
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

// log H(y) and log S(y) for evidence U(0, A) + v t, v ~ N(v, sv) untruncated,
// y >= A (separately: the band needs only one route).
inline double mtlba_logH(double y, double t, double A, double v, double sv)
{
  const double z1 = (y - t * v) / (t * sv);
  if (A <= 1e-10) return R::pnorm(z1, 0.0, 1.0, 1, 1);
  const double al = A / (t * sv);
  return mtlba_ldiff(mtlba_logG(z1), mtlba_logG(z1 - al)) - std::log(al);
}
inline double mtlba_logS(double y, double t, double A, double v, double sv)
{
  const double z1 = (y - t * v) / (t * sv);
  if (A <= 1e-10) return R::pnorm(z1, 0.0, 1.0, 0, 1);
  const double al = A / (t * sv);
  return mtlba_ldiff(mtlba_logG(-(z1 - al)), mtlba_logG(-z1)) - std::log(al);
}

// log P(lo_y < e < hi_y | v > 0) at decision time t; lowest: no threshold below.
// The route (S or H) is chosen by the smaller leading term, S(lo) vs H(hi).
inline double mtlba_logband(double t, double A, double v, double sv,
                            double lo_y, double hi_y, bool lowest)
{
  const double zv  = v / sv;
  const double lpv = R::pnorm(zv, 0.0, 1.0, 1, 1);
  const double hH  = mtlba_logH(hi_y, t, A, v, sv);
  const double lS  = lowest ? lpv : mtlba_logS(lo_y, t, A, v, sv);
  double lb;
  if (lS < hH) {
    lb = mtlba_ldiff(lS, mtlba_logS(hi_y, t, A, v, sv));
  } else {
    const double lH = lowest ? R::pnorm(-zv, 0.0, 1.0, 1, 1) : mtlba_logH(lo_y, t, A, v, sv);
    lb = mtlba_ldiff(hH, lH);
  }
  return lb - mtlba_pmax(lpv, std::log(1e-10));
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

// log kernel: accumulator in row rs stops at its threshold j while the one in
// row ro lies in band lo. -Inf for undefined values (non-finite parameters,
// dt <= 0).
inline double mtlba_logkernel(const MTLBACols& p, int K, double rt, int rs, int ro, int j, int lo)
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
  return mtlba_logf(dt, p.A[rs], Ts, p.v[rs], p.sv[rs]) +
         mtlba_logband(dt, p.A[ro], p.v[ro], p.sv[ro], ly, hi, lo == 0);
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

// Return codes of mtlba_trial_ll (the caller stops, outside any parallel region).
enum { MTLBA_OK = 0, MTLBA_RHO = 1, MTLBA_T0 = 2 };
inline const char* mtlba_error_message(int code)
{
  if (code == MTLBA_T0)
    return "MTLBA: censoring and truncation need t0 equal for both accumulators (t0 must not depend on lR or lM)";
  return "MTLBA: correlated drift rates (rho != 0) are not implemented yet";
}

// Trial log-likelihoods into ll_buf (indexed by trial), floored at min_ll.
// Observed trials: the response kernel; censored trials: the probability of
// stopping in the censored window(s); truncated trials: minus log P(LT < rt < UT).
// Returns MTLBA_RHO if a non-zero rho is met, MTLBA_T0 if a censored or
// truncated trial has unequal t0 (the caller stops, outside any parallel
// region).
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
  if (s.any_missing) cache.init(K);

  const int nt = (int)s.rt.size();
  for (int t = 0; t < nt; ++t) {
    const int r1 = s.row1[t];
    const int r2 = (r1 % 2 == 0) ? r1 + 1 : r1 - 1;
    if ((p.rho[r1] != 0 && !std::isnan(p.rho[r1])) || (p.rho[r2] != 0 && !std::isnan(p.rho[r2])))
      code = MTLBA_RHO;
    const int miss = s.miss[t];
    double ll = R_NegInf;
    if (miss == 0) {
      const int RR = s.RR[t];
      int a, j, lo;
      if (s.DK && RR == 1) {
        mtlba_lookup(1, 1, K, s.TC, a, j, lo);
        const double l1 = mtlba_logkernel(p, K, s.rt[t], a == 1 ? r1 : r2, a == 1 ? r2 : r1, j, lo);
        mtlba_lookup(2, 1, K, s.TC, a, j, lo);
        const double l2 = mtlba_logkernel(p, K, s.rt[t], a == 1 ? r1 : r2, a == 1 ? r2 : r1, j, lo);
        ll = mtlba_lse(l1, l2);
      } else {
        mtlba_lookup(s.R[t], RR, K, s.TC, a, j, lo);
        ll = mtlba_logkernel(p, K, s.rt[t], a == 1 ? r1 : r2, a == 1 ? r2 : r1, j, lo);
      }
    }
    if (miss != 0 || s.trunc[t]) {
      const double t0 = p.t0[r1];
      if (!std::isnan(t0) && !std::isnan(p.t0[r2]) && t0 != p.t0[r2]) {
        if (code == MTLBA_OK) code = MTLBA_T0;
        ll_buf[t] = min_ll;
        continue;
      }
      double PLT, SLT, PUT, SUT;
      cache.get(p, K, s.TC, r1, r2, s.LT[t], s.UT[t], PLT, SLT, PUT, SUT);
      if (miss != 0) {
        double P, S;
        if (miss == 1 || miss == 3) {
          mtlba_logstop(p, K, s.TC, s.LC[t] - t0, r1, r2, P, S);
          ll = mtlba_logwindow(PLT, SLT, P, S);
        }
        if (miss == 2 || miss == 3) {
          mtlba_logstop(p, K, s.TC, s.UC[t] - t0, r1, r2, P, S);
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
