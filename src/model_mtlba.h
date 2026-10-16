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
// Drift convention (c_options drifts_any): "any" (1, the default) conditions
// the bivariate normal drift pair on not both being negative, normaliser
// P(not both < 0), untruncated band probabilities for the other accumulator;
// "positive" (0) truncates both drifts to be positive, normaliser
// P(both > 0). A free lowest threshold (c_options free_lowest, balance of
// evidence only) is T_1 = T_2 exp(-c_{K-1}), possibly below A.
//
// Proportional parameterisation (c_options proportional = 1; absent: 0, the
// additive model above): the columns B, A and c1 .. c{K-1} hold b itself and
// proportions in (0, 1), T_{K-1} = b c1, A = T_{K-1} A (K = 1: A = b A), each
// lower threshold nested on the one above, T_j = A + (T_{j+1} - A) c_{K-j}, or
// T_1 = T_2 c_{K-1} when the lowest threshold is free. mtlba_trial_ll
// translates these rows exactly into the additive columns (natural A,
// B = b - A, criteria) before anything else reads them (mtlba_proportional),
// so the kernel below only ever sees the additive model.
//
// The band is formed on the log scale: under "positive" with the sentinels
// H(lowest) = Phi(-v/sv), S(lowest) = Phi(v/sv),
// band * Phi(v/sv) = S(lo) - S(hi) = H(hi) - H(lo) (a threshold below A takes
// the positive-drift forms Hp, Sp); under "any" H(lowest) = 0, S(lowest) = 1
// and the untruncated H, S for any threshold. Either way the route with the
// smaller leading term is taken; log G by the Mills ratio (z < 0) and
// Laplace's continued fraction (z < -5, 40 terms). Log tails use
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
//
// Confidence guessing (parameters gp, gB, pg1 .. pg{K-1}; see the R twin's
// "confidence guessing" block): an accumulator that is triggering ends the race
// at its guess threshold gb = A + B gB with a random rating. The accumulators
// trigger independently with their own gp (the default), or the whole trial is
// a guess with one gp (c_options guess_trial, BETE's mixture). A response is
// then a weighted sum of up to four kernels, the censoring/truncation terms a
// mixture over the trigger states with a triggering accumulator's thresholds
// capped at gb. Rows with gp = 0 take the plain kernel. Binomial guess ratings
// (c_options binomial_guess = 1; absent: 0, free pg): one column pq holds q,
// the guessed rating is 1 + Binomial(K - 1, q), and mtlba_trial_ll first
// fills scratch columns with the sequential pg_k = P(k) / P(>= k) that
// p.pg[] then points at (mtlba_binomial_pg).

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
  int  col_gp  = -1;       // confidence guessing: -1 if the model has none
  int  col_gB  = -1;
  bool guess_trial = false;    // one trigger per trial (else per accumulator)
  int  gb_constraint = 0;      // guess threshold moves to b with gp: 0 none, 1 own gp, 2 max gp
  bool any         = false;    // drifts "any" (not both negative) rather than "positive"
  bool free_lowest = false;    // lowest threshold T_1 = T_2 exp(-c_{K-1}) (may be below A)
  bool proportional = false;   // B, A, c columns are b and proportions (translated per call)
  bool binomial_guess = false; // guessed ratings 1 + Binomial(K - 1, q), q in column pq
  int  col_pq  = -1;
  std::vector<int>    col_c;   // c1 .. c{K-1}
  std::vector<int>    col_pg;  // pg1 .. pg{K-1}
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
  // confidence guessing (absent in models built before it existed: no guessing)
  auto find_col = [&](const std::string& x) {
    auto it = pt.name_to_base_idx.find(x);
    return it == pt.name_to_base_idx.end() ? -1 : (int)it->second;
  };
  s.col_gp = find_col("gp");
  s.col_gB = find_col("gB");
  for (int k = 1; k < s.K; ++k) {
    const int c = find_col("pg" + std::to_string(k));
    if (c < 0) break;
    s.col_pg.push_back(c);
  }
  s.col_pq = find_col("pq");

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
  for (int i = 0; i < on.size(); ++i) {
    if (on[i] == "guess_trial") s.guess_trial = o[i] != 0;
    if (on[i] == "gb_constraint") s.gb_constraint = (int)o[i];
    if (on[i] == "drifts_any") s.any = o[i] != 0;          // absent (older models): "positive"
    if (on[i] == "free_lowest") s.free_lowest = o[i] != 0;
    if (on[i] == "proportional") s.proportional = o[i] != 0;   // absent (older models): additive
    if (on[i] == "binomial_guess") s.binomial_guess = o[i] != 0;   // absent (older models): free pg
  }
  if (s.col_gp >= 0 && s.binomial_guess && (s.col_gB < 0 || s.col_pq < 0))
    Rcpp::stop("MTLBA: binomial guessing needs gp, gB and pq");
  if (s.col_gp >= 0 && !s.binomial_guess && (s.col_gB < 0 || (int)s.col_pg.size() != s.K - 1))
    Rcpp::stop("MTLBA: guessing needs gp, gB and pg1 .. pg%d", s.K - 1);
  if (s.free_lowest && s.TC)
    Rcpp::stop("MTLBA: a free lowest threshold needs the balance-of-evidence rule");
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

// log P(e <= y, v > 0) and log P(e > y, v > 0) for ANY y >= 0 (a threshold
// below A under "positive" drifts), as the R twin mtlba_logHSp: for y >= A
// they are H(y) - Phi(-zv) and S(y); for y < A, with z1 = (y - t v) / (t sv),
// al = A / (t sv), every term positive,
//   Hp(y) = [G(z1) - G(-zv)] / al - (y / A) Phi(-zv)
//   Sp(y) = ((A - y) / A) Phi(zv) + [G(zv) - G(-z1)] / al.
// Hp cancels as y -> 0: a narrow one, dz = y / (t sv) with dz max(1, |z1|)
// < 0.1, is (y / A) times the 4-node Gauss-Legendre mean over u in (0, 1) of
// P(-zv < Z < z1 - dz u) (positive terms), as the R twin.
inline double mtlba_logHp(double y, double t, double A, double v, double sv)
{
  const double zv = v / sv;
  if (!(y < A) || !(A > 1e-10))
    return mtlba_ldiff(mtlba_logH(y, t, A, v, sv), R::pnorm(-zv, 0.0, 1.0, 1, 1));
  const double al = A / (t * sv);
  const double z1 = (y - t * v) / (t * sv);
  const double dz = y / (t * sv);
  if (dz * mtlba_pmax(1.0, std::fabs(z1)) < 0.1) {
    static const double U4[4]  = {0.069431844202973714, 0.33000947820757187,
                                  0.66999052179242813, 0.93056815579702623};
    static const double LW4[4] = {-1.7491171784568635, -1.1206352922365659,
                                  -1.1206352922365659, -1.7491171784568635};
    double L[4], M = R_NegInf;
    for (int i = 0; i < 4; ++i) {
      L[i] = mtlnr_log_pnorm_interval(-zv, mtlba_pmax(z1 - dz * U4[i], -zv)) + LW4[i];   // clamped: y = 0
      if (std::isnan(L[i])) L[i] = R_NegInf;   // a one-ulp interval: zero (as the R twin)
      if (L[i] > M) M = L[i];
    }
    if (M == R_NegInf) return R_NegInf;
    double sum = 0.0;
    for (int i = 0; i < 4; ++i) sum += std::exp(L[i] - M);
    return std::log(y / A) + M + std::log(sum);
  }
  return mtlba_ldiff(mtlba_ldiff(mtlba_logG(z1), mtlba_logG(-zv)) - std::log(al),
                     std::log(y / A) + R::pnorm(-zv, 0.0, 1.0, 1, 1));
}
inline double mtlba_logSp(double y, double t, double A, double v, double sv)
{
  if (!(y < A) || !(A > 1e-10)) return mtlba_logS(y, t, A, v, sv);
  const double zv = v / sv;
  const double al = A / (t * sv);
  const double z1 = (y - t * v) / (t * sv);
  return mtlba_lse(std::log((A - y) / A) + R::pnorm(zv, 0.0, 1.0, 1, 1),
                   mtlba_ldiff(mtlba_logG(zv), mtlba_logG(-z1)) - std::log(al));
}

// log probability of the band (lo_y, hi_y) at decision time t; lowest: no
// threshold below. pos: log P(lo_y < e < hi_y, v > 0) (the Phi(v/sv)
// sentinels; a threshold below A takes Hp / Sp); !pos: log P(lo_y < e < hi_y)
// for the untruncated drift (sentinels H = 0, S = 1). The route (S or H) is
// chosen by the smaller leading term, S(lo) vs H(hi).
inline double mtlba_logband_un(double t, double A, double v, double sv,
                               double lo_y, double hi_y, bool lowest, bool pos = true)
{
  if (!pos) {
    const double hH = mtlba_logH(hi_y, t, A, v, sv);
    const double lS = lowest ? 0.0 : mtlba_logS(lo_y, t, A, v, sv);
    if (lS < hH) return mtlba_ldiff(lS, mtlba_logS(hi_y, t, A, v, sv));
    const double lH = lowest ? R_NegInf : mtlba_logH(lo_y, t, A, v, sv);
    return mtlba_ldiff(hH, lH);
  }
  const double zv  = v / sv;
  if (hi_y < A || (!lowest && lo_y < A)) {   // false for NaN: the usual route
    const double hH = mtlba_logHp(hi_y, t, A, v, sv);
    const double lS = lowest ? R::pnorm(zv, 0.0, 1.0, 1, 1) : mtlba_logSp(lo_y, t, A, v, sv);
    if (lS < hH) return mtlba_ldiff(lS, mtlba_logSp(hi_y, t, A, v, sv));
    const double lH = lowest ? R_NegInf : mtlba_logHp(lo_y, t, A, v, sv);
    return mtlba_ldiff(hH, lH);
  }
  const double hH  = mtlba_logH(hi_y, t, A, v, sv);
  const double lS  = lowest ? R::pnorm(zv, 0.0, 1.0, 1, 1) : mtlba_logS(lo_y, t, A, v, sv);
  if (lS < hH) return mtlba_ldiff(lS, mtlba_logS(hi_y, t, A, v, sv));
  const double lH = lowest ? R::pnorm(-zv, 0.0, 1.0, 1, 1) : mtlba_logH(lo_y, t, A, v, sv);
  return mtlba_ldiff(hH, lH);
}

// pos: log P(lo_y < e < hi_y | v > 0), the band over Phi(v/sv), floored at
// 1e-10; !pos: the untruncated band (nothing to normalise by).
inline double mtlba_logband(double t, double A, double v, double sv,
                            double lo_y, double hi_y, bool lowest, bool pos = true)
{
  if (!pos) return mtlba_logband_un(t, A, v, sv, lo_y, hi_y, lowest, false);
  return mtlba_logband_un(t, A, v, sv, lo_y, hi_y, lowest) -
         mtlba_pmax(R::pnorm(v / sv, 0.0, 1.0, 1, 1), std::log(1e-10));
}

// log LBA density of reaching threshold bs at t from an untruncated normal
// drift; pos divides by P(v > 0) (the positive-drift density, as LBA()).
inline double mtlba_logf(double t, double A, double bs, double v, double sv, bool pos = true)
{
  const double lden = pos ? mtlba_pmax(R::pnorm(v / sv, 0.0, 1.0, 1, 1), std::log(1e-10)) : 0.0;
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

// log of the correlated kernel times the normaliser (P(both > 0) under
// "positive", P(not both < 0) under "any"): the stopping accumulator (A_s,
// threshold bs, v_s, sv_s) reaches bs at dt while the other (A_o, v_o, sv_o)
// lies in its band (lo_y, hi_y, lowest). As the R twin mtlba_logkernel_rho
// with the default rule: log-scale integrand over the stopping drift x (the
// band at the conditional mean and sd, unnormalised), window as
// mtlba_quad_window, summed by log-sum-exp. A_s = 0: no integral.
inline double mtlba_logkernel_rho(double dt, double As, double bs, double vs, double svs,
                                  double Ao, double vo, double svo,
                                  double lo_y, double hi_y, bool lowest, double rho,
                                  bool pos = true)
{
  const double sc    = svo * std::sqrt(1.0 - rho * rho);
  const double slope = rho * svo / svs;
  auto lint = [&](double x) {
    const double mc = vo + slope * (x - vs);
    const double l = std::log(x) + R::dnorm(x, vs, svs, 1) +
                     mtlba_logband_un(dt, Ao, mc, sc, lo_y, hi_y, lowest, pos);
    return std::isnan(l) ? R_NegInf : l;
  };
  if (As <= 1e-10) {
    const double l = lint(bs / dt) - std::log(dt);
    return std::isnan(l) ? R_NegInf : l;
  }
  // window: the range of x clipped to within r of the interval (v_s, m); the
  // other drift's lower limit is the band's own under "any" (no truncation)
  const double lo_x = (bs - As) / dt;
  const double hi_x = bs / dt;
  const double yl = pos ? (lowest ? 0.0 : mtlba_pmax(0.0, (lo_y - Ao) / dt))
                        : (lowest ? R_NegInf : (lo_y - Ao) / dt);
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

// log P(not both drifts < 0) = log [P(v1 > 0) + P(v1 < 0, v2 > 0)], two
// positive terms (the second is the upper orthant of (-Z1, Z2), correlation -rho).
inline double mtlba_lognotneg(double v1, double sv1, double v2, double sv2, double rho)
{
  return std::log(mtlba_pmax(R::pnorm(v1 / sv1, 0.0, 1.0, 1, 0) +
                             mtlba_bvu(v1 / sv1, -v2 / sv2, -rho), 0.0));
}

// log normaliser of the drift convention
inline double mtlba_lognorm(double v1, double sv1, double v2, double sv2, double rho, bool pos)
{
  return pos ? mtlba_logppos(v1, sv1, v2, sv2, rho) : mtlba_lognotneg(v1, sv1, v2, sv2, rho);
}

// Correlated twin of mtlba_logstop for one parameter row and a finite bound
// y > 0, as the R twin mtlba_stop_rho1: one integral over accumulator 1's
// drift x > 0, composite Gauss-Legendre between the kinks of accumulator 1's
// crossing counts, all terms positive, over P(both > 0). T1, T2: thresholds
// 1 .. K of the two accumulators. x, w: work space. Under "any" (!pos) x runs
// over the whole line: for x > 0 the other drift is unrestricted, for x < 0
// it must be positive (a break at x = 0), over P(not both < 0).
inline void mtlba_stop_rho1(double y, double A1, double v1, double s1, const double* T1,
                            double A2, double v2, double s2, const double* T2,
                            double rho, int K, bool TC, double& lP, double& lS,
                            std::vector<double>& x, std::vector<double>& w, bool pos = true)
{
  lP = lS = R_NaN;
  bool fin = std::isfinite(y) && std::isfinite(A1) && std::isfinite(v1) && std::isfinite(s1) &&
             std::isfinite(A2) && std::isfinite(v2) && std::isfinite(s2) && std::isfinite(rho);
  for (int n = 0; fin && n < K; ++n) fin = std::isfinite(T1[n]) && std::isfinite(T2[n]);
  if (!fin || !(std::fabs(rho) < 1) || !(s1 > 0) || !(s2 > 0)) return;
  // window: v2 > 0 drags x towards m when v2's mean is negative
  const double m  = v1 + rho * s1 / s2 * ((v2 > 0 ? v2 : 0.0) - v2);
  const double hi = v1 > m ? v1 : m, lw = v1 < m ? v1 : m;
  const double d  = (pos && -hi > 0) ? -hi : 0.0;
  const double r  = std::sqrt(d * d + (MTLBA_STOP_HW * s1) * (MTLBA_STOP_HW * s1));
  const double a  = pos ? ((lw - r > 0) ? lw - r : 0.0) : lw - r;
  const double b  = hi + r;
  double brk[2 * MTLBA_MAX_K + 3];
  int nb = 0;
  brk[nb++] = a; brk[nb++] = b;
  for (int n = 0; n < K; ++n)
    for (double kink : {(T1[n] - A1) / y, T1[n] / y})
      if (kink > a && kink < b) brk[nb++] = kink;
  if (!pos && 0.0 > a && 0.0 < b) brk[nb++] = 0.0;
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
    // "positive": the other drift must be positive; "any": only when x < 0
    const bool pos2 = pos || xi < 0;
    double Sint = 0.0, Pint = (1.0 - Fle[K]) * (pos ? R::pnorm(mc / sc, 0.0, 1.0, 1, 0) : 1.0);
    for (int n1 = TC ? 0 : K - 1; n1 < K; ++n1) {
      // accumulator 1 at n1 crossings (BE: below b), accumulator 2 against T_{2,mm}
      const double wn = TC ? Fle[n1 + 1] - Fle[n1] : Fle[K];
      if (wn == 0) continue;
      const int mm = TC ? K - n1 : K;
      Sint += wn * std::exp(pos2 ? mtlba_logband_un(y, A2, mc, sc, 0.0, T2[mm - 1], true)
                                 : mtlba_logH(T2[mm - 1], y, A2, mc, sc));
      Pint += wn * std::exp(pos ? mtlba_logSp(T2[mm - 1], y, A2, mc, sc)
                                : mtlba_logS(T2[mm - 1], y, A2, mc, sc));
    }
    sumP += wi * Pint;
    sumS += wi * Sint;
  }
  const double lpp = mtlba_lognorm(v1, s1, v2, s2, rho, pos);
  lP = std::log(sumP) - lpp;
  lS = std::log(sumS) - lpp;
}

// Per-trial column pointers
struct MTLBACols {
  const double *v, *sv, *B, *A, *t0, *rho;
  const double *gp = nullptr, *gB = nullptr;   // confidence guessing (null: none)
  int gbc = 0;                                 // gb_constraint: 0 none, 1 own gp, 2 max gp
  bool pos = true;                             // drifts "positive" (else "any")
  bool free_lowest = false;                    // T_1 = T_2 exp(-c_{K-1})
  std::vector<const double*> c, pg;
};

// Threshold j (1 .. K; K = b) of the accumulator in row r. The R twin's
// arithmetic: T_j = A + (b - A) d_j, a free lowest threshold through its
// proportion d_1 = (T_2 exp(-c_{K-1}) - A) / (b - A).
inline double mtlba_threshold(const MTLBACols& p, int K, int r, int j)
{
  const double A = p.A[r];
  const double b = p.B[r] + A;
  if (j >= K) return b;
  if (j == 1 && p.free_lowest) {
    const double T2 = mtlba_threshold(p, K, r, 2);
    const double d1 = (T2 * std::exp(-p.c[K - 2][r]) - A) / (b - A);
    return A + (b - A) * d1;
  }
  double C = 0.0;
  for (int k = 0; k < K - j; ++k) C += p.c[k][r];
  return A + (b - A) * std::exp(-C);
}

// The log normaliser is the same for all trials of a design cell: a small
// cache, local to one likelihood call, keyed on the exact (v_s, sv_s, v_o,
// sv_o, rho); floored at log(1e-10).
struct MTLBAPposCache {
  static constexpr int MAX_ENTRIES = 16;
  int n = 0, next = 0;
  bool pos = true;
  double keys[MAX_ENTRIES][5], vals[MAX_ENTRIES];
  double get(double vs, double svs, double vo, double svo, double rho) {
    for (int e = 0; e < n; ++e) {
      const double* k = keys[e];
      if (k[0] == vs && k[1] == svs && k[2] == vo && k[3] == svo && k[4] == rho) return vals[e];
    }
    const double val = mtlba_pmax(mtlba_lognorm(vs, svs, vo, svo, rho, pos), std::log(1e-10));
    const int e = (n < MAX_ENTRIES) ? n++ : next++ % MAX_ENTRIES;
    keys[e][0] = vs; keys[e][1] = svs; keys[e][2] = vo; keys[e][3] = svo; keys[e][4] = rho;
    vals[e] = val;
    return val;
  }
};

// log kernel for explicit thresholds: the accumulator in row rs reaches Ts at
// rt - t0_s while the evidence of the one in row ro lies in (ly, hi) (lowest: no
// threshold below); rho the trial's drift-rate correlation (0: closed form).
// -Inf for undefined values (non-finite parameters, dt <= 0, an empty band, rho
// not inside (-1, 1)).
inline double mtlba_logkernel_y(const MTLBACols& p, double rt, int rs, int ro, double Ts,
                                double ly, double hi, bool lowest, double rho,
                                MTLBAPposCache& ppos)
{
  const double dt = rt - p.t0[rs];
  if (!(std::isfinite(dt) && dt > 0)) return R_NegInf;
  for (int r : {rs, ro})
    if (!std::isfinite(p.v[r]) || !std::isfinite(p.sv[r]) || !std::isfinite(p.A[r]) ||
        !std::isfinite(p.B[r] + p.A[r])) return R_NegInf;
  if (!std::isfinite(Ts) || !std::isfinite(hi)) return R_NegInf;
  if (!lowest && !(std::isfinite(ly) && hi > ly)) return R_NegInf;
  if (!(std::fabs(rho) < 1)) return R_NegInf;   // NaN too
  if (rho == 0) {
    const double l = mtlba_logf(dt, p.A[rs], Ts, p.v[rs], p.sv[rs], p.pos) +
                     mtlba_logband(dt, p.A[ro], p.v[ro], p.sv[ro], ly, hi, lowest, p.pos);
    // "any": the closed form is a product of untruncated terms over P(not both < 0)
    return p.pos ? l : l - ppos.get(p.v[rs], p.sv[rs], p.v[ro], p.sv[ro], 0.0);
  }
  return mtlba_logkernel_rho(dt, p.A[rs], Ts, p.v[rs], p.sv[rs], p.A[ro], p.v[ro], p.sv[ro],
                             ly, hi, lowest, rho, p.pos) -
         ppos.get(p.v[rs], p.sv[rs], p.v[ro], p.sv[ro], rho);
}

// log kernel: accumulator in row rs stops at its threshold j while the one in
// row ro lies in band lo.
inline double mtlba_logkernel(const MTLBACols& p, int K, double rt, int rs, int ro, int j, int lo,
                              double rho, MTLBAPposCache& ppos)
{
  return mtlba_logkernel_y(p, rt, rs, ro, mtlba_threshold(p, K, rs, j),
                           (lo == 0) ? 0.0 : mtlba_threshold(p, K, ro, lo),
                           mtlba_threshold(p, K, ro, lo + 1), lo == 0, rho, ppos);
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

// ---- confidence guessing ----

// Guess threshold of the accumulator in row r: gb = A + (b - A) gB, or with
// gb_constraint the proportion gB + g^2 (1 - gB), g the row's own gp or the
// larger of its trial's two (the other accumulator is the adjacent row; with a
// trial-level trigger the two rows hold the same gp).
inline double mtlba_gb(const MTLBACols& p, int r)
{
  const double A = p.A[r];
  double gB = p.gB[r];
  if (p.gbc) {
    const double g = (p.gbc == 2) ? mtlba_pmax(p.gp[r], p.gp[r ^ 1]) : p.gp[r];
    gB = gB + g * g * (1.0 - gB);
  }
  return A + ((p.B[r] + A) - A) * gB;
}

// log probability that a guess gives rating c: pg_c prod_{k<c} (1 - pg_k), the
// highest rating taking the rest.
inline double mtlba_lpg(const MTLBACols& p, int K, int r, int c)
{
  double l = 0.0;
  for (int k = 1; k < c; ++k) l += std::log1p(-p.pg[k - 1][r]);
  if (c < K) l += std::log(p.pg[c - 1][r]);
  return l;
}

// log(w1 exp(l1) + w2 exp(l2)) for weights >= 0 (a zero or NaN weight drops its term)
inline double mtlba_lmix(double w1, double l1, double w2, double l2)
{
  double a = (w1 > 0) ? std::log(w1) + l1 : R_NegInf;
  double b = (w2 > 0) ? std::log(w2) + l2 : R_NegInf;
  if (std::isnan(a)) a = R_NegInf;
  if (std::isnan(b)) b = R_NegInf;
  return mtlba_lse(a, b);
}

// log density of the response (R, RR) at rt, no don't-know merging; rows r1 / r2
// = accumulators 1 / 2. As the R twin mtlba_logresp.
inline double mtlba_logresp(const MTLBACols& p, int K, bool TC, bool gtrial, double rt,
                            int R, int RR, int r1, int r2, double rho, MTLBAPposCache& ppos)
{
  int a, j, lo;
  mtlba_lookup(R, RR, K, TC, a, j, lo);
  const int rs = (a == 1) ? r1 : r2, ro = (a == 1) ? r2 : r1;
  const double full = mtlba_logkernel(p, K, rt, rs, ro, j, lo, rho, ppos);
  if (!p.gp) return full;
  const double g1 = p.gp[r1], g2 = gtrial ? g1 : p.gp[r2];
  if (g1 == 0 && g2 == 0) return full;
  const int rw = (R == 1) ? r1 : r2, rv = (R == 1) ? r2 : r1;   // the choice guesses
  const double gb_w = mtlba_gb(p, rw), gb_v = mtlba_gb(p, rv);
  const double lpgc = mtlba_lpg(p, K, rw, RR);
  double ev, gu;
  if (gtrial) {
    ev = mtlba_lmix(1.0 - g1, full, 0.0, 0.0);
    gu = mtlba_lmix(g1, lpgc + mtlba_logkernel_y(p, rt, rw, rv, gb_w, 0.0, gb_v, true, rho, ppos),
                    0.0, 0.0);
  } else {
    const double gs = (a == 1) ? g1 : g2, go = (a == 1) ? g2 : g1;
    const double gw = (R == 1) ? g1 : g2, gv = (R == 1) ? g2 : g1;
    const double Ts = mtlba_threshold(p, K, rs, j);
    const double hi = mtlba_threshold(p, K, ro, lo + 1);
    const double ly = (lo == 0) ? 0.0 : mtlba_threshold(p, K, ro, lo);
    double cs = 1.0 - gs;
    if (TC && Ts < mtlba_gb(p, rs)) cs = 1.0;
    const double cap = mtlba_pmin(hi, mtlba_gb(p, ro));
    double tr = full;
    if (cap < hi)
      tr = (go > 0) ? mtlba_logkernel_y(p, rt, rs, ro, Ts, ly, cap, lo == 0, rho, ppos) : R_NegInf;
    ev = mtlba_lmix(cs * (1.0 - go), full, cs * go, tr);
    double U0;
    if (!TC) U0 = p.B[rv] + p.A[rv];
    else {
      int nw = 0;
      for (int jj = 1; jj <= K; ++jj) if (mtlba_threshold(p, K, rw, jj) < gb_w) ++nw;
      U0 = mtlba_threshold(p, K, rv, K - nw);
    }
    const double U1 = mtlba_pmin(U0, gb_v);
    const bool capped = U1 < U0;
    const double gk0 = (gw > 0 && (gv < 1 || !capped)) ?
      mtlba_logkernel_y(p, rt, rw, rv, gb_w, 0.0, U0, true, rho, ppos) : R_NegInf;
    double gk1 = gk0;
    if (capped)
      gk1 = (gw > 0 && gv > 0) ? mtlba_logkernel_y(p, rt, rw, rv, gb_w, 0.0, U1, true, rho, ppos)
                               : R_NegInf;
    gu = lpgc + mtlba_lmix(gw * (1.0 - gv), gk0, gw * gv, gk1);
  }
  const double out = mtlba_lse(ev, gu);
  return std::isnan(out) ? R_NegInf : out;
}

// log P(dt <= y) (lP) and log P(dt > y) (lS) for the race's stopping time at a
// finite decision-time bound y > 0, rho = 0; rows r1 / r2 = accumulators 1 / 2
// with thresholds T1, T2 (1 .. K). As the R twin mtlba_logstop: with n_a the
// thresholds crossed by accumulator a at y, F_a(n) = P(n_a <= n), Fc_a(n) =
// 1 - F_a(n), p_a(n) = P(n_a = n) (v_a > 0),
//   BE: S = F_1(K-1) F_2(K-1),   P = Fc_1(K-1) + F_1(K-1) Fc_2(K-1)
//   TC: S = sum_{n1<K} p_1(n1) F_2(K-1-n1),
//       P = Fc_1(K-1) + sum_{n1<K} p_1(n1) Fc_2(K-1-n1),
// all terms positive (log-sum-exp). NaN for undefined parameters.
// Under "any" drifts (!p.pos) the same positive-drift terms, unnormalised,
// are completed by the trials with one negative drift (N_a = P(v_a < 0)):
//   stopped  = [stopped_pos + N_1 S_2(T_{2,K}) + N_2 S_1(T_{1,K})] / P(not both < 0)
//   survivor = [survivor_pos + N_1 Hp_2(T_{2,K}) + N_2 Hp_1(T_{1,K})] / P(not both < 0).
inline void mtlba_logstop(const MTLBACols& p, int K, bool TC, double y, int r1, int r2,
                          const double* T1, const double* T2, double& lP, double& lS)
{
  double F[2][MTLBA_MAX_K], Fc[2][MTLBA_MAX_K], band[2][MTLBA_MAX_K], lPhi[2], lN[2];
  const int rows[2] = {r1, r2};
  const double* Ts[2] = {T1, T2};
  const bool pos = p.pos;
  for (int a = 0; a < 2; ++a) {
    const int r = rows[a];
    const double A = p.A[r], v = p.v[r], sv = p.sv[r];
    if (!std::isfinite(v) || !std::isfinite(sv) || !std::isfinite(A) ||
        !std::isfinite(p.B[r] + A)) { lP = lS = R_NaN; return; }
    lPhi[a] = R::pnorm(v / sv, 0.0, 1.0, 1, 1);
    lN[a]   = R::pnorm(-v / sv, 0.0, 1.0, 1, 1);
    const double lden = pos ? mtlba_pmax(lPhi[a], std::log(1e-10)) : 0.0;
    double lo_y = 0.0;
    for (int nn = 0; nn < K; ++nn) {
      const double hi = Ts[a][nn];
      if (!std::isfinite(hi)) { lP = lS = R_NaN; return; }
      F[a][nn]  = pos ? mtlba_logband(y, A, v, sv, 0.0, hi, true)
                      : mtlba_logband_un(y, A, v, sv, 0.0, hi, true);
      Fc[a][nn] = mtlba_logSp(hi, y, A, v, sv) - lden;
      band[a][nn] = (nn == 0) ? F[a][0] :
        (pos ? mtlba_logband(y, A, v, sv, lo_y, hi, false)
             : mtlba_logband_un(y, A, v, sv, lo_y, hi, false));
      lo_y = hi;
    }
  }
  const double l2 = pos ? 0.0 : lPhi[1];   // "any": accumulator 2's positivity explicit
  if (!TC) {
    lS = F[0][K - 1] + F[1][K - 1];
    lP = mtlba_lse(Fc[0][K - 1] + l2, F[0][K - 1] + Fc[1][K - 1]);
  } else {
    lS = band[0][0] + F[1][K - 1];
    lP = mtlba_lse(Fc[0][K - 1] + l2, band[0][0] + Fc[1][K - 1]);
    for (int n1 = 1; n1 < K; ++n1) {
      lS = mtlba_lse(lS, band[0][n1] + F[1][K - 1 - n1]);
      lP = mtlba_lse(lP, band[0][n1] + Fc[1][K - 1 - n1]);
    }
  }
  if (!pos) {
    const double lZ = mtlba_lse(lPhi[0], lN[0] + lPhi[1]);   // log P(not both < 0)
    lS = mtlba_lse(mtlba_lse(lS, lN[0] + F[1][K - 1]), lN[1] + F[0][K - 1]) - lZ;
    lP = mtlba_lse(mtlba_lse(lP, lN[0] + Fc[1][K - 1]), lN[1] + Fc[0][K - 1]) - lZ;
  }
}

// Stopping-time terms at decision-time bound y for any rho and with confidence
// guessing, as the R twins mtlba_logstop_any / mtlba_logstop_guess: rho = 0 the
// closed form, otherwise mtlba_stop_rho1; with guessing the mixture over the
// trigger states, a triggering accumulator's thresholds capped at its gb (the
// trial-level guess state is the BE race to the two guess thresholds). x, w:
// work space of the correlated integral.
inline void mtlba_stop_mix(const MTLBACols& p, int K, bool TC, bool gtrial, double y,
                           int r1, int r2, double rho, double& lP, double& lS,
                           std::vector<double>& x, std::vector<double>& w)
{
  if (std::isnan(y)) { lP = lS = R_NaN; return; }
  if (y == R_PosInf) { lP = 0.0; lS = R_NegInf; return; }
  if (!(y > 0)) { lP = R_NegInf; lS = 0.0; return; }
  double T[2][MTLBA_MAX_K], C[2][MTLBA_MAX_K];
  const int rows[2] = {r1, r2};
  for (int a = 0; a < 2; ++a)
    for (int j = 1; j <= K; ++j) T[a][j - 1] = mtlba_threshold(p, K, rows[a], j);
  auto state = [&](const double* T1, const double* T2, bool tc, double& P, double& S) {
    if (rho == 0) mtlba_logstop(p, K, tc, y, r1, r2, T1, T2, P, S);
    else mtlba_stop_rho1(y, p.A[r1], p.v[r1], p.sv[r1], T1, p.A[r2], p.v[r2], p.sv[r2], T2,
                         rho, K, tc, P, S, x, w, p.pos);   // NaN rho: NaN
  };
  const double g1 = p.gp ? p.gp[r1] : 0.0;
  const double g2 = p.gp ? (gtrial ? g1 : p.gp[r2]) : 0.0;
  if (g1 == 0 && g2 == 0) { state(T[0], T[1], TC, lP, lS); return; }
  for (int a = 0; a < 2; ++a) {
    const double gb = mtlba_gb(p, rows[a]);
    for (int j = 0; j < K; ++j) C[a][j] = mtlba_pmin(T[a][j], gb);
  }
  lP = lS = R_NegInf;
  auto add = [&](double wt, bool c1, bool c2, bool tc) {
    if (!(std::isnan(wt) || wt > 0)) return;
    double P, S;
    state(c1 ? C[0] : T[0], c2 ? C[1] : T[1], tc, P, S);
    lP = mtlba_lse(lP, std::log(wt) + P);
    lS = mtlba_lse(lS, std::log(wt) + S);
  };
  if (gtrial) {
    add(1.0 - g1, false, false, TC);
    add(g1, true, true, false);
  } else {
    add((1.0 - g1) * (1.0 - g2), false, false, TC);
    add(g1 * (1.0 - g2), true, false, TC);
    add((1.0 - g1) * g2, false, true, TC);
    add(g1 * g2, true, true, TC);
  }
}

// log P(a < dt <= b) from the (lP, lS) pairs at a and b: the route with the
// smaller leading term.
inline double mtlba_logwindow(double Pa, double Sa, double Pb, double Sb)
{
  if (Sa < Pb) return mtlba_ldiff(Sa, Sb);   // false for NaN
  return mtlba_ldiff(Pb, Pa);
}

// The stopping-time terms depend only on the trial's parameter rows and the
// bound, not on the rt, so trials of the same design cell share them. A small
// cache, local to one likelihood call, keyed on the exact values (the bound y,
// rho, then v, sv, A, B, gp, gB and the criteria of both rows); NaN never
// matches. The correlated terms are one integral per trigger state, so the
// cache is what makes them affordable.
struct MTLBAStopCache {
  static constexpr int MAX_ENTRIES = 128;
  int n_key = 0, n = 0, next = 0;
  std::vector<double> keys, vals;   // MAX_ENTRIES x n_key, MAX_ENTRIES x 2
  std::vector<double> key, x, w;
  void init(int K) {
    n_key = 2 + 2 * (6 + K - 1);
    keys.resize((size_t)MAX_ENTRIES * n_key);
    vals.resize((size_t)MAX_ENTRIES * 2);
    key.resize(n_key);
  }
  void get(const MTLBACols& p, int K, bool TC, bool gtrial, double y, int r1, int r2, double rho,
           double& lP, double& lS) {
    if (std::isnan(y)) { lP = lS = R_NaN; return; }
    if (y == R_PosInf) { lP = 0.0; lS = R_NegInf; return; }
    if (!(y > 0)) { lP = R_NegInf; lS = 0.0; return; }
    if (n_key == 0) init(K);
    int k = 0;
    key[k++] = y; key[k++] = rho;
    for (int r : {r1, r2}) {
      key[k++] = p.v[r]; key[k++] = p.sv[r]; key[k++] = p.A[r]; key[k++] = p.B[r];
      key[k++] = p.gp ? p.gp[r] : 0.0; key[k++] = p.gp ? p.gB[r] : 0.0;
      for (int c = 0; c < K - 1; ++c) key[k++] = p.c[c][r];
    }
    for (int e = 0; e < n; ++e) {
      const double* ke = keys.data() + (size_t)e * n_key;
      int i = 0;
      while (i < n_key && ke[i] == key[i]) ++i;
      if (i == n_key) { lP = vals[2 * e]; lS = vals[2 * e + 1]; return; }
    }
    mtlba_stop_mix(p, K, TC, gtrial, y, r1, r2, rho, lP, lS, x, w);
    const int e = (n < MAX_ENTRIES) ? n++ : next++ % MAX_ENTRIES;   // round-robin when full
    std::copy(key.begin(), key.end(), keys.begin() + (size_t)e * n_key);
    vals[2 * e] = lP; vals[2 * e + 1] = lS;
  }
};

// Proportional rows -> the additive model's columns, exactly (no root
// finding): with b = B, top threshold T = b c1 and A = T A (K = 1: A = b A),
// natural A, B = b - A, c1 = -log((T - A) / (b - A)) (K = 2 with a free lowest
// threshold: T_1 = b exp(-c1), so c1 = -log(c1_prop)) and c_k = -log(c_k_prop)
// for k >= 2 (T_j - A = (T_{j+1} - A) c_k above A; the free form
// T_1 = T_2 exp(-c_{K-1}) below it: the same identity). Proportions in (0, 1)
// give A < T < b; values out of bounds are voided by the caller.
inline void mtlba_proportional(MTLBACols& p, int K, int n, std::vector<double>& A,
                               std::vector<double>& B, std::vector<std::vector<double>>& c)
{
  A.resize(n); B.resize(n);
  c.assign(K > 1 ? K - 1 : 0, std::vector<double>(n));
  for (int r = 0; r < n; ++r) {
    const double b = p.B[r];
    if (K == 1) {
      A[r] = b * p.A[r];
    } else {
      const double T = b * p.c[0][r];
      A[r] = T * p.A[r];
      c[0][r] = (K == 2 && p.free_lowest) ? -std::log(p.c[0][r]) : -std::log((T - A[r]) / (b - A[r]));
      for (int k = 1; k < K - 1; ++k) c[k][r] = -std::log(p.c[k][r]);
    }
    B[r] = b - A[r];
  }
  p.A = A.data();
  p.B = B.data();
  for (int k = 0; k < K - 1; ++k) p.c[k] = c[k].data();
}

// Binomial guess ratings: q per row -> the sequential pg_k = P(k) / P(>= k),
// k = 1 .. K - 1, P(k) = choose(K - 1, k - 1) q^(k-1) (1 - q)^(K-k), on the log
// scale with the tails P(>= k) as log-sums of their terms (the R twin's
// arithmetic, mtlba_binomial_pg); pg_k = 1 where P(>= k) = 0.
inline void mtlba_binomial_pg(MTLBACols& p, int K, int n, const double* q,
                              std::vector<std::vector<double>>& pg)
{
  pg.assign(K - 1, std::vector<double>(n));
  std::vector<double> lch(K), lP(K);
  for (int k = 1; k <= K; ++k) {   // log choose(K - 1, k - 1), exact integer first
    double c = 1.0;
    for (int j = 1; j <= k - 1; ++j) c = c * (K - j) / j;
    lch[k - 1] = std::log(c);
  }
  for (int r = 0; r < n; ++r) {
    const double lq = std::log(q[r]), l1q = std::log1p(-q[r]);
    for (int k = 1; k <= K; ++k)
      lP[k - 1] = lch[k - 1] + (k > 1 ? (k - 1) * lq : 0.0) + (k < K ? (K - k) * l1q : 0.0);
    double tail = lP[K - 1];
    for (int k = K - 1; k >= 1; --k) {
      tail = mtlba_lse(lP[k - 1], tail);
      pg[k - 1][r] = (tail == R_NegInf) ? 1.0 : std::exp(lP[k - 1] - tail);
    }
  }
  p.pg.clear();
  for (int k = 0; k < K - 1; ++k) p.pg.push_back(pg[k].data());
}

// Return codes of mtlba_trial_ll (the caller stops, outside any parallel region).
enum { MTLBA_OK = 0, MTLBA_T0 = 2 };
inline const char* mtlba_error_message(int code)
{
  (void)code;
  return "MTLBA: censoring and truncation need t0 equal for both accumulators (t0 must not depend on lR or lM)";
}

// Trial log-likelihoods into ll_buf (indexed by trial), floored at min_ll.
// Observed trials: the response density (kernel, with guessing a mixture of
// kernels); censored trials: the probability of stopping in the censored
// window(s); truncated trials: minus log P(LT < rt < UT). rho is read from
// accumulator 1's row (a trial-level parameter). Proportional rows are first
// translated into the additive columns (mtlba_proportional). Returns MTLBA_T0
// if a censored or truncated trial has unequal t0 (the caller stops, outside
// any parallel region).
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
  if (s.col_gp >= 0) {
    p.gp = pt.base.colptr(s.col_gp);
    p.gB = pt.base.colptr(s.col_gB);
    if (!s.binomial_guess) for (int col : s.col_pg) p.pg.push_back(pt.base.colptr(col));
    p.gbc = s.gb_constraint;
  }
  p.pos = !s.any;
  p.free_lowest = s.free_lowest;
  const int K = s.K;
  // per-call scratch (thread-local: one call per particle) for the translated columns
  std::vector<double> propA, propB;
  std::vector<std::vector<double>> propc;
  if (s.proportional) mtlba_proportional(p, K, pt.base.nrow, propA, propB, propc);
  std::vector<std::vector<double>> binpg;
  if (s.binomial_guess && s.col_gp >= 0 && K > 1)
    mtlba_binomial_pg(p, K, pt.base.nrow, pt.base.colptr(s.col_pq), binpg);
  const bool gtrial = s.guess_trial;
  int code = MTLBA_OK;
  MTLBAStopCache cache;
  MTLBAPposCache ppos;
  ppos.pos = p.pos;

  const int nt = (int)s.rt.size();
  for (int t = 0; t < nt; ++t) {
    const int r1 = s.row1[t];
    const int r2 = (r1 % 2 == 0) ? r1 + 1 : r1 - 1;
    const double rho = p.rho[r1];
    const int miss = s.miss[t];
    double ll = R_NegInf;
    if (miss == 0) {
      const int RR = s.RR[t];
      if (s.DK && RR == 1)
        ll = mtlba_lse(mtlba_logresp(p, K, s.TC, gtrial, s.rt[t], 1, 1, r1, r2, rho, ppos),
                       mtlba_logresp(p, K, s.TC, gtrial, s.rt[t], 2, 1, r1, r2, rho, ppos));
      else
        ll = mtlba_logresp(p, K, s.TC, gtrial, s.rt[t], s.R[t], RR, r1, r2, rho, ppos);
    }
    if (miss != 0 || s.trunc[t]) {
      const double t0 = p.t0[r1];
      if (!std::isnan(t0) && !std::isnan(p.t0[r2]) && t0 != p.t0[r2]) {
        if (code == MTLBA_OK) code = MTLBA_T0;
        ll_buf[t] = min_ll;
        continue;
      }
      auto stop_at = [&](double y, double& P, double& S) {
        cache.get(p, K, s.TC, gtrial, y, r1, r2, rho, P, S);
      };
      double PLT, SLT, PUT, SUT;
      stop_at(s.LT[t] - t0, PLT, SLT);
      stop_at(s.UT[t] - t0, PUT, SUT);
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
