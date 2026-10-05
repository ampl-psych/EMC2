#ifndef MODEL_MTLNR_H
#define MODEL_MTLNR_H

// Correlated multiple-threshold log-normal race (MTLNR; Reynolds, Kvam, Osth &
// Heathcote, 2020). C++ twin of R/model_MTLNR.R; see there for the model.
//
// Joint per-trial kernel: each trial has two contiguous accumulator rows. With
// winner w, loser l, dt = rt - t0_w and rating RR = c (1 = lowest .. K),
//   ll = log dlnorm(dt; m_w, s_w) + log P(z_hi < Z < z_lo),
//   mc = m_l + rho * (s_l / s_w) * (log dt - m_w),  sc = s_l * sqrt((1-rho)(1+rho)),
//   z_hi = (log dt + C_{c-1} - mc) / sc,  z_lo = (log dt + C_c - mc) / sc,
// where C_k = c1 + ... + ck are the loser's natural-scale criteria (C_0 = 0,
// C_K = Inf), i.e. -log of the thresholds d_{K-k}. rho comes from the trial's
// first row, t0 from the winner's, the criteria from the loser's. Ttransform is
// not applied in C++, so the thresholds are built here from c1 .. c{K-1}.
//
// R::pnorm (Rmath pnorm_both) is used for the log tails: it is pure arithmetic
// with no R API call for a standard normal, so it is safe in the OpenMP loop
// (the pnorm parameter transform already calls it there).

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include "ParamTable.h"

struct MTLNRSpec {
  int K       = 1;
  int col_m   = -1;
  int col_s   = -1;
  int col_t0  = -1;
  int col_rho = -1;
  std::vector<int>    col_c;   // c1 .. c{K-1}
  std::vector<double> rt;      // per row
  std::vector<int>    RR;      // per row, 1 .. K
};

// Data and column lookups, once per call. Refuses what the likelihood does not
// (yet) define: other than two accumulators, missing/censored responses and
// truncation.
inline MTLNRSpec make_mtlnr_spec(const Rcpp::DataFrame& data, const ParamTable& pt, int n_lR)
{
  if (n_lR != 2) Rcpp::stop("MTLNR: needs exactly two accumulators per trial (found %d)", n_lR);
  MTLNRSpec s;
  s.col_m   = pt.base_index_for("m");
  s.col_s   = pt.base_index_for("s");
  s.col_t0  = pt.base_index_for("t0");
  s.col_rho = pt.base_index_for("rho");
  for (int k = 1; ; ++k) {
    auto it = pt.name_to_base_idx.find("c" + std::to_string(k));
    if (it == pt.name_to_base_idx.end()) break;
    s.col_c.push_back(it->second);
  }
  s.K = (int)s.col_c.size() + 1;

  const int n = data.nrow();
  if (n % 2 != 0) Rcpp::stop("MTLNR: odd number of data rows");
  Rcpp::CharacterVector nm = data.names();
  auto has = [&](const char* x) { for (int i = 0; i < nm.size(); ++i) if (nm[i] == x) return true; return false; };

  if (has("missingness")) {
    Rcpp::IntegerVector ms = data["missingness"];
    for (int i = 0; i < n; ++i)
      if (!Rcpp::IntegerVector::is_na(ms[i]))
        Rcpp::stop("MTLNR: missing or censored responses are not yet supported");
  }
  if (has("LT")) {
    Rcpp::NumericVector lt = Rcpp::as<Rcpp::NumericVector>(data["LT"]);
    for (int i = 0; i < n; ++i)
      if (lt[i] > 0 && std::isfinite(lt[i])) Rcpp::stop("MTLNR: truncation is not yet supported");
  }
  if (has("UT")) {
    Rcpp::NumericVector ut = Rcpp::as<Rcpp::NumericVector>(data["UT"]);
    for (int i = 0; i < n; ++i)
      if (std::isfinite(ut[i])) Rcpp::stop("MTLNR: truncation is not yet supported");
  }

  Rcpp::NumericVector rt = Rcpp::as<Rcpp::NumericVector>(data["rt"]);
  s.rt.assign(rt.begin(), rt.end());
  s.RR.assign(n, 1);
  if (has("RR")) {
    Rcpp::NumericVector rr = Rcpp::as<Rcpp::NumericVector>(data["RR"]);
    for (int i = 0; i < n; ++i) {
      const double v = rr[i];
      if (!(v >= 1 && v <= s.K && v == std::floor(v)))
        Rcpp::stop("MTLNR: RR must be an integer in 1 .. %d (row %d)", s.K, i + 1);
      s.RR[i] = (int)v;
    }
  } else if (s.K > 1) {
    Rcpp::stop("MTLNR: the data have no rating column RR");
  }
  return s;
}

// log(exp(la) - exp(lb)), la >= lb.
inline double mtlnr_log_diff_exp(double la, double lb)
{
  if (la == lb) return R_NegInf;
  const double x = lb - la;
  return la + (x > -M_LN2 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x)));
}

// log P(a < Z < b), a <= b, accurate in both tails and across 0. Across 0,
// Phi(b) - Phi(a) = (erf(b / sqrt 2) + erf(-a / sqrt 2)) / 2 (the R twin
// writes erf(|z| / sqrt 2) as pchisq(z^2, 1)).
inline double mtlnr_log_pnorm_interval(double a, double b)
{
  if (a >= 0)
    return mtlnr_log_diff_exp(R::pnorm(a, 0.0, 1.0, 0, 1), R::pnorm(b, 0.0, 1.0, 0, 1));
  if (b <= 0)
    return mtlnr_log_diff_exp(R::pnorm(b, 0.0, 1.0, 1, 1), R::pnorm(a, 0.0, 1.0, 1, 1));
  return std::log(0.5 * (std::erf(b * M_SQRT1_2) + std::erf(-a * M_SQRT1_2)));
}

// Trial log-likelihoods into ll_buf (indexed by trial) for the winner rows.
// Impossible or undefined trials get min_ll.
inline void mtlnr_trial_ll(const MTLNRSpec& s, const ParamTable& pt,
                           const std::vector<int>& idx_win, double min_ll,
                           double* __restrict__ ll_buf)
{
  const double* m   = pt.base.colptr(s.col_m);
  const double* sd  = pt.base.colptr(s.col_s);
  const double* t0  = pt.base.colptr(s.col_t0);
  const double* rho = pt.base.colptr(s.col_rho);
  const int Km1 = s.K - 1;
  std::vector<const double*> cc(Km1);
  for (int k = 0; k < Km1; ++k) cc[k] = pt.base.colptr(s.col_c[k]);
  const double LN_SQRT_2PI = 0.918938533204672741780329736406;

  const int n_win = (int)idx_win.size();
  for (int t = 0; t < n_win; ++t) {
    const int w    = idx_win[t];
    const int base = (w / 2) * 2;
    const int l    = (w == base) ? base + 1 : base;
    double& out = ll_buf[w / 2];

    const double dt = s.rt[w] - t0[w];
    const double r  = rho[base];
    const double sc = sd[l] * std::sqrt((1.0 - r) * (1.0 + r));
    if (!(dt > 0) || !std::isfinite(dt) || !(sc > 0)) { out = min_ll; continue; }

    const double ldt = std::log(dt);
    const double zw  = (ldt - m[w]) / sd[w];
    const double lpdf = -(LN_SQRT_2PI + 0.5 * zw * zw + std::log(dt * sd[w]));
    const double mc  = m[l] + r * (sd[l] / sd[w]) * (ldt - m[w]);

    // C_{c-1} and C_c from the loser's criteria
    const int c = s.RR[w];
    double C_lo = 0.0;
    for (int k = 0; k < c - 1; ++k) C_lo += cc[k][l];
    const double z_hi = (ldt + C_lo - mc) / sc;
    const double z_lo = (c == s.K) ? R_PosInf : (ldt + C_lo + cc[c - 1][l] - mc) / sc;

    const double ll = lpdf + mtlnr_log_pnorm_interval(z_hi, z_lo);
    out = (ll > min_ll) ? ll : min_ll;   // NaN -> min_ll
  }
}

#endif // MODEL_MTLNR_H
