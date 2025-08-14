#ifndef race_integrate_h
#define race_integrate_h

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>

using namespace Rcpp;
using namespace Numer;

class race_f : public Func {
private:
  NumericVector (*lpdf)(NumericVector, NumericMatrix, LogicalVector, double);
  NumericVector (*lccdf)(NumericVector, NumericMatrix, LogicalVector, double);
  const NumericMatrix pars;
  const LogicalVector winner;
  const double min_ll;
  const bool log_d;

public:
  race_f(
    NumericVector (*lpdf_)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericVector (*lccdf_)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericMatrix pars_,
    LogicalVector winner_,
    double min_ll_,
    bool log_d_ = false
  ) :
  lpdf(lpdf_),
  lccdf(lccdf_),
  pars(pars_),
  winner(winner_),
  min_ll(min_ll_),
  log_d(log_d_) {}

  double operator()(const double& x) const {
    int n_acc = pars.nrow();
    int n_winner = sum(winner);
    int n_loser = sum(!winner);
    // broadcast input response time x to repped vector, one per accumulator
    NumericVector t(n_acc, x);
    // initialise log likelihood
    double log_out = 0.;
    if (n_winner == 1) {
      // compute log density of winner accumulator finishing at time t
      NumericVector log_d = lpdf(t, pars, winner, min_ll);
      log_out += Rcpp::as<double>(log_d);
    }
    if (n_loser >= 1) {
      // (summed) log survival probability of losing accumulator(s)
      NumericVector log_s = lccdf(t, pars, !winner, min_ll);
      log_out += std::accumulate(log_s.begin(), log_s.end(), 0.);
    }
    // output is instantaneous (log) likelihood that winner accumulator finishes
    // at time t and the loser accumulator(s) have not yet finished by time t.
    return log_d ? log_out : std::exp(log_out);
  }
};


class race_f_jeroen : public Func {
private:
  Rcpp::NumericMatrix pars;
  Rcpp::LogicalVector winner;
  Rcpp::NumericVector (*dfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector);
  Rcpp::NumericVector (*pfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector);
  double min_ll;
  Rcpp::LogicalVector is_ok;

public:
  race_f_jeroen(Rcpp::NumericMatrix pars_,
        Rcpp::LogicalVector winner_,
        Rcpp::NumericVector (*dfun_)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
        Rcpp::NumericVector (*pfun_)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
        double min_ll_, Rcpp::LogicalVector is_ok_) :
  pars(pars_),
  winner(winner_),
  dfun(dfun_),
  pfun(pfun_),
  is_ok(is_ok_),
  min_ll(min_ll_) {}

  double operator()(const double& x) const
  {
    double accumulators = pars.nrow();
    NumericVector t(accumulators);
    t.fill(x);
    Rcpp::NumericVector d = dfun(t, pars, winner, exp(min_ll), is_ok);
//CHATGPT
//    double out = Rcpp::as<double>(d);
// ---- new guard: empty 'd' means no winner selected ----
    double out;
    if (d.size() == 0) {
        out = 1.0;                    // neutral multiplicative factor
    } else {
        out = Rcpp::as<double>(d);
    }
//CHATGPT

    if (accumulators > 1) {
      Rcpp::NumericVector p = 1 - pfun(t, pars, !winner, exp(min_ll), is_ok);
      double prod_p = std::accumulate(p.begin(), p.end(), 1.0,
                                      std::multiplies<double>());
      out *= prod_p;
    }
    return out;
  }
};


// --------------------------
// NB following code is currently unused; Frank originally developed this for
// the stop signal models, but we ended up using a more lightweight implementation,
// directly calling Numer::integrate instead of the following elaborate wrapper.
struct IntegrationResult {
  double value;
  double error_estimate;
  int error_code;
};

IntegrationResult my_integrate(
    const Func& race_pdf,
    double lower,
    double upper,
    bool check_err = true,
    double upper_retry = 10.,
    double fail_out = 0.,
    int max_subdivs = 100,
    double abs_tol = 1e-8,
    double rel_tol = 1e-6
) {

  // initial attempt at numerical integration
  double err_est;
  int err_code;
  double res = integrate(
    race_pdf, lower, upper,
    err_est, err_code,
    max_subdivs, abs_tol, rel_tol
  );

  if (!check_err) {
    return {res, err_est, err_code};
  }

  // integration could have failed for several reasons (7 different values for
  // err_code possible). If due to max. subdivisions reached (err_code == 1)
  // *AND* upper bound was Inf, try again with large finite upper bound.
  if (err_code == 1 && upper == R_PosInf) {
    double err_est_retry;
    int err_code_retry;
    double res_retry = integrate(
      race_pdf, lower, upper_retry,
      err_est_retry, err_code_retry,
      max_subdivs, abs_tol, rel_tol
    );
    // overwrite output elements
    res = res_retry;
    err_est = err_est_retry;
    err_code = err_code_retry;
  }

  // if integration (still) terminated abnormally, use fallback output
  if (err_code > 0) {
    res = fail_out;
  }

  return {res, err_est, err_code};
}
// --------------------------

//Jeroen
NumericVector f_integrate(Rcpp::NumericMatrix pars,
                          Rcpp::LogicalVector winner,
  Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
  Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                                double min_ll,
                                double lower,
                                double upper,
                                LogicalVector is_ok)
{
  race_f_jeroen f(pars, winner, dfun, pfun, min_ll, is_ok);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code);
  NumericVector out{res, err_est, (double) err_code};
  return out;
}

NumericVector f_integrate_slow(Rcpp::NumericMatrix pars,
                               Rcpp::LogicalVector winner,
  Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
  Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                          double min_ll,
                          double lower,
                          double upper,
                          LogicalVector is_ok)
{
  race_f_jeroen f(pars, winner, dfun, pfun, min_ll, is_ok);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code);
  if (err_code == 1 && upper == R_PosInf) {
    double err_est_hacky;
    int err_code_hacky;
    double res_hacky = integrate(f, lower, 10, err_est_hacky, err_code_hacky);
    NumericVector out{res_hacky, err_est_hacky, (double) err_code_hacky};
    return out;
  } else if (err_code > 2) {
    NumericVector out{min_ll, err_est, (double) err_code};
    return out;
  } else {
    NumericVector out{res, err_est, (double) err_code};
    return out;
  }
}

double pr_pt(Rcpp::NumericMatrix pars,
             Rcpp::LogicalVector winner,
  Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
  Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
             double min_ll,
             double LT,
             double UT,
             LogicalVector is_ok) {
  NumericVector pr = f_integrate(pars, winner, dfun, pfun, min_ll, 0, R_PosInf, is_ok);
  if ((pr[2] != 0) || traits::is_nan<REALSXP>(pr[0])) return NA_REAL;
  if (pr[0] == 0.0) return 0.0;

  NumericVector pt = f_integrate(pars, winner, dfun, pfun, min_ll, LT, UT, is_ok);
  if ((pt[2] != 0) || traits::is_nan<REALSXP>(pt[0])) return NA_REAL;

  double out = std::max(0.0, std::min(pr[0], 1.0)) / std::max(0.0, std::min(pt[0], 1.0));
  if (traits::is_infinite<REALSXP>(out)) return NA_REAL;

  return out;
}


double pLU(Rcpp::NumericMatrix pars,
           Rcpp::LogicalVector winner,
  Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
  Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
           double min_ll,
           double LT,
           double LC,
           double UC,
           double UT,
           LogicalVector is_ok) {
  NumericVector pL = f_integrate(pars, winner, dfun, pfun, min_ll, LT, LC, is_ok);
  if ((pL[2] != 0) || traits::is_nan<REALSXP>(pL[0])) {
    return NA_REAL;
  }
  NumericVector pU = f_integrate(pars, winner, dfun, pfun, min_ll, UC, UT, is_ok);
  if ((pU[2] != 0) || traits::is_nan<REALSXP>(pU[0])) {
    return NA_REAL;
  }
  double out = std::max(0.0, std::min(pL[0], 1.0)) + std::max(0.0, std::min(pU[0], 1.0));
  return out;
}

#endif
