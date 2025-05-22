#ifndef numerical_integration_h
#define numerical_integration_h

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

class race_f : public Func {
private:
  Rcpp::NumericMatrix pars;
  Rcpp::LogicalVector winner;
  Rcpp::NumericVector (*dfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double);
  Rcpp::NumericVector (*pfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double);
  double min_ll;

public:
  race_f(Rcpp::NumericMatrix pars_,
        Rcpp::LogicalVector winner_,
        Rcpp::NumericVector (*dfun_)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double),
        Rcpp::NumericVector (*pfun_)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double),
        double min_ll_) :
  pars(pars_),
  winner(winner_),
  dfun(dfun_),
  pfun(pfun_),
  min_ll(min_ll_) {}

  double operator()(const double& x) const
  {
    double accumulators = pars.nrow();
    NumericVector t(accumulators);
    t.fill(x);
    Rcpp::NumericVector d = dfun(t, pars, winner, exp(min_ll));
    double out = Rcpp::as<double>(d);

    if (accumulators > 1) {
      Rcpp::NumericVector p = 1 - pfun(t, pars, !winner, exp(min_ll));
      double prod_p = std::accumulate(p.begin(), p.end(), 1.0,
                                      std::multiplies<double>());
      out *= prod_p;
    }
    return out;
  }
};

NumericVector f_integrate(Rcpp::NumericMatrix pars,
                                Rcpp::LogicalVector winner,
                                Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
                                Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
                                double min_ll,
                                double lower,
                                double upper)
{
  race_f f(pars, winner, dfun, pfun, min_ll);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code);
  NumericVector out{res, err_est, (double) err_code};
  return out;
}

NumericVector f_integrate_slow(Rcpp::NumericMatrix pars,
                          Rcpp::LogicalVector winner,
                          Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
                          Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
                          double min_ll,
                          double lower,
                          double upper)
{
  race_f f(pars, winner, dfun, pfun, min_ll);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code);
  if (err_code == 1 && upper == R_PosInf) {
    double err_est_hacky;
    int err_code_hacky;
    double res_hacky = integrate(f, lower, 10, err_est_hacky, err_code_hacky);
    NumericVector out{res_hacky, err_est_hacky, (double) err_code_hacky};
    return out;
  } else if (err_code > 0) {
    NumericVector out{min_ll, err_est, (double) err_code};
    return out;
  } else {
    NumericVector out{res, err_est, (double) err_code};
    return out;
  }
}

double pr_pt(Rcpp::NumericMatrix pars,
             Rcpp::LogicalVector winner,
             Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
             Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
             double min_ll,
             double LT,
             double UT) {
  NumericVector pr = f_integrate(pars, winner, dfun, pfun, min_ll, 0, R_PosInf);
  if ((pr[2] != 0) || traits::is_nan<REALSXP>(pr[0])) return NA_REAL;
  if (pr[0] == 0.0) return 0.0;

  NumericVector pt = f_integrate(pars, winner, dfun, pfun, min_ll, LT, UT);
  if ((pt[2] != 0) || traits::is_nan<REALSXP>(pt[0])) return NA_REAL;

  double out = std::max(0.0, std::min(pr[0], 1.0)) / std::max(0.0, std::min(pt[0], 1.0));
  if (traits::is_infinite<REALSXP>(out)) return NA_REAL;

  return out;
}

double pLU(Rcpp::NumericMatrix pars,
           Rcpp::LogicalVector winner,
           Rcpp::NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
           Rcpp::NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
           double min_ll,
           double LT,
           double LC,
           double UC,
           double UT) {
  NumericVector pL = f_integrate(pars, winner, dfun, pfun, min_ll, LT, LC);
  if ((pL[2] != 0) | traits::is_nan<REALSXP>(pL[0])) {
    return NA_REAL;
  }
  NumericVector pU = f_integrate(pars, winner, dfun, pfun, min_ll, UC, UT);
  if ((pU[2] != 0) | traits::is_nan<REALSXP>(pU[0])) {
    return NA_REAL;
  }
  double out = std::max(0.0, std::min(pL[0], 1.0)) + std::max(0.0, std::min(pU[0], 1.0));
  return out;
}

#endif
