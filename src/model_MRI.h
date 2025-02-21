#ifndef mri_H
#define mri_H

#include <Rcpp.h>
using namespace Rcpp;

double c_log_likelihood_MRI(NumericMatrix pars, NumericVector y, LogicalVector is_ok,
                            int n, int m,
                            double min_ll = std::log(1e-10)) {
  NumericVector y_hat(n);
  double sum_yhat = 0.0;

  // Compute row sums of betas and accumulate for overall mean
  for (int i = 0; i < n; i++) {
    double s = 0.0;
    for (int j = 0; j < m - 1; j++) {
      s += pars(i, j);
    }
    y_hat[i] = s;
    sum_yhat += s;
  }

  double mean_y_hat = sum_yhat / n;

  // Center y_hat
  for (int i = 0; i < n; i++) {
    y_hat[i] -= mean_y_hat;
  }

  NumericVector ll(n);

  // Compute log likelihood using R::dnorm; each y[i] is assumed to be independent
  for (int i = 0; i < n; i++) {
    if(is_ok[i] == FALSE){
      ll[i] = R_NegInf;
    } else{
      // sigma is in the last column of pars for row i
      double sigma = pars(i, m - 1);
      ll[i] = R::dnorm(y[i], y_hat[i], sigma, true); // last argument: give_log = 1
    }
  }
  ll[is_na(ll)] = min_ll;
  ll[is_infinite(ll)] = min_ll;
  ll[ll < min_ll] = min_ll;
  return sum(ll);
}
//
double c_log_likelihood_MRI_white(NumericMatrix pars, NumericVector y, LogicalVector is_ok,
                                        int n, int m,
                                        double min_ll = std::log(1e-10)) {
  NumericVector y_hat(n);
  double sum_yhat = 0.0;
  // Compute row sums of the betas (assumed to be in the first m-2 columns)
  for (int i = 0; i < n; i++) {
    double s = 0.0;
    for (int j = 0; j < m - 2; j++) {
      s += pars(i, j);
    }
    y_hat[i] = s;
    sum_yhat += s;
  }

  double mean_y_hat = sum_yhat / n;

  // Center y_hat (i.e. omit the intercept by demeaning)
  for (int i = 0; i < n; i++) {
    y_hat[i] -= mean_y_hat;
  }

  NumericVector ll(n);

  // For the first observation: use the stationary variance directly.
  // sigma is assumed to be the stationary standard deviation (last column of pars).
  if (!is_ok[0]) {
    ll[0] = R_NegInf;
  } else {
    double sigma0 = pars(0, m - 1);
    ll[0] = R::dnorm(y[0], y_hat[0], sigma0, true);
  }

  // For t >= 2, use the AR(1) conditional likelihood.
  // The conditional variance is sigma^2 * (1 - rho^2) and
  // the conditional mean is y_hat[t] + rho*(y[t-1] - y_hat[t-1]).
  for (int i = 1; i < n; i++) {
    if (!is_ok[i]) {
      ll[i] = R_NegInf;
    } else {
      double sigma_i = pars(i, m - 1);
      double rho_i = pars(i, m-2);
      double cond_sd = sigma_i * std::sqrt(1 - rho_i * rho_i);
      double cond_mean = y_hat[i] + rho_i * (y[i - 1] - y_hat[i - 1]);
      ll[i] = R::dnorm(y[i], cond_mean, cond_sd, true);
    }
  }

  // Replace NA, -Inf, or very low likelihoods with min_ll.
  for (int i = 0; i < n; i++) {
    if (ISNAN(ll[i]) || !R_finite(ll[i]) || ll[i] < min_ll){
      ll[i] = min_ll;
    }
  }
  // Sum the log likelihoods over all observations.;
  return sum(ll);
}
//
//
NumericVector extract_y(DataFrame data) {
  CharacterVector names = data.names();
  int n = data.nrows();

  // Loop through the column names and return the first column that is not excluded
  for (int j = 0; j < names.size(); j++) {
    std::string nm = Rcpp::as<std::string>(names[j]);
    if (nm != "subjects" && nm != "run" && nm != "time" && nm != "trials") {
      return as<NumericVector>(data[nm]);
    }
  }

  // If no column is found, return an empty vector.
  return NumericVector(0);
}


#endif

