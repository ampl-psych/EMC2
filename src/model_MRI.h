#ifndef mri_H
#define mri_H

#include <Rcpp.h>
using namespace Rcpp;
//
// #include "mri.h"

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

  // LL calculation
  for (int i = 0; i < n; i++) {
    if(is_ok[i] == FALSE){
      ll[i] = R_NegInf;
    } else{
      // sigma is in the last column of pars for row i
      double sigma = pars(i, m - 1);
      ll[i] = R::dnorm(y[i], y_hat[i], sigma, true);
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

// // [[Rcpp::export]]
// double log_likelihood_double_gamma(NumericVector y,
//                                    NumericVector parameters,
//                                    double tr,
//                                    NumericVector frame_times,
//                                    List event_cache) {
//   // Define hyperparameters as constants.
//   int oversampling = 50;
//   double time_length = 32.0;
//   double onset = 0.0;
//   double undershoot = 12.0;
//   double dispersion = 0.9;
//   double u_dispersion = 0.9;
//   double ratio = 0.35;
//   double min_onset = -24.0;
//
//   // event_cache is a list with one element per condition.
//   int num_conditions = event_cache.size();
//
//   // Parameter vector format: [ beta_1, beta_2, ..., beta_m, free_delay, sigma ]
//   int total_params = parameters.size();
//   double free_delay = parameters[total_params - 2];
//   double sigma = parameters[total_params - 1];
//   NumericVector beta = parameters[Range(0, total_params - 3)];
//   if(total_params < 2)
//     stop("Parameter vector must contain beta weights, free delay, and sigma.");
//
//   double free_delay = parameters[total_params - 2];
//   double sigma = parameters[total_params - 1];
//   NumericVector beta = parameters[Range(0, total_params - 3)];
//   if(beta.size() != num_conditions)
//     stop("Number of beta weights must equal the number of conditions in event_cache.");
//
//   // Precompute common derived quantities.
//   int n_scanner = frame_times.size();  // number of scanner time points
//
//   // Recompute the current HRF kernel (for the "glover" model) using the free_delay.
//   // We use our updated build_hrf_kernel function.
//   NumericMatrix hkernel = build_hrf_kernel("glover", tr, oversampling, time_length, onset,
//                                            free_delay, undershoot, dispersion, u_dispersion, ratio);
//   // For "glover", the kernel is one-column. Reverse it.
//   NumericVector hkernel_col = hkernel(_, 0);
//   NumericVector hkernel_rev = reverse_vector(hkernel_col);
//
//   // All cached regressor FFTs should have the same length.
//   // Use the first condition to obtain the padded FFT length.
//   List first_cache = event_cache[0];
//   arma::cx_vec reg_fft_example = as<arma::cx_vec>(first_cache["regressor_fft"]);
//   int L = reg_fft_example.n_elem;
//
//   // Zero-pad the current HRF kernel to length L.
//   NumericVector kernel_padded(L);
//   int kernel_length = hkernel_rev.size();
//   for (int i = 0; i < kernel_length; i++) {
//     kernel_padded[i] = hkernel_rev[i];
//   }
//   for (int i = kernel_length; i < L; i++) {
//     kernel_padded[i] = 0.0;
//   }
//   arma::vec kernel_vec = as<arma::vec>(kernel_padded);
//   arma::cx_vec kernel_fft = arma::fft(conv_to<cx_vec>::from(kernel_vec));
//
//   // For each condition, use the cached FFT of the high-res event regressor.
//   // Multiply by the current HRF kernel FFT, perform inverse FFT, and downsample.
//   arma::mat X(n_scanner, num_conditions, fill::zeros);
//
//   for (int cond = 0; cond < num_conditions; cond++) {
//     List cond_cache = event_cache[cond];
//     NumericVector cond_hr_regressor = cond_cache["hr_regressor"];     // high-res regressor
//     NumericVector cond_hr_frame_times = cond_cache["hr_frame_times"];   // high-res time grid
//     arma::cx_vec regressor_fft = as<arma::cx_vec>(cond_cache["regressor_fft"]); // precomputed FFT
//
//     // Convolve: multiply cached FFT with current HRF kernel FFT.
//     arma::cx_vec conv_fft = regressor_fft % kernel_fft;
//     arma::cx_vec conv_ifft = arma::ifft(conv_fft);
//     arma::vec conv_full = arma::real(conv_ifft.head(cond_hr_regressor.size()));
//
//     // Downsample conv_full from cond_hr_frame_times to scanner frame_times.
//     // Wrap conv_full as a one-column NumericMatrix.
//     NumericVector conv_full_vec = wrap(conv_full);
//     NumericMatrix conv_full_mat(conv_full_vec.size(), 1);
//     for (int i = 0; i < conv_full_vec.size(); i++) {
//       conv_full_mat(i, 0) = conv_full_vec[i];
//     }
//     // Downsample using resample_matrix.
//     NumericMatrix ds = resample_matrix(conv_full_mat, cond_hr_frame_times, frame_times);
//     arma::vec ds_vec = as<arma::vec>(ds(_,0));
//     X.col(cond) = ds_vec;
//   }
//
//   // Compute the model prediction: X * beta.
//   arma::vec beta_vec = as<arma::vec>(beta);
//   arma::vec pred = X * beta_vec;
//
//   // Compute the log-likelihood: sum over time points of log-density.
//   arma::vec y_vec = as<arma::vec>(y);
//   int n = y_vec.n_elem;
//   double logLik = 0.0;
//   for (int i = 0; i < n; i++) {
//     logLik += R::dnorm(y_vec(i), pred[i], sigma, 1);
//   }
//
//   return logLik;
// }

#endif

