// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <set>
#include <vector>
#include <string>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat calculate_subject_means(const Rcpp::List& group_designs,
                                  const arma::colvec& params,
                                  const int n_subjects,
                                  const int n_pars) {

  // output matrix: rows = parameters, cols = subjects
  arma::mat subj_mu(n_pars, n_subjects, arma::fill::zeros);

  for (int s = 0; s < n_subjects; ++s) {
    int par_idx = 0;                       // reset for each subject
    for (int k = 0; k < n_pars; ++k) {

      // pull the k‑th design matrix and grab the s‑th row
      arma::mat   Xk   = group_designs[k];
      arma::rowvec xsk = Xk.row(s);
      int          p   = Xk.n_cols;

      // slice the relevant segment of params
      arma::colvec beta = params.subvec(par_idx, par_idx + p - 1);

      // x_sk %*% β_k
      subj_mu(k, s) = arma::as_scalar(xsk * beta);

      par_idx += p;
    }
  }
  return subj_mu;
}

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

static double const log2pi = std::log(2.0 * M_PI);

// ---- helper: one observation (row-vector) ----
inline double dmvnorm_row(const arma::rowvec& x,
                          const arma::rowvec& mean,
                          const arma::mat&    sigma,
                          const bool          logd = false)
{
  arma::mat  rooti = arma::inv( arma::trimatu( arma::chol(sigma) ) );
  arma::rowvec z   = (x - mean) * rooti;
  double logdens   = arma::sum( log(rooti.diag()) )
    - 0.5 * arma::dot(z, z)
    - 0.5 * x.n_elem * log2pi;
    return logd ? logdens : std::exp(logdens);
}

// ---- exported worker: n × p matrix in  ----
// [[Rcpp::export(name = "dmvnorm_cpp")]]
arma::vec dmvnorm_mat(const arma::mat&    x,
                      const arma::rowvec& mean,
                      const arma::mat&    sigma,
                      const bool          logd = false)
{
  const arma::uword n = x.n_rows;
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i)
    out(i) = dmvnorm_row(x.row(i), mean, sigma, /*logd=*/true);

  return logd ? out : arma::exp(out);
}


// [[Rcpp::export]]
arma::vec standard_subj_ll(const Rcpp::List& group_designs,
                           const arma::cube& theta_var,
                           const arma::mat&  theta_mu,
                           const arma::cube& alpha,
                           const int         n_subj)
{
  const int N = theta_var.n_slices;     // draws
  const int p = theta_mu.n_rows;        // parameters

  arma::vec lls(N, arma::fill::zeros);

  for (int i = 0; i < N; ++i) {

    const arma::mat  var_i   = theta_var.slice(i);   // p × p
    const arma::colvec reg_i = theta_mu.col(i);      // p × 1

    arma::mat subj_means = calculate_subject_means(
      group_designs, reg_i, n_subj, p); // p × n_subj
    const arma::mat alpha_i = alpha.slice(i);        // p × n_subj

    double ll_i = 0.0;
    for (int s = 0; s < n_subj; ++s) {
      arma::rowvec alpha_row = alpha_i.col(s).t();   // 1 × p
      ll_i += dmvnorm_row(alpha_row,
                          subj_means.col(s).t(),
                          var_i, true);
    }
    lls(i) = ll_i;
  }
  return lls;
}



// ---------------------------------------------------------
// Helper: next power of 2 (needed for FFT convolution)
// ---------------------------------------------------------
int next_power_of_two(int n) {
  int p = 1;
  while (p < n) p *= 2;
  return p;
}

// ---------------------------------------------------------
// Provided FFT-based convolution function.
// This exactly replicates R's convolve(x, rev(y), type="o").
// ---------------------------------------------------------
// [[Rcpp::export]]
arma::vec fft_convolve_equiv_cpp(const arma::vec &x, const arma::vec &y, bool conj_flag = true) {
  int n = x.n_elem;      // length(x)
  int ny = y.n_elem;     // length(y)
  int L = n + ny - 1;    // full convolution length

  // Pad x at the beginning with (ny - 1) zeros.
  arma::vec x_pad = join_vert(arma::zeros(ny - 1), x);
  // Pad y at the end with (n - 1) zeros.
  arma::vec y_pad = join_vert(y, arma::zeros(n - 1));

  // Zero-pad both to length M = next power of 2 >= L.
  int M = next_power_of_two(L);
  arma::vec X_input = join_vert(x_pad, arma::zeros(M - L));
  arma::vec Y_input = join_vert(y_pad, arma::zeros(M - L));

  // Compute FFTs.
  arma::cx_vec X_fft = arma::fft(conv_to<cx_vec>::from(X_input));
  arma::cx_vec Y_fft = arma::fft(conv_to<cx_vec>::from(Y_input));

  // Multiply elementwise.
  arma::cx_vec prod_fft;
  if (conj_flag)
    prod_fft = X_fft % conj(Y_fft);
  else
    prod_fft = X_fft % Y_fft;

  // Compute inverse FFT and normalize by L.
  arma::cx_vec conv_complex = arma::ifft(prod_fft);
  arma::vec conv_result = arma::real(conv_complex.head(L));
  return conv_result;
}

// ---------------------------------------------------------
// Helper: search_sorted.
// Returns the index (0-indexed) of the first element in x not less than v,
// mimicking R's findInterval(..., left.open=TRUE) (with adjustment).
// ---------------------------------------------------------
int search_sorted(const NumericVector &x, double v) {
  int n = x.size();
  int low = 0, high = n;
  while (low < high) {
    int mid = (low + high) / 2;
    if (x[mid] < v)
      low = mid + 1;
    else
      high = mid;
  }
  return low;
}

// ---------------------------------------------------------
// Helper: generate a sequence of n equally‐spaced numbers
// ---------------------------------------------------------
NumericVector seq_lin(double start, double end, int n) {
  NumericVector out(n);
  if(n == 1) {
    out[0] = start;
  } else {
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++) {
      out[i] = start + i * step;
    }
  }
  return out;
}

List sample_event_condition(NumericMatrix exp_condition, NumericVector frame_times,
                            int oversampling = 50, double min_onset = -24) {
  int n_frames = frame_times.size();
  double tmin = Rcpp::min(frame_times);
  double tmax = Rcpp::max(frame_times);
  // Compute number of high-resolution time points.
  double n_hr_calc = ((n_frames - 1) / (tmax - tmin)) *
    (tmax * (1 + 1.0 / (n_frames - 1)) - tmin - min_onset) * oversampling + 1;
  int n_hr = std::max(1, (int) std::round(n_hr_calc));

  NumericVector hr_frame_times = seq_lin(tmin + min_onset, tmax * (1 + 1.0 / (n_frames - 1)), n_hr);

  int n_events = exp_condition.nrow();
  NumericVector onsets(n_events), durations(n_events), values(n_events);
  for (int i = 0; i < n_events; i++) {
    onsets[i] = exp_condition(i, 0);
    durations[i] = exp_condition(i, 1);
    values[i] = exp_condition(i, 2);
  }

  int tmax_idx = hr_frame_times.size();
  NumericVector regressor(tmax_idx, 0.0);

  // Compute t_onset for each event.
  std::vector<int> t_onset(n_events);
  for (int i = 0; i < n_events; i++) {
    int idx = search_sorted(hr_frame_times, onsets[i]);
    if (idx > tmax_idx - 1) idx = tmax_idx - 1;
    t_onset[i] = idx;
    regressor[idx] += values[i];
  }

  // Compute t_offset for each event.
  std::vector<int> t_offset(n_events);
  for (int i = 0; i < n_events; i++) {
    int idx = search_sorted(hr_frame_times, onsets[i] + durations[i]);
    if (idx > tmax_idx - 1) idx = tmax_idx - 1;
    t_offset[i] = idx;
    if ((t_offset[i] < tmax_idx - 1) && (t_offset[i] == t_onset[i]))
      t_offset[i] = t_offset[i] + 1;
  }
  for (int i = 0; i < n_events; i++) {
    regressor[t_offset[i]] -= values[i];
  }

  // Compute the cumulative sum.
  for (int i = 1; i < regressor.size(); i++) {
    regressor[i] += regressor[i - 1];
  }

  return List::create(Named("regressor") = regressor,
                      Named("hr_frame_times") = hr_frame_times);
}


// ---------------------------------------------------------
// Helper: linear interpolation for a single value.
// (Implements rule=2: values outside the range get the boundary value.)
// ---------------------------------------------------------
double lininterp(const NumericVector &x, const NumericVector &y, double x0) {
  int n = x.size();
  if(x0 <= x[0]) return y[0];
  if(x0 >= x[n-1]) return y[n-1];
  int low = 0, high = n - 1;
  while(high - low > 1) {
    int mid = (low + high) / 2;
    if(x[mid] > x0)
      high = mid;
    else
      low = mid;
  }
  double t = (x0 - x[low]) / (x[high] - x[low]);
  return y[low] + t * (y[high] - y[low]);
}

// ---------------------------------------------------------
// Helper: resample (linearly interpolate) a vector.
// ---------------------------------------------------------
NumericVector resample_vector(const NumericVector &orig_x, const NumericVector &orig_y, const NumericVector &new_x) {
  int n_new = new_x.size();
  NumericVector result(n_new);
  for (int i = 0; i < n_new; i++) {
    result[i] = lininterp(orig_x, orig_y, new_x[i]);
  }
  return result;
}

// ---------------------------------------------------------
// Helper: resample a matrix column‐by‐column.
// ---------------------------------------------------------
NumericMatrix resample_matrix(const NumericMatrix &Y, const NumericVector &orig_x, const NumericVector &new_x) {
  int n_new = new_x.size();
  int ncol = Y.ncol();
  NumericMatrix result(n_new, ncol);
  for (int j = 0; j < ncol; j++) {
    NumericVector col = Y(_, j);
    NumericVector resampled = resample_vector(orig_x, col, new_x);
    for (int i = 0; i < n_new; i++) {
      result(i, j) = resampled[i];
    }
  }
  return result;
}

// ---------------------------------------------------------
// A generic resample_signal function that works for both vector and matrix input.
// ---------------------------------------------------------
SEXP resample_signal(SEXP signal, NumericVector orig_x, NumericVector new_x) {
  if (Rf_isMatrix(signal)) {
    NumericMatrix Y(signal);
    return resample_matrix(Y, orig_x, new_x);
  } else {
    NumericVector y(signal);
    return resample_vector(orig_x, y, new_x);
  }
}

// ---------------------------------------------------------
// resample_regressor: simply calls resample_signal.
// ---------------------------------------------------------
SEXP resample_regressor(SEXP hr_regressor, NumericVector hr_frame_times, NumericVector frame_times) {
  return resample_signal(hr_regressor, hr_frame_times, frame_times);
}

// =========================================================
// Updated HRF and Convolution Functions with Hyperparameters as Arguments
// =========================================================

// ----------------------------------------------
// Updated: compute_gamma_diff_hrf
// Now all HRF parameters are passed explicitly.
// ----------------------------------------------
// [[Rcpp::export]]
NumericVector compute_gamma_diff_hrf(double tr, int oversampling, double time_length, double onset,
                                     double delay, double undershoot, double dispersion,
                                     double u_dispersion, double ratio) {
  // 1) Determine time step dt and number of points
  double dt = tr / static_cast<double>(oversampling);
  int n_points = static_cast<int>(std::round(time_length / dt));
  if (n_points < 2) {
    n_points = 2;
  }

  // 2) Generate time_stamps from 0 .. (n_points-1)*dt
  Rcpp::NumericVector time_stamps(n_points);
  for (int i = 0; i < n_points; i++){
    time_stamps[i] = i * dt;
  }

  // 3) Subtract onset (same as Nilearn)
  for (int i = 0; i < n_points; i++){
    time_stamps[i] -= onset;
  }

  // 4) Nilearn’s gamma.pdf(...) calls actually do:
  //    shape = (delay / dispersion), loc = (dt / dispersion), scale = 1
  //    shape = (undershoot / u_dispersion), loc = (dt / u_dispersion), scale = 1
  //    We must manually shift t by loc and use scale=1 in R's dgamma().
  double loc_peak   = dt / dispersion;
  double loc_under  = dt / u_dispersion;

  Rcpp::NumericVector hrf(n_points, 0.0);

  for (int i = 0; i < n_points; i++){
    double tval = time_stamps[i];

    // Peak gamma
    double peak_val = 0.0;
    if (tval >= loc_peak) {
      // shape = (delay / dispersion), scale=1, argument=(tval - loc_peak)
      peak_val = R::dgamma(tval - loc_peak, delay / dispersion, /*scale=*/1.0, false);
    }

    // Undershoot gamma
    double under_val = 0.0;
    if (tval >= loc_under) {
      // shape = (undershoot / u_dispersion), scale=1, argument=(tval - loc_under)
      under_val = R::dgamma(tval - loc_under, undershoot / u_dispersion, /*scale=*/1.0, false);
    }

    hrf[i] = peak_val - ratio * under_val;
  }

  // 5) Normalize so the sum of the HRF = 1
  // double sum_hrf = std::accumulate(hrf.begin(), hrf.end(), 0.0);
  double max_hrf = max(hrf);
  if (max_hrf != 0.0) {
    for (int i = 0; i < n_points; i++){
      hrf[i] /= max_hrf;
    }
  }

  return hrf;
}

// ----------------------------------------------
// Updated: compute_glover_hrf
// Simply calls compute_gamma_diff_hrf with all parameters.
// ----------------------------------------------
// [[Rcpp::export]]
NumericVector compute_hrf(double tr, int oversampling, double time_length, double onset,
                                 double delay, double undershoot, double dispersion,
                                 double u_dispersion, double ratio) {
  return compute_gamma_diff_hrf(tr, oversampling, time_length, onset,
                                delay, undershoot, dispersion, u_dispersion, ratio);
}

// ----------------------------------------------
// Updated: compute_glover_time_derivative
// Now accepts additional HRF parameters and delta.
// ----------------------------------------------
// [[Rcpp::export]]
NumericVector compute_time_derivative(double tr, int oversampling, double time_length,
                                             double onset, double delay, double undershoot,
                                             double dispersion, double u_dispersion, double ratio,
                                             double delta = 0.1) {
  NumericVector hrf1 = compute_hrf(tr, oversampling, time_length, onset,
                                          delay, undershoot, dispersion, u_dispersion, ratio);
  NumericVector hrf2 = compute_hrf(tr, oversampling, time_length, onset + delta,
                                          delay, undershoot, dispersion, u_dispersion, ratio);
  int n = hrf1.size();
  NumericVector deriv(n);
  for (int i = 0; i < n; i++) {
    deriv[i] = (hrf1[i] - hrf2[i]) / delta;
  }
  return deriv;
}

NumericVector reverse_vector(const NumericVector &v) {
  int n = v.size();
  NumericVector rev(n);
  for (int i = 0; i < n; i++) {
    rev[i] = v[n - 1 - i];
  }
  return rev;
}


// ----------------------------------------------
// Updated: build_hrf_kernel
// Now requires HRF parameters as arguments.
// ----------------------------------------------
// [[Rcpp::export]]
NumericMatrix build_hrf_kernel(bool has_derivative, double tr, int oversampling,
                               double time_length, double onset, double delay,
                               double undershoot, double dispersion, double u_dispersion,
                               double ratio) {
  if (!has_derivative) {
    NumericVector hrf = compute_hrf(tr, oversampling, time_length, onset,
                                           delay, undershoot, dispersion, u_dispersion, ratio);
    NumericMatrix kernel(hrf.size(), 1);
    for (int i = 0; i < hrf.size(); i++) {
      kernel(i, 0) = hrf[i];
    }
    return kernel;
  } else {
    NumericVector hrf = compute_hrf(tr, oversampling, time_length, onset,
                                           delay, undershoot, dispersion, u_dispersion, ratio);
    NumericVector deriv = compute_time_derivative(tr, oversampling, time_length, onset,
                                                         delay, undershoot, dispersion, u_dispersion, ratio);
    int n = hrf.size();
    NumericMatrix kernel(n, 2);
    for (int i = 0; i < n; i++) {
      kernel(i, 0) = hrf[i];
      kernel(i, 1) = deriv[i];
    }
    return kernel;
  }
  return NumericMatrix(0);
}

// ----------------------------------------------
// Updated: compute_convolved_regressor
// Now accepts HRF parameters and passes them to build_hrf_kernel.
// ----------------------------------------------
List compute_convolved_regressor(NumericMatrix exp_condition, bool has_derivative,
                                 NumericVector frame_times, std::string con_id,
                                 int oversampling, double min_onset,
                                 double time_length, double onset, double delay,
                                 double undershoot, double dispersion, double u_dispersion,
                                 double ratio) {
  double tmin = Rcpp::min(frame_times);
  double tmax = Rcpp::max(frame_times);
  int n_frames = frame_times.size();
  double tr = (tmax - tmin) / (n_frames - 1);

  // High-resolution sampling.
  List cond_sample = sample_event_condition(exp_condition, frame_times, oversampling, min_onset);
  NumericVector hr_regressor = cond_sample["regressor"];
  NumericVector hr_frame_times = cond_sample["hr_frame_times"];

  // Build HRF kernel using updated parameters.
  NumericMatrix hkernel = build_hrf_kernel(has_derivative, tr, oversampling, time_length, onset,
                                           delay, undershoot, dispersion, u_dispersion, ratio);
  int n_basis = hkernel.ncol();
  int n_hr = hr_regressor.size();
  NumericMatrix conv_mat_full(n_hr, n_basis);
  for (int j = 0; j < n_basis; j++) {
    NumericVector hkernel_col = hkernel(_, j);
    NumericVector hkernel_rev = reverse_vector(hkernel_col);
    arma::vec hr_reg = as<arma::vec>(hr_regressor);
    arma::vec hrf_kernel = as<arma::vec>(hkernel_rev);
    arma::vec conv_full = fft_convolve_equiv_cpp(hr_reg, hrf_kernel, true);
    int L = hr_regressor.size();
    NumericVector conv_res(conv_full.memptr(), conv_full.memptr() + L);
    for (int i = 0; i < L; i++) {
      conv_mat_full(i, j) = conv_res[i];
    }
  }

  // Downsample to frame_times.
  NumericMatrix computed_regressors = resample_matrix(conv_mat_full, hr_frame_times, frame_times);

  // Gram-Schmidt orthogonalization.
  int ncols = computed_regressors.ncol();
  if (ncols > 1) {
    arma::mat X = as<arma::mat>(computed_regressors);
    int p = X.n_cols;
    for (int j = 1; j < p; j++) {
      for (int i = 0; i < j; i++) {
        double denom = dot(X.col(i), X.col(i));
        if (denom != 0) {
          double beta_val = dot(X.col(j), X.col(i)) / denom;
          X.col(j) -= beta_val * X.col(i);
        }
      }
    }
    computed_regressors = wrap(X);
  }

  // Set regressor names following the original behavior.
  CharacterVector reg_names;
  if (!has_derivative) {
    reg_names = CharacterVector::create(con_id);
  } else {
    reg_names = CharacterVector::create(con_id, con_id + std::string("_derivative"));
  }

  return List::create(Named("computed_regressors") = computed_regressors,
                      Named("regressor_names") = reg_names);
}

// ----------------------------------------------
// Updated: construct_design_matrix
// Now passes all HRF parameters as arguments.
// ----------------------------------------------
// [[Rcpp::export]]
DataFrame construct_design_matrix(NumericVector frame_times, DataFrame events,
                                  bool has_derivative, double min_onset, int oversampling,
                                  double time_length, double onset, double delay,
                                  double undershoot, double dispersion, double u_dispersion,
                                  double ratio, bool add_intercept) {
  // Extract event columns.
  NumericVector onset_vec = events["onset"];
  NumericVector duration = events["duration"];
  NumericVector modulation;
  if (events.containsElementNamed("modulation"))
    modulation = events["modulation"];
  else
    modulation = NumericVector(onset_vec.size(), 1.0);
  CharacterVector regressor_vec = events["regressor"];
  int n_events = regressor_vec.size();

  // Determine unique condition labels.
  std::set<std::string> cond_set;
  for (int i = 0; i < n_events; i++) {
    cond_set.insert(as<std::string>(regressor_vec[i]));
  }
  std::vector<std::string> conds(cond_set.begin(), cond_set.end());

  bool first = true;
  NumericMatrix regressor_matrix;
  std::vector<std::string> reg_names_all;
  for (size_t i = 0; i < conds.size(); i++) {
    std::string cond = conds[i];
    std::vector<int> indices;
    for (int j = 0; j < n_events; j++) {
      if (as<std::string>(regressor_vec[j]) == cond)
        indices.push_back(j);
    }
    int n_cond = indices.size();
    if (n_cond == 0) continue;
    // Build exp_condition: columns = onset, duration, modulation.
    NumericMatrix exp_condition(n_cond, 3);
    for (int j = 0; j < n_cond; j++) {
      int idx = indices[j];
      exp_condition(j, 0) = onset_vec[idx];
      exp_condition(j, 1) = duration[idx];
      exp_condition(j, 2) = modulation[idx];
    }
    // Compute convolved regressor for this condition.
    List out = compute_convolved_regressor(exp_condition, has_derivative, frame_times, cond, oversampling, min_onset,
                                           time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio);
    NumericMatrix reg = out["computed_regressors"];
    CharacterVector names = out["regressor_names"];

    // Append regressor columns.
    if (first) {
      regressor_matrix = reg;
      for (int k = 0; k < names.size(); k++) {
        reg_names_all.push_back(as<std::string>(names[k]));
      }
      first = false;
    } else {
      int n_rows = regressor_matrix.nrow();
      int n_old = regressor_matrix.ncol();
      int n_new = reg.ncol();
      NumericMatrix temp(n_rows, n_old + n_new);
      for (int r = 0; r < n_rows; r++) {
        for (int c = 0; c < n_old; c++) {
          temp(r, c) = regressor_matrix(r, c);
        }
      }
      for (int r = 0; r < n_rows; r++) {
        for (int c = 0; c < n_new; c++) {
          temp(r, n_old + c) = reg(r, c);
        }
      }
      regressor_matrix = temp;
      for (int k = 0; k < names.size(); k++) {
        reg_names_all.push_back(as<std::string>(names[k]));
      }
    }
  }

  // Create a data frame from the regressor matrix.
  DataFrame df = as<DataFrame>(wrap(regressor_matrix));
  df.attr("names") = reg_names_all;

  if (add_intercept) {
    int n_rows = frame_times.size();
    NumericVector intercept(n_rows, 1.0);
    df.push_back(intercept, "intercept");
  }

  // Set row names as the string conversion of frame_times.
  CharacterVector rn(frame_times.size());
  for (int i = 0; i < frame_times.size(); i++) {
    rn[i] = std::to_string(frame_times[i]);
  }
  df.attr("row.names") = rn;

  return df;
}



//
// // [[Rcpp::export]]
// NumericVector build_glover_hrf_kernel_numeric(double tr, int oversampling, double time_length, double onset,
//                                               double delay, double undershoot, double dispersion,
//                                               double u_dispersion, double ratio) {
//   // Compute the canonical Glover HRF using the provided hyperparameters.
//   NumericVector hrf = compute_hrf(tr, oversampling, time_length, onset,
//                                          delay, undershoot, dispersion, u_dispersion, ratio);
//   // For the Glover model, hrf is a single vector.
//   // Reverse the HRF kernel so that it can be used for convolution.
//   NumericVector hrf_rev = reverse_vector(hrf);
//   return hrf_rev;
// }
//
//
// // This function groups the events by condition and, for each condition,
// // computes the exp_condition matrix (onset, duration, modulation) and precomputes
// // the FFT of the zero‑padded high-resolution event regressor. It stores only the
// // TR, condition name, exp_condition, and regressor_fft in the cache.
//
// // [[Rcpp::export]]
// List build_event_cache_cpp(DataFrame events, NumericVector run_times,
//                            int oversampling = 50,
//                            double min_onset = -24.0,
//                            double time_length = 32.0,
//                            double onset = 0.0,
//                            double nominal_delay = 6.0,
//                            double undershoot = 12.0,
//                            double dispersion = 0.9,
//                            double u_dispersion = 0.9,
//                            double ratio = 0.35) {
//   // Compute TR from run_times.
//   double tr = run_times[1] - run_times[0];
//   int n_events = events.nrows();
//
//   // Get unique condition names from the "regressor" column.
//   CharacterVector regressor_col = events["regressor"];
//   std::set<std::string> cond_set;
//   for (int i = 0; i < n_events; i++) {
//     cond_set.insert(as<std::string>(regressor_col[i]));
//   }
//   std::vector<std::string> conditions(cond_set.begin(), cond_set.end());
//
//   List cache_list(conditions.size());
//
//   for (size_t j = 0; j < conditions.size(); j++) {
//     std::string cond = conditions[j];
//     // Subset events for this condition.
//     std::vector<int> idx;
//     for (int i = 0; i < n_events; i++) {
//       if(as<std::string>(regressor_col[i]) == cond)
//         idx.push_back(i);
//     }
//     int n_cond = idx.size();
//
//     // Build exp_condition matrix with columns: onset, duration, modulation.
//     NumericVector onset_vec = events["onset"];
//     NumericVector duration_vec = events["duration"];
//     NumericVector modulation;
//     if (events.containsElementNamed("modulation"))
//       modulation = events["modulation"];
//     else
//       modulation = NumericVector(n_events, 1.0);
//
//     NumericMatrix exp_condition(n_cond, 3);
//     for (int i = 0; i < n_cond; i++) {
//       int row = idx[i];
//       exp_condition(i, 0) = onset_vec[row];
//       exp_condition(i, 1) = duration_vec[row];
//       exp_condition(i, 2) = modulation[row];
//     }
//
//     // Compute high-resolution regressor and grid using sample_event_condition.
//     List sample_out = sample_event_condition(exp_condition, run_times, oversampling, min_onset);
//     NumericVector hr_regressor = sample_out["regressor"];
//     // M = length of hr_regressor.
//     int M = hr_regressor.size();
//
//     // Determine padded length for FFT using a nominal HRF kernel.
//     NumericVector hkernel_rev = build_glover_hrf_kernel_numeric(tr, oversampling, time_length, onset,
//                                                                 nominal_delay, undershoot, dispersion,
//                                                                 u_dispersion, ratio);
//     int kernel_length = hkernel_rev.size();
//     int padded_length = next_power_of_two(M + kernel_length - 1);
//
//     // Zero-pad hr_regressor using Armadillo's join_vert.
//     arma::vec hr_regressor_vec = as<arma::vec>(hr_regressor);
//     arma::vec zeros_pad = arma::zeros<arma::vec>(padded_length - M);
//     arma::vec padded_regressor = join_vert(hr_regressor_vec, zeros_pad);
//     arma::cx_vec regressor_fft = arma::fft(padded_regressor);
//
//     // Store TR, condition name, exp_condition, and regressor_fft.
//     List cond_cache = List::create(
//       Named("tr") = tr,
//       Named("cond_name") = cond,
//       Named("exp_condition") = exp_condition,
//       Named("regressor_fft") = regressor_fft
//     );
//
//     cache_list[j] = cond_cache;
//   }
//
//   return cache_list;
// }
//

// // [[Rcpp::export]]
// double log_likelihood_double_gamma(NumericVector y,
//                                    NumericVector parameters,
//                                    List event_cache) {
//   // hyperparameters
//   int oversampling = 50;
//   double time_length = 32.0;
//   double onset = 0.0;
//   double undershoot = 12.0;
//   double dispersion = 0.9;
//   double u_dispersion = 0.9;
//   double ratio = 0.35;
//   double min_onset = -24.0;
//
//   // TR stored in first condition
//   List firstCache = event_cache[0];
//   double tr = as<double>(firstCache["tr"]);
//
//   // Compute scanner frame_times (starting at 0).
//   int n_scanner = y.size();
//   NumericVector scanner_frame_times = seq_lin(0, tr * (n_scanner - 1), n_scanner);
//
//   // Number of conditions.
//   int num_conditions = event_cache.size();
//
//   // Betas first, then any additional parameters (sigma, delay, undershoot?, rho?)
//   int total_params = parameters.size();
//   double free_delay = parameters[total_params - 2];
//   double sigma = parameters[total_params - 1];
//   NumericVector beta = parameters[Range(0, total_params - 3)];
//   // Compute the current HRF kernel using free_delay.
//   NumericVector hkernel_rev = build_glover_hrf_kernel_numeric(tr, oversampling, time_length, onset,
//                                                               free_delay, undershoot, dispersion, u_dispersion, ratio);
//
//   // Retrieve padded length from the first condition's cached FFT.
//   arma::cx_vec first_fft = as<arma::cx_vec>(firstCache["regressor_fft"]);
//   int padded_length = first_fft.n_elem;
//   int kernel_length = as<arma::vec>(hkernel_rev).n_elem;
//
//   // Zero-pad HRF kernel using join_vert.
//   arma::vec padded_kernel = join_vert(as<arma::vec>(hkernel_rev), arma::zeros(padded_length - kernel_length));
//   arma::cx_vec kernel_fft = arma::fft(padded_kernel);
//
//   // Initialize prediction vector.
//   NumericVector pred(n_scanner);
//
//   // Loop over conditions.
//   for (int cond = 0; cond < num_conditions; cond++) {
//     List cond_cache = event_cache[cond];
//     // Retrieve the cached FFT.
//     arma::cx_vec cached_fft = as<arma::cx_vec>(cond_cache["regressor_fft"]);
//
//     // Recompute hr_regressor and hr_frame_times using the cached exp_condition.
//     NumericMatrix exp_condition = as<NumericMatrix>(cond_cache["exp_condition"]);
//     List sample_out = sample_event_condition(exp_condition, scanner_frame_times, oversampling, min_onset);
//     NumericVector hr_regressor = sample_out["regressor"];
//     NumericVector hr_frame_times = sample_out["hr_frame_times"];
//     int M = hr_regressor.size();
//
//     // Convolution: use cached FFT
//     arma::cx_vec conv_fft = cached_fft % kernel_fft;
//     arma::cx_vec conv_ifft = arma::ifft(conv_fft);
//     arma::vec conv_full = arma::real(conv_ifft.head(M));
//
//     // Downsample using resample_vector.
//     NumericVector ds_numeric = resample_vector(hr_frame_times, wrap(conv_full), scanner_frame_times);
//
//     // Prediction; sum of X_i * B_i
//     pred += beta[cond] * ds_numeric;
//   }
//
//   // LL
//   double logLik = 0.0;
//   for (int i = 0; i < n_scanner; i++) {
//     logLik += R::dnorm(y[i], pred[i], sigma, true);
//   }
//
//   return logLik;
// }
