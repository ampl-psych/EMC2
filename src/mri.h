#ifndef mri_dm
#define mri_dm

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <set>
#include <vector>
#include <string>
using namespace Rcpp;
using namespace arma;

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
  double dt = tr / oversampling;
  int n_points = std::max(1, (int) std::round(time_length / dt));
  NumericVector time_stamps = seq_lin(0, time_length, n_points);
  for (int i = 0; i < time_stamps.size(); i++) {
    time_stamps[i] = time_stamps[i] - onset;
  }
  int n = time_stamps.size();
  NumericVector peak_gamma(n), undershoot_gamma(n), hrf(n);
  for (int i = 0; i < n; i++) {
    double t_val = (time_stamps[i] - dt) / dispersion;
    peak_gamma[i] = R::dgamma(t_val, delay / dispersion, 1.0, 0);
    double t_val2 = (time_stamps[i] - dt) / u_dispersion;
    undershoot_gamma[i] = R::dgamma(t_val2, undershoot / u_dispersion, 1.0, 0);
    hrf[i] = peak_gamma[i] - ratio * undershoot_gamma[i];
  }
  // Normalize HRF by its maximum value.
  double max_val = hrf[0];
  for (int i = 1; i < n; i++) {
    if (hrf[i] > max_val)
      max_val = hrf[i];
  }
  for (int i = 0; i < n; i++) {
    hrf[i] /= max_val;
  }
  return hrf;
}

// ----------------------------------------------
// Updated: compute_glover_hrf
// Simply calls compute_gamma_diff_hrf with all parameters.
// ----------------------------------------------
// [[Rcpp::export]]
NumericVector compute_glover_hrf(double tr, int oversampling, double time_length, double onset,
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
NumericVector compute_glover_time_derivative(double tr, int oversampling, double time_length,
                                             double onset, double delay, double undershoot,
                                             double dispersion, double u_dispersion, double ratio,
                                             double delta = 0.1) {
  NumericVector hrf1 = compute_glover_hrf(tr, oversampling, time_length, onset,
                                          delay, undershoot, dispersion, u_dispersion, ratio);
  NumericVector hrf2 = compute_glover_hrf(tr, oversampling, time_length, onset + delta,
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
NumericMatrix build_hrf_kernel(std::string hrf_model, double tr, int oversampling,
                               double time_length, double onset, double delay,
                               double undershoot, double dispersion, double u_dispersion,
                               double ratio) {
  std::string model = hrf_model;
  std::transform(model.begin(), model.end(), model.begin(), ::tolower);
  if (model == "glover" || model == "glover_hrf") {
    NumericVector hrf = compute_glover_hrf(tr, oversampling, time_length, onset,
                                           delay, undershoot, dispersion, u_dispersion, ratio);
    NumericMatrix kernel(hrf.size(), 1);
    for (int i = 0; i < hrf.size(); i++) {
      kernel(i, 0) = hrf[i];
    }
    return kernel;
  } else if (model == "glover+derivative" || model == "glover + derivative") {
    NumericVector hrf = compute_glover_hrf(tr, oversampling, time_length, onset,
                                           delay, undershoot, dispersion, u_dispersion, ratio);
    NumericVector deriv = compute_glover_time_derivative(tr, oversampling, time_length, onset,
                                                         delay, undershoot, dispersion, u_dispersion, ratio);
    int n = hrf.size();
    NumericMatrix kernel(n, 2);
    for (int i = 0; i < n; i++) {
      kernel(i, 0) = hrf[i];
      kernel(i, 1) = deriv[i];
    }
    return kernel;
  } else {
    stop("Unsupported HRF model. Use 'glover' or 'glover + derivative'.");
  }
  return NumericMatrix(0);
}

// ----------------------------------------------
// Updated: compute_convolved_regressor
// Now accepts HRF parameters and passes them to build_hrf_kernel.
// ----------------------------------------------
// [[Rcpp::export]]
List compute_convolved_regressor(NumericMatrix exp_condition, std::string hrf_model,
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
  NumericMatrix hkernel = build_hrf_kernel(hrf_model, tr, oversampling, time_length, onset,
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
  std::string model = hrf_model;
  std::transform(model.begin(), model.end(), model.begin(), ::tolower);
  if (model == "glover" || model == "spm") {
    reg_names = CharacterVector::create(con_id);
  } else if (model == "glover+derivative" || model == "glover + derivative") {
    reg_names = CharacterVector::create(con_id, con_id + std::string("_derivative"));
  } else {
    reg_names = CharacterVector::create(con_id);
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
                                  std::string hrf_model, double min_onset, int oversampling,
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
    List out = compute_convolved_regressor(exp_condition, hrf_model, frame_times, cond, oversampling, min_onset,
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

#endif
