#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

namespace {

// Interpolate first hit on a circle along a straight segment p -> q
inline std::pair<double, double> interpolate_circle_hit_cpp(double px, double py,
                                                            double qx, double qy,
                                                            double a) {
  const double dx = qx - px;
  const double dy = qy - py;
  const double A = dx * dx + dy * dy;
  const double B = 2.0 * (px * dx + py * dy);
  const double C = (px * px + py * py) - (a * a);
  const double disc = B * B - 4.0 * A * C;

  bool ok = (A > 0.0) && R_finite(disc) && (disc >= 0.0);
  double s = NA_REAL;
  if (ok) {
    s = (-B + std::sqrt(std::max(disc, 0.0))) / (2.0 * A);
  }

  const double r0 = std::sqrt(px * px + py * py);
  const double r1 = std::sqrt(qx * qx + qy * qy);
  const double denom = (r1 - r0);
  const double s_rad = (denom != 0.0) ? (a - r0) / denom : 1.0;

  const bool bad = (!ok) || !R_finite(s) || (s < 0.0) || (s > 1.0);
  if (bad) s = s_rad;
  if (s < 0.0) s = 0.0;
  if (s > 1.0) s = 1.0;

  const double xi = px + s * dx;
  const double yi = py + s * dy;
  const double angle = std::atan2(yi, xi);
  return std::make_pair(s, angle);
}

inline int column_index(const StringVector& cn, const std::string& name) {
  for (int j = 0; j < cn.size(); ++j) {
    if (name == Rcpp::as<std::string>(cn[j])) return j;
  }
  return -1;
}

}

// [[Rcpp::export]]
DataFrame rCDM(NumericMatrix pars,
               Nullable<LogicalVector> ok = R_NilValue,
               double dt = 1e-5,
               int max_steps = 100000000) {
  RNGScope scope;

  const int n = pars.nrow();
  if (n == 0) {
    return DataFrame::create(Named("rt") = NumericVector(0),
                             Named("R")  = NumericVector(0));
  }

  LogicalVector ok_vec;
  if (ok.isNull()) {
    ok_vec = LogicalVector(n, true);
  } else {
    ok_vec = LogicalVector(ok);
    if (ok_vec.size() != n) stop("Length of 'ok' must match nrow(pars).");
  }

  StringVector cn = colnames(pars);
  if (cn.size() == 0) stop("pars must have column names.");

  const int idx_v     = column_index(cn, "v");
  const int idx_th    = column_index(cn, "theta");
  const int idx_a     = column_index(cn, "a");
  const int idx_t0    = column_index(cn, "t0");
  const int idx_sigma = column_index(cn, "sigma");
  const int idx_sv    = column_index(cn, "sv");
  if (idx_v < 0 || idx_th < 0 || idx_a < 0 || idx_t0 < 0 || idx_sigma < 0) {
    stop("pars must contain columns v, theta, a, t0, sigma.");
  }

  NumericVector v     = pars(_, idx_v);
  NumericVector theta = pars(_, idx_th);
  NumericVector a     = pars(_, idx_a);
  NumericVector t0    = pars(_, idx_t0);
  NumericVector sigma = pars(_, idx_sigma);
  NumericVector sv    = (idx_sv >= 0) ? pars(_, idx_sv) : NumericVector(n, 0.0);

  // map theta from (0,1) to [-pi, pi]
  const double two_pi = 6.28318530717958647692; // 2*pi
  for (int i = 0; i < n; ++i) {
    theta[i] = (theta[i] - 0.5) * two_pi;
  }

  std::vector<int> idx;
  idx.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (ok_vec[i]) idx.push_back(i);
  }
  NumericVector rt_out(n, NA_REAL);
  NumericVector R_out(n, NA_REAL);
  if (idx.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out);
  }

  std::vector<int> sel;
  sel.reserve(idx.size());
  for (int id : idx) {
    const double ai = a[id];
    const double t0i = t0[id];
    const double sigi = sigma[id];
    const bool valid = R_finite(ai) && ai > 0.0 &&
                       R_finite(t0i) && t0i >= 0.0 &&
                       R_finite(sigi) && sigi > 0.0;
    if (valid) sel.push_back(id);
  }
  if (sel.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out);
  }

  const int n_act = sel.size();
  std::vector<double> ax(n_act), t0_i(n_act), sig_i(n_act), sv_i(n_act);
  std::vector<double> mu_x(n_act), mu_y(n_act);
  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    ax[i]    = a[id];
    t0_i[i]  = t0[id];
    sig_i[i] = sigma[id];
    sv_i[i]  = sv[id];
    mu_x[i]  = v[id] * std::cos(theta[id]);
    mu_y[i]  = v[id] * std::sin(theta[id]);
  }

  // Between-trial drift variability
  {
    int k = 0;
    for (double s : sv_i) if (R_finite(s) && s > 0.0) ++k;
    if (k > 0) {
      NumericVector noise = rnorm(k * 2);
      int row = 0;
      for (int i = 0; i < n_act; ++i) {
        const double s = sv_i[i];
        if (R_finite(s) && s > 0.0) {
          mu_x[i] += noise[row] * s;
          mu_y[i] += noise[row + k] * s;
          ++row;
        }
      }
    }
  }

  // Initialize state
  std::vector<double> x0(n_act, 0.0), x1(n_act, 0.0), rt_vec(n_act, 0.0);
  std::vector<char> done(n_act, 0);
  std::vector<double> R_sel(n_act, NA_REAL), rt_sel(n_act, NA_REAL);

  std::vector<double> drift_x(n_act), drift_y(n_act), sqrt_dt_sigma(n_act), a_sq(n_act);
  for (int i = 0; i < n_act; ++i) {
    drift_x[i] = mu_x[i] * dt;
    drift_y[i] = mu_y[i] * dt;
    sqrt_dt_sigma[i] = std::sqrt(dt) * sig_i[i];
    a_sq[i] = ax[i] * ax[i];
  }

  int steps = 0;
  while (true) {
    bool all_done = true;
    for (int i = 0; i < n_act; ++i) if (!done[i]) { all_done = false; break; }
    if (all_done) break;
    ++steps;
    if (steps >= max_steps) stop("Maximum number of steps reached; consider increasing 'a' or 'dt'.");

    int k_active = 0;
    for (int i = 0; i < n_act; ++i) if (!done[i]) ++k_active;
    if (k_active == 0) continue;

    NumericVector noise = rnorm(k_active * 2);
    int row = 0;
    for (int i = 0; i < n_act; ++i) {
      if (done[i]) continue;
      const double px = x0[i];
      const double py = x1[i];
      const double sig_dt = sqrt_dt_sigma[i];
      const double nx = noise[row];
      const double ny = noise[row + k_active];
      const double qx = px + drift_x[i] + nx * sig_dt;
      const double qy = py + drift_y[i] + ny * sig_dt;

      x0[i] = qx;
      x1[i] = qy;
      rt_vec[i] += dt;

      const double r_sq = qx * qx + qy * qy;
      if (r_sq >= a_sq[i]) {
        const std::pair<double, double> hit = interpolate_circle_hit_cpp(px, py, qx, qy, ax[i]);
        R_sel[i]  = hit.second;
        rt_sel[i] = t0_i[i] + (rt_vec[i] - (1.0 - hit.first) * dt);
        done[i]   = 1;
      }
      ++row;
    }
  }

  // Map back to full length
  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    R_out[id]  = R_sel[i];
    rt_out[id] = rt_sel[i];
  }

  return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out);
}
