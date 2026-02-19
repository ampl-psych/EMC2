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

inline std::pair<double, std::vector<double>> interpolate_sphere_hit_cpp(const std::vector<double>& p,
                                                                          const std::vector<double>& q,
                                                                          double a) {
  const int d = p.size();
  std::vector<double> delta(d), x(d);
  double A = 0.0, B = 0.0, p2 = 0.0, q2 = 0.0;
  for (int k = 0; k < d; ++k) {
    delta[k] = q[k] - p[k];
    A += delta[k] * delta[k];
    B += p[k] * delta[k];
    p2 += p[k] * p[k];
    q2 += q[k] * q[k];
  }
  B *= 2.0;
  const double C = p2 - a * a;
  const double disc = B * B - 4.0 * A * C;

  bool ok = (A > 0.0) && R_finite(disc) && (disc >= 0.0);
  double s = NA_REAL;
  if (ok) {
    s = (-B + std::sqrt(std::max(disc, 0.0))) / (2.0 * A);
  }

  const double r0 = std::sqrt(p2);
  const double r1 = std::sqrt(q2);
  const double denom = (r1 - r0);
  const double s_rad = (denom != 0.0) ? (a - r0) / denom : 1.0;
  const bool bad = (!ok) || !R_finite(s) || (s < 0.0) || (s > 1.0);
  if (bad) s = s_rad;
  if (s < 0.0) s = 0.0;
  if (s > 1.0) s = 1.0;

  for (int k = 0; k < d; ++k) {
    x[k] = p[k] + s * delta[k];
  }
  return std::make_pair(s, x);
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
  const int idx_sigma = column_index(cn, "s");
  const int idx_sv    = column_index(cn, "sv");
  if (idx_v < 0 || idx_th < 0 || idx_a < 0 || idx_t0 < 0 || idx_sigma < 0) {
    stop("pars must contain columns v, theta, a, t0, s.");
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

// [[Rcpp::export]]
DataFrame rSDM(NumericMatrix pars,
               Nullable<LogicalVector> ok = R_NilValue,
               double dt = 1e-5,
               int max_steps = 100000000) {
  RNGScope scope;

  const int n = pars.nrow();
  if (n == 0) {
    return DataFrame::create(Named("rt") = NumericVector(0),
                             Named("R")  = NumericVector(0),
                             Named("R2") = NumericVector(0));
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

  const int idx_v      = column_index(cn, "v");
  const int idx_th1    = column_index(cn, "theta1");
  const int idx_th2    = column_index(cn, "theta2");
  const int idx_a      = column_index(cn, "a");
  const int idx_t0     = column_index(cn, "t0");
  const int idx_sigma  = column_index(cn, "s");
  const int idx_sv     = column_index(cn, "sv");
  if (idx_v < 0 || idx_th1 < 0 || idx_th2 < 0 || idx_a < 0 || idx_t0 < 0 || idx_sigma < 0) {
    stop("pars must contain columns v, theta1, theta2, a, t0, s.");
  }

  NumericVector v      = pars(_, idx_v);
  NumericVector theta1 = pars(_, idx_th1);
  NumericVector theta2 = pars(_, idx_th2);
  NumericVector a      = pars(_, idx_a);
  NumericVector t0     = pars(_, idx_t0);
  NumericVector sigma  = pars(_, idx_sigma);
  NumericVector sv     = (idx_sv >= 0) ? pars(_, idx_sv) : NumericVector(n, 0.0);

  const double pi = 3.14159265358979323846264338327950288;
  const double two_pi = 2.0 * pi;
  for (int i = 0; i < n; ++i) {
    theta1[i] = theta1[i] * pi;
    theta2[i] = (theta2[i] - 0.5) * two_pi;
  }

  std::vector<int> idx;
  idx.reserve(n);
  for (int i = 0; i < n; ++i) if (ok_vec[i]) idx.push_back(i);

  NumericVector rt_out(n, NA_REAL), R_out(n, NA_REAL), R2_out(n, NA_REAL);
  if (idx.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out);
  }

  std::vector<int> sel;
  sel.reserve(idx.size());
  for (int id : idx) {
    const bool valid = R_finite(a[id]) && a[id] > 0.0 &&
                       R_finite(t0[id]) && t0[id] >= 0.0 &&
                       R_finite(sigma[id]) && sigma[id] > 0.0;
    if (valid) sel.push_back(id);
  }
  if (sel.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out);
  }

  const int n_act = sel.size();
  std::vector<double> ax(n_act), t0_i(n_act), sig_i(n_act), sv_i(n_act);
  std::vector<double> mu0(n_act), mu1(n_act), mu2(n_act);
  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    const double st1 = std::sin(theta1[id]);
    ax[i] = a[id];
    t0_i[i] = t0[id];
    sig_i[i] = sigma[id];
    sv_i[i] = sv[id];
    mu0[i] = v[id] * std::cos(theta1[id]);
    mu1[i] = v[id] * st1 * std::cos(theta2[id]);
    mu2[i] = v[id] * st1 * std::sin(theta2[id]);
  }

  for (int i = 0; i < n_act; ++i) {
    const double s = sv_i[i];
    if (R_finite(s) && s > 0.0) {
      mu0[i] += R::rnorm(0.0, s);
      mu1[i] += R::rnorm(0.0, s);
      mu2[i] += R::rnorm(0.0, s);
    }
  }

  std::vector<double> x0(n_act, 0.0), x1(n_act, 0.0), x2(n_act, 0.0), rt_vec(n_act, 0.0);
  std::vector<char> done(n_act, 0);
  std::vector<double> R_sel(n_act, NA_REAL), R2_sel(n_act, NA_REAL), rt_sel(n_act, NA_REAL);
  std::vector<double> drift0(n_act), drift1(n_act), drift2(n_act), sqrt_dt_sigma(n_act), a_sq(n_act);
  for (int i = 0; i < n_act; ++i) {
    drift0[i] = mu0[i] * dt;
    drift1[i] = mu1[i] * dt;
    drift2[i] = mu2[i] * dt;
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

    NumericVector noise = rnorm(k_active * 3);
    int row = 0;
    for (int i = 0; i < n_act; ++i) {
      if (done[i]) continue;
      const double px0 = x0[i], px1 = x1[i], px2 = x2[i];
      const double sig_dt = sqrt_dt_sigma[i];
      const double q0 = px0 + drift0[i] + noise[row] * sig_dt;
      const double q1 = px1 + drift1[i] + noise[row + k_active] * sig_dt;
      const double q2 = px2 + drift2[i] + noise[row + 2 * k_active] * sig_dt;

      x0[i] = q0; x1[i] = q1; x2[i] = q2;
      rt_vec[i] += dt;

      const double r_sq = q0 * q0 + q1 * q1 + q2 * q2;
      if (r_sq >= a_sq[i]) {
        std::vector<double> p = {px0, px1, px2};
        std::vector<double> q = {q0, q1, q2};
        const std::pair<double, std::vector<double>> hit = interpolate_sphere_hit_cpp(p, q, ax[i]);
        const double hx0 = hit.second[0], hx1 = hit.second[1], hx2 = hit.second[2];
        R_sel[i] = std::atan2(std::sqrt(hx1 * hx1 + hx2 * hx2), hx0);
        R2_sel[i] = std::atan2(hx2, hx1);
        rt_sel[i] = t0_i[i] + (rt_vec[i] - (1.0 - hit.first) * dt);
        done[i] = 1;
      }
      ++row;
    }
  }

  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    R_out[id] = R_sel[i];
    R2_out[id] = R2_sel[i];
    rt_out[id] = rt_sel[i];
  }

  return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out);
}

// [[Rcpp::export]]
DataFrame rPSDM(NumericMatrix pars,
                Nullable<LogicalVector> ok = R_NilValue,
                double dt = 1e-5,
                int max_steps = 100000000) {
  RNGScope scope;

  const int n = pars.nrow();
  if (n == 0) {
    return DataFrame::create(Named("rt") = NumericVector(0),
                             Named("R") = NumericVector(0));
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

  const int idx_v = column_index(cn, "v");
  const int idx_th1 = column_index(cn, "theta1");
  const int idx_a = column_index(cn, "a");
  const int idx_t0 = column_index(cn, "t0");
  const int idx_sigma = column_index(cn, "s");
  const int idx_sv = column_index(cn, "sv");
  if (idx_v < 0 || idx_th1 < 0 || idx_a < 0 || idx_t0 < 0 || idx_sigma < 0) {
    stop("pars must contain columns v, theta1, a, t0, s.");
  }

  NumericVector v = pars(_, idx_v);
  NumericVector theta1 = pars(_, idx_th1);
  NumericVector a = pars(_, idx_a);
  NumericVector t0 = pars(_, idx_t0);
  NumericVector sigma = pars(_, idx_sigma);
  NumericVector sv = (idx_sv >= 0) ? pars(_, idx_sv) : NumericVector(n, 0.0);

  const double pi = 3.14159265358979323846264338327950288;
  const double phi = pi / 4.0;
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  for (int i = 0; i < n; ++i) {
    theta1[i] = theta1[i] * pi;
  }

  std::vector<int> idx;
  idx.reserve(n);
  for (int i = 0; i < n; ++i) if (ok_vec[i]) idx.push_back(i);

  NumericVector rt_out(n, NA_REAL), R_out(n, NA_REAL);
  if (idx.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out);
  }

  std::vector<int> sel;
  sel.reserve(idx.size());
  for (int id : idx) {
    const bool valid = R_finite(a[id]) && a[id] > 0.0 &&
      R_finite(t0[id]) && t0[id] >= 0.0 &&
      R_finite(sigma[id]) && sigma[id] > 0.0;
    if (valid) sel.push_back(id);
  }
  if (sel.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out);
  }

  const int n_act = sel.size();
  std::vector<double> ax(n_act), t0_i(n_act), sig_i(n_act), sv_i(n_act);
  std::vector<double> mu0(n_act), mu1(n_act), mu2(n_act);
  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    const double st1 = std::sin(theta1[id]);
    ax[i] = a[id];
    t0_i[i] = t0[id];
    sig_i[i] = sigma[id];
    sv_i[i] = sv[id];
    mu0[i] = v[id] * std::cos(theta1[id]);
    mu1[i] = v[id] * st1 * cphi;
    mu2[i] = v[id] * st1 * sphi;
  }

  for (int i = 0; i < n_act; ++i) {
    const double s = sv_i[i];
    if (R_finite(s) && s > 0.0) {
      mu0[i] += R::rnorm(0.0, s);
      mu1[i] += R::rnorm(0.0, s);
      mu2[i] += R::rnorm(0.0, s);
    }
  }

  std::vector<double> x0(n_act, 0.0), x1(n_act, 0.0), x2(n_act, 0.0), rt_vec(n_act, 0.0);
  std::vector<char> done(n_act, 0);
  std::vector<double> R_sel(n_act, NA_REAL), rt_sel(n_act, NA_REAL);
  std::vector<double> drift0(n_act), drift1(n_act), drift2(n_act), sqrt_dt_sigma(n_act), a_sq(n_act);
  for (int i = 0; i < n_act; ++i) {
    drift0[i] = mu0[i] * dt;
    drift1[i] = mu1[i] * dt;
    drift2[i] = mu2[i] * dt;
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

    NumericVector noise = rnorm(k_active * 3);
    int row = 0;
    for (int i = 0; i < n_act; ++i) {
      if (done[i]) continue;
      const double px0 = x0[i], px1 = x1[i], px2 = x2[i];
      const double sig_dt = sqrt_dt_sigma[i];
      const double q0 = px0 + drift0[i] + noise[row] * sig_dt;
      const double q1 = px1 + drift1[i] + noise[row + k_active] * sig_dt;
      const double q2 = px2 + drift2[i] + noise[row + 2 * k_active] * sig_dt;

      x0[i] = q0; x1[i] = q1; x2[i] = q2;
      rt_vec[i] += dt;

      const double r_sq = q0 * q0 + q1 * q1 + q2 * q2;
      if (r_sq >= a_sq[i]) {
        std::vector<double> p = {px0, px1, px2};
        std::vector<double> q = {q0, q1, q2};
        const std::pair<double, std::vector<double>> hit = interpolate_sphere_hit_cpp(p, q, ax[i]);
        const double hx0 = hit.second[0], hx1 = hit.second[1], hx2 = hit.second[2];
        R_sel[i] = std::atan2(std::sqrt(hx1 * hx1 + hx2 * hx2), hx0);
        rt_sel[i] = t0_i[i] + (rt_vec[i] - (1.0 - hit.first) * dt);
        done[i] = 1;
      }
      ++row;
    }
  }

  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    R_out[id] = R_sel[i];
    rt_out[id] = rt_sel[i];
  }

  return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out);
}

// [[Rcpp::export]]
DataFrame rPHSDM(NumericMatrix pars,
                 Nullable<LogicalVector> ok = R_NilValue,
                 double dt = 1e-5,
                 int max_steps = 100000000) {
  RNGScope scope;

  const int n = pars.nrow();
  if (n == 0) {
    return DataFrame::create(Named("rt") = NumericVector(0),
                             Named("R") = NumericVector(0),
                             Named("R2") = NumericVector(0));
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

  const int idx_v = column_index(cn, "v");
  const int idx_th1 = column_index(cn, "theta1");
  const int idx_th2 = column_index(cn, "theta2");
  const int idx_a = column_index(cn, "a");
  const int idx_t0 = column_index(cn, "t0");
  const int idx_sigma = column_index(cn, "s");
  const int idx_sv = column_index(cn, "sv");
  if (idx_v < 0 || idx_th1 < 0 || idx_th2 < 0 || idx_a < 0 || idx_t0 < 0 || idx_sigma < 0) {
    stop("pars must contain columns v, theta1, theta2, a, t0, s.");
  }

  NumericVector v = pars(_, idx_v);
  NumericVector theta1 = pars(_, idx_th1);
  NumericVector theta2 = pars(_, idx_th2);
  NumericVector a = pars(_, idx_a);
  NumericVector t0 = pars(_, idx_t0);
  NumericVector sigma = pars(_, idx_sigma);
  NumericVector sv = (idx_sv >= 0) ? pars(_, idx_sv) : NumericVector(n, 0.0);

  const double pi = 3.14159265358979323846264338327950288;
  const double phi = pi / 4.0;
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  for (int i = 0; i < n; ++i) {
    theta1[i] = theta1[i] * pi;
    theta2[i] = theta2[i] * pi;
  }

  std::vector<int> idx;
  idx.reserve(n);
  for (int i = 0; i < n; ++i) if (ok_vec[i]) idx.push_back(i);

  NumericVector rt_out(n, NA_REAL), R_out(n, NA_REAL), R2_out(n, NA_REAL);
  if (idx.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out);
  }

  std::vector<int> sel;
  sel.reserve(idx.size());
  for (int id : idx) {
    const bool valid = R_finite(a[id]) && a[id] > 0.0 &&
      R_finite(t0[id]) && t0[id] >= 0.0 &&
      R_finite(sigma[id]) && sigma[id] > 0.0;
    if (valid) sel.push_back(id);
  }
  if (sel.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out);
  }

  const int n_act = sel.size();
  std::vector<double> ax(n_act), t0_i(n_act), sig_i(n_act), sv_i(n_act);
  std::vector<double> mu0(n_act), mu1(n_act), mu2(n_act), mu3(n_act);
  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    const double st1 = std::sin(theta1[id]);
    const double st2 = std::sin(theta2[id]);
    ax[i] = a[id];
    t0_i[i] = t0[id];
    sig_i[i] = sigma[id];
    sv_i[i] = sv[id];
    mu0[i] = v[id] * std::cos(theta1[id]);
    mu1[i] = v[id] * st1 * std::cos(theta2[id]);
    mu2[i] = v[id] * st1 * st2 * cphi;
    mu3[i] = v[id] * st1 * st2 * sphi;
  }

  for (int i = 0; i < n_act; ++i) {
    const double s = sv_i[i];
    if (R_finite(s) && s > 0.0) {
      mu0[i] += R::rnorm(0.0, s);
      mu1[i] += R::rnorm(0.0, s);
      mu2[i] += R::rnorm(0.0, s);
      mu3[i] += R::rnorm(0.0, s);
    }
  }

  std::vector<double> x0(n_act, 0.0), x1(n_act, 0.0), x2(n_act, 0.0), x3(n_act, 0.0), rt_vec(n_act, 0.0);
  std::vector<char> done(n_act, 0);
  std::vector<double> R_sel(n_act, NA_REAL), R2_sel(n_act, NA_REAL), rt_sel(n_act, NA_REAL);
  std::vector<double> drift0(n_act), drift1(n_act), drift2(n_act), drift3(n_act), sqrt_dt_sigma(n_act), a_sq(n_act);
  for (int i = 0; i < n_act; ++i) {
    drift0[i] = mu0[i] * dt;
    drift1[i] = mu1[i] * dt;
    drift2[i] = mu2[i] * dt;
    drift3[i] = mu3[i] * dt;
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

    NumericVector noise = rnorm(k_active * 4);
    int row = 0;
    for (int i = 0; i < n_act; ++i) {
      if (done[i]) continue;
      const double px0 = x0[i], px1 = x1[i], px2 = x2[i], px3 = x3[i];
      const double sig_dt = sqrt_dt_sigma[i];
      const double q0 = px0 + drift0[i] + noise[row] * sig_dt;
      const double q1 = px1 + drift1[i] + noise[row + k_active] * sig_dt;
      const double q2 = px2 + drift2[i] + noise[row + 2 * k_active] * sig_dt;
      const double q3 = px3 + drift3[i] + noise[row + 3 * k_active] * sig_dt;

      x0[i] = q0; x1[i] = q1; x2[i] = q2; x3[i] = q3;
      rt_vec[i] += dt;

      const double r_sq = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
      if (r_sq >= a_sq[i]) {
        std::vector<double> p = {px0, px1, px2, px3};
        std::vector<double> q = {q0, q1, q2, q3};
        const std::pair<double, std::vector<double>> hit = interpolate_sphere_hit_cpp(p, q, ax[i]);
        const double hx0 = hit.second[0], hx1 = hit.second[1], hx2 = hit.second[2], hx3 = hit.second[3];
        R_sel[i] = std::atan2(std::sqrt(hx1 * hx1 + hx2 * hx2 + hx3 * hx3), hx0);
        R2_sel[i] = std::atan2(std::sqrt(hx2 * hx2 + hx3 * hx3), hx1);
        rt_sel[i] = t0_i[i] + (rt_vec[i] - (1.0 - hit.first) * dt);
        done[i] = 1;
      }
      ++row;
    }
  }

  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    R_out[id] = R_sel[i];
    R2_out[id] = R2_sel[i];
    rt_out[id] = rt_sel[i];
  }

  return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out);
}

// [[Rcpp::export]]
DataFrame rHSDM(NumericMatrix pars,
                Nullable<LogicalVector> ok = R_NilValue,
                double dt = 1e-5,
                int max_steps = 100000000) {
  RNGScope scope;

  const int n = pars.nrow();
  if (n == 0) {
    return DataFrame::create(Named("rt") = NumericVector(0),
                             Named("R")  = NumericVector(0),
                             Named("R2") = NumericVector(0),
                             Named("R3") = NumericVector(0));
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

  const int idx_v      = column_index(cn, "v");
  const int idx_th1    = column_index(cn, "theta1");
  const int idx_th2    = column_index(cn, "theta2");
  const int idx_th3    = column_index(cn, "theta3");
  const int idx_a      = column_index(cn, "a");
  const int idx_t0     = column_index(cn, "t0");
  const int idx_sigma  = column_index(cn, "s");
  const int idx_sv     = column_index(cn, "sv");
  if (idx_v < 0 || idx_th1 < 0 || idx_th2 < 0 || idx_th3 < 0 || idx_a < 0 || idx_t0 < 0 || idx_sigma < 0) {
    stop("pars must contain columns v, theta1, theta2, theta3, a, t0, s.");
  }

  NumericVector v      = pars(_, idx_v);
  NumericVector theta1 = pars(_, idx_th1);
  NumericVector theta2 = pars(_, idx_th2);
  NumericVector theta3 = pars(_, idx_th3);
  NumericVector a      = pars(_, idx_a);
  NumericVector t0     = pars(_, idx_t0);
  NumericVector sigma  = pars(_, idx_sigma);
  NumericVector sv     = (idx_sv >= 0) ? pars(_, idx_sv) : NumericVector(n, 0.0);

  const double pi = 3.14159265358979323846264338327950288;
  const double two_pi = 2.0 * pi;
  for (int i = 0; i < n; ++i) {
    theta1[i] = theta1[i] * pi;
    theta2[i] = theta2[i] * pi;
    theta3[i] = (theta3[i] - 0.5) * two_pi;
  }

  std::vector<int> idx;
  idx.reserve(n);
  for (int i = 0; i < n; ++i) if (ok_vec[i]) idx.push_back(i);

  NumericVector rt_out(n, NA_REAL), R_out(n, NA_REAL), R2_out(n, NA_REAL), R3_out(n, NA_REAL);
  if (idx.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out, Named("R3") = R3_out);
  }

  std::vector<int> sel;
  sel.reserve(idx.size());
  for (int id : idx) {
    const bool valid = R_finite(a[id]) && a[id] > 0.0 &&
                       R_finite(t0[id]) && t0[id] >= 0.0 &&
                       R_finite(sigma[id]) && sigma[id] > 0.0;
    if (valid) sel.push_back(id);
  }
  if (sel.empty()) {
    return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out, Named("R3") = R3_out);
  }

  const int n_act = sel.size();
  std::vector<double> ax(n_act), t0_i(n_act), sig_i(n_act), sv_i(n_act);
  std::vector<double> mu0(n_act), mu1(n_act), mu2(n_act), mu3(n_act);
  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    const double st1 = std::sin(theta1[id]);
    const double st2 = std::sin(theta2[id]);
    ax[i] = a[id];
    t0_i[i] = t0[id];
    sig_i[i] = sigma[id];
    sv_i[i] = sv[id];
    mu0[i] = v[id] * std::cos(theta1[id]);
    mu1[i] = v[id] * st1 * std::cos(theta2[id]);
    mu2[i] = v[id] * st1 * st2 * std::cos(theta3[id]);
    mu3[i] = v[id] * st1 * st2 * std::sin(theta3[id]);
  }

  for (int i = 0; i < n_act; ++i) {
    const double s = sv_i[i];
    if (R_finite(s) && s > 0.0) {
      mu0[i] += R::rnorm(0.0, s);
      mu1[i] += R::rnorm(0.0, s);
      mu2[i] += R::rnorm(0.0, s);
      mu3[i] += R::rnorm(0.0, s);
    }
  }

  std::vector<double> x0(n_act, 0.0), x1(n_act, 0.0), x2(n_act, 0.0), x3(n_act, 0.0), rt_vec(n_act, 0.0);
  std::vector<char> done(n_act, 0);
  std::vector<double> R_sel(n_act, NA_REAL), R2_sel(n_act, NA_REAL), R3_sel(n_act, NA_REAL), rt_sel(n_act, NA_REAL);
  std::vector<double> drift0(n_act), drift1(n_act), drift2(n_act), drift3(n_act), sqrt_dt_sigma(n_act), a_sq(n_act);
  for (int i = 0; i < n_act; ++i) {
    drift0[i] = mu0[i] * dt;
    drift1[i] = mu1[i] * dt;
    drift2[i] = mu2[i] * dt;
    drift3[i] = mu3[i] * dt;
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

    NumericVector noise = rnorm(k_active * 4);
    int row = 0;
    for (int i = 0; i < n_act; ++i) {
      if (done[i]) continue;
      const double px0 = x0[i], px1 = x1[i], px2 = x2[i], px3 = x3[i];
      const double sig_dt = sqrt_dt_sigma[i];
      const double q0 = px0 + drift0[i] + noise[row] * sig_dt;
      const double q1 = px1 + drift1[i] + noise[row + k_active] * sig_dt;
      const double q2 = px2 + drift2[i] + noise[row + 2 * k_active] * sig_dt;
      const double q3 = px3 + drift3[i] + noise[row + 3 * k_active] * sig_dt;

      x0[i] = q0; x1[i] = q1; x2[i] = q2; x3[i] = q3;
      rt_vec[i] += dt;

      const double r_sq = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
      if (r_sq >= a_sq[i]) {
        std::vector<double> p = {px0, px1, px2, px3};
        std::vector<double> q = {q0, q1, q2, q3};
        const std::pair<double, std::vector<double>> hit = interpolate_sphere_hit_cpp(p, q, ax[i]);
        const double hx0 = hit.second[0], hx1 = hit.second[1], hx2 = hit.second[2], hx3 = hit.second[3];
        R_sel[i] = std::atan2(std::sqrt(hx1 * hx1 + hx2 * hx2 + hx3 * hx3), hx0);
        R2_sel[i] = std::atan2(std::sqrt(hx2 * hx2 + hx3 * hx3), hx1);
        R3_sel[i] = std::atan2(hx3, hx2);
        rt_sel[i] = t0_i[i] + (rt_vec[i] - (1.0 - hit.first) * dt);
        done[i] = 1;
      }
      ++row;
    }
  }

  for (int i = 0; i < n_act; ++i) {
    const int id = sel[i];
    R_out[id] = R_sel[i];
    R2_out[id] = R2_sel[i];
    R3_out[id] = R3_sel[i];
    rt_out[id] = rt_sel[i];
  }

  return DataFrame::create(Named("rt") = rt_out, Named("R") = R_out, Named("R2") = R2_out, Named("R3") = R3_out);
}
