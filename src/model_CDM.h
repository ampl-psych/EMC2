#ifndef MODEL_CDM_H
#define MODEL_CDM_H

#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

namespace cdm_internal {

// Embedded Bessel J0 zeros and J1(zeros) values (same as R/model_CDM.R)
static const double ZEROS[] = {
  2.404825557695773, 5.520078110286311, 8.653727912911013, 11.791534439014281, 14.930917708487787,
  18.071063967910924, 21.21163662987926, 24.352471530749302, 27.493479132040253, 30.634606468431976,
  33.77582021357357, 36.917098353664045, 40.05842576462824, 43.19979171317673, 46.341188371661815,
  49.482609897397815, 52.624051841115, 55.76551075501998, 58.90698392608094, 62.04846919022717,
  65.18996480020687, 68.3314693298568, 71.47298160359374, 74.61450064370183, 77.75602563038805,
  80.89755587113763, 84.0390907769382, 87.18062984364116, 90.32217263721049, 93.46371878194478,
  96.60526795099626, 99.7468198586806, 102.8883742541948, 106.02993091645162, 109.17148964980538,
  112.3130502804949, 115.45461265366694, 118.59617663087253, 121.73774208795096, 124.87930891323295,
  128.02087700600833, 131.1624462752139, 134.30401663830546, 137.44558802028428, 140.58716035285428,
  143.72873357368974, 146.87030762579664, 150.01188245695477, 153.15345801922788, 156.29503426853353,
  159.43661116426316, 162.57818866894667, 165.71976674795502, 168.86134536923583, 172.0029245030782,
  175.14450412190274, 178.28608420007376, 181.42766471373105, 184.5692456406387, 187.71082696004936,
  190.85240865258152, 193.99399070010912, 197.1355730856614, 200.2771557933324, 203.41873880819864,
  206.56032211624446, 209.70190570429406, 212.8434895599495, 215.98507367153402, 219.12665802804057,
  222.2682426190843, 225.40982743485932, 228.5514124660988, 231.69299770403853, 234.83458314038324,
  237.97616876727565, 241.11775457726802, 244.2593405632957, 247.4009267186528, 250.54251303696995,
  253.6840995121931, 256.82568613856444, 259.9672729106045, 263.1088598230955, 266.2504468710659,
  269.39203404977604, 272.5336213547049, 275.67520878153744, 278.8167963261531, 281.9583839846149,
  285.09997175315954, 288.2415596281877, 291.3831476062552, 294.524735684065, 297.66632385845895,
  300.80791212641117, 303.9495004850206, 307.09108893150506, 310.232677463195, 313.37426607752786
};

static const double JVZ[] = {
  0.5191474972894667, -0.34026480655836816, 0.2714522999283819, -0.23245983136472478, 0.20654643307799606,
  -0.18772880304043937, 0.17326589422922986, -0.16170155068925005, 0.15218121377059451, -0.14416597768637318,
  0.13729694340850293, -0.13132462666866793, 0.1260694971272734, -0.12139862477175013, 0.11721119889066538,
  -0.1134291926164298, 0.10999114304627805, -0.10684788825471288, 0.10395957286936207, -0.10129349893394325,
  0.0988225538011999, -0.09652404046467991, 0.0943787939846764, -0.0923705048235533, 0.09048519416295771,
  -0.08871080244096977, 0.08703686332409762, -0.08545424291091487, 0.08395492928345759, -0.08253186130830982,
  0.08117878831953207, -0.07989015430874276, 0.07866100171930494, -0.07748689103965989, 0.07636383321829139,
  -0.07528823255205498, 0.0742568381822715, -0.07326670270620797, 0.0723151467023698, -0.071399728196232,
  0.07051821627333571, -0.06966856819003345, 0.06884890944684942, -0.068057516381685, 0.06729280091473994,
  -0.0665532971377112, 0.0658376494894271, -0.06514460230079303, 0.06447299052550669, -0.06382173150081485,
  0.06318981760571403, -0.06257630970331138, 0.06198033127024817, -0.06140106312970012, 0.060837738715961355,
  -0.06028963980834689, 0.05975609268041637, -0.05923646461756444, 0.05873016076204431, -0.058236621249651344,
  0.057755318606728306, -0.057285755379976114, 0.05682746197485655, -0.05637999468123289, 0.05594293386737833,
  -0.05551588232564261, 0.05509846375495183, -0.054690321366964154, 0.05429111660414672, -0.05390052795930567,
  0.05351824988721485, -0.053143991799969814, 0.05277747713856062, -0.05241844251392233, 0.052066636911401065,
  -0.051721820953175804, 0.05138376621371111, -0.051052254583793026, 0.05072707767912495, -0.05040803628983996,
  0.050094939867625976, -0.0497876060474634, 0.04948586020124814, -0.04918953502081812, 0.04889847012812101,
  -0.0486125117104598, 0.048331512178931856, -0.04805532984833819, 0.04778382863698613, -0.047516877784940494,
  0.04725435158939827, -0.04699612915597017, 0.04674209416475137, -0.04649213465015328, 0.046246142793549695,
  -0.04600401472786514, 0.045765650353301025, -0.045530953163457046, 0.045299830081161785, -0.04507219130337835
};

inline const std::vector<double>& weights() {
  static std::vector<double> w;
  if (w.empty()) {
    const size_t M = sizeof(ZEROS)/sizeof(ZEROS[0]);
    w.resize(M);
    for (size_t j = 0; j < M; ++j) w[j] = ZEROS[j] / JVZ[j];
  }
  return w;
}

inline const std::vector<double>& base_rates() {
  static std::vector<double> br;
  if (br.empty()) {
    const size_t M = sizeof(ZEROS)/sizeof(ZEROS[0]);
    br.resize(M);
    for (size_t j = 0; j < M; ++j) br[j] = 0.5 * ZEROS[j] * ZEROS[j];
  }
  return br;
}

inline double series_bessel_fpt_scalar(double t, double a, double sigma) {
  if (!(t > 0.0)) return 0.0;
  const std::vector<double>& w = weights();
  const std::vector<double>& br = base_rates();
  const double v_scale = (sigma * sigma) / (a * a);
  const double s = t * v_scale;
  double sum = 0.0;
  for (size_t j = 0; j < w.size(); ++j) {
    // Centered (e^{-λ s} - 1) to enforce f_T(0)=0 under truncation
    sum += w[j] * (std::exp(-br[j] * s) - 1.0);
  }
  return v_scale * sum;
}

inline double small_t_fpt_scalar(double t_scaled, double x_scaled) {
  // Uses first zero only via zeros[0]
  if (!(t_scaled > 0.0)) return 0.0;
  const double z1 = ZEROS[0];
  const double x = x_scaled;
  const double t = t_scaled;
  const double denom = std::sqrt(x + t) * std::pow(t, 1.5);
  if (!(denom > 0.0)) return 0.0;
  const double term1 = ((1.0 - x) * (1.0 + t)) / denom;
  const double term2 = std::exp(-0.5 * (1.0 - x) * (1.0 - x) / t - 0.5 * z1 * z1 * t);
  const double out = term1 * term2;
  if (!R_finite(out)) return 0.0;
  return out;
}

inline double percentile10(NumericVector v) {
  std::vector<double> vv;
  vv.reserve(v.size());
  for (int i = 0; i < v.size(); ++i) if (R_finite(v[i])) vv.push_back(v[i]);
  if (vv.empty()) return NA_REAL;
  std::sort(vv.begin(), vv.end());
  size_t idx = (size_t) std::floor(0.1 * (vv.size() - 1));
  return vv[idx];
}

} // namespace cdm_internal

// Returns log-density per trial (R_NegInf if invalid)
inline NumericVector c_dCDM(NumericVector rts, NumericVector Rs, NumericMatrix pars, LogicalVector is_ok) {
  using namespace cdm_internal;
  const int N = rts.size();
  NumericVector out(N);

  // Extract columns by index according to R model: v, theta, a, t0, sigma, sv
  NumericVector v   = pars(_, 0);
  NumericVector th  = pars(_, 1);
  NumericVector a   = pars(_, 2);
  NumericVector t0  = pars(_, 3);
  NumericVector sig = pars(_, 4);
  NumericVector sv  = pars(_, 5);

  // map theta from (0,1) to [-pi, pi]
  const double PI = 3.14159265358979323846264338327950288;
  for (int i = 0; i < N; ++i) {
    th[i] = (th[i] - 0.5) * 2.0 * PI;
  }

  // Compute tt and FPT series
  NumericVector tt(N);
  for (int i = 0; i < N; ++i) {
    double d = rts[i] - t0[i];
    tt[i] = (R_finite(d) && d > 0.0) ? d : 0.0;
  }

  // Base FPT
  NumericVector fpt(N);
  for (int i = 0; i < N; ++i) {
    fpt[i] = series_bessel_fpt_scalar(tt[i], a[i], sig[i]);
  }

  // Small-tt safeguard: replace lowest 10% tt by small-t approximation (scaled by v)
  double t_min = percentile10(tt);
  if (R_finite(t_min)) {
    for (int i = 0; i < N; ++i) {
      if (tt[i] < t_min) {
        const double ai = a[i];
        const double sigi = sig[i];
        const double vi = (sigi * sigi) / (ai * ai);
        const double a_sq = ai * ai;
        const double fsm = vi * small_t_fpt_scalar(tt[i] * vi, 0.001 / a_sq);
        fpt[i] = fsm;
      }
    }
  }

  // Joint density computation
  for (int i = 0; i < N; ++i) {
    if (!is_ok[i]) { out[i] = R_NegInf; continue; }

    const double tti = tt[i];
    if (!(tti > 0.0)) { out[i] = R_NegInf; continue; }

    const double ai = a[i];
    const double vi = v[i];
    const double thi = th[i];
    const double svi = sv[i];
    const double cosR = std::cos(Rs[i]);
    const double sinR = std::sin(Rs[i]);
    const double x0 = ai * cosR;
    const double x1 = ai * sinR;

    const double mu_x = vi * std::cos(thi);
    const double mu_y = vi * std::sin(thi);
    const double mu_sq = vi * vi;

    // Guard rail: kappa = a * |mu| / sigma^2 < 5
    const double kappa = (ai * std::sqrt(mu_sq)) / (sig[i] * sig[i]);
    if (!(kappa < 4.0)) { out[i] = R_NegInf; continue; }

    double dens = 0.0;
    if (!(svi > 0.0)) {
      // Constant drift
      const double term1 = ai * (mu_x * cosR + mu_y * sinR);
      const double term2 = 0.5 * mu_sq * tti;
      dens = std::exp(term1 - term2) * fpt[i];
    } else {
      // Integrated over Gaussian drift variability
      const double eta2 = svi * svi;
      const double denom = eta2 * tti + 1.0;
      const double fixed = 1.0 / std::sqrt(denom);

      const double exp0 = -0.5 * (mu_x * mu_x) / eta2 + 0.5 * (x0 * eta2 + mu_x) * (x0 * eta2 + mu_x) / (eta2 * denom);
      const double exp1 = -0.5 * (mu_y * mu_y) / eta2 + 0.5 * (x1 * eta2 + mu_y) * (x1 * eta2 + mu_y) / (eta2 * denom);
      dens = (fixed * std::exp(exp0)) * (fixed * std::exp(exp1)) * fpt[i];
    }

    // Normalize by angle support 2*pi and guard
    dens *= (1.0 / (2.0 * PI));
    if (!(dens > 0.0) || !R_finite(dens)) {
      out[i] = R_NegInf;
    } else {
      out[i] = std::log(dens);
    }
  }

  return out;
}


#endif // MODEL_CDM_H


