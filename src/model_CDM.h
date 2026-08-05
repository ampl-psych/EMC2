#ifndef MODEL_CDM_H
#define MODEL_CDM_H

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cfloat>

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

// Zeros of J1 and J2(zeros) for HSDM FPT expansion
static const double ZEROS_HSDM[] = {
  3.831705970207512, 7.015586669815619, 10.173468135062722, 13.323691936314223, 16.470630050877634,
  19.615858510468243, 22.760084380592772, 25.903672087618382, 29.046828534916855, 32.189679910974405,
  35.332307550083868, 38.474766234771614, 41.617094212814450, 44.759318997652819, 47.901460887185436,
  51.043535183571500, 54.185553641061318, 57.327525437901009, 60.469457845347492, 63.611356698481231,
  66.753226734098490, 69.895071837495777, 73.036895225573829, 76.178699584641464, 79.320487175476302,
  82.462259914373561, 85.604019436350228, 88.745767144926305, 91.887504251694992, 95.029231808044699,
  98.170950730790778, 101.312661823038724, 104.454365791282754, 107.596063259509165, 110.737754780899209,
  113.879440847594992, 117.021121898892432, 120.162798328149009, 123.304470488635715, 126.446138698516592,
  129.587803245103999, 132.729464388509626, 135.871122364788988, 139.012777388659714, 142.154429655859019,
  145.296079345195920, 148.437726620342232, 151.579371631401415, 154.721014516285948, 157.862655401930311,
  161.004294405362003, 164.145931634649628, 167.287567189744095, 170.429201163226651, 173.570833640975934,
  176.712464702763782, 179.854094422788393, 182.995722870152974, 186.137350109295539, 189.278976200376036,
  192.420601199625736, 195.562225159662603, 198.703848129777072, 201.845470156190913, 204.987091282292369,
  208.128711548850021, 211.270330994207796, 214.411949654461921, 217.553567563624142, 220.695184753769325,
  223.836801255171679, 226.978417096429439, 230.120032304579070, 233.261646905200593, 236.403260922514278,
  239.544874379469832, 242.686487297828677, 245.828099698239782, 248.969711600309921, 252.111323022668557,
  255.252933983028100, 258.394544498239497, 261.536154584344047, 264.677764256621458, 267.819373529634561,
  270.960982417270714, 274.102590932780686, 277.244199088814526, 280.385806897455552, 283.527414370251392,
  286.669021518243426, 289.810628351994410, 292.952234881613890, 296.093841116782471, 299.235447066774100,
  302.377052740477495, 305.518658146415589, 308.660263292764398, 311.801868187370417, 314.943472837767160
};

static const double JVZ_HSDM[] = {
  0.402759395702553, -0.300115752526133, 0.249704877057843, -0.218359407247873, 0.196465371468657,
  -0.180063375344316, 0.167184600473818, -0.156724986252852, 0.148011109972778, -0.140605798183982,
  0.134211240310001, -0.128616622072070, 0.123667960769837, -0.119249812010690, 0.115273694120168,
  -0.111670496859211, 0.108385348943683, -0.105374055395235, 0.102600567103397, -0.100035146811523,
  0.097653015783173, -0.095433339020535, 0.093358453290455, -0.091413272155921, 0.089584821964856,
  -0.087861876039410, 0.086234663413288, -0.084694634803424, 0.083234272982226, -0.081846937926486,
  0.080526739448403, -0.079268431724519, 0.078067325407949, -0.076919213961391, 0.075820311569167,
  -0.074767200537075, 0.073756786512857, -0.072786260189239, 0.071853064408847, -0.070954865793037,
  0.070089530177261, -0.069255101263766, 0.068449782005188, -0.067671918315546, 0.066919984772397,
  -0.066192572028753, 0.065488375698257, -0.064806186514098, 0.064144881592667, -0.063503416658322,
  0.062880819106749, -0.062276181802089, 0.061688657517819, -0.061117453943897, 0.060561829193241,
  -0.060021087749576, 0.059494576806311, -0.058981682952622, 0.058481829168482, -0.057994472095174,
  0.057519099551918, -0.057055228272818, 0.056602401841390, -0.056160188802597, 0.055728180934662,
  -0.055305991664887, 0.054893254615537, -0.054489622267340, 0.054094764729526, -0.053708368606511,
  0.053330135952378, -0.052959783305230, 0.052597040794298, -0.052241651313420, 0.051893369755131,
  -0.051551962300183, 0.051217205757823, -0.050888886952581, 0.050566802153761, -0.050250756544151,
  0.049940563724815, -0.049636045253094, 0.049337030211228, -0.049043354803226, 0.048754861977815,
  -0.048471401075509, 0.048192827497991, -0.047919002398152, 0.047649792389297, -0.047385069272114,
  0.047124709778146, -0.046868595328605, 0.046616611807442, -0.046368649347702, 0.046124602130248,
  -0.045884368194008, 0.045647849256998, -0.045414950547369, 0.045185580643861, -0.044959651325018
};

inline const std::vector<double>& weights() {
  static std::vector<double> w;
  if (w.empty()) {
    const size_t M = sizeof(ZEROS)/sizeof(ZEROS[0]);
    w.resize(M);
    for (size_t j = 0; j < M; ++j) {
      w[j] = ZEROS[j] / JVZ[j];
    }
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

inline const std::vector<double>& weights_hsdm() {
  static std::vector<double> w;
  if (w.empty()) {
    const size_t M = sizeof(ZEROS_HSDM)/sizeof(ZEROS_HSDM[0]);
    w.resize(M);
    for (size_t j = 0; j < M; ++j) {
      w[j] = (ZEROS_HSDM[j] * ZEROS_HSDM[j]) / JVZ_HSDM[j];
    }
  }
  return w;
}

inline const std::vector<double>& base_rates_hsdm() {
  static std::vector<double> br;
  if (br.empty()) {
    const size_t M = sizeof(ZEROS_HSDM)/sizeof(ZEROS_HSDM[0]);
    br.resize(M);
    for (size_t j = 0; j < M; ++j) br[j] = 0.5 * ZEROS_HSDM[j] * ZEROS_HSDM[j];
  }
  return br;
}

inline double series_bessel_fpt_scalar(double t, double a, double sigma) {
  if (!(t >= 0.0)) return 0.0;
  const std::vector<double>& w = weights();
  const std::vector<double>& br = base_rates();
  const double v_scale = (sigma * sigma) / (a * a);
  const double s = t * v_scale;
  double sum = 0.0;
  for (size_t j = 0; j < w.size(); ++j) {
    sum += w[j] * std::exp(-br[j] * s);
  }
  return v_scale * sum;
}

inline double series_sdm_fpt_scalar(double t, double a, double sigma, int max_n = 500) {
  if (!(t >= 0.0)) return 0.0;
  const double v_scale = (sigma * sigma) / (a * a);
  const double s = t * v_scale;
  const double PI = 3.14159265358979323846264338327950288;
  const double pi2 = PI * PI;
  double sum = 0.0;
  for (int n = 1; n <= max_n; ++n) {
    const double nn = static_cast<double>(n) * static_cast<double>(n);
    const double sign = (n % 2 == 0) ? -1.0 : 1.0;
    const double lambda = 0.5 * pi2 * nn;
    const double w = sign * pi2 * nn;
    sum += w * std::exp(-lambda * s);
  }
  return v_scale * sum;
}

inline double series_hsdm_fpt_scalar(double t, double a, double sigma) {
  if (!(t >= 0.0)) return 0.0;
  const std::vector<double>& w = weights_hsdm();
  const std::vector<double>& br = base_rates_hsdm();
  const double v_scale = (sigma * sigma) / (a * a);
  const double s = t * v_scale;
  double sum = 0.0;
  for (size_t j = 0; j < w.size(); ++j) {
    sum += w[j] * std::exp(-br[j] * s);
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
  const double term1 = ((1.0 - x) * (1.0 + t) * (1.0 + t)) / denom;
  const double term2 = std::exp(-0.5 * (1.0 - x) * (1.0 - x) / t - 0.5 * z1 * z1 * t);
  const double out = term1 * term2;
  if (!R_finite(out)) return 0.0;
  return out;
}

inline double small_t_fpt_sdm_scalar(double t_scaled, double x_scaled) {
  if (!(t_scaled > 0.0)) return 0.0;
  const double x = x_scaled;
  const double t = t_scaled;
  const double PI = 3.14159265358979323846264338327950288;
  const double denom = (x + t) * std::pow(t, 1.5);
  if (!(denom > 0.0)) return 0.0;
  const double term1 = ((1.0 - x) * std::pow(1.0 + t, 2.5)) / denom;
  const double term2 = std::exp(-0.5 * (1.0 - x) * (1.0 - x) / t - 0.5 * PI * PI * t);
  const double out = term1 * term2;
  if (!R_finite(out)) return 0.0;
  return out;
}

inline double small_t_fpt_hsdm_scalar(double t_scaled, double x_scaled) {
  if (!(t_scaled > 0.0)) return 0.0;
  const double z1 = ZEROS_HSDM[0];
  const double x = x_scaled;
  const double t = t_scaled;
  const double denom = (x + t) * std::sqrt(x + t) * std::pow(t, 1.5);
  if (!(denom > 0.0)) return 0.0;
  const double term1 = ((1.0 - x) * std::pow(1.0 + t, 3.0)) / denom;
  const double term2 = std::exp(-0.5 * (1.0 - x) * (1.0 - x) / t - 0.5 * z1 * z1 * t);
  const double out = term1 * term2;
  if (!R_finite(out)) return 0.0;
  return out;
}

inline double log_density_core_scalar(double log_fpt, double tt, double a, double v, double sig2,
                                      double sv, double cos_term, int dim, double log_surface) {
  const double v2 = v * v;
  if (sv <= 0.0) {
    return log_fpt - 0.5 * (v2 * tt) / sig2 + (a * v * cos_term) / sig2 - log_surface;
  }
  const double sv2 = sv * sv;
  const double tt_v = std::max(tt, DBL_EPSILON);
  const double D = sig2 + sv2 * tt_v;
  const double log_det = 0.5 * static_cast<double>(dim) * (std::log(sig2) - std::log(D));
  const double log_rad = 0.5 * (a * a * sv2) / (sig2 * D);
  const double log_ang = (a * v * cos_term) / D;
  const double log_tim = -0.5 * (v2 * tt_v) / D;
  return log_fpt + log_det + log_rad + log_ang + log_tim - log_surface;
}

inline double log_i0_stable_scalar(double x) {
  const double ax = std::fabs(x);
  double i0_scaled = R::bessel_i(ax, 0.0, 2.0);
  if (!R_finite(i0_scaled) || i0_scaled <= 0.0) i0_scaled = DBL_MIN;
  return ax + std::log(i0_scaled);
}

} // namespace cdm_internal

// Returns log-density per trial (R_NegInf if invalid)
inline NumericVector c_dCDM(NumericVector rts, NumericVector Rs, NumericMatrix pars, LogicalVector is_ok) {
  using namespace cdm_internal;
  const int N = rts.size();
  NumericVector out(N);

  // Extract columns by index according to R model: v, theta, a, t0, s, sv
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

  const double s0 = 0.002;
  const double s1 = 0.02;
  const double two_pi = 2.0 * PI;

  // Joint density computation
  for (int i = 0; i < N; ++i) {
    if (!is_ok[i]) { out[i] = R_NegInf; continue; }

    const double ai = a[i];
    const double sigi = sig[i];
    const double vi = v[i];
    const double thi = th[i];
    const double svi = sv[i];
    const double rti = rts[i];
    const double t0i = t0[i];

    if (!R_finite(ai) || !R_finite(sigi) || ai <= 0.0 || sigi <= 0.0 ||
        !R_finite(vi) || !R_finite(thi) || !R_finite(rti) || !R_finite(t0i)) {
      out[i] = R_NegInf; continue;
    }

    const double tt = std::max(0.0, rti - t0i);
    const double sig2 = sigi * sigi;
    const double v2 = vi * vi;

    const double cosR = std::cos(Rs[i]);
    const double sinR = std::sin(Rs[i]);
    const double cosT = std::cos(thi);
    const double sinT = std::sin(thi);
    const double cosRth = cosR * cosT + sinR * sinT;

    // Zero-drift FPT with small-t blend (matches R implementation)
    const double v_scale = sig2 / (ai * ai);
    const double s = tt * v_scale;
    const double fpt_series = series_bessel_fpt_scalar(tt, ai, sigi);
    const double fpt_small = v_scale * small_t_fpt_scalar(s, 1e-8 / (ai * ai));
    const double w = std::min(std::max((s - s0) / (s1 - s0), 0.0), 1.0);
    double fpt = (1.0 - w) * fpt_small + w * fpt_series;
    if (!R_finite(fpt) || fpt < 0.0) fpt = 0.0;
    const double log_fpt = std::log(std::max(fpt, DBL_MIN));

    double log_raw = R_NegInf;
    if (svi <= 0.0) {
      // Constant drift
      const double log_time = -0.5 * (v2 * tt) / sig2;
      const double log_ang = (ai * vi * cosRth) / sig2;
      log_raw = log_fpt + log_time + log_ang - std::log(two_pi);
    } else {
      // Integrated Gaussian drift variability
      const double sv2 = svi * svi;
      const double tt_v = std::max(tt, DBL_EPSILON);
      const double D = sig2 + sv2 * tt_v;
      const double log_det = -std::log(D) + std::log(sig2);
      const double log_rad = 0.5 * (ai * ai * sv2) / (sig2 * D);
      const double log_ang = (ai * vi * cosRth) / D;
      const double log_tim = -0.5 * (v2 * tt_v) / D;
      log_raw = log_fpt + log_det + log_rad + log_ang + log_tim - std::log(two_pi);
    }

    out[i] = log_raw;
  }

  return out;
}

// Returns log-density per trial for SDM (R, R2)
inline NumericVector c_dSDM(NumericVector rts, NumericVector Rs, NumericVector R2s,
                            NumericMatrix pars, LogicalVector is_ok) {
  using namespace cdm_internal;
  const int N = rts.size();
  NumericVector out(N);

  NumericVector v   = pars(_, 0);
  NumericVector th1 = pars(_, 1);
  NumericVector th2 = pars(_, 2);
  NumericVector a   = pars(_, 3);
  NumericVector t0  = pars(_, 4);
  NumericVector sig = pars(_, 5);
  NumericVector sv  = (pars.ncol() > 6) ? pars(_, 6) : NumericVector(N, 0.0);

  const double PI = 3.14159265358979323846264338327950288;
  const double two_pi = 2.0 * PI;
  for (int i = 0; i < N; ++i) {
    th1[i] = th1[i] * PI;
    th2[i] = (th2[i] - 0.5) * two_pi;
  }

  const double s0 = 0.002;
  const double s1 = 0.02;
  const double log_surface = std::log(4.0 * PI);

  for (int i = 0; i < N; ++i) {
    if (!is_ok[i]) { out[i] = R_NegInf; continue; }
    const double ai = a[i];
    const double sigi = sig[i];
    const double vi = v[i];
    const double t1 = th1[i];
    const double t2 = th2[i];
    const double svi = sv[i];
    const double rti = rts[i];
    const double t0i = t0[i];
    const double r1 = Rs[i];
    const double r2 = R2s[i];

    if (!R_finite(ai) || !R_finite(sigi) || ai <= 0.0 || sigi <= 0.0 ||
        !R_finite(vi) || !R_finite(t1) || !R_finite(t2) || !R_finite(rti) || !R_finite(t0i) ||
        !R_finite(r1) || !R_finite(r2)) {
      out[i] = R_NegInf; continue;
    }

    const double tt = std::max(0.0, rti - t0i);
    const double sig2 = sigi * sigi;
    const double v_scale = sig2 / (ai * ai);
    const double s = tt * v_scale;
    const double fpt_series = series_sdm_fpt_scalar(tt, ai, sigi);
    const double fpt_small = v_scale * small_t_fpt_sdm_scalar(s, 1e-8 / (ai * ai));
    const double w = std::min(std::max((s - s0) / (s1 - s0), 0.0), 1.0);
    double fpt = (1.0 - w) * fpt_small + w * fpt_series;
    if (!R_finite(fpt) || fpt < 0.0) fpt = 0.0;
    const double log_fpt = std::log(std::max(fpt, DBL_MIN));

    const double sr = std::sin(r1);
    const double cr = std::cos(r1);
    const double x0 = cr;
    const double x1 = sr * std::cos(r2);
    const double x2 = sr * std::sin(r2);

    const double st1 = std::sin(t1);
    const double m0 = std::cos(t1);
    const double m1 = st1 * std::cos(t2);
    const double m2 = st1 * std::sin(t2);
    const double cos_term = x0 * m0 + x1 * m1 + x2 * m2;

    out[i] = log_density_core_scalar(log_fpt, tt, ai, vi, sig2, svi, cos_term, 3, log_surface);
  }

  return out;
}

// Returns log-density per trial for HSDM (R, R2, R3)
inline NumericVector c_dHSDM(NumericVector rts, NumericVector Rs, NumericVector R2s, NumericVector R3s,
                             NumericMatrix pars, LogicalVector is_ok) {
  using namespace cdm_internal;
  const int N = rts.size();
  NumericVector out(N);

  NumericVector v   = pars(_, 0);
  NumericVector th1 = pars(_, 1);
  NumericVector th2 = pars(_, 2);
  NumericVector th3 = pars(_, 3);
  NumericVector a   = pars(_, 4);
  NumericVector t0  = pars(_, 5);
  NumericVector sig = pars(_, 6);
  NumericVector sv  = (pars.ncol() > 7) ? pars(_, 7) : NumericVector(N, 0.0);

  const double PI = 3.14159265358979323846264338327950288;
  const double two_pi = 2.0 * PI;
  for (int i = 0; i < N; ++i) {
    th1[i] = th1[i] * PI;
    th2[i] = th2[i] * PI;
    th3[i] = (th3[i] - 0.5) * two_pi;
  }

  const double s0 = 0.002;
  const double s1 = 0.02;
  const double log_surface = std::log(2.0 * PI * PI);

  for (int i = 0; i < N; ++i) {
    if (!is_ok[i]) { out[i] = R_NegInf; continue; }
    const double ai = a[i];
    const double sigi = sig[i];
    const double vi = v[i];
    const double t1 = th1[i];
    const double t2 = th2[i];
    const double t3 = th3[i];
    const double svi = sv[i];
    const double rti = rts[i];
    const double t0i = t0[i];
    const double r1 = Rs[i];
    const double r2 = R2s[i];
    const double r3 = R3s[i];

    if (!R_finite(ai) || !R_finite(sigi) || ai <= 0.0 || sigi <= 0.0 ||
        !R_finite(vi) || !R_finite(t1) || !R_finite(t2) || !R_finite(t3) ||
        !R_finite(rti) || !R_finite(t0i) || !R_finite(r1) || !R_finite(r2) || !R_finite(r3)) {
      out[i] = R_NegInf; continue;
    }

    const double tt = std::max(0.0, rti - t0i);
    const double sig2 = sigi * sigi;
    const double v_scale = sig2 / (ai * ai);
    const double s = tt * v_scale;
    const double fpt_series = series_hsdm_fpt_scalar(tt, ai, sigi);
    const double fpt_small = v_scale * small_t_fpt_hsdm_scalar(s, 1e-8 / (ai * ai));
    const double w = std::min(std::max((s - s0) / (s1 - s0), 0.0), 1.0);
    double fpt = (1.0 - w) * fpt_small + w * fpt_series;
    if (!R_finite(fpt) || fpt < 0.0) fpt = 0.0;
    const double log_fpt = std::log(std::max(fpt, DBL_MIN));

    const double sr1 = std::sin(r1);
    const double cr1 = std::cos(r1);
    const double sr2 = std::sin(r2);
    const double cr2 = std::cos(r2);
    const double x0 = cr1;
    const double x1 = sr1 * cr2;
    const double x2 = sr1 * sr2 * std::cos(r3);
    const double x3 = sr1 * sr2 * std::sin(r3);

    const double st1 = std::sin(t1);
    const double st2 = std::sin(t2);
    const double m0 = std::cos(t1);
    const double m1 = st1 * std::cos(t2);
    const double m2 = st1 * st2 * std::cos(t3);
    const double m3 = st1 * st2 * std::sin(t3);
    const double cos_term = x0 * m0 + x1 * m1 + x2 * m2 + x3 * m3;

    out[i] = log_density_core_scalar(log_fpt, tt, ai, vi, sig2, svi, cos_term, 4, log_surface);
  }

  return out;
}

// Returns log-density per trial for PSDM (R)
inline NumericVector c_dPSDM(NumericVector rts, NumericVector Rs, NumericMatrix pars, LogicalVector is_ok) {
  using namespace cdm_internal;
  const int N = rts.size();
  NumericVector out(N);

  NumericVector v   = pars(_, 0);
  NumericVector th1 = pars(_, 1);
  NumericVector a   = pars(_, 2);
  NumericVector t0  = pars(_, 3);
  NumericVector sig = pars(_, 4);
  NumericVector sv  = (pars.ncol() > 5) ? pars(_, 5) : NumericVector(N, 0.0);

  const double PI = 3.14159265358979323846264338327950288;
  for (int i = 0; i < N; ++i) {
    th1[i] = th1[i] * PI;
  }

  const double s0 = 0.002;
  const double s1 = 0.02;
  const double log_two = std::log(2.0);

  for (int i = 0; i < N; ++i) {
    if (!is_ok[i]) { out[i] = R_NegInf; continue; }
    const double ai = a[i];
    const double sigi = sig[i];
    const double vi = v[i];
    const double t1 = th1[i];
    const double svi = sv[i];
    const double rti = rts[i];
    const double t0i = t0[i];
    const double r1 = Rs[i];

    if (!R_finite(ai) || !R_finite(sigi) || ai <= 0.0 || sigi <= 0.0 ||
        !R_finite(vi) || !R_finite(t1) || !R_finite(rti) || !R_finite(t0i) || !R_finite(r1)) {
      out[i] = R_NegInf; continue;
    }

    const double tt = std::max(0.0, rti - t0i);
    const double sig2 = sigi * sigi;
    const double v2 = vi * vi;
    const double v_scale = sig2 / (ai * ai);
    const double s = tt * v_scale;
    const double fpt_series = series_sdm_fpt_scalar(tt, ai, sigi);
    const double fpt_small = v_scale * small_t_fpt_sdm_scalar(s, 1e-8 / (ai * ai));
    const double w = std::min(std::max((s - s0) / (s1 - s0), 0.0), 1.0);
    double fpt = (1.0 - w) * fpt_small + w * fpt_series;
    if (!R_finite(fpt) || fpt < 0.0) fpt = 0.0;
    const double log_fpt = std::log(std::max(fpt, DBL_MIN));

    const double A = std::cos(t1) * std::cos(r1);
    const double B = std::sin(t1) * std::sin(r1);

    double log_raw = R_NegInf;
    if (svi <= 0.0) {
      const double k = (ai * vi) / sig2;
      const double log_base = log_fpt - 0.5 * (v2 * tt) / sig2;
      log_raw = log_base + k * A + log_i0_stable_scalar(k * B) - log_two;
    } else {
      const double sv2 = svi * svi;
      const double tt_v = std::max(tt, DBL_EPSILON);
      const double D = sig2 + sv2 * tt_v;
      const double k = (ai * vi) / D;
      const double log_base = log_fpt +
        0.5 * 3.0 * (std::log(sig2) - std::log(D)) +
        0.5 * (ai * ai * sv2) / (sig2 * D) -
        0.5 * (v2 * tt_v) / D;
      log_raw = log_base + k * A + log_i0_stable_scalar(k * B) - log_two;
    }
    out[i] = log_raw;
  }

  return out;
}

// Returns log-density per trial for PHSDM (R, R2)
inline NumericVector c_dPHSDM(NumericVector rts, NumericVector Rs, NumericVector R2s,
                              NumericMatrix pars, LogicalVector is_ok) {
  using namespace cdm_internal;
  const int N = rts.size();
  NumericVector out(N);

  NumericVector v   = pars(_, 0);
  NumericVector th1 = pars(_, 1);
  NumericVector th2 = pars(_, 2);
  NumericVector a   = pars(_, 3);
  NumericVector t0  = pars(_, 4);
  NumericVector sig = pars(_, 5);
  NumericVector sv  = (pars.ncol() > 6) ? pars(_, 6) : NumericVector(N, 0.0);

  const double PI = 3.14159265358979323846264338327950288;
  for (int i = 0; i < N; ++i) {
    th1[i] = th1[i] * PI;
    th2[i] = th2[i] * PI;
  }

  const double s0 = 0.002;
  const double s1 = 0.02;
  const double log_pi = std::log(PI);

  for (int i = 0; i < N; ++i) {
    if (!is_ok[i]) { out[i] = R_NegInf; continue; }
    const double ai = a[i];
    const double sigi = sig[i];
    const double vi = v[i];
    const double t1 = th1[i];
    const double t2 = th2[i];
    const double svi = sv[i];
    const double rti = rts[i];
    const double t0i = t0[i];
    const double r1 = Rs[i];
    const double r2 = R2s[i];

    if (!R_finite(ai) || !R_finite(sigi) || ai <= 0.0 || sigi <= 0.0 ||
        !R_finite(vi) || !R_finite(t1) || !R_finite(t2) || !R_finite(rti) || !R_finite(t0i) ||
        !R_finite(r1) || !R_finite(r2)) {
      out[i] = R_NegInf; continue;
    }

    const double tt = std::max(0.0, rti - t0i);
    const double sig2 = sigi * sigi;
    const double v2 = vi * vi;
    const double v_scale = sig2 / (ai * ai);
    const double s = tt * v_scale;
    const double fpt_series = series_hsdm_fpt_scalar(tt, ai, sigi);
    const double fpt_small = v_scale * small_t_fpt_hsdm_scalar(s, 1e-8 / (ai * ai));
    const double w = std::min(std::max((s - s0) / (s1 - s0), 0.0), 1.0);
    double fpt = (1.0 - w) * fpt_small + w * fpt_series;
    if (!R_finite(fpt) || fpt < 0.0) fpt = 0.0;
    const double log_fpt = std::log(std::max(fpt, DBL_MIN));

    const double st1 = std::sin(t1);
    const double ct1 = std::cos(t1);
    const double st2 = std::sin(t2);
    const double ct2 = std::cos(t2);
    const double sr1 = std::sin(r1);
    const double cr1 = std::cos(r1);
    const double sr2 = std::sin(r2);
    const double cr2 = std::cos(r2);
    const double A = ct1 * cr1 + st1 * ct2 * sr1 * cr2;
    const double B = st1 * st2 * sr1 * sr2;

    double log_raw = R_NegInf;
    if (svi <= 0.0) {
      const double k = (ai * vi) / sig2;
      const double log_base = log_fpt - 0.5 * (v2 * tt) / sig2;
      log_raw = log_base + k * A + log_i0_stable_scalar(k * B) - log_pi;
    } else {
      const double sv2 = svi * svi;
      const double tt_v = std::max(tt, DBL_EPSILON);
      const double D = sig2 + sv2 * tt_v;
      const double k = (ai * vi) / D;
      const double log_base = log_fpt +
        0.5 * 4.0 * (std::log(sig2) - std::log(D)) +
        0.5 * (ai * ai * sv2) / (sig2 * D) -
        0.5 * (v2 * tt_v) / D;
      log_raw = log_base + k * A + log_i0_stable_scalar(k * B) - log_pi;
    }
    out[i] = log_raw;
  }

  return out;
}


#endif // MODEL_CDM_H
