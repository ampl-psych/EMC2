// K-category joint neural likelihood (kind "flow_cat"): the flow_joint design
// (flow_ddm.cpp) generalised from two responses to K categories, with the
// non-decision time applied OUTSIDE the network by an exact uniform
// convolution and with censoring at per-trial bounds. Written by the rtcon2
// project (stage R2) for RTCON2 (Ratcliff & Starns, 2013), where the category
// is the (choice, rating) pair k = rating_fold(R, RR). Stage R3 adds the `mix`
// block: a latent trial strength w ~ N(mu, s^2) integrated by Gauss-Hermite
// quadrature outside the network, the network conditioning on the race's own
// parameters (the drifts a v_k(w), sb, B_1..B_K) only.
//
// Model. For a design row with context theta (the network's n_ctx inputs on
// the sampled scale) the network gives
//   log P(k | theta), k = 1..K      a K-way log-softmax head (classifier_mlp,
//                                   n_ctx inputs, K outputs)
//   F_k(tau), f_k(tau)              the flow's defective CDF / density of the
//                                   DECISION TIME tau for category k: a
//                                   rational-quadratic spline on log tau whose
//                                   knots come from flow_mlp on [theta, one-hot(k)]
//                                   (n_ctx + K inputs), p_k Phi(z_k(log tau)).
// With rt = tau + t_er, t_er ~ U(t0 - st0/2, t0 + st0/2) independent of tau
// (the card's `conv`: mean = t0, range = st0, EMC2 parameters that are not
// network inputs):
//   observed (rt, k):   L = [F_k(rt - t0 + st0/2) - F_k(rt - t0 - st0/2)] / st0
//                       (st0 = 0: f_k(rt - t0))
//   lower censored:     L = P(rt < LC) = (1/st0) int sum_k F_k(LC - t0 - u) du over
//                       u in [-st0/2, st0/2]  (16-node Gauss-Legendre; one node at st0 = 0)
//   upper censored:     L = 1 - the same at UC
// F_k(tau) = 0 for tau <= 0. Out-of-box rows: log L = -Inf (rejection sentinel).
//
// Mix (card `mix`, stage R3). The registration-facing inputs are the model's
// natural parameters per design row (context_names, n_src of them: for RTCON2
// mu, s, C1, g2..g5, a, sb, B1..B6). Per row and quadrature node m,
//   w_m = mu + sqrt(2) s x_m,  weight omega_m (gh_w, sum 1),
//   c = cumsum(C1, g2..g5) with the simulator's middle-out clamp,
//   v_1 = Phi(c_1 - w), v_k = Phi(c_k - w) - Phi(c_{k-1} - w), v_K = 1 - Phi(c_{K-1} - w),
//   network inputs (mix.net_names): a v_1..a v_K, then the direct sources
//   (sb, B_1..B_K), each on its own transform (mix.net_transforms) and then the
//   scaler; one head call and K flow conditionings per (row, node), and per
//   trial the mixture over nodes of the quantities above:
//   P(k) = sum_m omega_m p_m(k), f(rt, k) = sum_m omega_m p_m(k) f_m(rt - t_er | k),
//   censored masses likewise. The in-box check is on the natural sources
//   against bounds_natural (the training box of the sources).
//
// API:
//   cat_build(fl)                               -> external pointer
//   cat_eval_trials_cpp(ptr, Theta, rt, t0, st0, k, miss, LC, UC)
//        trial-wise rows (Theta n x n_ctx on the sampled scale, or for a mix
//        card n x n_src natural sources; k 1..K, NA on censored rows; miss 0
//        observed / 1 lower / 2 upper) -> log_pdf, and for observed rows the
//        convolved defective CDF of category k at rt
//   native branch: nle_cat_native() (nle_native.h, model_NN.h)

#include "nle_flow.h"
#include "nle_native.h"
using namespace Rcpp;

struct MixSpec {
  bool on = false;
  int M = 1;                                  // quadrature nodes (1: no mixture)
  std::vector<double> gh_x, gh_logw;          // nodes, log weights (weights sum to 1)
  double clamp = 0.25;                        // minimum criterion gap (the simulator's mind)
  int n_src = 0;                              // registration-facing natural inputs per row
  std::vector<double> src_lower, src_upper;   // their training box (natural scale, closed)
  std::vector<int> net_tf;                    // transform code per network input (0 identity, 1 log, 2 probit, 3 zerobox)
  std::vector<double> net_tfc;                // zerobox constant per network input (NaN otherwise)
  int i_mu = -1, i_s = -1, i_a = -1;          // source columns
  std::vector<int> i_c;                       // source columns of C1, g2.. (K - 1 of them)
  std::vector<int> direct_src;                // network inputs K.. : source column each
};

struct CatModel {
  int n_ctx = 0, K = 0;
  nle::Spline sp;
  std::vector<double> scaler_mean, scaler_scale;
  std::vector<double> lower, upper;       // bounding box (sampled scale) of the network inputs (plain cards)
  nle::Mlp flow;                          // input n_ctx + K (one-hot category)
  nle::Mlp clf;                           // input n_ctx, output K logits
  MixSpec mix;
};

static const int GL_N = 16;
static const double GL_X[GL_N] = {-0.98940093499164994, -0.9445750230732326, -0.86563120238783164, -0.75540440835500311, -0.61787624440264377, -0.45801677765722737, -0.28160355077925892, -0.095012509837637441, 0.095012509837637441, 0.28160355077925892, 0.45801677765722737, 0.61787624440264377, 0.75540440835500311, 0.86563120238783164, 0.9445750230732326, 0.98940093499164994};
static const double GL_W[GL_N] = {0.027152459411754055, 0.062253523938647609, 0.095158511682492994, 0.12462897125553395, 0.14959598881657665, 0.16915651939500262, 0.18260341504492361, 0.18945061045506859, 0.18945061045506859, 0.18260341504492361, 0.16915651939500262, 0.14959598881657665, 0.12462897125553395, 0.095158511682492994, 0.062253523938647609, 0.027152459411754055};

// log(Phi(z_hi) - Phi(z_lo)) for z_hi >= z_lo, accurate in both tails; -Inf when empty.
static inline double log_phi_diff(double z_hi, double z_lo) {
  if (!(z_hi > z_lo)) return R_NegInf;
  double v;
  if (z_lo > 0.0) v = R::pnorm(z_lo, 0.0, 1.0, 0, 0) - R::pnorm(z_hi, 0.0, 1.0, 0, 0);
  else            v = R::pnorm(z_hi, 0.0, 1.0, 1, 0) - R::pnorm(z_lo, 0.0, 1.0, 1, 0);
  return v > 0.0 ? std::log(v) : R_NegInf;
}

// log(sum_m exp(a_m)) over M values (-Inf allowed).
static inline double logsumexp_n(const double* a, int M) {
  double mx = R_NegInf;
  for (int m = 0; m < M; ++m) mx = std::fmax(mx, a[m]);
  if (!std::isfinite(mx)) return mx;
  double s = 0.0;
  for (int m = 0; m < M; ++m) s += std::exp(a[m] - mx);
  return mx + std::log(s);
}

// The mix's network inputs (natural scale, n_ctx values) of source row src at node w.
static inline void mix_net_natural(const CatModel& m, const double* src, double w, double* out) {
  const MixSpec& mx = m.mix; const int K = m.K;
  double c[32];
  c[0] = src[mx.i_c[0]];
  for (int i = 1; i < K - 1; ++i) c[i] = c[i - 1] + src[mx.i_c[i]];
  // the simulator's middle-out clamp (perturb_criteria with csd = 0)
  const int nc = K - 1, middle = (nc - 1) / 2;
  for (int i = 0; i < (nc - 1) / 2; ++i) {
    c[middle + i + 1] = std::fmax(c[middle + i + 1], c[middle + i] + mx.clamp);
    c[middle - i - 1] = std::fmin(c[middle - i - 1], c[middle - i] - mx.clamp);
  }
  if (nc % 2 == 0) c[nc - 1] = std::fmax(c[nc - 1], c[nc - 2] + mx.clamp);
  const double a = src[mx.i_a];
  double cum = 0.0;
  for (int i = 0; i < K - 1; ++i) {
    const double F = R::pnorm(c[i] - w, 0.0, 1.0, 1, 0);
    out[i] = a * (F - cum); cum = F;
  }
  out[K - 1] = a * (1.0 - cum);
  for (size_t j = 0; j < mx.direct_src.size(); ++j) out[K + j] = src[mx.direct_src[j]];
}

// Conditioning: Xnet holds the network-scale (pre-scaler) inputs of every (in-box
// row, node) as columns, column ri * M + mm. Returns the head's log-softmax per
// column (lpk, K per column) and the flow's knots per (column, category) in kb
// (knot column = col * K + q).
static void cat_condition(const CatModel& m, const arma::mat& Xnet, std::vector<double>& lpk, nle::KnotBlock& kb) {
  const int nc = m.n_ctx, K = m.K, NB = m.sp.num_bins, C = (int)Xnet.n_cols;
  lpk.assign((size_t)C * K, 0.0);
  for (int c0 = 0; c0 < C; c0 += nle::NLE_BLOCK) {
    const int nb = std::min(nle::NLE_BLOCK, C - c0);
    arma::mat X(nc, nb);
    for (int c = 0; c < nb; ++c)
      for (int j = 0; j < nc; ++j) X(j, c) = (Xnet(j, c0 + c) - m.scaler_mean[j]) / m.scaler_scale[j];
    const arma::mat logit = nle::mlp_forward_batch(m.clf, X);
    for (int c = 0; c < nb; ++c) {
      double mxv = logit(0, c);
      for (int q = 1; q < K; ++q) mxv = std::fmax(mxv, logit(q, c));
      double s = 0.0;
      for (int q = 0; q < K; ++q) s += std::exp(logit(q, c) - mxv);
      const double lse = mxv + std::log(s);
      for (int q = 0; q < K; ++q) lpk[(size_t)(c0 + c) * K + q] = logit(q, c) - lse;
    }
  }
  const int F = C * K;
  kb.reset(NB, std::max(F, 1));
  for (int c0 = 0; c0 < F; c0 += nle::NLE_BLOCK) {
    const int nb = std::min(nle::NLE_BLOCK, F - c0);
    arma::mat X(nc + K, nb, arma::fill::zeros);
    for (int c = 0; c < nb; ++c) {
      const int f = c0 + c, col = f / K, q = f % K;
      for (int j = 0; j < nc; ++j) X(j, c) = (Xnet(j, col) - m.scaler_mean[j]) / m.scaler_scale[j];
      X(nc + q, c) = 1.0;
    }
    const arma::mat raw = nle::mlp_forward_batch(m.flow, X);
    for (int c = 0; c < nb; ++c) nle::build_knots(m.sp, raw.colptr(c), kb, c0 + c);
  }
}

// Trials t = 0..n-1 with distinct parameter row uid[t] in 0..U-1; inb_id[u] is the
// row's index among the in-box rows (-1 out of box); M nodes with log weights
// logw (M = 1, logw = {0} without a mixture). lpk / kb from cat_condition.
// log_pdf (n) always; cdf (n, may be null): the convolved defective CDF of
// category k at rt for observed rows (NaN for censored rows).
static void cat_eval_trials(const CatModel& m, const std::vector<int>& inb_id, int M, const std::vector<double>& logw,
                            const std::vector<double>& lpk, const nle::KnotBlock& kb, const int* uid,
                            const double* rt, const double* t0, const double* st0, const int* k,
                            const int* miss, const double* LC, const double* UC, int n,
                            double* log_pdf, double* cdf) {
  const int K = m.K, NB = m.sp.num_bins;
  std::vector<double> acc(M);
  auto z_of = [&](int col, double tau, double& logdet) -> double {
    if (!(tau > 0.0)) { logdet = 0.0; return R_NegInf; }
    double z;
    nle::rqs_inverse_one(kb.xc(col), kb.yc(col), kb.dc(col), NB, std::log(tau), z, logdet);
    return z;
  };
  auto Phi = [](double z) -> double { return std::isfinite(z) ? R::pnorm(z, 0.0, 1.0, 1, 0) : (z < 0 ? 0.0 : 1.0); };
  // total CDF at decision time tau (window node) for row ri: sum_m w_m sum_q p_mq Phi(z_mq(tau))
  auto Ftot = [&](int ri, double tau) -> double {
    if (!(tau > 0.0)) return 0.0;
    double F = 0.0, ld;
    for (int mm = 0; mm < M; ++mm) {
      const int col = ri * M + mm; double Fm = 0.0;
      for (int q = 0; q < K; ++q) Fm += std::exp(lpk[(size_t)col * K + q]) * Phi(z_of(col * K + q, tau, ld));
      F += std::exp(logw[mm]) * Fm;
    }
    return F;
  };
  for (int t = 0; t < n; ++t) {
    const int u = uid[t];
    if (cdf) cdf[t] = NA_REAL;
    if (inb_id[u] < 0) { log_pdf[t] = R_NegInf; continue; }
    const int ri = inb_id[u];
    const double T0 = t0[t], ST = st0[t];
    if (!(T0 >= 0.0) || !(ST >= 0.0)) { log_pdf[t] = R_NegInf; continue; }
    double ld;
    if (miss[t] == 0) {
      const int q = k[t];
      if (q < 0 || q >= K) { log_pdf[t] = R_NegInf; continue; }
      const double tau = rt[t] - T0;
      if (ST == 0.0) {
        if (!(tau > 0.0)) { log_pdf[t] = R_NegInf; if (cdf) cdf[t] = 0.0; continue; }
        double cF = 0.0;
        for (int mm = 0; mm < M; ++mm) {
          const int col = ri * M + mm; const double lp = lpk[(size_t)col * K + q];
          const double z = z_of(col * K + q, tau, ld);
          acc[mm] = logw[mm] + lp + R::dnorm(z, 0.0, 1.0, 1) + ld - std::log(tau);
          if (cdf) cF += std::exp(logw[mm] + lp) * Phi(z);
        }
        log_pdf[t] = logsumexp_n(acc.data(), M);
        if (cdf) cdf[t] = cF;
      } else {
        double cF = 0.0;
        for (int mm = 0; mm < M; ++mm) {
          const int col = ri * M + mm; const double lp = lpk[(size_t)col * K + q];
          const double z_hi = z_of(col * K + q, tau + 0.5 * ST, ld), z_lo = z_of(col * K + q, tau - 0.5 * ST, ld);
          acc[mm] = logw[mm] + lp + log_phi_diff(z_hi, z_lo) - std::log(ST);
          if (cdf) {
            double a2 = 0.0;
            for (int g = 0; g < GL_N; ++g) a2 += 0.5 * GL_W[g] * Phi(z_of(col * K + q, tau - 0.5 * ST * GL_X[g], ld));
            cF += std::exp(logw[mm] + lp) * a2;
          }
        }
        log_pdf[t] = logsumexp_n(acc.data(), M);
        if (cdf) cdf[t] = cF;
      }
    } else if (miss[t] == 1 || miss[t] == 2) {
      const double bound = (miss[t] == 1 ? LC[t] : UC[t]) - T0;
      double mass = 0.0;
      const int ng = (ST == 0.0) ? 1 : GL_N;
      for (int g = 0; g < ng; ++g) {
        const double tau = (ST == 0.0) ? bound : bound - 0.5 * ST * GL_X[g];
        const double w = (ST == 0.0) ? 1.0 : 0.5 * GL_W[g];
        mass += w * Ftot(ri, tau);
      }
      if (miss[t] == 2) mass = 1.0 - mass;
      log_pdf[t] = mass > 0.0 ? std::log(mass) : R_NegInf;
    } else log_pdf[t] = R_NegInf;                   // code 3 (both) and others: not supported
  }
}

// Plain card: Theta_u (U rows x n_ctx, row-major) on the sampled scale; in-box
// rows conditioned once (M = 1).
static void cat_eval_core(const CatModel& m, const double* Theta_u, int U, const int* uid,
                          const double* rt, const double* t0, const double* st0, const int* k,
                          const int* miss, const double* LC, const double* UC, int n,
                          double* log_pdf, double* cdf) {
  const int nc = m.n_ctx;
  std::vector<int> inb_id(U, -1); int Uin = 0;
  for (int u = 0; u < U; ++u)
    if (nle::in_box(m.lower, m.upper, Theta_u + (size_t)u * nc)) inb_id[u] = Uin++;
  arma::mat Xnet(nc, std::max(Uin, 1));
  for (int u = 0; u < U; ++u)
    if (inb_id[u] >= 0) for (int j = 0; j < nc; ++j) Xnet(j, inb_id[u]) = Theta_u[(size_t)u * nc + j];
  std::vector<double> lpk; nle::KnotBlock kb;
  if (Uin > 0) cat_condition(m, Xnet, lpk, kb);
  const std::vector<double> logw(1, 0.0);
  cat_eval_trials(m, inb_id, 1, logw, lpk, kb, uid, rt, t0, st0, k, miss, LC, UC, n, log_pdf, cdf);
}

// Mix card: Src_u (U rows x n_src, row-major) natural-scale sources; in-box rows
// conditioned at every node.
static void cat_eval_core_mix(const CatModel& m, const double* Src_u, int U, const int* uid,
                              const double* rt, const double* t0, const double* st0, const int* k,
                              const int* miss, const double* LC, const double* UC, int n,
                              double* log_pdf, double* cdf) {
  const MixSpec& mx = m.mix; const int nc = m.n_ctx, ns = mx.n_src, M = mx.M;
  std::vector<int> inb_id(U, -1); int Uin = 0;
  for (int u = 0; u < U; ++u)
    if (nle::in_box(mx.src_lower, mx.src_upper, Src_u + (size_t)u * ns)) inb_id[u] = Uin++;
  arma::mat Xnet(nc, std::max(Uin * M, 1));
  std::vector<double> nat(nc);
  for (int u = 0; u < U; ++u) {
    if (inb_id[u] < 0) continue;
    const double* src = Src_u + (size_t)u * ns;
    const double mu = src[mx.i_mu], s = src[mx.i_s];
    for (int mm = 0; mm < M; ++mm) {
      const double w = mu + M_SQRT2 * s * mx.gh_x[mm];
      mix_net_natural(m, src, w, nat.data());
      const int col = inb_id[u] * M + mm;
      for (int j = 0; j < nc; ++j) Xnet(j, col) = nle::to_net_scale(nat[j], mx.net_tf[j], mx.net_tfc[j]);
    }
  }
  std::vector<double> lpk; nle::KnotBlock kb;
  if (Uin > 0) cat_condition(m, Xnet, lpk, kb);
  cat_eval_trials(m, inb_id, M, mx.gh_logw, lpk, kb, uid, rt, t0, st0, k, miss, LC, UC, n, log_pdf, cdf);
}

static int tf_code(const std::string& s) {
  if (s == "identity") return 0;
  if (s == "log") return 1;
  if (s == "probit") return 2;
  if (s == "zerobox") return 3;
  stop("cat_build: unsupported transform '%s' in mix.net_transforms.", s.c_str());
}

static void load_cat_from_list(List fl, CatModel& m) {
  nle::load_spline(fl["spline"], m.sp, "cat_build");
  List scaler = fl["scaler"];
  m.scaler_mean = as<std::vector<double>>(scaler["mean"]);
  m.scaler_scale = as<std::vector<double>>(scaler["scale"]);
  m.n_ctx = (int)m.scaler_mean.size();
  if (!fl.containsElementNamed("n_categories")) stop("cat_build: the card has no n_categories.");
  m.K = as<int>(fl["n_categories"]);
  if (m.K < 2) stop("cat_build: n_categories must be >= 2.");
  if (m.K > 30) stop("cat_build: at most 30 categories.");
  nle::load_mlp(fl["flow_mlp"], m.flow, "cat_build");
  nle::load_mlp(fl["classifier_mlp"], m.clf, "cat_build");
  if (m.flow.n_in() != m.n_ctx + m.K)
    stop("cat_build: flow input dim must be n_params + n_categories (one-hot category).");
  if (m.flow.n_out() != 3 * m.sp.num_bins + 1)
    stop("cat_build: flow output dim must be 3 * num_bins + 1.");
  if (m.clf.n_in() != m.n_ctx) stop("cat_build: classifier input dim must be n_params.");
  if (m.clf.n_out() != m.K) stop("cat_build: classifier must have n_categories outputs (the logits).");
  if (fl.containsElementNamed("mix") && !Rf_isNull(fl["mix"])) {
    List mx = fl["mix"]; MixSpec& s = m.mix; s.on = true;
    s.M = as<int>(mx["nodes"]);
    s.gh_x = as<std::vector<double>>(mx["gh_x"]);
    std::vector<double> w = as<std::vector<double>>(mx["gh_w"]);
    if ((int)s.gh_x.size() != s.M || (int)w.size() != s.M || s.M < 1) stop("cat_build: mix.gh_x / gh_w need `nodes` entries.");
    double sw = 0.0; for (double v : w) sw += v;
    if (std::fabs(sw - 1.0) > 1e-8) stop("cat_build: mix.gh_w must sum to 1.");
    s.gh_logw.resize(s.M); for (int i = 0; i < s.M; ++i) s.gh_logw[i] = std::log(w[i]);
    s.clamp = mx.containsElementNamed("clamp") ? as<double>(mx["clamp"]) : 0.25;
    if (!(as<std::string>(mx["map"]) == "rtcon2_drift")) stop("cat_build: unknown mix.map (only rtcon2_drift).");
    CharacterVector src = fl["context_names"]; s.n_src = src.size();
    auto src_idx = [&](const std::string& nm) -> int {
      for (int j = 0; j < s.n_src; ++j) if (as<std::string>(src[j]) == nm) return j;
      stop("cat_build: mix names source '%s', which is not in context_names.", nm.c_str());
    };
    s.i_mu = src_idx(as<std::string>(mx["mean"])); s.i_s = src_idx(as<std::string>(mx["sd"]));
    CharacterVector mi = mx["map_inputs"];
    if (mi.size() != m.K) stop("cat_build: mix.map_inputs must name C1, the K - 2 gaps and a (K entries).");
    for (int j = 0; j < m.K - 1; ++j) s.i_c.push_back(src_idx(as<std::string>(mi[j])));
    s.i_a = src_idx(as<std::string>(mi[m.K - 1]));
    CharacterVector nn = mx["net_names"];
    if (nn.size() != m.n_ctx) stop("cat_build: mix.net_names must have one entry per network input (the scaler's length).");
    List ntf = mx["net_transforms"]; List nzc = mx.containsElementNamed("net_zerobox_c") ? List(mx["net_zerobox_c"]) : List();
    List direct = mx["direct"]; CharacterVector dnames = direct.names();
    for (int j = 0; j < m.n_ctx; ++j) {
      const std::string nm = as<std::string>(nn[j]);
      s.net_tf.push_back(tf_code(as<std::string>(ntf[nm])));
      s.net_tfc.push_back(s.net_tf.back() == 3 ? as<double>(nzc[nm]) : NAN);
      if (j >= m.K) s.direct_src.push_back(src_idx(as<std::string>(direct[nm])));
    }
    if ((int)s.direct_src.size() != m.n_ctx - m.K) stop("cat_build: mix.direct must map every network input after the K drifts.");
    List bn = fl["bounds_natural"];
    s.src_lower = as<std::vector<double>>(bn["lower"]); s.src_upper = as<std::vector<double>>(bn["upper"]);
    if ((int)s.src_lower.size() != s.n_src || (int)s.src_upper.size() != s.n_src) stop("cat_build: bounds_natural length != context_names.");
  } else {
    List bounds = fl["bounds_sampled"];
    m.lower = as<std::vector<double>>(bounds["lower"]);
    m.upper = as<std::vector<double>>(bounds["upper"]);
    if ((int)m.lower.size() != m.n_ctx || (int)m.upper.size() != m.n_ctx)
      stop("cat_build: bounds_sampled length != n_params.");
  }
}

// [[Rcpp::export]]
SEXP cat_build(List fl) {
  CatModel* m = new CatModel();
  try { load_cat_from_list(fl, *m); } catch (...) { delete m; throw; }
  XPtr<CatModel> ptr(m, true);
  ptr.attr("class") = "cat_flow_model";
  ptr.attr("n_categories") = m->K;
  ptr.attr("mix") = m->mix.on;
  return ptr;
}

// [[Rcpp::export]]
bool cat_ptr_valid(SEXP ptr_) {
  return TYPEOF(ptr_) == EXTPTRSXP && R_ExternalPtrAddr(ptr_) != nullptr;
}

// Trial-wise rows. theta: n x n_ctx (sampled scale), or for a mix card n x n_src
// (natural-scale sources); rt, t0, st0, LC, UC: n (s); k: n in 1..K (NA on
// censored rows); miss: n in {0, 1, 2}.
// [[Rcpp::export]]
List cat_eval_trials_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt, NumericVector t0,
                         NumericVector st0, IntegerVector k, IntegerVector miss,
                         NumericVector LC, NumericVector UC) {
  XPtr<CatModel> ptr(ptr_);
  const CatModel& m = *ptr;
  const int nc = m.mix.on ? m.mix.n_src : m.n_ctx, n = rt.size();
  if (theta.nrow() != n || theta.ncol() != nc) stop("theta must be length(rt) x %d.", nc);
  if (t0.size() != n || st0.size() != n || k.size() != n || miss.size() != n || LC.size() != n || UC.size() != n)
    stop("rt, t0, st0, k, miss, LC and UC must have equal length.");
  std::vector<int> k0(n), ms(n);
  for (int t = 0; t < n; ++t) {
    ms[t] = IntegerVector::is_na(miss[t]) ? 0 : miss[t];
    k0[t] = IntegerVector::is_na(k[t]) ? -1 : k[t] - 1;
    if (ms[t] == 0 && (k0[t] < 0 || k0[t] >= m.K)) stop("k must be in 1..%d on observed rows.", m.K);
  }
  NumericVector log_pdf(n), cdf(n);
  const nle::RowIndex ix = nle::index_rows(theta.begin(), n, nc);
  const int U = ix.U();
  std::vector<double> Theta_u((size_t)U * nc);
  for (int u = 0; u < U; ++u)
    for (int j = 0; j < nc; ++j) Theta_u[(size_t)u * nc + j] = theta(ix.first[u], j);
  if (m.mix.on)
    cat_eval_core_mix(m, Theta_u.data(), U, ix.uid.data(), rt.begin(), t0.begin(), st0.begin(), k0.data(),
                      ms.data(), LC.begin(), UC.begin(), n, log_pdf.begin(), cdf.begin());
  else
    cat_eval_core(m, Theta_u.data(), U, ix.uid.data(), rt.begin(), t0.begin(), st0.begin(), k0.data(),
                  ms.data(), LC.begin(), UC.begin(), n, log_pdf.begin(), cdf.begin());
  return List::create(_["log_pdf"] = log_pdf, _["cdf"] = cdf);
}

// ---------------------------------------------------------------------------
// Native calc_ll branch (nle_native.h)
// ---------------------------------------------------------------------------

const CatModel* nle_cat_handle(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || !Rf_inherits(ptr, "cat_flow_model") || R_ExternalPtrAddr(ptr) == nullptr)
    stop("Neural likelihood: the evaluator is not a live K-category flow.");
  return static_cast<const CatModel*>(R_ExternalPtrAddr(ptr));
}

// the inputs the parameter table must supply per row: the sources of a mix card, else the network inputs
int nle_cat_n_ctx(const CatModel* m) { return m->mix.on ? m->mix.n_src : m->n_ctx; }
int nle_cat_n_categories(const CatModel* m) { return m->K; }

void nle_cat_native(const CatModel* m, const double* th, int n, const int* tf, const double* tfc,
                    const double* rt, const double* t0, const double* st0, const int* k,
                    const int* miss, const double* LC, const double* UC, double* log_pdf) {
  if (m->mix.on) {
    // natural-scale sources, distinct rows once; the registration's transforms are not applied
    // (the mix applies its own to the derived network inputs)
    const int ns = m->mix.n_src;
    const nle::RowIndex ix = nle::index_rows(th, n, ns);
    const int U = ix.U();
    std::vector<double> Src_u((size_t)U * ns);
    for (int u = 0; u < U; ++u)
      for (int j = 0; j < ns; ++j) Src_u[(size_t)u * ns + j] = th[ix.first[u] + (size_t)j * n];
    cat_eval_core_mix(*m, Src_u.data(), U, ix.uid.data(), rt, t0, st0, k, miss, LC, UC, n, log_pdf, nullptr);
    return;
  }
  std::vector<double> Theta_u;
  const nle::RowIndex ix = nle::index_net_rows(th, n, m->n_ctx, tf, tfc, Theta_u);
  cat_eval_core(*m, Theta_u.data(), ix.U(), ix.uid.data(), rt, t0, st0, k, miss, LC, UC, n, log_pdf, nullptr);
}
