// Plain-MLP likelihood networks: a single MLP whose one output is the joint
// log density log p(rt, R | theta) of a two-response model. Two card kinds
// share this evaluator (R/nn_register.R):
//   "regression_joint"  direct-regression nets of the NLE project (GELU(tanh)
//                       hidden layers, standardised inputs, -Inf outside the box)
//   "mlp_joint"         likelihood approximation networks (LANs, e.g. HSSM's
//                       tanh MLPs converted by inst/scripts/onnx_to_card.py):
//                       raw inputs, out-of-box rows return the card's log floor
//
// The weights are loaded once (mlp_lik_build -> external pointer). Each input
// of the network is a parameter (entered on the sampled scale), rt, log rt or
// the response; `input_layout` says which, in the network's order. Rows whose
// parameters leave the training box are not evaluated (they get `oob`).
// Evaluation is batched: NLE_BLOCK rows per matrix product, no state, so the
// core is reentrant (native calc_ll branch, nle_native.h / model_NN.h).
//
//   mlp_lik_build(card, lower_s, upper_s, oob) -> external pointer
//   mlp_lik_eval_cpp(ptr, theta, rt, R)        -> log density per row (theta:
//                                                 n x n_ctx, sampled scale)

#include "nle_flow.h"
#include "nle_native.h"
using namespace Rcpp;

enum { IN_RT = -1, IN_LOG_RT = -2, IN_R = -3 };

struct MlpLik {
  nle::Mlp mlp;
  int n_ctx = 0;
  std::vector<int> src;                    // per network input: context index, or IN_*
  std::vector<double> in_mean, in_scale;   // empty = no input standardisation
  double resp[2] = {1.0, 2.0};             // network input coding responses 1, 2
  bool out_scaled = false;
  double out_mean = 0.0, out_scale = 1.0;  // output = raw * scale + mean
  std::vector<double> lower, upper;        // training box, sampled scale
  double oob = R_NegInf;                   // value for rows outside the box
};

// Th: n x n_ctx column-major on the sampled scale. Rows with NaN or out-of-box
// parameters, or a time <= 0 where the network takes one, get m.oob.
static void mlp_lik_core(const MlpLik& m, const double* Th, int n, const double* tn,
                         const int* R, double* out) {
  const int n_in = m.mlp.n_in();
  std::vector<int> rows;
  std::vector<double> th(m.n_ctx);
  rows.reserve(n);
  for (int t = 0; t < n; ++t) {
    for (int j = 0; j < m.n_ctx; ++j) th[j] = Th[t + (size_t)j * n];
    if (nle::in_box(m.lower, m.upper, th.data()) && tn[t] > 0.0) rows.push_back(t);
    else out[t] = m.oob;
  }
  const int B = nle::NLE_BLOCK;
  arma::mat X;
  for (size_t s = 0; s < rows.size(); s += B) {
    const int nb = (int)std::min<size_t>(B, rows.size() - s);
    X.set_size(n_in, nb);
    for (int c = 0; c < nb; ++c) {
      const int t = rows[s + c];
      double* x = X.colptr(c);
      for (int i = 0; i < n_in; ++i) {
        const int q = m.src[i];
        double v;
        if (q >= 0) v = Th[t + (size_t)q * n];
        else if (q == IN_RT) v = tn[t];
        else if (q == IN_LOG_RT) v = std::log(tn[t]);
        else v = m.resp[R[t] - 1];
        x[i] = m.in_mean.empty() ? v : (v - m.in_mean[i]) / m.in_scale[i];
      }
    }
    const arma::mat Y = nle::mlp_forward_batch(m.mlp, X);
    for (int c = 0; c < nb; ++c) {
      const double y = Y(0, c);
      out[rows[s + c]] = m.out_scaled ? y * m.out_scale + m.out_mean : y;
    }
  }
}

// TRUE if the external pointer still holds a live network.
// [[Rcpp::export]]
bool mlp_lik_valid(SEXP ptr_) {
  return TYPEOF(ptr_) == EXTPTRSXP && Rf_inherits(ptr_, "mlp_lik") &&
         R_ExternalPtrAddr(ptr_) != nullptr;
}

// card: a normalised card member (`mlp`, `context_names`, `input_layout`,
// optional `scaler`, `response_values`, `output_scaler`). lower_s/upper_s: the
// training box on the sampled scale in context_names order. oob: the value of
// rows outside it (-Inf, or the card's log floor).
// [[Rcpp::export]]
SEXP mlp_lik_build(List card, NumericVector lower_s, NumericVector upper_s, double oob) {
  MlpLik* m = new MlpLik();
  try {
    nle::load_mlp(card["mlp"], m->mlp, "mlp_lik_build", true);
    if (m->mlp.n_out() != 1)
      stop("mlp_lik_build: the network must have a single output (the joint log density).");
    CharacterVector ctx = card["context_names"];
    CharacterVector lay = card["input_layout"];
    m->n_ctx = ctx.size();
    if (lay.size() != m->mlp.n_in())
      stop("mlp_lik_build: input_layout has %d entries; the network has %d inputs.",
           (int)lay.size(), m->mlp.n_in());
    if (lower_s.size() != m->n_ctx || upper_s.size() != m->n_ctx)
      stop("mlp_lik_build: the box needs one bound per context name (%d).", m->n_ctx);
    for (int i = 0; i < lay.size(); ++i) {
      const std::string nm = as<std::string>(lay[i]);
      int q = IN_R;
      if (nm == "rt") q = IN_RT;
      else if (nm == "log_rt") q = IN_LOG_RT;
      else if (nm != "R") {
        q = -100;
        for (int j = 0; j < ctx.size(); ++j) if (as<std::string>(ctx[j]) == nm) q = j;
        if (q == -100) stop("mlp_lik_build: input_layout entry '%s' is not a context name.", nm.c_str());
      }
      m->src.push_back(q);
    }
    if (card.containsElementNamed("scaler") && !Rf_isNull(card["scaler"])) {
      List sc = card["scaler"];
      m->in_mean = as<std::vector<double>>(sc["mean"]);
      m->in_scale = as<std::vector<double>>(sc["scale"]);
      if ((int)m->in_mean.size() != m->mlp.n_in() || (int)m->in_scale.size() != m->mlp.n_in())
        stop("mlp_lik_build: scaler needs one mean and scale per network input.");
    }
    if (card.containsElementNamed("response_values") && !Rf_isNull(card["response_values"])) {
      NumericVector rv = card["response_values"];
      if (rv.size() != 2) stop("mlp_lik_build: response_values must have two entries.");
      m->resp[0] = rv[0]; m->resp[1] = rv[1];
    }
    if (card.containsElementNamed("output_scaler") && !Rf_isNull(card["output_scaler"])) {
      List os = card["output_scaler"];
      m->out_scaled = true;
      m->out_mean = as<double>(os["mean"]);
      m->out_scale = as<double>(os["scale"]);
    }
    m->lower.assign(lower_s.begin(), lower_s.end());
    m->upper.assign(upper_s.begin(), upper_s.end());
    m->oob = oob;
  } catch (...) { delete m; throw; }
  XPtr<MlpLik> ptr(m, true);
  ptr.attr("class") = "mlp_lik";
  return ptr;
}

// Joint log density per row. theta: n x n_ctx on the sampled scale; rt > 0;
// R in {1, 2}.
// [[Rcpp::export]]
NumericVector mlp_lik_eval_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt, IntegerVector R) {
  if (!mlp_lik_valid(ptr_)) stop("The neural-likelihood evaluator is not live; rebuild it.");
  XPtr<MlpLik> ptr(ptr_);
  const int n = rt.size();
  if (theta.nrow() != n || theta.ncol() != ptr->n_ctx) stop("theta must be length(rt) x %d.", ptr->n_ctx);
  if (R.size() != n) stop("rt and R must have equal length.");
  for (int t = 0; t < n; ++t)
    if (R[t] != 1 && R[t] != 2) stop("R must be 1 or 2.");
  NumericVector out(n);
  mlp_lik_core(*ptr, theta.begin(), n, rt.begin(), R.begin(), out.begin());
  return out;
}

// ---------------------------------------------------------------------------
// Native calc_ll branch (nle_native.h)
// ---------------------------------------------------------------------------

const MlpLik* nle_mlp_handle(SEXP ptr) {
  if (!mlp_lik_valid(ptr)) stop("Neural likelihood: the evaluator is not a live MLP likelihood.");
  return static_cast<const MlpLik*>(R_ExternalPtrAddr(ptr));
}

int nle_mlp_n_ctx(const MlpLik* m) { return m->n_ctx; }

void nle_mlp_native(const MlpLik* m, const double* th, int n, const int* tf,
                    const double* tn, const int* R, double* log_pdf) {
  std::vector<double> Th((size_t)n * m->n_ctx);
  for (int j = 0; j < m->n_ctx; ++j)
    for (int t = 0; t < n; ++t)
      Th[t + (size_t)j * n] = nle::to_net_scale(th[t + (size_t)j * n], tf[j]);
  mlp_lik_core(*m, Th.data(), n, tn, R, log_pdf);
}
