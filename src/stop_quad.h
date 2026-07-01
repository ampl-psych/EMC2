#ifndef stop_quad_h
#define stop_quad_h

// ---------------------------------------------------------------------------
// cens replacement for Zach's gl_quad.h.
//
// Zach's branch integrated the stop-success integral with a self-contained
// fixed Gauss-Legendre rule (gl_quad.h) for the "gl" route and a bundled GSL
// (qags/qagiu) for the "integrate" route. The cens branch already ships its own
// adaptive cubature (src/hcubature.cpp) exposed through src/gauss.h, with a 1D
// Gauss-Kronrod panel (gauss_kronrod, the xd7/wd7/gwd7 7/15 rule) and an
// adaptive driver (hcubature). Per the SS-port brief we DROP both GL and GSL and
// re-back every stop-success route on cens's own integrator, while keeping the
// same "fixed-point" strategy: a deterministic finite window plus a fixed,
// controllable number of quadrature panels.
//
// Routes (same option semantics as Zach):
//   * "gl"        -> ss_fixed_integrate: deterministic composite Gauss-Kronrod
//                    over [lo, ub], n_nodes panels. Fast, fixed cost. (The label
//                    "gl" is kept for API/test back-compat; the engine is G-K.)
//   * "integrate" -> ss_adaptive_integrate: cens hcubature over [lo, ub],
//                    tolerance-driven. Robust fallback.
//   * "analytic"  -> ss_exg_analytic.h closed form (unchanged; not here).
//   * "auto"      -> analytic where it applies, else "gl" with an auto
//                    panel-density bump (gl_auto_nodes).
//
// The finite-window heuristic (the "fixed point strategy" to preserve) is
// [muS - K_SIGMA*sigS - K_TAU*tauS , muS + K_SIGMA*sigS + K_TAU*tauS] clipped to
// the stop lower bound; the model headers build lo/ub from these constants and
// the R routes (stop_success_gl.R) mirror them. Keep all three in sync.
// ---------------------------------------------------------------------------

#include <cmath>
#include <cstddef>
#include "gauss.h"   // gauss_kronrod, hcubature, struct one_d

constexpr double SS_WINDOW_K_SIGMA = 8.0;
constexpr double SS_WINDOW_K_TAU   = 16.0;

// Integrand signature: cens's hcubature/gauss_kronrod callback form.
//   int f(unsigned dim, const double* x, void* params, unsigned fdim, double* retval)
// (single evaluation point in x[0]; writes the integrand value to retval[0]).
using ss_integrand_fn = int (*)(unsigned, const double*, void*, unsigned, double*);

// ---------------------------------------------------------------------------
// Process-global stop_method configuration for the LIVE stop-success path.
// Set from R via emc2_set_stop_method() (see model_SS_EXG.h). SSEXG()/SSRDEX()
// store stop_method/stop_n_nodes in the model list and calc_ll_manager() calls
// the setter once per likelihood call before entering C++. Process-global (NOT
// thread_local): set once per call, read many times; forked workers inherit it,
// PSOCK workers each run calc_ll_manager() themselves.
// NB: emc2_get_stop_method() indexes its method_names array by these values —
// keep the two in sync if methods are added or reordered.
// ---------------------------------------------------------------------------
enum StopMethod {
  STOP_METHOD_AUTO      = 0,   // DEFAULT
  STOP_METHOD_INTEGRATE = 1,   // adaptive hcubature over the window
  STOP_METHOD_GL        = 2,   // fixed composite Gauss-Kronrod, n_nodes panels
  STOP_METHOD_ANALYTIC  = 3    // EXG n_go==1 closed form (GL fallback otherwise)
};

struct StopMethodConfig {
  int method  = STOP_METHOD_AUTO;
  int n_nodes = 64;
};

inline StopMethodConfig& stop_method_config() {
  static StopMethodConfig cfg;
  return cfg;
}

// "auto" panel-density bump (mirror of Zach's gl_auto_nodes). A fixed panel
// count under-resolves a sharp (near-Gaussian, small-tauS) stop density when the
// window spans many stop-sigma widths; the Gaussian core width ~sigS does NOT
// shrink as tauS -> 0, so hold a minimum node density per sigS across the window.
//   want  = GL_NODES_PER_SIG * (ub - lo) / sigS
//   n_eff = clamp(roundup(want, GL_NODE_STEP), n_nodes, GL_MAX_NODES)
// n_nodes is a floor (never reduced); quantising up to a multiple keeps the
// distinct panel counts few. Mirrored in R by gl_auto_nodes_R().
constexpr double GL_NODES_PER_SIG = 6.0;
constexpr int    GL_NODE_STEP     = 32;
constexpr int    GL_MAX_NODES     = 256;

inline int gl_auto_nodes(int n_nodes, double lo, double ub, double sigS) {
  if (!(sigS > 0.0) || !(ub > lo)) return n_nodes;
  const double want = GL_NODES_PER_SIG * (ub - lo) / sigS;
  if (!(want > n_nodes)) return n_nodes;
  int n_eff = static_cast<int>(std::ceil(want / GL_NODE_STEP)) * GL_NODE_STEP;
  if (n_eff < n_nodes)      n_eff = n_nodes;
  if (n_eff > GL_MAX_NODES) n_eff = GL_MAX_NODES;
  return n_eff;
}

// Each cens gauss_kronrod panel spends a 7/15-point Gauss-Kronrod rule (15
// integrand evals). Map Zach's GL "node count" to a comparable number of G-K
// panels so n_nodes keeps tracking cost/resolution and the convergence tests
// (more nodes -> more accuracy) still hold.
inline int ss_panels_for_nodes(int n_nodes) {
  int p = (n_nodes + 14) / 15;
  return p < 1 ? 1 : p;
}

// --- "gl" route: deterministic composite Gauss-Kronrod over [lo, ub] ---------
inline double ss_fixed_integrate(ss_integrand_fn f, void* params,
                                 double lo, double ub, int n_nodes) {
  if (!(ub > lo)) return 0.0;
  const int n_panels = ss_panels_for_nodes(n_nodes);
  const double h = (ub - lo) / n_panels;
  double total = 0.0;
  for (int i = 0; i < n_panels; ++i) {
    one_d out;
    const double a = lo + i * h;
    const double b = (i == n_panels - 1) ? ub : (a + h);
    gauss_kronrod(a, b, out, params, f);
    total += out.result;
  }
  return total;
}

// --- "integrate" route: cens adaptive hcubature over [lo, ub] -----------------
inline double ss_adaptive_integrate(ss_integrand_fn f, void* params,
                                    double lo, double ub,
                                    double abs_tol, double rel_tol,
                                    std::size_t max_eval) {
  if (!(ub > lo)) return 0.0;
  double a = lo, b = ub, val = 0.0, err = 0.0;
  hcubature(f, params, 1, &a, &b, max_eval, abs_tol, rel_tol, &val, &err);
  return val;
}

#endif // stop_quad_h
