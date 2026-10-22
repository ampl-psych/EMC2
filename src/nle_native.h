// Entry points of the neural-likelihood evaluators (flow_ddm.cpp,
// flow_race.cpp, mlp_lik.cpp, flow_cat.cpp) for the native calc_ll branch (model_NN.h, particle_ll.cpp).
// Kept free of Armadillo so that translation units built on Rcpp.h can
// include it.
//
// Resolve an evaluator on the main thread (nle_*_handle stop()s if the
// external pointer is not a live evaluator of that kind); the *_native
// functions are reentrant (no R API calls) and may run inside OpenMP threads.

#ifndef EMC2_NLE_NATIVE_H
#define EMC2_NLE_NATIVE_H

#include <Rinternals.h>

struct DdmEnsemble;
struct FlowModel;
struct MlpLik;
struct CatModel;

const DdmEnsemble* nle_ddm_handle(SEXP ptr);
const FlowModel*   nle_flow_handle(SEXP ptr);
const MlpLik*      nle_mlp_handle(SEXP ptr);
const CatModel*    nle_cat_handle(SEXP ptr);
int nle_ddm_n_ctx(const DdmEnsemble* e);
int nle_flow_n_ctx(const FlowModel* f);
int nle_mlp_n_ctx(const MlpLik* m);
int nle_cat_n_ctx(const CatModel* m);
int nle_cat_n_categories(const CatModel* m);

// th: m x n_ctx natural-scale network inputs (column-major, the network's
// order); tf: per input 0 identity, 1 log, 2 probit, 3 zerobox; tfc: per input
// the zerobox constant c (asinh(x / c); NaN for the other codes); tn: the network's time
// (> 0). Distinct rows are conditioned once. Out-of-box rows: log pdf -Inf,
// log survivor 0.

// Joint DDM flow: log p(tn, R | theta), R in {1, 2}.
void nle_ddm_native(const DdmEnsemble* e, const double* th, int m, const int* tf,
                    const double* tfc, const double* tn, const int* R, double* log_pdf);

// The same, with the defective CDF P(time <= tn, R | theta) as well (scaled DDM cards).
void nle_ddm_native_full(const DdmEnsemble* e, const double* th, int m, const int* tf,
                         const double* tfc, const double* tn, const int* R, double* log_pdf,
                         double* cdf);

// Plain-MLP likelihood (mlp_lik.cpp; regression_joint and mlp_joint cards):
// log p(tn, R | theta) per row, R in {1, 2}. Rows outside the training box
// get the network's out-of-box value (-Inf, or the card's log floor).
void nle_mlp_native(const MlpLik* m, const double* th, int n, const int* tf,
                    const double* tfc, const double* tn, const int* R, double* log_pdf);

// Race flow: per row, want = 1 -> log pdf into log_pdf, want = 2 -> log
// survivor into log_sf (the other output is left untouched).
void nle_flow_native(const FlowModel* f, const double* th, int m, const int* tf,
                     const double* tfc, const double* tn, const unsigned char* want,
                     double* log_pdf, double* log_sf);

// K-category joint flow with the uniform non-decision convolution and censoring
// (flow_cat.cpp): per row rt, t0, st0 (s), k in 0..K-1 (observed rows), miss
// 0 observed / 1 lower-censored at LC / 2 upper-censored at UC. th: the network
// inputs as above, or for a card with a `mix` block the natural-scale sources
// (the transform codes are then not applied).
void nle_cat_native(const CatModel* m, const double* th, int n, const int* tf, const double* tfc,
                    const double* rt, const double* t0, const double* st0, const int* k,
                    const int* miss, const double* LC, const double* UC, double* log_pdf);

#endif
