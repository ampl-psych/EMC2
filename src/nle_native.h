// Entry points of the neural-likelihood evaluators (flow_ddm.cpp,
// flow_race.cpp) for the native calc_ll branch (model_NN.h, particle_ll.cpp).
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

const DdmEnsemble* nle_ddm_handle(SEXP ptr);
const FlowModel*   nle_flow_handle(SEXP ptr);
int nle_ddm_n_ctx(const DdmEnsemble* e);
int nle_flow_n_ctx(const FlowModel* f);

// th: m x n_ctx natural-scale network inputs (column-major, the network's
// order); tf: per input 0 identity, 1 log, 2 probit; tn: the network's time
// (> 0). Distinct rows are conditioned once. Out-of-box rows: log pdf -Inf,
// log survivor 0.

// Joint DDM flow: log p(tn, R | theta), R in {1, 2}.
void nle_ddm_native(const DdmEnsemble* e, const double* th, int m, const int* tf,
                    const double* tn, const int* R, double* log_pdf);

// Race flow: per row, want = 1 -> log pdf into log_pdf, want = 2 -> log
// survivor into log_sf (the other output is left untouched).
void nle_flow_native(const FlowModel* f, const double* th, int m, const int* tf,
                     const double* tn, const unsigned char* want,
                     double* log_pdf, double* log_sf);

#endif
