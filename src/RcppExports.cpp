// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dlba
NumericVector dlba(NumericVector t, NumericVector A, NumericVector b, NumericVector v, NumericVector sv, bool posdrift);
RcppExport SEXP _EMC2_dlba(SEXP tSEXP, SEXP ASEXP, SEXP bSEXP, SEXP vSEXP, SEXP svSEXP, SEXP posdriftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sv(svSEXP);
    Rcpp::traits::input_parameter< bool >::type posdrift(posdriftSEXP);
    rcpp_result_gen = Rcpp::wrap(dlba(t, A, b, v, sv, posdrift));
    return rcpp_result_gen;
END_RCPP
}
// plba
NumericVector plba(NumericVector t, NumericVector A, NumericVector b, NumericVector v, NumericVector sv, bool posdrift);
RcppExport SEXP _EMC2_plba(SEXP tSEXP, SEXP ASEXP, SEXP bSEXP, SEXP vSEXP, SEXP svSEXP, SEXP posdriftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sv(svSEXP);
    Rcpp::traits::input_parameter< bool >::type posdrift(posdriftSEXP);
    rcpp_result_gen = Rcpp::wrap(plba(t, A, b, v, sv, posdrift));
    return rcpp_result_gen;
END_RCPP
}
// dWald
NumericVector dWald(NumericVector t, NumericVector v, NumericVector B, NumericVector A, NumericVector t0);
RcppExport SEXP _EMC2_dWald(SEXP tSEXP, SEXP vSEXP, SEXP BSEXP, SEXP ASEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(dWald(t, v, B, A, t0));
    return rcpp_result_gen;
END_RCPP
}
// pWald
NumericVector pWald(NumericVector t, NumericVector v, NumericVector B, NumericVector A, NumericVector t0);
RcppExport SEXP _EMC2_pWald(SEXP tSEXP, SEXP vSEXP, SEXP BSEXP, SEXP ASEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(pWald(t, v, B, A, t0));
    return rcpp_result_gen;
END_RCPP
}
// fft_convolve_equiv_cpp
arma::vec fft_convolve_equiv_cpp(const arma::vec& x, const arma::vec& y, bool conj_flag);
RcppExport SEXP _EMC2_fft_convolve_equiv_cpp(SEXP xSEXP, SEXP ySEXP, SEXP conj_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type conj_flag(conj_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(fft_convolve_equiv_cpp(x, y, conj_flag));
    return rcpp_result_gen;
END_RCPP
}
// compute_gamma_diff_hrf
NumericVector compute_gamma_diff_hrf(double tr, int oversampling, double time_length, double onset, double delay, double undershoot, double dispersion, double u_dispersion, double ratio);
RcppExport SEXP _EMC2_compute_gamma_diff_hrf(SEXP trSEXP, SEXP oversamplingSEXP, SEXP time_lengthSEXP, SEXP onsetSEXP, SEXP delaySEXP, SEXP undershootSEXP, SEXP dispersionSEXP, SEXP u_dispersionSEXP, SEXP ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< int >::type oversampling(oversamplingSEXP);
    Rcpp::traits::input_parameter< double >::type time_length(time_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type onset(onsetSEXP);
    Rcpp::traits::input_parameter< double >::type delay(delaySEXP);
    Rcpp::traits::input_parameter< double >::type undershoot(undershootSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type u_dispersion(u_dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_gamma_diff_hrf(tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio));
    return rcpp_result_gen;
END_RCPP
}
// compute_hrf
NumericVector compute_hrf(double tr, int oversampling, double time_length, double onset, double delay, double undershoot, double dispersion, double u_dispersion, double ratio);
RcppExport SEXP _EMC2_compute_hrf(SEXP trSEXP, SEXP oversamplingSEXP, SEXP time_lengthSEXP, SEXP onsetSEXP, SEXP delaySEXP, SEXP undershootSEXP, SEXP dispersionSEXP, SEXP u_dispersionSEXP, SEXP ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< int >::type oversampling(oversamplingSEXP);
    Rcpp::traits::input_parameter< double >::type time_length(time_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type onset(onsetSEXP);
    Rcpp::traits::input_parameter< double >::type delay(delaySEXP);
    Rcpp::traits::input_parameter< double >::type undershoot(undershootSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type u_dispersion(u_dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_hrf(tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio));
    return rcpp_result_gen;
END_RCPP
}
// compute_time_derivative
NumericVector compute_time_derivative(double tr, int oversampling, double time_length, double onset, double delay, double undershoot, double dispersion, double u_dispersion, double ratio, double delta);
RcppExport SEXP _EMC2_compute_time_derivative(SEXP trSEXP, SEXP oversamplingSEXP, SEXP time_lengthSEXP, SEXP onsetSEXP, SEXP delaySEXP, SEXP undershootSEXP, SEXP dispersionSEXP, SEXP u_dispersionSEXP, SEXP ratioSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< int >::type oversampling(oversamplingSEXP);
    Rcpp::traits::input_parameter< double >::type time_length(time_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type onset(onsetSEXP);
    Rcpp::traits::input_parameter< double >::type delay(delaySEXP);
    Rcpp::traits::input_parameter< double >::type undershoot(undershootSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type u_dispersion(u_dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_time_derivative(tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio, delta));
    return rcpp_result_gen;
END_RCPP
}
// build_hrf_kernel
NumericMatrix build_hrf_kernel(bool has_derivative, double tr, int oversampling, double time_length, double onset, double delay, double undershoot, double dispersion, double u_dispersion, double ratio);
RcppExport SEXP _EMC2_build_hrf_kernel(SEXP has_derivativeSEXP, SEXP trSEXP, SEXP oversamplingSEXP, SEXP time_lengthSEXP, SEXP onsetSEXP, SEXP delaySEXP, SEXP undershootSEXP, SEXP dispersionSEXP, SEXP u_dispersionSEXP, SEXP ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type has_derivative(has_derivativeSEXP);
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< int >::type oversampling(oversamplingSEXP);
    Rcpp::traits::input_parameter< double >::type time_length(time_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type onset(onsetSEXP);
    Rcpp::traits::input_parameter< double >::type delay(delaySEXP);
    Rcpp::traits::input_parameter< double >::type undershoot(undershootSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type u_dispersion(u_dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(build_hrf_kernel(has_derivative, tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio));
    return rcpp_result_gen;
END_RCPP
}
// construct_design_matrix
DataFrame construct_design_matrix(NumericVector frame_times, DataFrame events, bool has_derivative, double min_onset, int oversampling, double time_length, double onset, double delay, double undershoot, double dispersion, double u_dispersion, double ratio, bool add_intercept);
RcppExport SEXP _EMC2_construct_design_matrix(SEXP frame_timesSEXP, SEXP eventsSEXP, SEXP has_derivativeSEXP, SEXP min_onsetSEXP, SEXP oversamplingSEXP, SEXP time_lengthSEXP, SEXP onsetSEXP, SEXP delaySEXP, SEXP undershootSEXP, SEXP dispersionSEXP, SEXP u_dispersionSEXP, SEXP ratioSEXP, SEXP add_interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type frame_times(frame_timesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< bool >::type has_derivative(has_derivativeSEXP);
    Rcpp::traits::input_parameter< double >::type min_onset(min_onsetSEXP);
    Rcpp::traits::input_parameter< int >::type oversampling(oversamplingSEXP);
    Rcpp::traits::input_parameter< double >::type time_length(time_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type onset(onsetSEXP);
    Rcpp::traits::input_parameter< double >::type delay(delaySEXP);
    Rcpp::traits::input_parameter< double >::type undershoot(undershootSEXP);
    Rcpp::traits::input_parameter< double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type u_dispersion(u_dispersionSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    Rcpp::traits::input_parameter< bool >::type add_intercept(add_interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_design_matrix(frame_times, events, has_derivative, min_onset, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio, add_intercept));
    return rcpp_result_gen;
END_RCPP
}
// calc_ll
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants, List designs, String type, List bounds, List transforms, List pretransforms, CharacterVector p_types, double min_ll, List trend);
RcppExport SEXP _EMC2_calc_ll(SEXP p_matrixSEXP, SEXP dataSEXP, SEXP constantsSEXP, SEXP designsSEXP, SEXP typeSEXP, SEXP boundsSEXP, SEXP transformsSEXP, SEXP pretransformsSEXP, SEXP p_typesSEXP, SEXP min_llSEXP, SEXP trendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p_matrix(p_matrixSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type constants(constantsSEXP);
    Rcpp::traits::input_parameter< List >::type designs(designsSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< List >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< List >::type transforms(transformsSEXP);
    Rcpp::traits::input_parameter< List >::type pretransforms(pretransformsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type p_types(p_typesSEXP);
    Rcpp::traits::input_parameter< double >::type min_ll(min_llSEXP);
    Rcpp::traits::input_parameter< List >::type trend(trendSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_ll(p_matrix, data, constants, designs, type, bounds, transforms, pretransforms, p_types, min_ll, trend));
    return rcpp_result_gen;
END_RCPP
}
// c_add_charvectors
CharacterVector c_add_charvectors(CharacterVector x, CharacterVector y);
RcppExport SEXP _EMC2_c_add_charvectors(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_add_charvectors(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EMC2_dlba", (DL_FUNC) &_EMC2_dlba, 6},
    {"_EMC2_plba", (DL_FUNC) &_EMC2_plba, 6},
    {"_EMC2_dWald", (DL_FUNC) &_EMC2_dWald, 5},
    {"_EMC2_pWald", (DL_FUNC) &_EMC2_pWald, 5},
    {"_EMC2_fft_convolve_equiv_cpp", (DL_FUNC) &_EMC2_fft_convolve_equiv_cpp, 3},
    {"_EMC2_compute_gamma_diff_hrf", (DL_FUNC) &_EMC2_compute_gamma_diff_hrf, 9},
    {"_EMC2_compute_hrf", (DL_FUNC) &_EMC2_compute_hrf, 9},
    {"_EMC2_compute_time_derivative", (DL_FUNC) &_EMC2_compute_time_derivative, 10},
    {"_EMC2_build_hrf_kernel", (DL_FUNC) &_EMC2_build_hrf_kernel, 10},
    {"_EMC2_construct_design_matrix", (DL_FUNC) &_EMC2_construct_design_matrix, 13},
    {"_EMC2_calc_ll", (DL_FUNC) &_EMC2_calc_ll, 11},
    {"_EMC2_c_add_charvectors", (DL_FUNC) &_EMC2_c_add_charvectors, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_EMC2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
