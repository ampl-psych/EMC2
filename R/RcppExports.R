# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

sp_new <- function(iter, lambda_varimax, q, p, dim_all_c, all_c, lambda_hat, st, cost_matrix, perm) {
    .Call(`_EMC2_sp_new`, iter, lambda_varimax, q, p, dim_all_c, all_c, lambda_hat, st, cost_matrix, perm)
}

calculate_subject_means <- function(group_designs, params, n_subjects, n_pars) {
    .Call(`_EMC2_calculate_subject_means`, group_designs, params, n_subjects, n_pars)
}

dlba <- function(t, A, b, v, sv, posdrift = TRUE) {
    .Call(`_EMC2_dlba`, t, A, b, v, sv, posdrift)
}

plba <- function(t, A, b, v, sv, posdrift = TRUE) {
    .Call(`_EMC2_plba`, t, A, b, v, sv, posdrift)
}

dWald <- function(t, v, B, A, t0) {
    .Call(`_EMC2_dWald`, t, v, B, A, t0)
}

pWald <- function(t, v, B, A, t0) {
    .Call(`_EMC2_pWald`, t, v, B, A, t0)
}

pEXG <- function(q, mu = 5., sigma = 1., tau = 1., lower_tail = TRUE, log_p = FALSE) {
    .Call(`_EMC2_pEXG`, q, mu, sigma, tau, lower_tail, log_p)
}

dEXG <- function(x, mu = 5., sigma = 1., tau = 1., log_d = FALSE) {
    .Call(`_EMC2_dEXG`, x, mu, sigma, tau, log_d)
}

dEXGrace <- function(dt, mu, sigma, tau) {
    .Call(`_EMC2_dEXGrace`, dt, mu, sigma, tau)
}

stopfn_exg <- function(t, mu, sigma, tau, SSD) {
    .Call(`_EMC2_stopfn_exg`, t, mu, sigma, tau, SSD)
}

pEXG_RDEX <- function(q, mu = 5., sigma = 1., tau = 1., lower_tail = TRUE, log_p = FALSE) {
    .Call(`_EMC2_pEXG_RDEX`, q, mu, sigma, tau, lower_tail, log_p)
}

dEXG_RDEX <- function(x, mu = 5., sigma = 1., tau = 1., log_d = FALSE) {
    .Call(`_EMC2_dEXG_RDEX`, x, mu, sigma, tau, log_d)
}

pigt0_RDEX <- function(t, k = 1., l = 1.) {
    .Call(`_EMC2_pigt0_RDEX`, t, k, l)
}

digt0_RDEX <- function(t, k = 1., l = 1.) {
    .Call(`_EMC2_digt0_RDEX`, t, k, l)
}

pigt_RDEX <- function(t, k = 1, l = 1, a = .1, threshold = 1e-10) {
    .Call(`_EMC2_pigt_RDEX`, t, k, l, a, threshold)
}

digt_RDEX <- function(t, k = 1., l = 1., a = .1, threshold = 1e-10) {
    .Call(`_EMC2_digt_RDEX`, t, k, l, a, threshold)
}

dWald_RDEX <- function(t, v, B, A, t0) {
    .Call(`_EMC2_dWald_RDEX`, t, v, B, A, t0)
}

pWald_RDEX <- function(t, v, B, A, t0) {
    .Call(`_EMC2_pWald_RDEX`, t, v, B, A, t0)
}

dRDEXrace <- function(dt, mu, sigma, tau, v, B, A, t0, exgWinner = TRUE) {
    .Call(`_EMC2_dRDEXrace`, dt, mu, sigma, tau, v, B, A, t0, exgWinner)
}

stopfn_rdex <- function(t, n_acc, mu, sigma, tau, v, B, A, t0, SSD) {
    .Call(`_EMC2_stopfn_rdex`, t, n_acc, mu, sigma, tau, v, B, A, t0, SSD)
}

fft_convolve_equiv_cpp <- function(x, y, conj_flag = TRUE) {
    .Call(`_EMC2_fft_convolve_equiv_cpp`, x, y, conj_flag)
}

compute_gamma_diff_hrf <- function(tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio) {
    .Call(`_EMC2_compute_gamma_diff_hrf`, tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio)
}

compute_hrf <- function(tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio) {
    .Call(`_EMC2_compute_hrf`, tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio)
}

compute_time_derivative <- function(tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio, delta = 0.1) {
    .Call(`_EMC2_compute_time_derivative`, tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio, delta)
}

build_hrf_kernel <- function(has_derivative, tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio) {
    .Call(`_EMC2_build_hrf_kernel`, has_derivative, tr, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio)
}

construct_design_matrix <- function(frame_times, events, has_derivative, min_onset, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio, add_intercept) {
    .Call(`_EMC2_construct_design_matrix`, frame_times, events, has_derivative, min_onset, oversampling, time_length, onset, delay, undershoot, dispersion, u_dispersion, ratio, add_intercept)
}

calc_ll <- function(p_matrix, data, constants, designs, type, bounds, transforms, pretransforms, p_types, min_ll, trend) {
    .Call(`_EMC2_calc_ll`, p_matrix, data, constants, designs, type, bounds, transforms, pretransforms, p_types, min_ll, trend)
}

c_add_charvectors <- function(x, y) {
    .Call(`_EMC2_c_add_charvectors`, x, y)
}

