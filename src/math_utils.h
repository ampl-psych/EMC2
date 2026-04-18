#pragma once

// Platform-agnostic vectorised math operations.
// Implementations are in math_utils_apple.cpp (Apple/Accelerate)
// and math_utils_generic.cpp (Linux/auto-vectorisation).

// --- exp ---

// In-place: x[i] = exp(x[i])
void vec_exp(double* x, int n);

// Out-of-place: dst[i] = exp(src[i])
void vec_exp(double* dst, const double* src, int n);

// In-place with additive offset: x[i] = offset + exp(x[i])
void vec_exp_offset(double* x, int n, double offset);

// --- log ---

// In-place: x[i] = log(x[i])
void vec_log(double* x, int n);

// Out-of-place: dst[i] = log(src[i])
void vec_log(double* dst, const double* src, int n);

// --- sqrt ---

// In-place: x[i] = sqrt(x[i])
void vec_sqrt(double* x, int n);

// Out-of-place: dst[i] = sqrt(src[i])
void vec_sqrt(double* dst, const double* src, int n);
