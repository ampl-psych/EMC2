#ifndef __APPLE__
// Generic implementation of math_utils.h.
// Simple loops that auto-vectorise via libmvec (_ZGVdN4v_exp etc.)
// under -O3 -march=native on Linux/GCC.

#include <cmath>
#include "math_utils.h"

// --- exp ---

void vec_exp(double* x, int n)
{
  for (int i = 0; i < n; ++i)
    x[i] = std::exp(x[i]);
}

void vec_exp(double* dst, const double* src, int n)
{
  for (int i = 0; i < n; ++i)
    dst[i] = std::exp(src[i]);
}

void vec_exp_offset(double* x, int n, double offset)
{
  for (int i = 0; i < n; ++i)
    x[i] = offset + std::exp(x[i]);
}

// --- log ---

void vec_log(double* x, int n)
{
  for (int i = 0; i < n; ++i)
    x[i] = std::log(x[i]);
}

void vec_log(double* dst, const double* src, int n)
{
  for (int i = 0; i < n; ++i)
    dst[i] = std::log(src[i]);
}

// --- sqrt ---

void vec_sqrt(double* x, int n)
{
  for (int i = 0; i < n; ++i)
    x[i] = std::sqrt(x[i]);
}

void vec_sqrt(double* dst, const double* src, int n)
{
  for (int i = 0; i < n; ++i)
    dst[i] = std::sqrt(src[i]);
}

#endif
