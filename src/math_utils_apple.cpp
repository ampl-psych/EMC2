#ifdef __APPLE__
// Apple/Accelerate implementation of math_utils.h.
// vvexp / vvlog / vvsqrt are hand-tuned for NEON/AMX on Apple Silicon.
//
// Accelerate must be included before any R/Rcpp headers to avoid the
// COMPLEX typedef collision between vecLib (DSPComplex) and R's COMPLEX() macro.
#include <Accelerate/Accelerate.h>
#undef COMPLEX  // vecLib defines COMPLEX as DSPComplex; R's macro takes priority

#include "math_utils.h"

// --- exp ---

void vec_exp(double* x, int n)
{
  vvexp(x, x, &n);
}

void vec_exp(double* dst, const double* src, int n)
{
  vvexp(dst, src, &n);
}

void vec_exp_offset(double* x, int n, double offset)
{
  vvexp(x, x, &n);
  if (offset != 0.0)
    vDSP_vsaddD(x, 1, &offset, x, 1, (vDSP_Length)n);
}

// --- log ---

void vec_log(double* x, int n)
{
  vvlog(x, x, &n);
}

void vec_log(double* dst, const double* src, int n)
{
  vvlog(dst, src, &n);
}

// --- sqrt ---

void vec_sqrt(double* x, int n)
{
  vvsqrt(x, x, &n);
}

void vec_sqrt(double* dst, const double* src, int n)
{
  vvsqrt(dst, src, &n);
}

#endif // __APPLE__
