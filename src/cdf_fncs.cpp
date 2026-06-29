
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "cdf_fncs.h"
#include "tools.h"
#include "ddm_functions_inline.h"

/* DISTRIBUTION */

/* P term in the distribution function */
double logP(int pm, double a, double v, double w) {
    return ddm_logP(pm, a, v, w);
}

/* calculate number of terms needed for short t */
double Ks(double t, double v, double a, double w, double eps)
{
    return ddm_Ks(t, v, a, w, eps);
}

/* calculate number of terms needed for large t */
double Kl(double t, double v, double a, double w, double err) {
    return ddm_Kl(t, v, a, w, err);
}

/* calculate terms of the sum for short t */
double logFs(double t, double v, double a, double w, int K)
{
    return ddm_logFs_faithful(t, v, a, w, K);
}

/* calculate terms of the sum for large t */
double logFl(double q, double v, double a, double w, int K)
{
    return ddm_logFl(q, v, a, w, K);
}

/* calculate distribution */
double pwiener(double q, double a, double v, double w, double err, int K, int epsFLAG) {
    return pwiener_inline(q, a, v, w, err, K, epsFLAG);
}
