
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "tools.h"
#include "pdf_fncs.h"
#include "fncs_seven.h"
#include "gauss.h"
#include <thread>


/* DENSITY */

/* dependencies */

/* integrand density */
int int_ddiff(unsigned dim, const double *x, void *p, unsigned fdim, double *retval) {
  my_params *params = static_cast<my_params*>(p);
  double t = (params->t);
  int low_or_up = (params->low_or_up);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st = (params->st);
  double errorW = (params->errorW);
  int K = (params->K);
  int epsFLAG = (params->epsFLAG);
  // double *val_ptr = (params->val_ptr);

  // usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
  //double temp = sv ? pow(x[0], 2) : 0;
  //double y = sv ? x[0] / (1 - temp) : 0;
  //double nu = sv ? v + sv * y : v;
  //double omega = sv ? (sw ? w + sw * (x[1] - 0.5) : w) : (sw ? w + sw * (x[0] - 0.5) : w);
  //double tau = sv ? ( sw ? (st ? t0 + st * x[2] : t0) : (st ? t0 + st * x[1] : t0) ) : ( sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0) );
  // usually: 0  = omega (w), 1 = tau (t0), depending on whether sv, sw, or st = 0
  double omega = sw ? w + sw * (x[0] - 0.5) : w;
  double tau = sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0);

  if (t - tau <= 0) retval[0] = 0.0;
  else {
    double ldW = dwiener(low_or_up * (t - tau), a, v, omega, sv, errorW, K, epsFLAG);

    double temp2 = 0;
    //if (sv) temp2 = -0.5 * pow(y, 2) - M_LN_SQRT_PI - 0.5 * M_LN2 + log1p(temp) - 2 * log1p(-temp);

    double integrand = exp(ldW + temp2);

    retval[0] = integrand;
  }
  return 0;
}

/* calculate density for 7-param diffusion */
void ddiff(int choice, double t, int low_or_up, double a, double v, double t0, double w, double sw, double sv, double st, double myerr, int K, int epsFLAG, int Neval, double *derivF, double *Rerr) {

  //double result;
  //double error;

  double value;
  // double valueln;

  double *val_ptr = &value;
  double errorW = myerr*0.1;

  my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st, errorW, K, epsFLAG, val_ptr};

  int dim = (sw!=0)+(st!=0);

  double *xmin = (double*)R_Calloc(dim, double);
  double *xmax = (double*)R_Calloc(dim, double);

  // 0  = s (v); 1 = u (w), 2 = v (w)
  //	if(sv) {
  //		xmin[0] = -1; xmax[0] = 1;
  //		for (int i = 1; i < dim; i++) {
  //			xmin[i] = 0;
  //			xmax[i] = 1;
  //		}
  //	} else {
  for (int i = 0; i < dim; i++) {
    xmin[i] = 0;
    xmax[i] = 1;
  }
  //	}
  if (st) xmax[dim-1] = fmin(1.0, (t-t0)/st);

  double reltol = 0.0;

  double abstol = myerr*0.9;

  double val, err;

  int Meval = Neval;

  //	printf("%u-dim integral, tolerance = %g\n", dim, tol);
  hcubature(int_ddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err);

  //if(err > abstol) Rprintf("absolute error not achieved: %g < %g\n", abstol, err);

  R_Free(xmin); R_Free(xmax);
  *derivF = val;
  if (*Rerr < err+errorW) *Rerr = err+errorW;

}
