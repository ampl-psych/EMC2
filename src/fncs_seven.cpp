
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

/* integrand d/da */
int int_daddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

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
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

		double temp2 = 0;
		//if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dadwiener(low_or_up * (t-tau), a, v, omega, sv, ldW, val_ptr, errorW, K, epsFLAG);
		double dda = val_ptr[0];

		double integrand = dda * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dv */
int int_dvddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

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
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

		double temp2 = 0;
		//if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dvdwiener(low_or_up * (t-tau), a, v, omega, sv, ldW, val_ptr);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dt0 */
int int_dt0ddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

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
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

		double temp2 = 0;
		//if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		dtdwiener(t-tau, a, -low_or_up*v, wn, sv, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = -val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dz */
int int_dwddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

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
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

		double temp2 = 0;
		//if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dwdwiener(low_or_up * (t-tau), a, v, omega, sv, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dsz */
int int_dswddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	//double temp = sv ? pow(x[0], 2) : 0;
	//double y = sv ? x[0] / (1 - temp) : 0;
	//double nu = sv ? v + sv * y : v;
	double omega = w + sw * (x[0] - 0.5);
	double tau = st ? t0 + st * x[1] : t0;
	double temp_sw = x[0]-0.5;

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

		double temp2 = 0;
		//if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		dwdwiener(low_or_up * (t-tau), a, v, omega, sv, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = temp_sw * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dsv */
int int_dsvddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	//double temp = pow(x[0], 2);
	//double y = x[0] / (1 - temp);
	//double nu = v + sv * y;
	double y = 1;
	double omega = sw ? w + sw * (x[0] - 0.5) : w;
	double tau = sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0);

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

	  double temp2 = 0;
		//temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		//dvdwiener(low_or_up * (t-tau), a, v, omega, sv, ldW, val_ptr);
		dsvdwiener(low_or_up * (t-tau), a, v, omega, sv, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = y * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dst0 */
int int_dst0ddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

	// usually: 0  = s (v); 1 = u (w), 2 = v (t), depending on whether sv, sw, or st = 0
	//double temp = sv ? pow(x[0], 2) : 0;
	//double y = sv ? x[0] / (1 - temp) : 0;
	//double nu = sv ? v + sv * y : v;
	double omega = sw ? w + sw * (x[0] - 0.5) : w;
	double tau = sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0);
	double temp_st0 = sw ? -x[1] : -x[0];

	if (t - tau <= 0) retval[0] = 0.0;
	else {
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

		double temp2 = 0;
		//if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		dtdwiener(t-tau, a, -low_or_up*v, wn, sv, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = temp_st0 * val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}

/* integrand d/dt */
int int_dtddiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
	double *val_ptr = (params->val_ptr);

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
		double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

		double temp2 = 0;
		//if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

		double wn = omega;
		if (low_or_up==1) wn = 1-omega;

		dtdwiener(t-tau, a, -low_or_up*v, wn, sv, ldW, val_ptr, errorW, K, epsFLAG);

		double integrand = val_ptr[0] * exp(temp2);

		retval[0] = integrand;
	}
	return 0;
}


/* integrand d/dst0 ___CDF___ */
int int_dst0pdiff(unsigned dim, const double* x, void* p, unsigned fdim, double* retval) {
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
  double omega = sw ? w + sw * (x[0] - 0.5) : w;
  double tau = sw ? (st ? t0 + st * x[1] : t0) : (st ? t0 + st * x[0] : t0);
  double temp_st0 = sw ? -x[1] : -x[0];

  if (t - tau <= 0) retval[0] = 0.0;
  else {
    double ldW = dwiener(low_or_up * (t-tau), a, v, omega, sv, errorW, K, epsFLAG);

    double temp2 = 0;
    //if (sv) temp2 = - 0.5*pow(y, 2) - M_LN_SQRT_PI - 0.5*M_LN2 + log1p(temp) - 2*log1p(-temp);

    double integrand = temp_st0 * exp(ldW + temp2);

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
	double errorW1 = sw, errorW2 = sw, abstol1 = sw, abstol2 = sw;
	if (!myerr) myerr = 1.e-12;
	if (choice == 5 || choice == 7) {
	  if (st) {
	    errorW1 *= myerr*0.05;
	    errorW2 *= myerr*0.05;
	    abstol1 *= myerr*0.9*2/3;
	    abstol2 *= myerr*0.9*1/3;
	  } else{
	    errorW1 *= myerr*0.05;
	    errorW2 *= myerr*0.05;
	    abstol1 *= myerr*0.9;
	  }
	}
	double errorW = myerr*0.1;

	my_params params1; my_params params2; my_params params3;
	if (choice == 7) {
	  params1 = {t, low_or_up, a, v, t0, w, sw, sv, st, errorW1, K, epsFLAG, val_ptr};
	  params2 = {t, low_or_up, a, v, t0+st, w, sw, sv, 0, errorW2, K, epsFLAG, val_ptr};
	}
	if (choice == 5) {
	  params1 = {t, low_or_up, a, v, t0, w, sw, sv, st, errorW1, K, epsFLAG, val_ptr};
	  params2 = {t, low_or_up, a, v, t0, w+0.5*sw, 0, sv, st, errorW2/2, K, epsFLAG, val_ptr};
	  params3 = {t, low_or_up, a, v, t0, w-0.5*sw, 0, sv, st, errorW2/2, K, epsFLAG, val_ptr};
	}
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
	switch (choice) {
		case 0: hcubature(int_ddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
	  case 1: hcubature(int_daddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
		case 2: hcubature(int_dvddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
		case 3: hcubature(int_dt0ddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
		case 4: hcubature(int_dwddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
		case 5: //hcubature(int_dswddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
		  double tmp5, tmp5err;
		  if (st) {
		    hcubature(int_ddiff, &params1, dim, xmin, xmax, Meval, abstol1, reltol, &val, &err);
		    tmp5 = -val;
		    tmp5err = err;
		    hcubature(int_ddiff, &params2, 1, xmin, xmax, Meval, abstol2/2, reltol, &val, &err);
		    tmp5 += 0.5*val;
		    tmp5err += err;
		    hcubature(int_ddiff, &params3, 1, xmin, xmax, Meval, abstol2/2, reltol, &val, &err);
		    val = 1/sw*(tmp5 + 0.5*val);
		    err += tmp5err;
		    errorW = errorW1 + errorW2;
		  } else {
		    hcubature(int_ddiff, &params1, dim, xmin, xmax, Meval, abstol1, reltol, &val, &err);
		    tmp5 = val;
		    val = 1/sw*(-tmp5 +
		      0.5*(exp(dwiener(low_or_up * (t-t0), a, v, w+0.5*sw, sv, errorW2/2, K, epsFLAG)) +
		           exp(dwiener(low_or_up * (t-t0), a, v, w-0.5*sw, sv, errorW2/2, K, epsFLAG))));
		    errorW = errorW1 + errorW2;
		  }
		  break;
		case 6: hcubature(int_dsvddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
		case 7: //hcubature(int_dst0ddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
		  double tmp7, tmp7err;
		  if (sw) {
		    hcubature(int_ddiff, &params1, dim, xmin, xmax, Meval, abstol1, reltol, &val, &err);
		    tmp7 = val;
		    tmp7err = err;
		    hcubature(int_ddiff, &params2, 1, xmin, xmax, Meval, abstol2, reltol, &val, &err);
		    val = 1/st*(val - tmp7);
		    err += tmp7err;
		    errorW = errorW1 + errorW2;
		  } else {
		    hcubature(int_ddiff, &params1, dim, xmin, xmax, Meval, abstol1, reltol, &val, &err);
		    tmp7 = val;
		    val = 1/st*(-tmp7 + exp(dwiener(low_or_up * (t-(t0+st)), a, v, w, sv, errorW2, K, epsFLAG)));
		    errorW = errorW1 + errorW2;
		  }
		  break;
		case 8: hcubature(int_dtddiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
	  case 9: hcubature(int_dst0pdiff, &params, dim, xmin, xmax, Meval, abstol, reltol, &val, &err); break;
	}
	//if(err > abstol) Rprintf("absolute error not achieved: %g < %g\n", abstol, err);

	R_Free(xmin); R_Free(xmax);
	*derivF = val;
	if (*Rerr < err+errorW) *Rerr = err+errorW;

}
