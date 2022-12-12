#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

/* Robust normal distribution for LBA */

double pnormP(double q, double mean = 0.0, double sd = 1.0,
                    bool lower = true, bool log = false, bool robust = false){
	if (robust == true){
		if (q < -7.) {
			return 0.;
		} else if (q > 7.){
			return 1.;
		}
	}
	return R::pnorm(q, mean, sd, lower, log);
}

double dnormP(double x, double mean = 0.0, double sd = 1.0,
                    bool log = false, bool robust = false){
	if (robust == true){
		if (x < -7.) {
			return 0.;
		} else if (x > 7.){
			return 1.;
		}
	}
	return R::dnorm(x, mean, sd, log);
}


/* LBA Functions */

double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true, bool robust = false){
	double denom = 1.;
	if (posdrift) {
		denom = pnormP(v / sv, 0., 1., true, false, robust);
		if (denom < 1e-10)
			denom = 1e-10;
	}

	double cdf;

	if (A > 1e-10){
		double zs = t * sv;
		double cmz = b - t * v;
		double xx = cmz - A;
		double cz = cmz / zs;
		double cz_max = xx / zs;
		cdf = (1. + (zs * (dnormP(cz_max, 0., 1., false, robust) - dnormP(cz, 0., 1., false, robust))
                 + xx * pnormP(cz_max, 0., 1., true, false, robust) - cmz * pnormP(cz, 0., 1., true, false, robust))/A) / denom;
	} else {
		cdf = pnormP(b / t, v, sv, false, false, robust) / denom;
	}

	if (cdf < 0.) {
		return 0.;
	} else if (cdf > 1.){
		return 1.;
	}
	return cdf;
}

double dlba_norm(double t, double A,double b, double v, double sv,
                 bool posdrift = true, bool robust = false){
	double denom = 1.;
	if (posdrift) {
		denom = pnormP(v / sv, 0., 1., true, false, robust);
		if (denom < 1e-10)
			denom = 1e-10;
	}

	double pdf;

	if (A > 1e-10){
		double zs = t * sv;
		double cmz = b - t * v;;
		double cz = cmz / zs;
		double cz_max = (cmz - A) / zs;
		pdf = (v * (pnormP(cz, 0., 1., true, false, robust) - pnormP(cz_max, 0., 1., true, false, robust)) +
			sv * (dnormP(cz_max, 0., 1., false, robust) - dnormP(cz, 0., 1., false, robust))) / (A * denom);
	} else {
		pdf = dnormP(b / t, v, sv, false, robust) * b / (t * t * denom);
	}

	if (pdf < 0.) {
		return 0.;
	}
	return pdf;
}


// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
  NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
  bool posdrift = true, bool robust = false)

{
	int n = t.size();
	NumericVector pdf(n);

	for (int i = 0; i < n; i++){
	  pdf[i] = dlba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift, robust);
	}
	return pdf;
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
  NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
  bool posdrift = true, bool robust = false)

{
	int n = t.size();
	NumericVector cdf(n);

	for (int i = 0; i < n; i++){
		cdf[i] = plba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift, robust);
	}
	return cdf;
}

