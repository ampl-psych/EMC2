
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "tools.h"
#include "pdf_fncs.h"
#include <cmath>


/* DENSITY */

/* calculate number of terms needed for short t */
double ks(double t, double w, double eps) {
	double K1 = (sqrt(2.0 * t) + w) / 2.0;
	double u_eps = fmin(-1.0, M_LN2 + M_LNPI + 2.0 * std::log(t) + 2.0 * (eps));
	double	arg = -t * (u_eps - sqrt(-2.0 * u_eps - 2.0));
	double 	K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
	return ceil(fmax(K1, K2));
}

/* calculate number of terms needed for large t */
double kl(double q, double v, double w, double err) {
	double K1 = 1.0 / (M_PI * sqrt(q)), K2=0.0;
	double temp = -2.0 * (std::log(M_PI * q) + err);
	if (temp>=0) K2 = sqrt(temp/(pow(M_PI, 2) * q));
	return ceil(fmax(K1,K2));
}

/* calculate terms of the sum for short t */
double logfs(double t, double w, int K) {
	if (w == 0) return R_NegInf;
	double	fplus = R_NegInf, fminus = R_NegInf, twot = 2.0 * t;
	if (K > 0)
		for (int k = K; k >= 1; k--) {
			double temp1 = w + 2.0 * k, temp2 = w - 2.0 * k;

			fplus = logsum(std::log(temp1) - pow(temp1, 2) / twot, fplus);
			fminus = logsum(std::log(-temp2) - pow(temp2, 2) / twot, fminus);
		}
	fplus = logsum(std::log(w) - pow(w, 2) / twot, fplus);
	return  -0.5 * M_LN2 - M_LN_SQRT_PI - 1.5 * std::log(t) + logdiff(fplus, fminus);
}

double logfs_linear(double t, double w, int K) {
  if (w == 0.0) return R_NegInf;
  if (K < 0) K = 0;

  const double twot     = 2.0 * t;
  const double inv_twot = 1.0 / twot;

  double sum_plus  = 0.0;
  double sum_minus = 0.0;

  // Main sum over k
  if (K > 0) {
    #pragma omp simd reduction(+:sum_plus,sum_minus)
    for (int k = 1; k <= K; ++k) {
      const double temp1 = w + 2.0 * k;
      const double temp2 = w - 2.0 * k;

      const double e1 = std::exp(-temp1 * temp1 * inv_twot);
      const double e2 = std::exp(-temp2 * temp2 * inv_twot);

      sum_plus  += temp1  * e1;
      sum_minus += (-temp2) * e2;  // -temp2 > 0
    }
  }

  // k = 0 term (the w term)
  {
    const double e0 = std::exp(-w * w * inv_twot);
    sum_plus += w * e0;
  }

  // sum_minus may be zero if K == 0
  const double S = sum_plus - sum_minus;
  if (!(S > 0.0)) { // guard against <=0 due to numerical issues
    return R_NegInf;
  }

  return -0.5 * M_LN2 - M_LN_SQRT_PI - 1.5 * std::log(t) + std::log(S);
}



/* calculate terms of the sum for large t */
double logfl(double q, double v, double w, int K) {
	if (w == 0) return R_NegInf;
	double fplus = R_NegInf, fminus = R_NegInf;
	double halfq = q / 2.0;
	for (int k = K; k >= 1; k--) {
		double temp = k * M_PI;
		double check = sin(temp * w);
		if (check > 0) fplus = logsum(std::log(static_cast<double>(k)) - pow(temp, 2) * halfq + std::log(check), fplus);
		else fminus = logsum(std::log(static_cast<double>(k)) - pow(temp, 2) * halfq + std::log(-check), fminus);
	}
	return	logdiff(fplus, fminus) + M_LNPI;
}

double logfl_linear(double q, double v, double w, int K) {
  (void)v;  // v unused
  if (w == 0.0) return R_NegInf;
  if (K <= 0)   return R_NegInf;

  const double halfq = q / 2.0;

  double sum_plus  = 0.0;
  double sum_minus = 0.0;

  #pragma omp simd reduction(+:sum_plus,sum_minus)
  for (int k = 1; k <= K; ++k) {
    const double temp  = k * M_PI;
    const double s     = std::sin(temp * w);  // can be +/-/0
    const double e     = std::exp(- (temp * temp) * halfq);

    // signed contribution
    const double base  = k * s * e;  // might be pos/neg/zero

    // add positive and negative parts separately, branchless
    const double pos   = (base > 0.0) ? base  : 0.0;
    const double neg   = (base < 0.0) ? -base : 0.0;

    sum_plus  += pos;
    sum_minus += neg;
  }

  const double S = sum_plus - sum_minus;
  if (!(S > 0.0)) return R_NegInf;

  return std::log(S) + M_LNPI;
}


/* calculate density */
double dwiener(double q, double a, double vn, double wn, double sv, double err, int K, int epsFLAG) {
	if (q == 0.0) {
		return R_NegInf;
	}
	double kll, kss, ans, v, w;
	if(!epsFLAG && K==0) {
		err = -27.63102;  // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102;  // exp(err) = 1.e-12
	else if(epsFLAG) err = std::log(err);

	if (q >= 0) {
		w = 1.0 - wn;
		v = -vn;
	}
	else {
		q = fabs(q);
		w = wn;
		v = vn;
	}

	double q_asq = q / pow(a, 2);
	ans = 0.0;

	/* calculate the number of terms needed for short t*/
	double eta_sqr = pow(sv, 2);
	double temp = 1 + eta_sqr * q;
	double lg1 = (eta_sqr * pow(a * w, 2) - 2 * a * v * w - pow(v, 2) * q) / 2.0 / temp - 2 *std::log(a) - 0.5 *std::log(temp);
	//double lg1 = (-v * a * w - (pow(v, 2)) * q / 2.0) - 2.0*std::log(a);
	double es = (err - lg1);
	kss = ks(q_asq, w, es);
	/* calculate the number of terms needed for large t*/
	double el = es;
	kll = kl(q_asq, v, w, el);

	// if small t is better
	if (2 * kss <= kll) {
		if((epsFLAG && kss<K) || !epsFLAG) kss = K;
		ans = lg1 + logfs_linear(q_asq, w, static_cast<int>(kss));
	}
	// if large t is better
	else {
		if((epsFLAG && kll<K) || !epsFLAG) kll = K;
		ans = lg1 + logfl_linear(q_asq, v, w, static_cast<int>(kll));
	}

	return ans;
}
