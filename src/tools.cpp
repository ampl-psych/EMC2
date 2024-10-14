
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann
#define R_NO_REMAP
//#include <cmath>
#include <mutex>
#include "tools.h"
#include <Rinternals.h>

std::mutex mtx_samp, mtx_RCUI;


double logsum(double xa, double xb) {
	double temp;
	if (xa <= -INFINITY) return xb;
	if (xb <= -INFINITY) return xa;
	if (xa > xb) temp = xa + log1p(exp(xb - xa));
	else temp = xb + log1p(exp(xa - xb));
	return temp;
}

double logdiff(double xa, double xb) {
	double result;
	if (xb <= -INFINITY) return(xa);
	if (xa <= -INFINITY) return(xb);
	if (xa > xb) result = (xa + log1p(-exp(xb - xa)));
	else
		if (xb > xa) result = (xb + log1p(-exp(xa - xb)));
		else result = -INFINITY;
	return result;
}

double rexp(double x) {
	double result;
	if (x <= 700.0)
		result = exp(x);
	else
		result = exp(700.00);
	return result;
}

double lognormal(double x) {
	return	-0.5*x*x - M_LN_SQRT_PI - 0.5*M_LN2;
}

/*
The original function lnnorm was written by
	Jean Marie Linhart
	StataCorp LP
	jlinhart@stata.com
	January 4, 2008
and later modified
*/
double lnnorm(double z)
{
        int lower ;

        double z2, y, s, p1, q1, p2, q2, t, a1, a2 ;
        double n, m ;

	if (z==0.0e0) return(std::log(0.50e0)) ;

	if (z > LNNORM_MAX_X) return(0.0e0);
        if (z <= LNNORM_MIN_X) return(-0.5e0*z*z);

        if (z<0.0e0) {
                z= -z ;
                lower=1 ;
        }
        else    lower=0 ;
        //upper = !lower ;

        z2 = z*z ;

        y = exp(-0.5*z2) / SQR2PI ;
        n = y/z ;

        if (!( ( z>2.00e0))) {
                z *= y ;
                s=z ;
                t=0.0e0 ;

                for (n=3.0e0;s!=t;n+=2.0e0) {
                        t=s ;
                        z *= ((z2)/ (n)) ;
                        s+=z ;
                }
                if (lower) return(std::log(0.50e0-s)) ;
                return(std::log(0.50e0+s)) ;
        }

        a1=2.0e0 ;
        a2=0.0e0 ;
        n=z2+3.0e0 ;
        p1=1.0e0 ;
        q1=z ;
        p2=(n-1.0e0) ;
        q2=n*z ;
        m = p1/q1 ;
        t = p2/q2 ;
        s = m ;

        for (n+=4.0e0; m!=t && s!=t; n+=4.0e0) {
                a1 -= 8.0 ;
                a2 += a1 ;
                s = a2*p1 + n*p2 ;
                p1=p2 ;
                p2=s ;
                s = a2*q1 + n*q2 ;
                q1=q2 ;
                q2=s ;
                s=m ;
                m=t ;
                if (q2>1.0e30) {
                        p1 /= 1.0e30 ;
                        p2 /= 1.0e30 ;
                        q1 /= 1.0e30 ;
                        q2 /= 1.0e30 ;
                }
 		            t = p2/q2;
        }
        t = lower ? std::log(t) - 0.5*z2 - std::log(SQR2PI) : log1p(-y*t);
        return(t) ;
}

double logMill(double x) {
	double m;
	if (x > 1.0e5) return -std::log(x);
	m = lnnorm(-x) - lognormal(x);
	return m;
}

double phi1(double x, double y, double v) {
	return exp(2 * v*y) - exp(2 * v*x);
}

double lower_bound_var(double a, double vn, double wn) {
	double z = a * wn, phiza = phi1(z, a, vn);
	double temp = (-2 * a * phi1(0, z, vn) * (2 * vn * a * phi1(z, 2 * a, vn) + phi1(0, a, vn) * phiza)) * exp(2 * vn * a) / (pow(vn, 3) * pow(phi1(0, a, vn) * phiza, 2));
	temp += (4 * vn * z * (2 * a - z) * exp(2 * vn * (z + a)) + z * phi1(2 * z, 2 * a, vn)) / pow(vn, 3) / pow(phiza, 2);
	if (temp < 0) {
		Rprintf("! %20g%20g%20g%20g\n", a, vn, wn, temp);
		temp = 0.1;
	}
	return temp;
}

double coth(double x) {
	double temp = (exp(2 * x) + 1.0) / expm1(2*x);
	return temp;
}

double lower_bound_time(double a, double vn, double wn) {
	double temp, amw = a * (1 - wn);
	if (fabs(vn) < 1e-5) {
		temp = (pow(a, 2) - pow(amw, 2)) / 3.0;

	}
	else {
		temp = a  * coth(a*vn) - amw  * coth(vn*amw);
		temp /= vn;
	}
	return temp;
}

double exp_mean(int pm, double a, double v, double w) {
	if (pm == 1) { v = -v; w = 1 - w; }
	return lower_bound_time(a, v, w);

}

double oneuni() {
	double u;
	std::lock_guard<std::mutex> guard(mtx_samp);
	do {
		GetRNGstate();
		u = unif_rand();
		PutRNGstate();
	} while (u <= 0 || u >= 1);
	return u;

}

double oneuniL() {
	double u;
	std::lock_guard<std::mutex> guard(mtx_samp);
	do {
		GetRNGstate();
		u = unif_rand();
		PutRNGstate();
	} while (u < 0 || u >= 1);
	return u;

}

void R_CheckUserInterruptGuarded() {
  std::lock_guard<std::mutex> guard(mtx_RCUI);
  R_CheckUserInterrupt();
}


/* ----------------------------------- */





/* The following functions ...

rat_eval
small
intermediate
tail
gsl_cdf_ugaussian_Pinv
cheb_eval_e
erfc_xlt1_data
erfc_xlt1_cs
erfc_x15_data
erfc_x15_cs
erfc_x510_data
erfc_x510_cs
erfc8_sum
erfc8
gsl_sf_erfc_e
gsl_sf_erfc

... are copied from the GNU scientific library version 2.6 ...

... and the following functions ...

gsl_ran_gaussian
gsl_ran_gaussian_tail
gsl_ran_ugaussian_tail
gsl_error

... are also copied from the GNU scientific library version 2.6 but modified, too

Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover. GPL 3
Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003 Gerard Jungman. GPL 3
Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough. GPL 3


*/

double rat_eval(const double a[], const size_t na,
          const double b[], const size_t nb, const double x) {
  size_t i, j;
  double u, v, r;
  u = a[na - 1];
  for (i = na - 1; i > 0; i--)
    {
      u = x * u + a[i - 1];
    }
  v = b[nb - 1];
  for (j = nb - 1; j > 0; j--)
    {
      v = x * v + b[j - 1];
    }
  r = u / v;
  return r;
}

double small(double q) {
  const double a[8] = { 3.387132872796366608, 133.14166789178437745,
    1971.5909503065514427, 13731.693765509461125,
    45921.953931549871457, 67265.770927008700853,
    33430.575583588128105, 2509.0809287301226727
  };
  const double b[8] = { 1.0, 42.313330701600911252,
    687.1870074920579083, 5394.1960214247511077,
    21213.794301586595867, 39307.89580009271061,
    28729.085735721942674, 5226.495278852854561
  };
  double r = 0.180625 - q * q;
  double x = q * rat_eval (a, 8, b, 8, r);
  return x;
}

double intermediate(double r) {
  const double a[] = { 1.42343711074968357734, 4.6303378461565452959,
    5.7694972214606914055, 3.64784832476320460504,
    1.27045825245236838258, 0.24178072517745061177,
    0.0227238449892691845833, 7.7454501427834140764e-4
  };
  const double b[] = { 1.0, 2.05319162663775882187,
    1.6763848301838038494, 0.68976733498510000455,
    0.14810397642748007459, 0.0151986665636164571966,
    5.475938084995344946e-4, 1.05075007164441684324e-9
  };
  double x = rat_eval (a, 8, b, 8, (r - 1.6));
  return x;
}

double tail(double r) {
  const double a[] = { 6.6579046435011037772, 5.4637849111641143699,
    1.7848265399172913358, 0.29656057182850489123,
    0.026532189526576123093, 0.0012426609473880784386,
    2.71155556874348757815e-5, 2.01033439929228813265e-7
  };
  const double b[] = { 1.0, 0.59983220655588793769,
    0.13692988092273580531, 0.0148753612908506148525,
    7.868691311456132591e-4, 1.8463183175100546818e-5,
    1.4215117583164458887e-7, 2.04426310338993978564e-15
  };
  double x = rat_eval (a, 8, b, 8, (r - 5.0));
  return x;
}

double gsl_cdf_ugaussian_Pinv(const double P) {
  double r, x, pp;
  double dP = P - 0.5;
  if (P == 1.0)
    {
      return INFINITY;
    }
  else if (P == 0.0)
    {
      return -INFINITY;
    }

  if (fabs (dP) <= 0.425)
    {
      x = small (dP);

      return x;
    }

  pp = (P < 0.5) ? P : 1.0 - P;
  r = sqrt (-log (pp));
  if (r <= 5.0)
    {
      x = intermediate (r);
    }
  else
    {
      x = tail (r);
    }

  if (P < 0.5)
    {
      return -x;
    }
  else
    {
      return x;
    }
}

static inline int cheb_eval_e(const cheb_series * cs, const double x, gsl_sf_result * result);

static inline int cheb_eval_e(const cheb_series * cs,
            const double x,
            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  {
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}

/* modified gsl_ran_gaussian and renamed */
double onenorm() {
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      //x = -1 + 2 * gsl_rng_uniform_pos (r);
			x = -1 + 2 * oneuni();
      //y = -1 + 2 * gsl_rng_uniform_pos (r);
			y = -1 + 2 * oneuni();

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
	// return sigma * y * sqrt (-2.0 * log (r2) / r2);
	return 1.0 * y * sqrt (-2.0 * log (r2) / r2);
}


/* modified gsl_ran_gaussian_tail */
double gsl_ran_gaussian_tail(const double a, const double sigma) {
  /* Returns a gaussian random variable larger than a
   * This implementation does one-sided upper-tailed deviates.
   */

  double s = a / sigma;

  if (s < 1)
    {
      /* For small s, use a direct rejection method. The limit s < 1
         can be adjusted to optimise the overall efficiency */

      double x;

      do
        {
          //x = gsl_ran_gaussian (r, 1.0);
					x = onenorm();
        }
      while (x < s);
      return x * sigma;
    }
  else
    {
      /* Use the "supertail" deviates from the last two steps
       * of Marsaglia's rectangle-wedge-tail method, as described
       * in Knuth, v2, 3rd ed, pp 123-128.  (See also exercise 11, p139,
       * and the solution, p586.)
       */

      double u, v, x;

      do
        {
          u = oneuniL(); //gsl_rng_uniform (r);
          do
            {
              v = oneuniL(); //gsl_rng_uniform (r);
            }
          while (v == 0.0);
          x = sqrt (s * s - 2 * log (v));
        }
      while (x * u > s);
      return x * sigma;
    }
}

/* modified gsl_ran_ugaussian_tail */
double gsl_ran_ugaussian_tail(const double a) {
  return gsl_ran_gaussian_tail (a, 1.0) ;
}

static double erfc_xlt1_data[20] = {
  1.06073416421769980345174155056,
 -0.42582445804381043569204735291,
  0.04955262679620434040357683080,
  0.00449293488768382749558001242,
 -0.00129194104658496953494224761,
 -0.00001836389292149396270416979,
  0.00002211114704099526291538556,
 -5.23337485234257134673693179020e-7,
 -2.78184788833537885382530989578e-7,
  1.41158092748813114560316684249e-8,
  2.72571296330561699984539141865e-9,
 -2.06343904872070629406401492476e-10,
 -2.14273991996785367924201401812e-11,
  2.22990255539358204580285098119e-12,
  1.36250074650698280575807934155e-13,
 -1.95144010922293091898995913038e-14,
 -6.85627169231704599442806370690e-16,
  1.44506492869699938239521607493e-16,
  2.45935306460536488037576200030e-18,
 -9.29599561220523396007359328540e-19
};
static cheb_series erfc_xlt1_cs = {
  erfc_xlt1_data,
  19,
  -1, 1,
  12
};

static double erfc_x15_data[25] = {
  0.44045832024338111077637466616,
 -0.143958836762168335790826895326,
  0.044786499817939267247056666937,
 -0.013343124200271211203618353102,
  0.003824682739750469767692372556,
 -0.001058699227195126547306482530,
  0.000283859419210073742736310108,
 -0.000073906170662206760483959432,
  0.000018725312521489179015872934,
 -4.62530981164919445131297264430e-6,
  1.11558657244432857487884006422e-6,
 -2.63098662650834130067808832725e-7,
  6.07462122724551777372119408710e-8,
 -1.37460865539865444777251011793e-8,
  3.05157051905475145520096717210e-9,
 -6.65174789720310713757307724790e-10,
  1.42483346273207784489792999706e-10,
 -3.00141127395323902092018744545e-11,
  6.22171792645348091472914001250e-12,
 -1.26994639225668496876152836555e-12,
  2.55385883033257575402681845385e-13,
 -5.06258237507038698392265499770e-14,
  9.89705409478327321641264227110e-15,
 -1.90685978789192181051961024995e-15,
  3.50826648032737849245113757340e-16
};
static cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  16
};

static double erfc_x510_data[20] = {
  1.11684990123545698684297865808,
  0.003736240359381998520654927536,
 -0.000916623948045470238763619870,
  0.000199094325044940833965078819,
 -0.000040276384918650072591781859,
  7.76515264697061049477127605790e-6,
 -1.44464794206689070402099225301e-6,
  2.61311930343463958393485241947e-7,
 -4.61833026634844152345304095560e-8,
  8.00253111512943601598732144340e-9,
 -1.36291114862793031395712122089e-9,
  2.28570483090160869607683087722e-10,
 -3.78022521563251805044056974560e-11,
  6.17253683874528285729910462130e-12,
 -9.96019290955316888445830597430e-13,
  1.58953143706980770269506726000e-13,
 -2.51045971047162509999527428316e-14,
  3.92607828989125810013581287560e-15,
 -6.07970619384160374392535453420e-16,
  9.12600607264794717315507477670e-17
};
static cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  12
};

gsl_error_handler_t * gsl_error_handler = NULL;

/* modified gsl_error */
void gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_error_handler)
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return ;
    }

  Rprintf ("ERROR");

  //fflush (stdout);
  Rprintf ("Default GSL error handler invoked.\n");
  //fflush (stderr);

  //abort ();
}

static double erfc8_sum(double x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */

  static double P[] = {
      2.97886562639399288862,
      7.409740605964741794425,
      6.1602098531096305440906,
      5.019049726784267463450058,
      1.275366644729965952479585264,
      0.5641895835477550741253201704
  };
  static double Q[] = {
      3.3690752069827527677,
      9.608965327192787870698,
      17.08144074746600431571095,
      12.0489519278551290360340491,
      9.396034016235054150430579648,
      2.260528520767326969591866945,
      1.0
  };
  double num=0.0, den=0.0;
  int i;

  num = P[5];
  for (i=4; i>=0; --i) {
      num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
      den = x*den + Q[i];
  }

  return num/den;
}

inline static double erfc8(double x)
{
  double e;
  e = erfc8_sum(x);
  e *= exp(-x*x);
  return e;
}

int gsl_sf_erfc_e(double x, gsl_sf_result * result) {
  const double ax = fabs(x);
  double e_val, e_err;

  /* CHECK_POINTER(result) */

  if(ax <= 1.0) {
    double t = 2.0*ax - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&erfc_xlt1_cs, t, &c);
    e_val = c.val;
    e_err = c.err;
  }
  else if(ax <= 5.0) {
    double ex2 = exp(-x*x);
    double t = 0.5*(ax-3.0);
    gsl_sf_result c;
    cheb_eval_e(&erfc_x15_cs, t, &c);
    e_val = ex2 * c.val;
    e_err = ex2 * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON);
  }
  else if(ax < 10.0) {
    double exterm = exp(-x*x) / ax;
    double t = (2.0*ax - 15.0)/5.0;
    gsl_sf_result c;
    cheb_eval_e(&erfc_x510_cs, t, &c);
    e_val = exterm * c.val;
    e_err = exterm * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON + GSL_DBL_EPSILON);
  }
  else {
    e_val = erfc8(ax);
    e_err = (x*x + 1.0) * GSL_DBL_EPSILON * fabs(e_val);
  }

  if(x < 0.0) {
    result->val  = 2.0 - e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    result->val  = e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}

double gsl_sf_erfc(double x)
{
  EVAL_RESULT(gsl_sf_erfc_e(x, &result));
}


/* ----------------------------------- */


void check_interruption(void *ptr) {
  R_CheckUserInterrupt();
}

int is_interruption() {
  return !(R_ToplevelExec(check_interruption, nullptr));
}
