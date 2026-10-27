// Normal probabilities and scaled tails. The rational coefficients are from
// R's pnorm, based on Cody's Algorithm 715.
// Copyright (C) 2000-2024 The R Core Team, 2003 The R Foundation,
// 1998 Ross Ihaka. Distributed under GPL version 2 or later, without
// any warranty; see https://www.R-project.org/Licenses/.
#ifndef PNORM_UTILS_H
#define PNORM_UTILS_H

#include <cmath>
#include <algorithm>
#include "nan_check.h"
#include "r_constants.h"

namespace pnorm_detail {
constexpr double INV_SQRT_2PI = 0.39894228040143267794;
constexpr double LOG_SQRT_2PI = 0.91893853320467274178;

// Phi(x) - 1/2, for |x| <= 0.67448975.
inline double central_normal(double x) {
  const double u = x*x;
  return x * ((((0.065682337918207449113*u + 2.2352520354606839287)*u
      + 161.02823106855587881)*u + 1067.6894854603709582)*u
      + 18154.981253343561249) /
      ((((u + 47.20258190468824187)*u + 976.09855173777669322)*u
      + 10260.932208618978205)*u + 45507.789335026729956);
}

// H(x) = exp(x*x/2) Phi(-x). Its derivative is evaluated from the
// rational functions, avoiding cancellation in x*H(x)-phi(0).
inline double normal_tail_factor(double x, double* derivative = nullptr) {
  if (is_pos_inf(x)) { if (derivative) *derivative = -0.0; return 0.0; }
  if (x <= 0.67448975) {
    const double h = std::exp(0.5*x*x) * (0.5 - central_normal(x));
    if (derivative) *derivative = x*h - INV_SQRT_2PI;
    return h;
  }
  if (x <= 5.6568542494923801952) {
    constexpr double c[] = {1.0765576773720192317e-8, 0.39894151208813466764,
      8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
      2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124,
      9842.7148383839780218};
    constexpr double d[] = {1.0, 22.266688044328115691, 235.38790178262499861,
      1519.377599407554805, 6485.558298266760755, 18615.571640885098091,
      34900.952721145977266, 38912.003286093271411, 19685.429676859990727};
    double n=c[0], den=d[0], dn=0.0, dd=0.0;
    for (int i=1; i<9; ++i) {
      if (derivative) { dn=dn*x+n; dd=dd*x+den; }
      n=n*x+c[i]; den=den*x+d[i];
    }
    const double h=n/den;
    if (derivative) *derivative=(dn-h*dd)/den;
    return h;
  }
  const double inv=1.0/x, u=inv*inv;
  constexpr double p[] = {0.02307344176494017303, 0.21589853405795699,
    0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
    2.9112874951168792e-5};
  constexpr double q[] = {1.0, 1.28426009614491121, 0.468238212480865118,
    0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5};
  double n=p[0], den=q[0], dn=0.0, dd=0.0;
  for (int i=1; i<6; ++i) {
    if (derivative) { dn=dn*u+n; dd=dd*u+den; }
    n=n*u+p[i]; den=den*u+q[i];
  }
  const double r=n/den;
  if (derivative) *derivative=-u*(INV_SQRT_2PI-3*u*r-2*u*u*(dn-r*dd)/den);
  return (INV_SQRT_2PI-u*r)*inv;
}

// Fixed Gauss-Legendre integration on [0,1], used only for differences
// whose endpoints are too close for subtraction. No allocation or recursion.
template<class F> inline double integrate_unit(F f) {
  constexpr double x[] = {0.095012509837637440185,0.28160355077925891323,
    0.45801677765722738634,0.61787624440264374845,0.75540440835500303390,
    0.86563120238783174388,0.94457502307323257608,0.98940093499164993260};
  constexpr double w[] = {0.18945061045506849629,0.18260341504492358887,
    0.16915651939500253819,0.14959598881657673208,0.12462897125553387205,
    0.09515851168249278481,0.06225352393864789286,0.02715245941175409485};
  double sum=0;
  for (int i=0;i<8;++i) sum += w[i]*(f(0.5-0.5*x[i])+f(0.5+0.5*x[i]));
  return 0.5*sum;
}

// H(x)-H(x+width), x >= -0.67448975 and width >= 0.
inline double normal_tail_difference(double x, double width) {
  if (width < 0.1*(1.0+std::max(0.0,x))) {
    return width*integrate_unit([&](double u) {
      double derivative;
      normal_tail_factor(x+width*u,&derivative);
      return -derivative;
    });
  }
  return normal_tail_factor(x)-normal_tail_factor(x+width);
}

inline double gaussian_exponential(double x) {
  const double e=(-0.5*x)*x;
  if (e < -746.0) return 0.0;
  return std::exp(e) * (1.0+std::fma(-0.5*x,x,-e));
}

// Reuse exp(-x*x/2) when a distribution also needs the normal density.
inline double normal_probability(double x,double gaussian,bool lower) {
  if (std::fabs(x)<=0.67448975)
    return 0.5+(lower ? central_normal(x) : -central_normal(x));
  const double tail=gaussian*normal_tail_factor(std::fabs(x));
  return (x<0)==lower ? tail : 1.0-tail;
}

inline double pnorm_std(double x, bool lower, bool logp) {
  if (is_nan(x)) return x;
  if (is_inf(x)) {
    const bool one=(x>0)==lower;
    return logp ? (one ? 0.0 : neg_inf()) : double(one);
  }
  const double y=std::fabs(x);
  if (y <= 0.67448975) {
    const double p=0.5+(lower ? central_normal(x) : -central_normal(x));
    return logp ? std::log(p) : p;
  }
  const double h=normal_tail_factor(y);
  const bool small=(x<0)==lower;
  if (logp && small) return (-0.5*y)*y+std::log(h);
  const double tail=gaussian_exponential(y)*h;
  return logp ? std::log1p(-tail) : (small ? tail : 1.0-tail);
}

inline double normal_short_interval(double lo,double width,bool logp) {
  const double centre=lo+0.5*width;
  // Centred Gaussian integral. The branch bound keeps |centre*width|
  // and width below 0.1; terms after degree eight are below roundoff.
  const double z=centre*width, u=width*width;
  double hm2=1.0,hm1=z,coefficient=0.25,correction=0.0;
  for (int n=2;n<=8;++n) {
    const double h=z*hm1-(n-1)*u*hm2;
    coefficient/=2*(n+1);
    if (n%2==0) correction+=coefficient*h;
    hm2=hm1;hm1=h;
  }
  if (!logp) return gaussian_exponential(centre)*INV_SQRT_2PI*width*(1+correction);
  return (-0.5*centre)*centre-LOG_SQRT_2PI+std::log(width)+std::log1p(correction);
}
inline double normal_interval(double lo, double hi, bool logp=false) {
  if (is_nan(lo) || is_nan(hi) || lo>hi) return na_real();
  if (lo==hi) return logp ? neg_inf() : 0.0;
  if (is_neg_inf(lo)) return pnorm_std(hi,true,logp);
  if (is_pos_inf(hi)) return pnorm_std(lo,false,logp);
  const double width=hi-lo, centre=lo+0.5*width;
  if (width*(1.0+std::fabs(centre)) < 0.1)
    return normal_short_interval(lo,width,logp);
  const bool lower=lo+hi<0.0;
  if (!logp)
    return pnorm_std(lower ? hi : lo,lower,false)-pnorm_std(lower ? lo : hi,lower,false);
  const double a=pnorm_std(lower ? hi : lo,lower,true);
  const double b=pnorm_std(lower ? lo : hi,lower,true);
  if (is_neg_inf(a)) return a;
  return a+std::log(-std::expm1(b-a));
}
// Preserve a known width when forming hi-lo would lose significant digits.
inline double normal_interval_width(double lo,double width,bool logp=false) {
  if (width>=0 && width*(1.0+std::fabs(lo+0.5*width))<0.1)
    return normal_short_interval(lo,width,logp);
  return normal_interval(lo,lo+width,logp);
}

namespace as241 {

// Coefficients for the central region (|p - 0.5| <= 0.425)
constexpr double A0 =  3.3871328727963666080e0;
constexpr double A1 =  1.3314166789178437745e+2;
constexpr double A2 =  1.9715909503065514427e+3;
constexpr double A3 =  1.3731693765509461125e+4;
constexpr double A4 =  4.5921953931549871457e+4;
constexpr double A5 =  6.7265770927008700853e+4;
constexpr double A6 =  3.3430575583588128105e+4;
constexpr double A7 =  2.5090809287301226727e+3;

constexpr double B1 =  4.2313330701600911252e+1;
constexpr double B2 =  6.8718700749205790830e+2;
constexpr double B3 =  5.3941960214247511077e+3;
constexpr double B4 =  2.1213794301586595867e+4;
constexpr double B5 =  3.9307895800092710610e+4;
constexpr double B6 =  2.8729085735721942674e+4;
constexpr double B7 =  5.2264952788528545610e+3;

// Coefficients for the tail regions (0.425 < |p - 0.5| < 0.5)
constexpr double C0 =  1.42343711074721209650e0;
constexpr double C1 =  4.63033784615654529590e0;
constexpr double C2 =  5.76949722146864628717e0;
constexpr double C3 =  3.64784832476320460504e0;
constexpr double C4 =  1.27045825245236838258e0;
constexpr double C5 =  2.41780725177450611770e-1;
constexpr double C6 =  2.27001535109994502416e-2;
constexpr double C7 =  7.74545433931955401503e-4;

constexpr double D1 =  2.05319162663775882187e0;
constexpr double D2 =  1.67638483950684182780e0;
constexpr double D3 =  6.89767334985100004550e-1;
constexpr double D4 =  1.48103976427480074590e-1;
constexpr double D5 =  1.51986665636164571966e-2;
constexpr double D6 =  5.47593808499534494600e-4;
constexpr double D7 =  1.05075007164441684324e-9;

// Coefficients for sqrt(-log(min(p, 1-p))) > 5.
constexpr double E0 =  6.65790464350110377720e0;
constexpr double E1 =  5.46378491116411436990e0;
constexpr double E2 =  1.78482653991729133580e0;
constexpr double E3 =  2.96560571828504891230e-1;
constexpr double E4 =  2.65321895265761230930e-2;
constexpr double E5 =  1.24266094738807843860e-3;
constexpr double E6 =  2.71155556874348757815e-5;
constexpr double E7 =  2.01033439929228813265e-7;

constexpr double F1 =  5.99832206555887937690e-1;
constexpr double F2 =  1.36929880922735805310e-1;
constexpr double F3 =  1.48753612908506508940e-2;
constexpr double F4 =  7.86869131145613259100e-4;
constexpr double F5 =  1.84631831751005468180e-5;
constexpr double F6 =  1.42151175831644588870e-7;
constexpr double F7 =  2.04426310338993978564e-15;

constexpr double SPLIT1 = 0.425;
constexpr double SPLIT2 = 5.0;
constexpr double CONST1 = 0.180625;   // SPLIT1^2
constexpr double CONST2 = 1.6;

[[gnu::always_inline]] inline double qnorm(double p)
{
  if (p <= 0.0) return -INFINITY;
  if (p >= 1.0) return  INFINITY;

  const double q = p - 0.5;

  if (std::fabs(q) <= SPLIT1) {
    // Central region
    const double r = CONST1 - q * q;
    return q * (((((((A7*r+A6)*r+A5)*r+A4)*r+A3)*r+A2)*r+A1)*r+A0) /
      (((((((B7*r+B6)*r+B5)*r+B4)*r+B3)*r+B2)*r+B1)*r+1.0);
  }

  // Tail regions — work with the smaller of p and 1-p
  double r = std::sqrt(-std::log(q < 0.0 ? p : 1.0 - p));
  double x;

  if (r <= SPLIT2) {
    r -= CONST2;
    x = (((((((C7*r+C6)*r+C5)*r+C4)*r+C3)*r+C2)*r+C1)*r+C0) /
      (((((((D7*r+D6)*r+D5)*r+D4)*r+D3)*r+D2)*r+D1)*r+1.0);
  } else {
    r -= SPLIT2;
    x = (((((((E7*r+E6)*r+E5)*r+E4)*r+E3)*r+E2)*r+E1)*r+E0) /
      (((((((F7*r+F6)*r+F5)*r+F4)*r+F3)*r+F2)*r+F1)*r+1.0);
  }

  return q < 0.0 ? -x : x;
}

} // namespace as241

inline double cpp_qnorm_std(double p)
{
  return as241::qnorm(p);
}


inline double fast_dnorm_std(double x) { return INV_SQRT_2PI*gaussian_exponential(x); }
inline double fast_dnorm(double x,double mean,double sd) { return fast_dnorm_std((x-mean)/sd)/sd; }
inline double fast_log_dnorm(double x,double mean,double sd) {
  const double z=(x-mean)/sd;
  return (-0.5*z)*z-LOG_SQRT_2PI-std::log(sd);
}
inline double fast_dlnorm(double t,double m,double s) {
  return fast_dnorm_std((std::log(t)-m)/s)/(t*s);
}
} // namespace pnorm_detail

#define PNORM_STD(x, lower, logp) pnorm_detail::pnorm_std((x),(lower),(logp))
#define QNORM_STD(p) pnorm_detail::cpp_qnorm_std((p))
#define PLNORM(t,m,s) PNORM_STD((std::log(t)-(m))/(s),true,false)
#define DLNORM(t,m,s) pnorm_detail::fast_dlnorm((t),(m),(s))
#define DNORM_STD(x) pnorm_detail::fast_dnorm_std((x))
#define DNORM(x,mean,sd) pnorm_detail::fast_dnorm((x),(mean),(sd))
#define LOG_DNORM(x,mean,sd) pnorm_detail::fast_log_dnorm((x),(mean),(sd))
#endif
