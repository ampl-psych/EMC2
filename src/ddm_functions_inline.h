#ifndef DDM_FUNCTIONS_INLINE_H
#define DDM_FUNCTIONS_INLINE_H

#include <cmath>
#include <algorithm>
#include "utility_functions.h"
#include "tools.h"

/* DENSITY HELPERS */

inline double ddm_ks(double t, double w, double eps) {
    double K1 = (std::sqrt(2.0 * t) + w) / 2.0;
    double u_eps = std::fmin(-1.0, M_LN2 + M_LNPI + 2.0 * std::log(t) + 2.0 * (eps));
    double arg = -t * (u_eps - std::sqrt(-2.0 * u_eps - 2.0));
    double K2 = (arg > 0) ? 0.5 * (std::sqrt(arg) - w) : K1;
    return std::ceil(std::fmax(K1, K2));
}

inline double ddm_kl(double q, double v, double w, double err) {
    double K1 = 1.0 / (M_PI * std::sqrt(q)), K2 = 0.0;
    double temp = -2.0 * (std::log(M_PI * q) + err);
    if (temp >= 0) K2 = std::sqrt(temp / (M_PI * M_PI * q));
    return std::ceil(std::fmax(K1, K2));
}

inline double ddm_logfs_linear(double t, double w, int K) {
    if (w == 0.0) return R_NegInf;
    if (K < 0) K = 0;

    const double inv_twot = 1.0 / (2.0 * t);
    double sum_plus = 0.0;
    double sum_minus = 0.0;

    if (K > 0) {
        #pragma omp simd reduction(+:sum_plus,sum_minus)
        for (int k = 1; k <= K; ++k) {
            const double temp1 = w + 2.0 * k;
            const double temp2 = w - 2.0 * k;
            sum_plus += temp1 * std::exp(-temp1 * temp1 * inv_twot);
            sum_minus += (-temp2) * std::exp(-temp2 * temp2 * inv_twot);
        }
    }

    sum_plus += w * std::exp(-w * w * inv_twot);
    const double sum = sum_plus - sum_minus;
    if (!(sum > 0.0)) return R_NegInf;

    return -0.5 * M_LN2 - M_LN_SQRT_PI - 1.5 * std::log(t) + std::log(sum);
}

inline double ddm_logfl_linear(double q, double v, double w, int K) {
    (void)v;
    if (w == 0.0) return R_NegInf;
    if (K <= 0) return R_NegInf;

    const double halfq = q / 2.0;
    double sum_plus = 0.0;
    double sum_minus = 0.0;

    #pragma omp simd reduction(+:sum_plus,sum_minus)
    for (int k = 1; k <= K; ++k) {
        const double temp = k * M_PI;
        const double base = k * std::sin(temp * w) * std::exp(-(temp * temp) * halfq);
        sum_plus += (base > 0.0) ? base : 0.0;
        sum_minus += (base < 0.0) ? -base : 0.0;
    }

    const double sum = sum_plus - sum_minus;
    if (!(sum > 0.0)) return R_NegInf;

    return std::log(sum) + M_LNPI;
}

/* DISTRIBUTION HELPERS */

inline double ddm_logP(int pm, double a, double v, double w) {
    const double em1 = 1.0 - 1.0e-6;
    if (pm == 1) { v = -v; w = 1.0 - w; }
    if (std::abs(v) == 0.0) return std::log1p(-w);
    
    double prob;
    double e = (-2.0 * v * a * (1.0 - w));
    if (e < 0) {
        double tt = std::exp(e);
        if (tt >= em1) return std::log1p(-w);
        tt = std::log1p(-tt) - logdiff(2 * v * a * w, e);
        prob = tt;
    }
    else {
        double tt = std::exp(-e);
        if (tt >= em1) return std::log1p(-w);
        tt = std::log1p(-tt) - std::log1p(-std::exp(2 * v * a));
        prob = tt;
    }
    return prob;
}

inline double ddm_Ks(double t, double v, double a, double w, double eps) {
    double K1 = 0.5 * (std::abs(v) / a * t - w);
    double arg = std::fmax(0, std::fmin(1, std::exp(v * a * w + v * v * t / 2 + (eps)) / 2));
    double K2 = (arg == 0) ? INFINITY : (arg == 1) ? -INFINITY : -std::sqrt(t) / 2 / a * gsl_cdf_ugaussian_Pinv(arg);
    return std::ceil(std::fmax(K1, K1 + K2));
}

inline double ddm_Kl(double t, double v, double a, double w, double err) {
    double api = a / M_PI, vsq = v * v;
    double sqrtL1 = std::sqrt(1 / t) * api;
    double sqrtL2 = std::sqrt(std::fmax(1.0, -2 / t * api * api * (err + std::log(M_PI * t / 2 * (vsq + (M_PI / a) * (M_PI / a))) + v * a * w + vsq * t / 2)));
    return std::ceil(std::fmax(sqrtL1, sqrtL2));
}

inline double ddm_logFs(double t, double v, double a, double w, int K) {
    double fplus = R_NegInf, fminus = R_NegInf;
    double sqt = std::sqrt(t), temp = -v * a * w - v * v * t / 2;
    double vt = v * t;

    for (int k = K; k >= 0; k--) {
        double rj = a * (2 * k + w);
        double dj = lognormal(rj / sqt);
        double pos1 = dj + logMill((rj - vt) / sqt);
        double pos2 = dj + logMill((rj + vt) / sqt);
        fplus = logsum(logsum(pos1, pos2), fplus);
        rj = a * (2.0 * k + 2.0 - w);
        dj = lognormal(rj / sqt);
        double neg1 = dj + logMill((rj - vt) / sqt);
        double neg2 = dj + logMill((rj + vt) / sqt);
        fminus = logsum(logsum(neg1, neg2), fminus);
    }

    return logsum(fplus, -fminus) + temp; // Original used logdiff, but that expects log(A) and log(B)
    // Actually, looking at original logFs: return logdiff(fplus, fminus)+temp;
    // My manual translation was a bit loose. Let's stick to the original names and logic.
}

// Re-translating more faithfully:
inline double ddm_logFs_faithful(double t, double v, double a, double w, int K) {
    double fplus = R_NegInf, fminus = R_NegInf;
    double sqt = std::sqrt(t), temp = -v * a * w - v * v * t / 2;
    double vt = v * t;

    for (int k = K; k >= 0; k--) {
        double rj = a * (2 * k + w);
        double dj = lognormal(rj / sqt);
        double pos1 = dj + logMill((rj - vt) / sqt);
        double pos2 = dj + logMill((rj + vt) / sqt);
        fplus = logsum(logsum(pos1, pos2), fplus);
        rj = a * (2.0 * k + 2.0 - w);
        dj = lognormal(rj / sqt);
        double neg1 = dj + logMill((rj - vt) / sqt);
        double neg2 = dj + logMill((rj + vt) / sqt);
        fminus = logsum(logsum(neg1, neg2), fminus);
    }

    return logdiff(fplus, fminus) + temp;
}

inline double ddm_logFl(double q, double v, double a, double w, int K) {
    double fplus = R_NegInf, fminus = R_NegInf;
    double la = std::log(a), lv = std::log(std::abs(v));
    for (int k = K; k >= 1; k--) {
        double temp0 = std::log(k * 1.0), temp1 = k * M_PI, temp2 = temp1 * w;
        double check = std::sin(temp2);
        if (check > 0) {
            double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * (temp1 / a) * (temp1 / a) * q + std::log(check);
            fplus = logsum(temp, fplus);
        }
        else if (check < 0) {
            double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * (temp1 / a) * (temp1 / a) * q + std::log(-check);
            fminus = logsum(temp, fminus);
        }
    }
    double F = logdiff(fplus, fminus);
    return (F - v * a * w - 0.5 * v * v * q);
}

/* MAIN INLINE FUNCTIONS */

inline double dwiener_inline(double q, double a, double vn, double wn, double sv, double err, int K, int epsFLAG) {
    if (q == 0.0) return R_NegInf;
    double kll, kss, ans, v, w;
    if (!epsFLAG && K == 0) {
        err = -27.63102;  // exp(err) = 1.e-12
        epsFLAG = 1;
    }
    else if (!epsFLAG && K > 0) err = -27.63102;
    else if (epsFLAG) err = std::log(err);

    if (q >= 0) { w = 1.0 - wn; v = -vn; }
    else { q = std::abs(q); w = wn; v = vn; }

    double q_asq = q / (a * a);
    double eta_sqr = sv * sv;
    double temp = 1 + eta_sqr * q;
    double lg1 = (eta_sqr * (a * w) * (a * w) - 2 * a * v * w - v * v * q) / 2.0 / temp - 2 * std::log(a) - 0.5 * std::log(temp);
    double es = (err - lg1);
    kss = ddm_ks(q_asq, w, es);
    kll = ddm_kl(q_asq, v, w, es);

    if (2 * kss <= kll) {
        if ((epsFLAG && kss < K) || !epsFLAG) kss = K;
        ans = lg1 + ddm_logfs_linear(q_asq, w, static_cast<int>(kss));
    }
    else {
        if ((epsFLAG && kll < K) || !epsFLAG) kll = K;
        ans = lg1 + ddm_logfl_linear(q_asq, v, w, static_cast<int>(kll));
    }
    return ans;
}

inline double pwiener_inline(double q, double a, double v, double w, double err, int K, int epsFLAG) {
    double Kll, Kss, ans;
    if (!epsFLAG && K == 0) {
        err = -27.63102; // exp(err) = 1.e-12
        epsFLAG = 1;
    }
    else if (!epsFLAG && K > 0) err = -27.63102;
    else if (epsFLAG) err = std::log(err);

    if (emc2_isinf(q)) return ddm_logP(0, a, v, w);

    Kss = ddm_Ks(q, v, a, w, err);
    Kll = ddm_Kl(q, v, a, w, err);
    double lg = M_LN2 + M_LNPI - 2.0 * std::log(a);

    if (3 * Kss < Kll) {
        if ((epsFLAG && Kss < K) || !epsFLAG) Kss = K;
        ans = ddm_logFs_faithful(q, v, a, w, static_cast<int>(Kss));
    }
    else {
        if ((epsFLAG && Kll < K) || !epsFLAG) Kll = K;
        ans = logdiff(ddm_logP(0, a, v, w), lg + ddm_logFl(q, v, a, w, static_cast<int>(Kll)));
    }
    return ans;
}

#endif
