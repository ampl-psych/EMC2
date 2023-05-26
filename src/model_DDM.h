/* RFastDM.cpp - Main source file for the RCpp implementation of fast-dm
 *
 * Copyright (C) 2006  Jochen Voss, Andreas Voss.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301 USA.
 */

#ifndef DDM_h
#define DDM_h

#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <cfloat>
#include <iostream>
#include <sstream>

// Include all .hpp files
// Note: This is bad organisation, but Rcpp (and especially RStudio's "sourceCpp" don't seem to handle
//       projects with multiple .cpp files)


// Used by both PDF and CDF (also includes tuning params)
#include "Parameters_DDM.h"
#include "utility_functions.h"


// While not enforced, this is the global parameters Singleton
//   To be created and freed in the d_, p_ and r_ calls in RFastDM.cpp
Parameters *g_Params;

#define BOUNDARY_LOWER 0
#define BOUNDARY_UPPER 1

using namespace Rcpp;



#define EPSILON 1e-6

// Forward declarations
double g_minus (double t);
double g_plus  (double t);

static double integral_t0_g_minus (double t, Parameters *params);
static double integral_z_g_minus  (double t, Parameters *params);
static double integral_v_g_minus  (double t, double zr, Parameters *params);

static double g_minus_no_var     (double t, double a, double zr, double v);
static double g_minus_small_time (double t, double zr, int N);
static double g_minus_large_time (double t, double zr, int N);

// TODO: Make sure these function names are accurate
static double integrate_z_over_t  (Parameters *params, double a, double b, double step_width);
static double integrate_v_over_zr (Parameters *params, double a, double b, double t, double step_width);


// Main calls
// NumericVector density (NumericVector rts, int boundary)
// {
//   int length = rts.length();
//   NumericVector out(length);
//
//   if (boundary == 1) {// Calc upper
//     for (int i = 0; i < length; i++) {
//       out[i] =  g_plus(rts[i]);
//     }
//   }
//   else { for (int i = 0; i < length; i++) {// Calc lower
//     out[i] = -g_minus(rts[i]);
//     }
//   }
//   return out;
// }

double g_minus(double t)
{
  return integral_t0_g_minus (t - g_Params->t0 - 0.5*g_Params->d, g_Params);
}

double g_plus(double t)
{
  // Make a copy so we don't disturb our params
  // (?TODO: we could optimise the object creation out and just set them back after the call)
  Parameters new_params(*g_Params);
  new_params.zr = 1 - g_Params->zr;
  new_params.v = -g_Params->v;

  return integral_t0_g_minus (t - new_params.t0 + 0.5*new_params.d, &new_params);
}



static double integral_t0_g_minus (double t, Parameters *params)
{
  double res;

  if (params->st0 < params->TUNE_ST0_EPSILON) // 170501   was == 0)
  {
    res = integral_z_g_minus (t, params);
  }
  else
  {
    res = integrate_z_over_t(params,
                             t - .5*params->st0,
                             t + .5*params->st0, params->TUNE_INT_T0) / params->st0;
  }

  return res;
}


static double integral_z_g_minus (double t, Parameters *params)
{
  double res;

  if (t <= 0) return 0;

  if (params->szr < params->TUNE_SZ_EPSILON) // 170501   was == 0)
  {
    res = integral_v_g_minus (t, params->zr, params);
  }
  else
  {
    res = integrate_v_over_zr(params, params->zr - .5*params->szr, params->zr + .5*params->szr,
                              t, params->TUNE_INT_Z) / params->szr;
  }
  return res;
}


static double integral_v_g_minus (double t, double zr, Parameters *params)
{
  double a = params->a;
  double v = params->v;
  double sv = params->sv;

  int N_small, N_large;
  double simple, factor, eps;

  double ta = t/(a*a);

  factor = 1 / (a*a * sqrt(t * sv*sv + 1)) * exp(-0.5 * (v*v*t + 2*v*a*zr - a*a * zr*zr * sv*sv) / (t*sv*sv+1));

  if (std::isinf(factor))
  {
    return 0;
  }

  eps = EPSILON / factor;

  if (params->sv == 0)
  {
    return g_minus_no_var(t, a, zr, v);
  }

  N_large = (int)ceil(1 / (M_PI*sqrt(t)));
  if (M_PI*ta*eps < 1)
  {
    N_large = std::max(N_large, (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
  }

  if (2*sqrt(2*M_PI*ta)*eps < 1)
  {
    N_small = (int)ceil(fmax(sqrt(ta) + 1, 2 + sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
  }
  else
  {
    N_small = 2;
  }

  if (N_small < N_large)
  {
    simple = g_minus_small_time(t/(a*a), zr, N_small);
  }
  else
  {
    simple = g_minus_large_time(t/(a*a), zr, N_large);
  }
  return factor * simple;
}


static double g_minus_no_var(double t, double a, double zr, double v)
{
  int N_small, N_large;
  double simple, factor, eps;
  double ta = t/(a*a);

  factor = exp(-a*zr*v - 0.5*v*v*t) / (a*a);
  if (std::isinf(factor)) { return 0; }

  eps = EPSILON / factor;

  N_large = (int)ceil (1/ (M_PI*sqrt(t)));
  if (M_PI*ta*eps < 1)
  {
    N_large = std::max(N_large, (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
  }

  if (2*sqrt(2*M_PI*ta)*eps < 1)
  {
    N_small = (int)ceil(fmax(sqrt(ta) + 1, 2 + sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
  }
  else
  {
    N_small = 2;
  }

  if (N_small < N_large)
  {
    simple = g_minus_small_time(t/(a*a), zr, N_small);
  }
  else
  {
    simple = g_minus_large_time(t/(a*a), zr, N_large);
  }
  return factor * simple;
}


static double g_minus_small_time(double t, double zr, int N)
{
  int i;
  double sum = 0;
  double d;

  for(i = -N/2; i <= N/2; i++)
  {
    d = 2*i + zr;
    sum += exp(-d*d / (2*t)) * d;
  }

  return sum / sqrt(2*M_PI*t*t*t);
}

static double g_minus_large_time(double t, double zr, int N)
{
  int i;
  double sum = 0;
  double d;

  for(i = 1; i <= N; i++)
  {
    d = i * M_PI;
    sum += exp(-0.5 * d*d * t) * sin(d*zr) * i;
  }

  return sum * M_PI;
}

// CONVERSION NOTE: Simplest way to deal with the integrate function is to remove
//                  the clever recursiveness and instead (ugh) duplicate code
static double integrate_z_over_t (Parameters *params, double a, double b, double step_width)
{
  double width = b-a;
  double tmp_N = width / step_width;
  if (std::isnan(tmp_N)) tmp_N = 20;
  int N = std::max(4, static_cast<int>(tmp_N));
  double step = std::max(width / N, EPSILON);
  double x;
  double result = 0;

  for(x = a+0.5*step; x < b; x += step)
  {
    result += step * integral_z_g_minus(x, params);
  }
  return result;
}

static double integrate_v_over_zr (Parameters *params, double a, double b, double t, double step_width)
{
  double width = b-a;
  double tmp_N = width / step_width;
  if (std::isnan(tmp_N)) tmp_N = 20;
  int N = std::max(4, static_cast<int>(tmp_N));
  double step = std::max(width / N, EPSILON);
  double x;
  double result = 0;

  for(x = a+0.5*step; x < b; x += step)
  {
    result += step * integral_v_g_minus (t, x, params);
  }
  return result;
}

// R-callable PDF for fastdm - pass boundary to retrieve (1 = lower, 2 = upper)
// [[Rcpp::export]]
NumericVector d_DDM_c (NumericVector rts, NumericVector R, List group_idx, NumericMatrix pars, double precision=3, bool stop_on_error=true)
{
  //0 = "v", 1 = "a", 2= "sv", 3 = "t0", 4 = "st0", 5 = "s", 6 = "Z", 7 = "SZ", 8 = "DP"),
  int length = rts.length();
  NumericVector out(length, 0.0);
  NumericVector pars_tmp(8);
  for(int k = 0; k < group_idx.length(); k ++){
    NumericVector total_idx = group_idx[k];
    int first_idx = total_idx[0] - 1;
    if((pars(first_idx,0) > 20 || pars(first_idx,1) > 10 || pars(first_idx,2) > 10 || pars(first_idx, 7) > .999 || pars(first_idx,4) > .2)){
      continue;
    } else if((pars(first_idx,2) != 0 && pars(first_idx,2) < 0.001)){
      continue;
    } else if(pars(first_idx,7) != 0 && pars(first_idx,7) < 0.001){
      continue;
    } else{
      // now we create a new vector with:  a/s, v/s, t0, d, sz/a, sv/s, st0, z/a
      pars_tmp[0] = pars(first_idx,1)/pars(first_idx,5); //a/s
      pars_tmp[1] = pars(first_idx,0)/pars(first_idx,5); //v/s
      pars_tmp[2] = pars(first_idx,3) + pars(first_idx,4)/2; //t0 = t0 + st0/2
      pars_tmp[3] = pars(first_idx,3) *(2*pars(first_idx,8)-1); //d = t0 * 2*DP -1
      if(pars(first_idx,6) < (1 - pars(first_idx,6))){
        pars_tmp[4] = 2*pars(first_idx,7) * pars(first_idx,6); //sz = 2*SZ * Z
      } else{
        pars_tmp[4] = 2*pars(first_idx,7) * (1 - pars(first_idx,6));
      }
      pars_tmp[5] = pars(first_idx,2)/pars(first_idx,5); //sv/s
      pars_tmp[6] = pars(first_idx,4); // st0
      pars_tmp[7] = pars(first_idx,6);
      // z Note that for z and sz the EMC R code multiplies by a, but the R ddiffusion divides by a, so avoiding both
      g_Params = new Parameters (pars_tmp, precision);
      bool is_bad = is_true(any(is_infinite(pars_tmp) == TRUE));
      if(!is_bad){
        if(R[first_idx] == 2){
          for (int i = 0; i < total_idx.length(); i++) {
            int idx = total_idx[i] - 1;
            out[idx] =  g_plus(rts[idx]); //upper
          }
        } else{
          for (int i = 0; i < total_idx.length(); i++) {
            int idx = total_idx[i] -1;
            out[idx] = -g_minus(rts[idx]); //lower
          }
        }
      }
    }
    delete g_Params;
  }
  return abs(out);
}


NumericMatrix Ntransform_DDM(NumericMatrix x) {
  NumericMatrix out(clone(x));
  CharacterVector is_log = {"a","sv","t0","st0","s"};
  CharacterVector is_probit = {"Z","SZ","DP"};
  LogicalVector col_idx_log = contains_multiple(colnames(x), is_log);
  LogicalVector col_idx_probit = contains_multiple(colnames(x), is_probit);
  for(int i = 0; i < x.ncol(); i ++){
    if(col_idx_log[i] == TRUE){
      out (_, i) = exp(out(_, i));
    };
    if(col_idx_probit[i] == TRUE){
      out (_, i) = pnorm_multiple(out(_, i));
    };
  };
  return(out);
}

NumericVector transform_DDM(NumericVector x){
  return(x);
}



#endif


