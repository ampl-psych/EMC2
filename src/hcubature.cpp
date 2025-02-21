// /*
//  * For the code in hcubature.cpp., the package JuliaMath/HCubature.jl written
//  * in Julia by Steven G. Johnson served as a template. It comes with the
//  * following MIT "Expat" license:
//  *
//  * Copyright (c) 2017: Steven G. Johnson.
//  *
//  * Permission is hereby granted, free of charge, to any person obtaining a copy
//  * of this software and associated documentation files (the "Software"), to deal
//  * in the Software without restriction, including without limitation the rights
//  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  * copies of the Software, and to permit persons to whom the Software is
//  * furnished to do so, subject to the following conditions:
//  *
//  * The above copyright notice and this permission notice shall be included in
//  * all copies or substantial portions of the Software.
//  *
//  */
//
// /*
//  * The JuliaMath/HCubature.jl was translated to C/C++ by
//  * Prof. Dr. Karl Christoph Klauer
//  */
//
//
// #include "gauss.h"
// #include "tools.h"
// #include <set>
//
// int choose(int n, int k) {
//   if (k > n) return 0;
//   if (k * 2 > n) k = n - k;
//   if (k == 0) return 1;
//
//   int result = n;
//   for (int i = 2; i <= k; ++i) {
//     result *= (n - i + 1);
//     result /= i;
//   }
//   return result;
// }
//
// /** [combination c n p x]
//  * get the [x]th lexicographically ordered set of [p] elements in [n]
//  * output is in [c], and should be sizeof(int)*[p]
//  * "Algorithm 515: Generation of a Vector from the Lexicographical Index"; Buckles, B. P., and Lybanon, M. ACM Transactions on Mathematical Software, Vol. 3, No. 2, June 1977.
//  * User lucaroni from https://stackoverflow.com/questions/561/how-to-use-combinations-of-sets-as-test-data#794
//  */
//
// void combination(int* c, int n, int p, int x) {
//   int i, r, k = 0;
//   for (i = 0; i < p - 1; i++) {
//     c[i] = (i != 0) ? c[i - 1] : 0;
//     do {
//       c[i]++;
//       r = choose(n - c[i], p - (i + 1));
//       k = k + r;
//     } while (k < x);
//     k = k - r;
//   }
//   if (p > 1) c[p - 1] = c[p - 2] + x - k; else c[0] = x;
// }
//
// void combos(int k, double lambda, int n, vector<vector<double>>& p) {
//   int* c = (int*)malloc(k * sizeof(int));
//
//   for (int i = 1; i != choose(n, k) + 1; i++) {
//     vector<double> temp;
//     for (int x = 0; x != n; x++) temp.push_back(0.0);
//     combination(c, n, k, i);
//     for (int j = 0; j != k; j++) temp[c[j] - 1] = lambda;
//     p.push_back(temp);
//   }
//   free(c);
// }
//
// void increment(vector<bool>& index, int k, double lambda, int n, int* c, vector<double>& temp) {
//   // temp size n, all elements initially zero
//   if (index.size() == 0) {
//     index.push_back(false);
//     for (int j = 0; j != k; j++) temp[c[j] - 1] = lambda;
//     return;
//   }
//   int first_zero = 0;
//   while ((first_zero < index.size()) && index[first_zero]) first_zero++;
//   if (first_zero == index.size()) {
//     index.flip();
//     for (int j = 0; j != index.size(); j++) temp[c[j] - 1] *= -1;
//     index.push_back(true);
//     temp[c[index.size() - 1] - 1] = -lambda;
//   }
//   else {
//     for (int i = 0; i != first_zero + 1; i++) {
//       index[i] = !index[i];
//       temp[c[i] - 1] *= -1;
//     }
//   }
// }
//
// void signcombos(int k, double lambda, int n, vector<vector<double>>& p) {
//   int* c = (int*)malloc(k * sizeof(int));
//   for (int i = 1; i != choose(n, k) + 1; i++) {
//     vector<double> temp;
//     for (int x = 0; x != n; x++) temp.push_back(0.0);
//     combination(c, n, k, i);
//     vector<bool> index; index.clear();
//     for (int j = 0; j != pow(2, k); j++) {
//       increment(index, k, lambda, n, c, temp);
//       p.push_back(temp);
//     }
//   }
//   free(c);
// }
//
// void gauss_kronrod(double a, double b, one_d& out, void* pars, int integrand(unsigned dim, const double* x, void* p,
//                                                                              unsigned fdim, double* retval)) {
//   double c = 0.5 * (a + b);
//   double delta = 0.5 * (b - a);
//   double f0;
//   integrand(1, &c, pars, 1, &f0);
//   double I = f0 * wd7[7], Idash = f0 * gwd7[3];
//   for (int i = 0; i != 7; i++) {
//     double deltax = delta * xd7[i], cp = c + deltax, cm = c - deltax;
//     double  fx;
//     integrand(1, &cp, pars, 1, &fx);
//     double temp;
//     integrand(1, &cm, pars, 1, &temp);
//     fx += temp;
//     I += fx * wd7[i];
//     if (i % 2 == 1) Idash += fx * gwd7[i / 2];
//   }
//   double V = fabs(delta);
//   I *= V;
//   Idash *= V;
//   out.result = I;
//   out.err = fabs(I - Idash);
// }
//
// void make_GenzMalik(int n, GenzMalik& g) {
//   double l4 = sqrt(9 * 1.0 / 10);
//   double l2 = sqrt(9 * 1.0 / 70);
//   double l3 = l4;
//   double l5 = sqrt(9 * 1.0 / 19);
//
//   int twopn = pow(2, n);
//
//   g.w[0] = twopn * ((12824 - 9120 * n + 400 * n * n) * 1.0 / 19683);
//   g.w[1] = twopn * (980.0 / 6561);
//   g.w[2] = twopn * ((1820 - 400 * n) * 1.0 / 19683);
//   g.w[3] = twopn * (200.0 / 19683);
//   g.w[4] = 6859.0 / 19683;
//   g.wd[3] = twopn * (25.0 / 729);
//   g.wd[2] = twopn * ((265 - 100 * n) * 1.0 / 1458);
//   g.wd[1] = twopn * (245.0 / 486);
//   g.wd[0] = twopn * ((729 - 950 * n + 50 * n * n) * 1.0 / 729);
//
//   combos(1, l2, n, g.p[0]);
//   combos(1, l3, n, g.p[1]);
//   signcombos(2, l4, n, g.p[2]);
//   signcombos(n, l5, n, g.p[3]);
// }
//
// void clean_GenzMalik(GenzMalik& g) {
//   for (int j = 0; j != 4; j++)
//     for (int i = 0; i != g.p[j].size(); i++) g.p[j][i].clear();
// }
//
// void integrate_GenzMalik(GenzMalik g, int n, const double* a, const double* b, one_d& out, void* pars, int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval)) {
//   double* c = (double*)malloc(n * sizeof(double));
//   double* deltac = (double*)malloc(n * sizeof(double));
//
//   for (int i = 0; i != n; i++) c[i] = (a[i] + b[i]) / 2;
//   for (int i = 0; i != n; i++) deltac[i] = fabs(b[i] - a[i]) / 2;
//   double v = 1.0;
//   for (int i = 0; i != n; i++) v *= deltac[i];
//
//   if (v == 0.0) {
//     out.err = 0.0;
//     out.result = 0.0;
//     out.kdivide = 0;
//     return;
//   }
//
//   double f1;
//   integrand(n, c, pars, 1, &f1);
//   double f2 = 0.0, f3 = 0.0;
//   double twelvef1 = 12 * f1;
//
//   double maxdivdiff = 0.0;
//   double* divdiff = (double*)malloc(n * sizeof(double));
//   double* p2 = (double*)malloc(n * sizeof(double));
//   double* p3 = (double*)malloc(n * sizeof(double));
//   double* cc = (double*)malloc(n * sizeof(double));
//
//   for (int i = 0; i != n; i++) {
//
//     for (int j = 0; j != n; j++) p2[j] = deltac[j] * g.p[0][i][j];
//
//     for (int j = 0; j != n; j++) cc[j] = c[j] + p2[j];
//     double f2i;
//     integrand(n, cc, pars, 1, &f2i);
//     for (int j = 0; j != n; j++) cc[j] = c[j] - p2[j];
//     double temp;
//     integrand(n, cc, pars, 1, &temp);
//     f2i += temp;
//
//
//     for (int j = 0; j != n; j++) p3[j] = deltac[j] * g.p[1][i][j];
//     for (int j = 0; j != n; j++) cc[j] = c[j] + p3[j];
//     double f3i;
//     integrand(n, cc, pars, 1, &f3i);
//     for (int j = 0; j != n; j++) cc[j] = c[j] - p3[j];
//     integrand(n, cc, pars, 1, &temp);
//     f3i += temp;
//     f2 += f2i;
//     f3 += f3i;
//     divdiff[i] = fabs(f3i + twelvef1 - 7 * f2i);
//
//   }
//   free(p2); free(p3);
//   double* p4 = (double*)malloc(n * sizeof(double));
//   double f4 = 0.0;
//   for (int i = 0; i != g.p[2].size(); i++) {
//
//     for (int j = 0; j != n; j++) p4[j] = deltac[j] * g.p[2][i][j];
//     for (int j = 0; j != n; j++) cc[j] = c[j] + p4[j];
//     double temp;
//     integrand(n, cc, pars, 1, &temp);
//     f4 += temp;
//   }
//   free(p4);
//   double f5 = 0.0;
//   double* p5 = (double*)malloc(n * sizeof(double));
//   for (int i = 0; i != g.p[3].size(); i++) {
//
//     for (int j = 0; j != n; j++) p5[j] = deltac[j] * g.p[3][i][j];
//
//     for (int j = 0; j != n; j++) cc[j] = c[j] + p5[j];
//     double temp;
//     integrand(n, cc, pars, 1, &temp);
//     f5 += temp;
//   }
//   free(p5); free(cc);
//   double I = v * (g.w[0] * f1 + g.w[1] * f2 + g.w[2] * f3 + g.w[3] * f4 + g.w[4] * f5);
//   double Idash = v * (g.wd[0] * f1 + g.wd[1] * f2 + g.wd[2] * f3 + g.wd[3] * f4);
//   double E = fabs(I - Idash);
// #
//   int kdivide = 0;
//   double deltaf = E / (pow(10, n) * v);
//   for (int i = 0; i != n; i++) {
//     double delta = divdiff[i] - maxdivdiff;
//     if (delta > deltaf) {
//       kdivide = i;
//       maxdivdiff = divdiff[i];
//     }
//     else if ((fabs(delta) <= deltaf) && (deltac[i] > deltac[kdivide])) kdivide = i;
//   }
//   out.result = I;
//   out.err = E;
//   out.kdivide = kdivide;
//   free(c); free(deltac); free(divdiff);
// }
//
// class Box {
// public:
//   Box(double* a, double* b, double I, double err, int kdivide) : a(a), b(b), I(I), E(err), kdiv(kdivide) {};
//   bool operator<(const Box& box) const { return E > box.E; }
//   double* a;
//   double* b;
//   double I;
//   double E;
//   int kdiv;
// };
//
// class Box make_box(int n, const double* a, const double* b, one_d out) {
//   double* ac = (double*)malloc(n * sizeof(double));
//   double* bc = (double*)malloc(n * sizeof(double));
//   memcpy(ac, a, n * sizeof(double));
//   memcpy(bc, b, n * sizeof(double));
//   Box box(ac, bc, out.result, out.err, out.kdivide);
//   return box;
// }
//
// void delete_box(Box& box) {
//   if (box.a) free(box.a);
//   if (box.b) free(box.b);
// }
//
// int hcubature(int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval), void* pars, unsigned n, const double* a, const double* b,
//               size_t maxEval, double reqAbsError, double reqRelError, double* val, double* err) {
//
//   one_d out;
//   GenzMalik g;
//
//   if (n == 1) gauss_kronrod(a[0], b[0], out, pars, integrand);
//   else {
//     make_GenzMalik(n, g);
//     integrate_GenzMalik(g, n, a, b, out, pars, integrand);
//   }
//   int numevals = (n == 1) ? 15 : 1 + 4 * n + 2 * n * (n - 1) + pow(2, n);
//   int evals_per_box = numevals;
//   int kdiv = out.kdivide;
//   err[0] = out.err;
//   val[0] = out.result;
//   // convergence test
//   if ((err[0] <= std::max(reqRelError * fabs(val[0]), reqAbsError)) || ((maxEval!=0) && (numevals >= maxEval))) {
//     //        std::cout << numevals << std::endl;
//     return 0;
//   }
//
//   std::multiset<Box> ms;
//   ms.insert(make_box(n, a, b, out));
//
//   while (true) {
//     std::multiset<Box>::iterator it;
//     it = ms.begin();
//     Box box = *it;
//     ms.erase(it);
//     // split along dimension kdiv
//     double w = (box.b[box.kdiv] - box.a[box.kdiv]) / 2;
//     double* ma = (double*)malloc(n * sizeof(double));
//     memcpy(ma, box.a, n * sizeof(double));
//     ma[box.kdiv] += w;
//     double* mb = (double*)malloc(n * sizeof(double));
//     memcpy(mb, box.b, n * sizeof(double));
//     mb[box.kdiv] -= w;
//
//     if (n == 1) gauss_kronrod(ma[0], box.b[0], out, pars, integrand);
//     else {
//       integrate_GenzMalik(g, n, ma, box.b, out, pars, integrand);
//     }
//     Box box1 = make_box(n, ma, box.b, out);
//     ms.insert(box1);
//
//     if (n == 1) gauss_kronrod(box.a[0], mb[0], out, pars, integrand);
//     else {
//       integrate_GenzMalik(g, n, box.a, mb, out, pars, integrand);
//     }
//     Box box2 = make_box(n, box.a, mb, out);
//     ms.insert(box2);
//     val[0] += box1.I + box2.I - box.I;
//     err[0] += box1.E + box2.E - box.E;
//     numevals += 2 * evals_per_box;
//     delete_box(box);
//     free(ma); free(mb);
//     if (((err[0] <= std::max(reqRelError * fabs(val[0]), reqAbsError)) || ((maxEval != 0) && (numevals >= maxEval))) || !(std::isfinite(val[0])) ) {
//       break;
//     }
//   }
//   val[0] = 0.0;
//   err[0] = 0.0;
//
//   for (std::multiset<Box>::reverse_iterator rit = ms.rbegin(); rit != ms.rend(); rit++) {
//     val[0] += (*rit).I;
//     err[0] += (*rit).E;
//   }
//
//   for (std::multiset<Box>::iterator it = ms.begin(); it != ms.end(); it++) {
//     Box box = *it;
//     delete_box((box));
//   }
//   clean_GenzMalik(g);
//   return 0;
// }


#include "gauss.h"
#include <queue>
#include <cmath>
#include <cstring>

//multivariate integration routines

int choose(int n, int k) {
  if (k > n) return 0;
  if (k * 2 > n) k = n - k;
  if (k == 0) return 1;

  int result = n;
  for (int i = 2; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

/** [combination c n p x]
 * get the [x]th lexicographically ordered set of [p] elements in [n]
 * output is in [c], and should be sizeof(int)*[p]
 * "Algorithm 515: Generation of a Vector from the Lexicographical Index"; Buckles, B. P., and Lybanon, M. ACM Transactions on Mathematical Software, Vol. 3, No. 2, June 1977.
 * User lucaroni from https://stackoverflow.com/questions/561/how-to-use-combinations-of-sets-as-test-data#794
 */

void combination(int* c, int n, int p, int x) {
  int i, r, k = 0;
  for (i = 0; i < p - 1; i++) {
    c[i] = (i != 0) ? c[i - 1] : 0;
    do {
      c[i]++;
      r = choose(n - c[i], p - (i + 1));
      k = k + r;
    } while (k < x);
    k = k - r;
  }
  if (p > 1) c[p - 1] = c[p - 2] + x - k; else c[0] = x;
}

void combos(int k, double lambda, int n, vector<vector<double>>& p) {
  int* c = (int*)malloc(k * sizeof(int));
  int cnk = choose(n, k) + 1;
  for (int i = 1; i != cnk; i++) {
    vector<double> temp(n, 0.0);
    combination(c, n, k, i);
    for (int j = 0; j != k; j++) temp[c[j] - 1] = lambda;
    p.push_back(temp);
  }
  free(c);
}

void increment(vector<bool>& index, int k, double lambda, int n, int* c, vector<double>& temp) {
  // temp size n, all elements initially zero
  if (index.size() == 0) {
    index.push_back(false);
    for (int j = 0; j != k; j++) temp[c[j] - 1] = lambda;
    return;
  }
  int first_zero = 0;
  while ((first_zero < static_cast<int>(index.size())) && index[first_zero]) first_zero++;
  if (first_zero == static_cast<int>(index.size())) {
    index.flip();
    for (int j = 0; j != static_cast<int>(index.size()); j++) temp[c[j] - 1] *= -1;
    index.push_back(true);
    temp[c[index.size() - 1] - 1] = -lambda;
  }
  else {
    int fzp1 = first_zero + 1;
    for (int i = 0; i != fzp1; i++) {
      index[i] = !index[i];
      temp[c[i] - 1] *= -1;
    }
  }
}

void signcombos(int k, double lambda, int n, vector<vector<double>>& p) {
  int* c = (int*)malloc(k * sizeof(int));
  int cnk = choose(n, k) + 1;
  for (int i = 1; i != cnk; i++) {
    vector<double> temp(n, 0.0);
    combination(c, n, k, i);
    vector<bool> index; index.clear();
    int p2k = pow(2, k);
    for (int j = 0; j != p2k; j++) {
      increment(index, k, lambda, n, c, temp);
      p.push_back(temp);
    }
  }
  free(c);
}

void gauss_kronrod(double a, double b, one_d& out, void* pars, int integrand(unsigned dim, const double* x, void* p,
                                                                             unsigned fdim, double* retval)) {
  double c = 0.5 * (a + b);
  double delta = 0.5 * (b - a);
  double f0;
  integrand(1, &c, pars, 1, &f0);
  double I = f0 * wd7[7], Idash = f0 * gwd7[3];
  for (int i = 0; i != 7; i++) {
    double deltax = delta * xd7[i], cp = c + deltax, cm = c - deltax;
    double  fx;
    integrand(1, &cp, pars, 1, &fx);
    double temp;
    integrand(1, &cm, pars, 1, &temp);
    fx += temp;
    I += fx * wd7[i];
    if (i % 2 == 1) Idash += fx * gwd7[i / 2];
  }
  double V = fabs(delta);
  I *= V;
  Idash *= V;
  out.result = I;
  out.err = fabs(I - Idash);
}

void make_GenzMalik(int n, GenzMalik& g) {
  double l4 = sqrt(9 * 1.0 / 10);
  double l2 = sqrt(9 * 1.0 / 70);
  double l3 = l4;
  double l5 = sqrt(9 * 1.0 / 19);

  int twopn = pow(2, n);

  g.w[0] = twopn * ((12824 - 9120 * n + 400 * n * n) * 1.0 / 19683);
  g.w[1] = twopn * (980.0 / 6561);
  g.w[2] = twopn * ((1820 - 400 * n) * 1.0 / 19683);
  g.w[3] = twopn * (200.0 / 19683);
  g.w[4] = 6859.0 / 19683;
  g.wd[3] = twopn * (25.0 / 729);
  g.wd[2] = twopn * ((265 - 100 * n) * 1.0 / 1458);
  g.wd[1] = twopn * (245.0 / 486);
  g.wd[0] = twopn * ((729 - 950 * n + 50 * n * n) * 1.0 / 729);

  combos(1, l2, n, g.p[0]);
  combos(1, l3, n, g.p[1]);
  signcombos(2, l4, n, g.p[2]);
  signcombos(n, l5, n, g.p[3]);
}

void clean_GenzMalik(GenzMalik& g) {
  for (int j = 0; j != 4; j++) {
    int gpjs = g.p[j].size();
    for (int i = 0; i != gpjs; i++) g.p[j][i].clear();
  }
}

void integrate_GenzMalik(GenzMalik g, int n, const double* a, const double* b, one_d& out, void* pars, int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval)) {
  double* c = (double*)malloc(n * sizeof(double));
  double* deltac = (double*)malloc(n * sizeof(double));

  for (int i = 0; i != n; i++) c[i] = (a[i] + b[i]) / 2;
  for (int i = 0; i != n; i++) deltac[i] = fabs(b[i] - a[i]) / 2;
  double v = 1.0;
  for (int i = 0; i != n; i++) v *= deltac[i];

  if (v == 0.0) {
    out.err = 0.0;
    out.result = 0.0;
    out.kdivide = 0;
    return;
  }

  double f1;
  integrand(n, c, pars, 1, &f1);
  double f2 = 0.0, f3 = 0.0;
  double twelvef1 = 12 * f1;

  double maxdivdiff = 0.0;
  double* divdiff = (double*)malloc(n * sizeof(double));
  double* p2 = (double*)malloc(n * sizeof(double));
  double* p3 = (double*)malloc(n * sizeof(double));
  double* cc = (double*)malloc(n * sizeof(double));

  for (int i = 0; i != n; i++) {

    for (int j = 0; j != n; j++) p2[j] = deltac[j] * g.p[0][i][j];

    for (int j = 0; j != n; j++) cc[j] = c[j] + p2[j];
    double f2i;
    integrand(n, cc, pars, 1, &f2i);
    for (int j = 0; j != n; j++) cc[j] = c[j] - p2[j];
    double temp;
    integrand(n, cc, pars, 1, &temp);
    f2i += temp;


    for (int j = 0; j != n; j++) p3[j] = deltac[j] * g.p[1][i][j];
    for (int j = 0; j != n; j++) cc[j] = c[j] + p3[j];
    double f3i;
    integrand(n, cc, pars, 1, &f3i);
    for (int j = 0; j != n; j++) cc[j] = c[j] - p3[j];
    integrand(n, cc, pars, 1, &temp);
    f3i += temp;
    f2 += f2i;
    f3 += f3i;
    divdiff[i] = fabs(f3i + twelvef1 - 7 * f2i);

  }
  free(p2); free(p3);
  double* p4 = (double*)malloc(n * sizeof(double));
  double f4 = 0.0;
  int gp2s = g.p[2].size(), gp3s = g.p[3].size();
  for (int i = 0; i != gp2s; i++) {

    for (int j = 0; j != n; j++) p4[j] = deltac[j] * g.p[2][i][j];
    for (int j = 0; j != n; j++) cc[j] = c[j] + p4[j];
    double temp;
    integrand(n, cc, pars, 1, &temp);
    f4 += temp;
  }
  free(p4);
  double f5 = 0.0;
  double* p5 = (double*)malloc(n * sizeof(double));
  for (int i = 0; i != gp3s; i++) {

    for (int j = 0; j != n; j++) p5[j] = deltac[j] * g.p[3][i][j];

    for (int j = 0; j != n; j++) cc[j] = c[j] + p5[j];
    double temp;
    integrand(n, cc, pars, 1, &temp);
    f5 += temp;
  }
  free(p5); free(cc);
  double I = v * (g.w[0] * f1 + g.w[1] * f2 + g.w[2] * f3 + g.w[3] * f4 + g.w[4] * f5);
  double Idash = v * (g.wd[0] * f1 + g.wd[1] * f2 + g.wd[2] * f3 + g.wd[3] * f4);
  double E = fabs(I - Idash);
#
  int kdivide = 0;
  double deltaf = E / (pow(10, n) * v);
  for (int i = 0; i != n; i++) {
    double delta = divdiff[i] - maxdivdiff;
    if (delta > deltaf) {
      kdivide = i;
      maxdivdiff = divdiff[i];
    }
    else if ((fabs(delta) <= deltaf) && (deltac[i] > deltac[kdivide])) kdivide = i;
  }
  out.result = I;
  out.err = E;
  out.kdivide = kdivide;
  free(c); free(deltac); free(divdiff);
}

class Box {
public:
  Box(double* a, double* b, double I, double err, int kdivide) : a(a), b(b), I(I), E(err), kdiv(kdivide) {};
  bool operator<(const Box& box) const { return E < box.E; }
  double* a;
  double* b;
  double I;
  double E;
  int kdiv;
};

class Box make_box(int n, const double* a, const double* b, one_d out) {
  double* ac = (double*)malloc(n * sizeof(double));
  double* bc = (double*)malloc(n * sizeof(double));
  memcpy(ac, a, n * sizeof(double));
  memcpy(bc, b, n * sizeof(double));
  Box box(ac, bc, out.result, out.err, out.kdivide);
  return box;
}

void delete_box(Box& box) {
  if (box.a) free(box.a);
  if (box.b) free(box.b);
}

int hcubature(int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval), void* pars, unsigned n, const double* a, const double* b,
              size_t maxEval, double reqAbsError, double reqRelError, double* val, double* err) {

  one_d out;
  GenzMalik g;

  if (n == 1) gauss_kronrod(a[0], b[0], out, pars, integrand);
  else {
    make_GenzMalik(n, g);
    integrate_GenzMalik(g, n, a, b, out, pars, integrand);
  }
  int numevals = (n == 1) ? 15 : 1 + 4 * n + 2 * n * (n - 1) + pow(2, n);
  int evals_per_box = numevals;
  err[0] = out.err;
  val[0] = out.result;
  // convergence test
  if ((err[0] <= std::max(reqRelError * fabs(val[0]), reqAbsError)) || ((maxEval!=0) && (numevals >= static_cast<int>(maxEval)))) {
    //        std::cout << numevals << std::endl;
    return 0;
  }

  std::priority_queue<Box> ms;
  ms.push(make_box(n, a, b, out));

  while (true) {
    Box box = ms.top();
    ms.pop();
    // split along dimension kdiv
    double w = (box.b[box.kdiv] - box.a[box.kdiv]) / 2;
    double* ma = (double*)malloc(n * sizeof(double));
    memcpy(ma, box.a, n * sizeof(double));
    ma[box.kdiv] += w;
    double* mb = (double*)malloc(n * sizeof(double));
    memcpy(mb, box.b, n * sizeof(double));
    mb[box.kdiv] -= w;

    if (n == 1) gauss_kronrod(ma[0], box.b[0], out, pars, integrand);
    else {
      integrate_GenzMalik(g, n, ma, box.b, out, pars, integrand);
    }
    Box box1 = make_box(n, ma, box.b, out);
    ms.push(box1);

    if (n == 1) gauss_kronrod(box.a[0], mb[0], out, pars, integrand);
    else {
      integrate_GenzMalik(g, n, box.a, mb, out, pars, integrand);
    }
    Box box2 = make_box(n, box.a, mb, out);
    ms.push(box2);
    val[0] += box1.I + box2.I - box.I;
    err[0] += box1.E + box2.E - box.E;
    numevals += 2 * evals_per_box;
    delete_box(box);
    free(ma); free(mb);
    if (((err[0] <= std::max(reqRelError * fabs(val[0]), reqAbsError)) || ((maxEval != 0) && (numevals >= static_cast<int>(maxEval)))) || !(std::isfinite(val[0])) ) {
      break;
    }
  }
  val[0] = 0.0;
  err[0] = 0.0;

  while (!ms.empty()) {
    Box box = ms.top();
    val[0] += box.I;
    err[0] += box.E;
    delete_box(box);
    ms.pop();
  }
  clean_GenzMalik(g);
  return 0;
}
