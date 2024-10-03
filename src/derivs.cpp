
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include <cmath>
#include <thread>
#include "pdf_fncs.h"
#include "fncs_seven.h"
#include "tools.h"


/* PDF of Wiener diffusion */
void PDF(double *t, double *a, double *v, double *w, double *sv, double eps, int *resp, int K, int N, int epsFLAG, double *Rpdf, double *Rlogpdf) {
    for(int i = 0; i < N; i++) {
      if (i % 1024 == 0) R_CheckUserInterrupt();
      double pm = (resp[i]==1) ? 1.0 : -1.0;
      double ld = dwiener(t[i]*pm, a[i], v[i], w[i], sv[i], eps, K, epsFLAG);
      Rlogpdf[i] = ld;
      Rpdf[i] = exp(ld);
    }
}

/* PDF of 7-param diffusion */
void PDF7(int choice, double *t, int *resp, double *a, double *v, double *t0, double *w, double *sw, double *sv, double *st, double err, int K, int N, int epsFLAG, double *Rval, double *Rlogval, double *Rerr, int Neval) {
  for(int i = 0; i < N; i++) {
    if (i % 1024 == 0) R_CheckUserInterrupt();
    double low_or_up = (resp[i]==1) ? 1.0 : -1.0;
    Rerr[i] = 0;
    ddiff(choice, t[i], low_or_up, a[i], v[i], t0[i], w[i], sw[i], sv[i], st[i], err, K, epsFLAG, Neval, &Rval[i], &Rerr[i]);
    if (choice == 0) {
      Rlogval[i] = std::log(Rval[i]);
    }
  }
}
