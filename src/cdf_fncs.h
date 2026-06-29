
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#ifndef CDF_FNCS_H
#define CDF_FNCS_H



/* dependencies for DISTRIBUTION */
double logP(int, double, double, double);
double Ks(double, double, double, double, double);
double Kl(double, double, double, double, double);
double logFs(double, double, double, double, int);
double logFl(double, double, double, double, int);
/*-----------------------------------------------*/

/* DISTRIBUTION */
double pwiener(double, double, double, double, double, int, int);
/*-----------------------------------------------*/

#endif
