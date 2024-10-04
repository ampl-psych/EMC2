
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#ifndef PDF_FNCS_H
#define PDF_FNCS_H



/* dependencies for DENSITY */
double ks(double, double, double);
double kl(double, double, double, double);
double logfs(double, double, int);
double logfl(double, double, double, int);
/*-----------------------------------------------*/

/* dependencies for d/dt and d/da DENSITY */
double dtaks(double, double, double);
double dtakl(double, double, double, double);
void logdtfs(double, double, int, double &, int &);
void logdtfl(double, double, int, double &, int &);
/*-----------------------------------------------*/

/* dependencies for d/dw DENSITY */
double dwks(double, double, double);
double dwkl(double, double, double);
void logdwfs(double, double, int, double &, int &);
void logdwfl(double, double, double, int, double &, int &);
/*-----------------------------------------------*/

/* dependencies for d/dt d/da d/dv d/dw DENSITY */
void dxks(double, double, double, double, double, double &, double &);
void dxkl(double, double, double, double, double, double &, double &);
void logdxfs(double, double, int, int, double &, double &, int &, int &);
void logdxfl(double, double, int, int, double &, double &, int &, int &);
/*-----------------------------------------------*/

/* DENSITY */
double dwiener(double, double, double, double, double, double, int, int);
/*-----------------------------------------------*/

/* d/dt DENSITY */
void dtdwiener(double, double, double, double, double, double, double *, double, int, int);
/*-----------------------------------------------*/

/* d/da DENSITY */
void dadwiener(double, double, double, double, double, double, double *, double, int, int);
/*-----------------------------------------------*/

/* d/dv DENSITY */
void dvdwiener(double, double, double, double, double, double, double *);
/*-----------------------------------------------*/

/* d/dw DENSITY */
void dwdwiener(double, double, double, double, double, double, double *, double, int, int);
/*-----------------------------------------------*/

/* d/dsv DENSITY */
void dsvdwiener(double, double, double, double, double, double, double *, double, int, int);
/*-----------------------------------------------*/

/* d/dt d/da d/dv d/dw DENSITY */
void dxdwiener(double, double, double, double, double, double, int, int, double *, double *, double *);
/*-----------------------------------------------*/



#endif
