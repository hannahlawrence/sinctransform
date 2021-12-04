#ifndef __sinctransform__hpp__included__
#define __sinctransform__hpp__included__

int sinc1d(int ifl,int numlocs,double *a1,double *klocs, complex<double> *q,double tol,complex<double> *res,int quad);
int sincsq1d(int ifl,int numlocs,double *a1,double *klocs, complex<double> *q,double tol,complex<double> *res,int quad);

int sinc2d(int ifl,int numlocs,double *a1,double *a2,double *klocs_d1,double *klocs_d2,complex<double> *q,double tol, complex<double> *res,int quad);
int sincsq2d(int ifl,int numlocs,double *a1,double *a2,double *klocs_d1,double *klocs_d2,complex<double> *q,double tol,complex<double> *res,int quad);

int sinc3d(int ifl,int numlocs,double *a1,double *a2,double *a3,double *klocs_d1,double *klocs_d2,double *klocs_d3,complex<double> *q,double tol, complex<double> *res,int quad);
int sincsq3d(int ifl,int numlocs,double *a1,double *a2,double *a3,double *klocs_d1,double *klocs_d2,double *klocs_d3,complex<double> *q,double tol,complex<double> *res,int quad);

#endif