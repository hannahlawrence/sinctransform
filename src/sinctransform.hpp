#ifndef __sinctransform__hpp__included__
#define __sinctransform__hpp__included__

int sinc1d(int ifl,int numlocs,double *klocs, complex<double> *q,double tol,complex<double> *res);
int sincsq1d(int ifl,int numlocs,double *klocs, complex<double> *q,double tol,complex<double> *res);

int sinc2d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,complex<double> *q,double tol, complex<double> *res);
int sincsq2d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,complex<double> *q,double tol,complex<double> *res);

int sinc3d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3,complex<double> *q,double tol, complex<double> *res);
int sincsq3d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3,complex<double> *q,double tol,complex<double> *res);

#endif