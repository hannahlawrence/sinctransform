#ifndef __sinctransform__hpp__included__
#define __sinctransform__hpp__included__

int sinc1d(int ifl,int numlocs,double *klocs, double *q,double tol,double *res);
int sincsq1d(int ifl,int numlocs,double *klocs, double *q,double tol, double *res);

int sinc2d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *q,double tol, double *res);
int sincsq2d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *q,double tol,double *res);

int sinc3d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3,double *q,double tol, double *res);
int sincsq3d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3,double *q,double tol,double *res);

#endif