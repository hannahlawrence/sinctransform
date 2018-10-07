#ifndef __sincdirect__hpp__included__
#define __sincdirect__hpp__included__

void directsinc1d(int ifl,int numlocs,int numeval,double *a1,double *klocs, complex<double>* q,complex<double> *ans);
void directsincsq1d(int ifl, int numlocs,int numeval,double *a1,double *klocs, complex<double>* q,complex<double> *ans);
void directsinc2d(int ifl,int numlocs,int numeval,double *a1,double *a2,double *klocs_d1,double *klocs_d2, complex<double>* q, complex<double> *ans);
void directsincsq2d(int ifl, int numlocs,int numeval,double *a1,double *a2,double *klocs_d1,double *klocs_d2, complex<double>* q, complex<double> *ans);
void directsinc3d(int ifl,int numlocs,int numeval,double *a1,double *a2,double *a3,double *klocs_d1,double *klocs_d2,double *klocs_d3, complex<double>* q,complex<double> *ans);
void directsincsq3d(int ifl, int numlocs,int numeval,double *a1,double *a2,double *a3,double *klocs_d1,double *klocs_d2,double *klocs_d3, complex<double>* q,complex<double> *ans);

#endif