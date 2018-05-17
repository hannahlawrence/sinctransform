#ifndef __sincdirect__hpp__included__
#define __sincdirect__hpp__included__

void directsinc1d(int ifl,int numlocs,double *klocs, complex<double>* q,complex<double> *ans, double pr);
void directsincsq1d(int ifl, int numlocs,double *klocs, complex<double>* q,complex<double> *ans, double pr);
void directsinc2d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2, complex<double>* q, complex<double> *ans, double pr);
void directsincsq2d(int ifl, int numlocs,double *klocs_d1,double *klocs_d2, complex<double>* q, complex<double> *ans, double pr);
void directsinc3d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3, complex<double>* q,complex<double> *ans, double pr);
void directsincsq3d(int ifl, int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3, complex<double>* q,complex<double> *ans, double pr);

#endif