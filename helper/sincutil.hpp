#ifndef __sincutil__hpp__included__
#define __sincutil__hpp__included__

#include <iostream>
#include <math.h>
#include <complex>

using namespace std;

int equals(double x,double y,double eps);
void printarr_double(double* arr,int len);
void printarr_cdouble(complex<double>* arr,int len);
double geterr(double* a, double* b, int n);
void randarr(double lb,double ub,int n,double* arr);
void randcarr(double lb,double ub,int n,complex<double>* arr);
double getcerr(complex<double>* x, complex<double>* y, int n);

#endif