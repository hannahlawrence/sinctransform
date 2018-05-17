#include <iostream>
#include <math.h>
#include <complex>
#include "sincutil.hpp"

int equals(double x,double y,double eps) // Test whether two numbers are less than eps apart (for testing equality of doubles)
{
	double diff=x>y ? x-y : y-x;
	if(diff<eps)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void printarr_double(double* arr,int len) // Print a double array of length len in one line
{
	for(int a=0;a<len;a++)
	{
		cout<<arr[a]<<" ";
	}
	cout<<"\n";
} 
void printarr_cdouble(complex<double>* arr,int len) // Print a complex double array of length len in one line
{
	for(int a=0;a<len;a++)
	{
		cout<<arr[a]<<" ";
	}
	cout<<"\n";
}

double geterr(double* x, double* y, int n) // Calculate sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2 + ... + (x[n]-y[n])^2)
{
	double err=0;
	for(int a=0;a<n;a++)
	{
		err=err+pow(x[a]-y[a],2);
	}
	err=sqrt(err);
	return err;
}

double getcerr(complex<double>* x, complex<double>* y, int n) // Calculate sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2 + ... + (x[n]-y[n])^2)
{
	double err=0;
	complex<double> idiff;
	for(int a=0;a<n;a++)
	{
		idiff=x[a]-y[a];

		err=err+pow(real(idiff),2)+pow(imag(idiff),2);
	}
	err=sqrt(err);
	return err;
}

void randarr(double lb,double ub,int n,double* arr) // Populate arr, of length n, with random values between lb and ub
{
	for(int a=0;a<n;a++)
	{
		arr[a]=rand()/(double)RAND_MAX; // double between 0 and 1
		arr[a]=(ub-lb)*arr[a] + lb;
	}
}

void randcarr(double lb,double ub,int n,complex<double>* arr) // Populate arr, of length n, with random values between lb and ub
{
	double tempa;
	double tempb;
	for(int a=0;a<n;a++)
	{
		tempa=rand()/(double)RAND_MAX;
		tempa=(ub-lb)*tempa + lb;

		tempb=rand()/(double)RAND_MAX;
		tempb=(ub-lb)*tempb + lb;
		arr[a]=complex<double>(tempa,tempb);
	}
}