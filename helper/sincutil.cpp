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
	cout.precision(17);
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

double geterr(double* x, double* y, int n) // Calculate ||x-y||_2 / ||y||_2 over first n elements
{
	double err=0;
	double temp;
	for(int a=0;a<n;a++)
	{
		temp=x[a]-y[a];
		err=err+(temp*temp);
	}
	err=sqrt(err);
	double ynorm=0;
	for(int a=0;a<n;a++)
	{
		ynorm=ynorm+(y[a]*y[a]);
	}
	ynorm=sqrt(ynorm);
	return err/ynorm;
}

double getcerr(complex<double>* x, complex<double>* y, int n) // Calculate ||x-y||_2 / ||y||_2 over first n elements
{
	double err=0;
	double temp1,temp2;
	complex<double> idiff;
	for(int a=0;a<n;a++)
	{
		idiff=x[a]-y[a];
		temp1=real(idiff);
		temp2=imag(idiff);
		err=err+(temp1*temp1)+(temp2*temp2);
	}
	err=sqrt(err);
	double ynorm=0;
	for(int a=0;a<n;a++)
	{
		ynorm=ynorm+(real(y[a])*real(y[a]))+(imag(y[a])*imag(y[a]));
	}
	ynorm=sqrt(ynorm);
	return err/ynorm;
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