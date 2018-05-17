#include <iostream>
#include "sincutil.hpp"
#include "sinctransform.hpp"
#include "directsinc.hpp"
#include <complex>

int main() // A simple example of usage for sinc1d and sincsq1d
{
	// Set precision
	double pr=1e-6;

	// Set sinc convention
	int ifl=1;

	// Create input arrays
	int numlocs=5;
	double klocs[5]={0.1, 0.3, 0.5, 0.7, 0.9};
	complex<double> q[5]={0.2+1i, 0.4, 0.6-0.4i, 0.8+0.3i, 1.0-0.7i};

	// Test sinc1d and print the output
	complex<double>* myout_sinc1d=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	int s_err=sinc1d(ifl,numlocs,klocs,q,pr,myout_sinc1d);
	cout<<"---Sinc1d Output---\n";
	printarr_cdouble(myout_sinc1d,numlocs);

	// Test sincsq1d and print the output
	complex<double>* myout_sincsq1d=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	s_err=sincsq1d(ifl,numlocs,klocs,q,pr,myout_sincsq1d); 
	cout<<"---Sincsq1d Output---\n";
	printarr_cdouble(myout_sincsq1d,numlocs);
}