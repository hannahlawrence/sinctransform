#include <iostream>
#include "sincutil.hpp"
#include "sinctransform.hpp"

int main() // A simple example of usage for sinc2d and sincsq2d
{
	// Set precision
	double pr=1e-14;

	// Set sinc convention
	int ifl=0;

	// Create input arrays
	int numlocs=5;
	double klocs_d1[5]={0.1, 0.3, 0.5, 0.7, 0.9}; // X-coordinates of k locations
	double klocs_d2[5]={0.2, 3, .1, 2, 10}; // Y-coordinates of k locations
	complex<double> q[5]={0.2+1i, 0.4, 0.6-0.4i, 0.8+0.3i, 1.0-0.7i};

	// Test sinc1d and print the output
	complex<double>* myout_sinc2d=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	int s_err=sinc2d(ifl,numlocs,klocs_d1,klocs_d2,q,pr,myout_sinc2d); 
	cout<<"---Sinc2d Output---\n";
	printarr_cdouble(myout_sinc2d,numlocs);

	// Test sincsq1d and print the output
	complex<double>* myout_sincsq2d=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	s_err=sincsq2d(ifl,numlocs,klocs_d1,klocs_d2,q,pr,myout_sincsq2d);
	cout<<"---Sincsq2d Output---\n";
	printarr_cdouble(myout_sincsq2d,numlocs);
	
	return s_err;
}
