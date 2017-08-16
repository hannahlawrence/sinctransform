#include <iostream>
#include "sincutil.hpp"
#include "sinctransform.hpp"

int main() // A simple example of usage for sinc3d and sincsq3d
{
	// Set precision
	double pr=1e-6;

	// Set sinc convention
	int ifl=1;

	// Create input arrays
	int numlocs=5;
	double klocs_d1[5]={0.1, 0.3, 0.5, 0.7, 0.9}; // X-coordinates of k locations
	double klocs_d2[5]={1, .3, 2, .4, 5}; // Y-coordinates of k locations
	double klocs_d3[5]={2, 3, 5, 7, 11}; // Z-coordinates of k locations
	double q[5]={0.2, 0.4, 0.6, 0.8, 1.0};

	// Test sinc1d and print the output
	double* myout_sinc3d=(double*)malloc(sizeof(double)*numlocs);
	int s_err=sinc3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,pr,myout_sinc3d); 
	cout<<"---Sinc3d Output---\n";
	printarr_double(myout_sinc3d,numlocs);

	// Test sincsq1d and print the output
	double* myout_sincsq3d=(double*)malloc(sizeof(double)*numlocs);
	s_err=sincsq3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,pr,myout_sincsq3d);
	cout<<"---Sincsq3d Output---\n";
	printarr_double(myout_sincsq3d,numlocs);

}