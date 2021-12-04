#include <iostream>
#include "sincutil.hpp"
#include "sinctransform.hpp"
#include "directsinc.hpp"
#include <vector>

int main() // A simple example of usage for sinc3d and sincsq3d
{
	// Set precision
	double pr=1e-8;

	// Set sinc convention
	int ifl=1;

	// Create input arrays
	int numlocs=5;
	double klocs_d1[5]={0.1, 0.3, 0.5, 0.7, 0.9}; // X-coordinates of k locations
	double klocs_d2[5]={1, .3, 2, .4, 5}; // Y-coordinates of k locations
	double klocs_d3[5]={2, 3, 5, 7, 11}; // Z-coordinates of k locations
	double a1[5]={0.2,0.4,0.6,0.8,1.0};
	double a2[5]={3,1,2,0,1};
	double a3[5]={1,2,3,4,5};
	complex<double> q[5]={{0.2,1}, {0.4,0}, {0.6,-0.4}, {0.8,0.3}, {1.0,-0.7}};

	// Test sinc3d and print the output
	std::vector<complex <double>> myout_sinc3d(numlocs);
	int s_err=sinc3d(ifl,numlocs,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q,pr,myout_sinc3d.data(),1); 
	cout<<"---Sinc3d Output---\n";
	printarr_cdouble(myout_sinc3d.data(),numlocs);

	std::vector<complex <double>> corr(numlocs);
	directsinc3d(ifl,numlocs,numlocs,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q,corr.data()); 
	cout<<"Correct: \n";
	printarr_cdouble(corr.data(),numlocs);
	double err=getcerr(myout_sinc3d.data(),corr.data(),numlocs);
	cout<<"Error: "<<err<<"\n";

	// Test sincsq3d and print the output
	std::vector<complex <double>> myout_sincsq3d(numlocs);
	s_err=sincsq3d(ifl,numlocs,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q,pr,myout_sincsq3d.data(),1);
	cout<<"---Sincsq3d Output---\n";
	printarr_cdouble(myout_sincsq3d.data(),numlocs);

	directsincsq3d(ifl,numlocs,numlocs,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q,corr.data()); 
	cout<<"Correct: \n";
	printarr_cdouble(corr.data(),numlocs);
	err=getcerr(myout_sincsq3d.data(),corr.data(),numlocs);
	cout<<"Error: "<<err<<"\n";

	return s_err;
}
