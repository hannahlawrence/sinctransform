#include <iostream>
#include "sincutil.hpp"
#include "sinctransform.hpp"
#include "directsinc.hpp"
#include <vector>

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
	double a1[5]={0.2,0.4,0.6,0.8,1.0};
	double a2[5]={3,1,2,0,1};
	complex<double> q[5]={{0.2,1}, {0.4,0}, {0.6,-0.4}, {0.8,0.3}, {1.0,-0.7}};
	//complex<double> q[5]={0.2+1i, 0.4, 0.6-0.4i, 0.8+0.3i, 1.0-0.7i};

	// Test sinc2d and print the output
	std::vector<complex <double>> myout_sinc2d(numlocs);
	int s_err=sinc2d(ifl,numlocs,a1,a2,klocs_d1,klocs_d2,q,pr,myout_sinc2d.data(),1); 
	cout<<"---Sinc2d Output---\n";
	printarr_cdouble(myout_sinc2d.data(),numlocs);

	std::vector<complex <double>> corr(numlocs);
	directsinc2d(ifl,numlocs,numlocs,a1,a2,klocs_d1,klocs_d2,q,corr.data()); 
	cout<<"Correct: \n";
	printarr_cdouble(corr.data(),numlocs);
	double err=getcerr(myout_sinc2d.data(),corr.data(),numlocs);
	cout<<"Error: "<<err<<"\n";

	// Test sincsq2d and print the output
	std::vector<complex <double>> myout_sincsq2d(numlocs);
	s_err=sincsq2d(ifl,numlocs,a1,a2,klocs_d1,klocs_d2,q,pr,myout_sincsq2d.data(),1);
	cout<<"---Sincsq2d Output---\n";
	printarr_cdouble(myout_sincsq2d.data(),numlocs);
	
	std::vector<complex <double>> corr2(numlocs);
	directsincsq2d(ifl,numlocs,numlocs,a1,a2,klocs_d1,klocs_d2,q,corr2.data()); 
	cout<<"Correct: \n";
	printarr_cdouble(corr2.data(),numlocs);
	err=getcerr(myout_sincsq2d.data(),corr2.data(),numlocs);
	cout<<"Error: "<<err<<"\n";

	return s_err;
}
