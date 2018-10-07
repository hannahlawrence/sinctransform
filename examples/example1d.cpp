#include <iostream>
#include "sincutil.hpp"
#include "sinctransform.hpp"
#include "directsinc.hpp"
#include <vector>
//#include <complex>

int main() // A simple example of usage for sinc1d and sincsq1d
{
	// Set precision
	double pr=1e-12;

	// Set sinc convention
	int ifl=1;

	// Create input arrays
	int numlocs=5;
	double klocs[5]={100, 0.3, 15, 0.7, 0.9};
	double a1[5]={0.2,0.4,6,0.8,1.0};
	complex<double> q[5]={{2,1}, {0.4,0}, {0.6,-0.4}, {0.8,0.3}, {1.0,-0.7}};

	// Test sinc1d and print the output
	std::vector<complex<double>> myout_sinc1d(numlocs);
	int s_err=sinc1d(ifl,numlocs,a1,klocs,q,pr,myout_sinc1d.data(),1);
	cout<<"---Sinc1d Output---\n";
	printarr_cdouble(myout_sinc1d.data(),numlocs);

	std::vector<complex <double>> corr(numlocs);
	directsinc1d(ifl,numlocs,numlocs,a1,klocs,q,corr.data()); 
	cout<<"Correct: \n";
	printarr_cdouble(corr.data(),numlocs);
	double err=getcerr(myout_sinc1d.data(),corr.data(),numlocs);
	cout<<"Error: "<<err<<"\n";

	// Test sincsq1d and print the output
	std::vector<complex<double>> myout_sincsq1d(numlocs);
	s_err=sincsq1d(ifl,numlocs,a1,klocs,q,pr,myout_sincsq1d.data(),1); 
	cout<<"---Sincsq1d Output---\n";
	printarr_cdouble(myout_sincsq1d.data(),numlocs);
	directsincsq1d(ifl,numlocs,numlocs,a1,klocs,q,corr.data()); 
	cout<<"Correct: \n";
	printarr_cdouble(corr.data(),numlocs);
	err=getcerr(myout_sincsq1d.data(),corr.data(),numlocs);
	cout<<"Error: "<<err<<"\n";
	return s_err;
}
