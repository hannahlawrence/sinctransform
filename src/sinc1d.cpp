#include <iostream>
#include <iomanip>
#include <complex>
#include <math.h> 
#include <ctime> 
#include "finufft.h" // for non-uniform Fourier transform
#include "fastgl.hpp" //for Legendre-Gauss weights
#include "sincutil.hpp"
#include "sinctransform.hpp"

using namespace std;

// Contains sinc1d and sincsq1d

int sinc1d(int ifl,int numlocs,double *klocs_, complex<double> *q,double tol, complex<double> *res)
{
	/*  
	Computes res[j] = sum sinc(klocs_[k]-klocs_[j]) * q[j]
	             	   k

	Inputs:
		ifl = sinc convention
			0: sinc(x) = sin(x)/x
			1: sinc(x)=sin(pi*x)/(pi*x)
		klocs_ = (real) sample locations
		q = sample strengths
		tol = requested precision

	Returns:
		0: if success
		Other: error code, as returned by finufft (see finufft documentation)
	*/
	double pi=4*atan(1); 
	double lowesttol=1e-16; //best precision that finufft can take
	double newtol=(tol/1000 > lowesttol) ? tol/1000 : lowesttol;

	double *klocs=(double*)malloc(sizeof(double)*numlocs);

	// find rkmax
	double rkmax=0;
	double temp;
	complex<double>* qc=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	for(int a=0;a<numlocs;a++)
	{
		temp=abs(klocs_[a]);
		if (temp>rkmax)
		{
			rkmax=temp;
		}
		//qc[a]=complex<double> (q[a],0); //make q complex for finufft1d3
		qc[a]=q[a];
		klocs[a]=klocs_[a]; // to ensure that program does not alter input
	}

	// adjust based on ifl convention
	if (ifl==1)
	{
		rkmax=pi*rkmax;
		for(int a=0;a<numlocs;a++)
		{
			klocs[a]=klocs[a]*pi;
		}
	}

	// set nx
	int rsamp=2;
	int nx=ceil(rsamp*round(rkmax+3));

	// get lgwt weights
	double* xx=(double*)malloc(sizeof(double)*nx);
	double* ww=(double*)malloc(sizeof(double)*nx);
	for(int a=1;a<=nx;a++)
	{
		fastgl::QuadPair p = fastgl::GLPair((size_t)nx,(size_t)a);
		xx[a-1]=p.x();
		ww[a-1]=p.weight;
	}

	// call 1 to finufft; h_at_xx will be complex
	nufft_opts opts; finufft_default_opts(opts);	
	complex<double>* h_at_xx=(complex<double>*)malloc(sizeof(complex<double>)*nx);

	int ier1=finufft1d3(numlocs,klocs,qc,-1,newtol,nx,xx,h_at_xx,opts);
	if (ier1 != 0)
	{
		cout<<"Error: call 1 to finufft1d3 failed: "<<ier1<<"\n";
		return ier1;
	}

	// call 2 to finufft
	// source strengths will need to be complex!
	complex<double>* weighted=(complex<double>*)malloc(sizeof(complex<double>)*nx);
	for(int a=0;a<nx;a++)
	{
		weighted[a]=0.5*h_at_xx[a]*ww[a];
	}
	complex<double>* wtrans=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	int ier2=finufft1d3(nx,xx,weighted,1,newtol,numlocs,klocs,wtrans,opts);
	if (ier2 != 0)
	{
		cout<<"Error: call 2 to finufft1d3 failed: "<<ier2<<"\n";
		return ier2;
	}

	//make wtrans real (double array)

	for(int a=0;a<numlocs;a++)
	{
		res[a]=wtrans[a];
	}

	free(wtrans);
	free(weighted);
	free(h_at_xx);
	free(qc);
	free(xx);
	free(ww);
	return 0;
}

int sincsq1d(int ifl,int numlocs,double *klocs_, complex<double> *q,double tol, complex<double> *res)
{
	/*  
	Computes res[j] = sum sinc^2(klocs_[k]-klocs_[j]) * q[j]
	             	   k

	Inputs:
		ifl = sinc convention
			0: sinc(x) = sin(x)/x
			1: sinc(x)=sin(pi*x)/(pi*x)
		klocs_ = (real) sample locations
		q = sample strengths
		tol = requested precision

	Returns:
		0: if success
		Other: error code, as returned by finufft (see finufft documentation)
	*/
	double pi=4*atan(1);
	double lowesttol=1e-16; //hopefully 1e-16, which is best precision that finufft can take?
	double newtol=(tol/1000 > lowesttol) ? tol/1000 : lowesttol;

	double *klocs=(double*)malloc(sizeof(double)*numlocs);
	
	// find rkmax
	double rkmax=0;
	double temp;
	complex<double>* qc=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	for(int a=0;a<numlocs;a++)
	{
		temp=abs(klocs_[a]);
		if (temp>rkmax)
		{
			rkmax=temp;
		}
		//qc[a]=complex<double> (q[a],0); //make q complex for finufft1d3
		qc[a]=q[a];
		klocs[a]=klocs_[a]; // to ensure that program does not alter input
	}

	// adjust based on ifl convention
	if (ifl==1)
	{
		rkmax=pi*rkmax;
		for(int a=0;a<numlocs;a++)
		{
			klocs[a]=klocs[a]*pi;
		}
	}

	// set nx
	int rsamp=2;
	int nx=ceil(rsamp*round(rkmax+3));
	
	// get lgwt weights
	double* xx=(double*)malloc(sizeof(double)*2*nx);
	double* ww=(double*)malloc(sizeof(double)*2*nx);
	double tempx,tempw;
	for(int a=1;a<=nx;a++)
	{
		fastgl::QuadPair p = fastgl::GLPair((size_t)nx,(size_t)a);
		tempx=p.x();
		tempw=p.weight;
		xx[a-1]=tempx-1;
		xx[a-1+nx]=tempx+1;
		ww[a-1]=tempw;
		ww[a+nx-1]=tempw;
	}

	// call 1 to finufft; h_at_xx will be complex
	nufft_opts opts; finufft_default_opts(opts);	
	complex<double>* h_at_xx=(complex<double>*)malloc(sizeof(complex<double>)*2*nx);
	int ier1=finufft1d3(numlocs,klocs,qc,-1,newtol,2*nx,xx,h_at_xx,opts);
	if (ier1 != 0)
	{
		cout<<"Issue: call 1 to finufft1d3 failed\n";
		return ier1;
	}

	// call 2 to finufft
	// source strengths will need to be complex!
	complex<double>* weighted=(complex<double>*)malloc(sizeof(complex<double>)*2*nx);
	for(int a=0;a<2*nx;a++)
	{
		weighted[a]=0.25*h_at_xx[a]*ww[a]*(2-abs(xx[a]));
	}
	complex<double>* wtrans=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	int ier2=finufft1d3(2*nx,xx,weighted,1,newtol,numlocs,klocs,wtrans,opts);
	if (ier2 != 0)
	{
		cout<<"Issue: call 2 to finufft1d3 failed\n";
		return ier2;
	}

	//make wtrans real (double array)
	for(int a=0;a<numlocs;a++)
	{
		res[a]=wtrans[a];
	}

	free(weighted);
	free(wtrans);
	free(qc);
	free(h_at_xx);
	free(xx);
	free(ww);
	return 0;
}

