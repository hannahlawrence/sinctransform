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

// Contains sinc2d and sincsq2d

int sinc2d(int ifl,int numlocs,double *klocs_d1_,double *klocs_d2_,double *q,double tol,double *res)
{
	/*  
	Computes res[j] = sum sinc(klocs_d1_[k]-klocs_d1_[j]) * sinc(klocs_d2_[k]-klocs_d2_[j]) * q[j]
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
	double lowesttol=1e-16; // best precision that finufft can take
	double newtol=(tol/1000 > lowesttol) ? tol/1000 : lowesttol;
	float rkmaxx=0;
	float rkmaxy=0;
	double* klocs_d1=(double*)malloc(sizeof(double)*numlocs);
	double* klocs_d2=(double*)malloc(sizeof(double)*numlocs);
	complex<double>* qc=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	for(int a=0;a<numlocs;a++)
	{
		double a_d1=abs(klocs_d1_[a]); 
		double a_d2=abs(klocs_d2_[a]);
		if(a_d1>rkmaxx)
		{
			rkmaxx=a_d1;
		}
		if(a_d2>rkmaxy)
		{
			rkmaxy=a_d2;
		}
		qc[a]=complex<double> (q[a],0); //make q complex for finufft1d3
		klocs_d1[a]=klocs_d1_[a]; // to ensure that program does not alter input
		klocs_d2[a]=klocs_d2_[a];
	}
	if (ifl==1)
	{
		rkmaxx=pi*rkmaxx;
		rkmaxy=pi*rkmaxy;
		for(int a=0;a<numlocs;a++)
		{
			klocs_d1[a]=klocs_d1[a]*pi;
			klocs_d2[a]=klocs_d2[a]*pi;
		}
	}
	int rsamp=2;
	int nx=ceil(rsamp*round(rkmaxx+3));
	int ny=ceil(rsamp*round(rkmaxy+3));

	//create grid
	double* allxx=(double*)malloc(sizeof(double)*nx*ny);
	double* allyy=(double*)malloc(sizeof(double)*nx*ny);
	double* allww=(double*)malloc(sizeof(double)*nx*ny);
	int counter=0;
	double tempx,tempwx,tempy,tempwy;
	for(int a=0;a<ny;a++)
	{
		fastgl::QuadPair py = fastgl::GLPair((size_t)ny,(size_t)(a+1));
		tempy=py.x();
		tempwy=py.weight;
		for(int b=0;b<nx;b++)
		{
			fastgl::QuadPair px = fastgl::GLPair((size_t)nx,(size_t)(b+1));
			tempx=px.x();
			tempwx=px.weight;
			allxx[counter]=tempx;
			allyy[counter]=tempy;
			allww[counter]=tempwx*tempwy;
			counter++;
		}
	}

	nufft_opts opts;
	//h_at_xx will be complex
	complex<double>* h_at_xxyy=(complex<double>*)malloc(sizeof(complex<double>)*nx*ny);
	int ier1=finufft2d3(numlocs,klocs_d1,klocs_d2,qc,-1,newtol,nx*ny,allxx,allyy,h_at_xxyy,opts);
	if (ier1 != 0)
	{
		cout<<"Error: call 1 to finufft2d3 failed: "<<ier1<<"\n";
		return ier1;
	}
	complex<double>* weighted=(complex<double>*)malloc(sizeof(complex<double>)*nx*ny);
	for(int a=0;a<(nx*ny);a++)
	{
		weighted[a]=0.25*h_at_xxyy[a]*allww[a];
	}

	complex<double>* wtrans=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	int ier2=finufft2d3(nx*ny,allxx,allyy,weighted,1,newtol,numlocs,klocs_d1,klocs_d2,wtrans,opts);
	if (ier2 != 0)
	{
		cout<<"Error: call 2 to finufft2d3 failed: "<<ier2<<"\n";
		return ier2;
	}

	//make wtrans real (double array)
	for(int a=0;a<numlocs;a++)
	{
		res[a]=real(wtrans[a]);
	}

	free(wtrans);
	free(weighted);
	free(h_at_xxyy);
	free(qc);
	free(allxx);
	free(allyy);
	free(allww);
	return 0;
}

int sincsq2d(int ifl,int numlocs,double *klocs_d1_,double *klocs_d2_,double *q,double tol, double *res)
{
	/*  
	Computes res[j] = sum sinc^2(klocs_d1_[k]-klocs_d1_[j]) * sinc^2(klocs_d2_[k]-klocs_d2_[j]) * q[j]
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
	double lowesttol=1e-16; // best precision that finufft can take
	double newtol=(tol/1000 > lowesttol) ? tol/1000 : lowesttol;
	float rkmaxx=0;
	float rkmaxy=0;
	double* klocs_d1=(double*)malloc(sizeof(double)*numlocs);
	double* klocs_d2=(double*)malloc(sizeof(double)*numlocs);
	complex<double>* qc=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	for(int a=0;a<numlocs;a++)
	{
		double a_d1=(klocs_d1_[a] > 0) ? klocs_d1_[a] : -1*klocs_d1_[a];
		double a_d2=(klocs_d2_[a] > 0) ? klocs_d2_[a] : -1*klocs_d2_[a];
		if(a_d1>rkmaxx)
		{
			rkmaxx=a_d1;
		}
		if(a_d2>rkmaxy)
		{
			rkmaxy=a_d2;
		}
		qc[a]=complex<double> (q[a],0); //make q complex for finufft1d3
		klocs_d1[a]=klocs_d1_[a]; // to ensure that program does not alter input
		klocs_d2[a]=klocs_d2_[a];
	}
	if (ifl==1)
	{
		rkmaxx=pi*rkmaxx;
		rkmaxy=pi*rkmaxy;
		for(int a=0;a<numlocs;a++)
		{
			klocs_d1[a]=klocs_d1[a]*pi;
			klocs_d2[a]=klocs_d2[a]*pi;
		}
	}
	int rsamp=2;
	int nx=ceil(rsamp*round(rkmaxx+3));
	int ny=ceil(rsamp*round(rkmaxy+3));

	double* xx=(double*)malloc(sizeof(double)*nx*2);
	double* wwx=(double*)malloc(sizeof(double)*nx*2);
	double tempx;
	double tempw;
	for(int a=1;a<=nx;a++)
	{
		fastgl::QuadPair p = fastgl::GLPair((size_t)nx,(size_t)a);
		tempx=p.x();
		tempw=p.weight;
		xx[a-1]=tempx-1;
		wwx[a-1]=tempw*(2-abs(xx[a-1]));
		xx[a+nx-1]=tempx+1;
		wwx[a+nx-1]=tempw*abs(2-abs(xx[a+nx-1]));

	}

	double* yy=(double*)malloc(sizeof(double)*ny*2);
	double* wwy=(double*)malloc(sizeof(double)*ny*2);
	for(int a=1;a<=ny;a++)
	{
		fastgl::QuadPair p = fastgl::GLPair((size_t)ny,(size_t)a);
		tempx=p.x();
		tempw=p.weight;
		yy[a-1]=tempx-1;
		wwy[a-1]=tempw*(2-abs(yy[a-1]));
		yy[a+ny-1]=tempx+1;
		wwy[a+ny-1]=tempw*(2-abs(yy[a+ny-1]));
	}

	//create grid
	double* allxx=(double*)malloc(sizeof(double)*4*nx*ny);
	double* allyy=(double*)malloc(sizeof(double)*4*nx*ny);
	double* allww=(double*)malloc(sizeof(double)*4*nx*ny);
	int counter=0;

	for(int a=0;a<(2*ny);a++)
	{
		for(int b=0;b<(2*nx);b++)
		{
			allxx[counter]=xx[b];
			allyy[counter]=yy[a];
			allww[counter]=wwx[b]*wwy[a];
			counter++;
		}
	}
	
	nufft_opts opts;
	complex<double>* h_at_xxyy=(complex<double>*)malloc(sizeof(complex<double>)*4*nx*ny);
	int ier1=finufft2d3(numlocs,klocs_d1,klocs_d2,qc,-1,newtol,4*nx*ny,allxx,allyy,h_at_xxyy,opts);
	if (ier1 != 0)
	{
		cout<<"Error: call 1 to finufft2d3 failed: "<<ier1<<"\n";
		return ier1;
	}
	complex<double>* weighted=(complex<double>*)malloc(sizeof(complex<double>)*4*nx*ny);
	for(int a=0;a<(4*nx*ny);a++)
	{
		weighted[a]=(1.0/16)*h_at_xxyy[a]*allww[a];
	}

	complex<double>* wtrans=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	int ier2=finufft2d3(4*nx*ny,allxx,allyy,weighted,1,newtol,numlocs,klocs_d1,klocs_d2,wtrans,opts);
	if (ier2 != 0)
	{
		cout<<"Error: call 2 to finufft2d3 failed: "<<ier2<<"\n";
		return ier2;
	}

	//make wtrans real (double array)
	for(int a=0;a<numlocs;a++)
	{
		res[a]=real(wtrans[a]);
	}
	free(wtrans);
	free(weighted);
	free(h_at_xxyy);
	free(qc);
	free(allxx);
	free(allyy);
	free(allww);
	free(xx);
	free(yy);
	free(wwx);
	free(wwy);
	return 0;
}

