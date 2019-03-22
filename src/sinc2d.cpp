#include <iostream>
#include <iomanip>
#include <complex>
#include <math.h> 
#include <ctime> 
#include "finufft.h" // for non-uniform Fourier transform
#include "fastgl.hpp" //for Legendre-Gauss weights
#include "sincutil.hpp"
#include "sinctransform.hpp"
#include "constants.hpp"
#include <vector>

using namespace std;

// Contains sinc2d and sincsq2d

int sinc2d(int ifl,int numlocs,double *a1_,double *a2_,double *klocs_d1_,double *klocs_d2_,complex<double> *q,double tol,complex<double> *res,int quad)
{
	/*  
	Computes res[j] = sum sinc(a1_[k]-klocs_d1_[j]) * sinc(a2_[k]-klocs_d2_[j]) * q[j]
	             	   k

	Inputs:
		ifl = sinc convention
			0: sinc(x) = sin(x)/x
			1: sinc(x)=sin(pi*x)/(pi*x)
		a1_ = (real) evaluation locations
		a2_ = (real) evaluation locations
		klocs_ = (real) sample locations
		q = sample strengths
		tol = requested precision
		quad = quadrature mode; 0 for legendre, 1 for corrected trapezoidal

	Returns:
		0: if success
		Other: error code, as returned by finufft (see finufft documentation)
	*/

	double pi=4*atan(1);
	double lowesttol=1e-15; // best precision that finufft can take
	double newtol=(tol > lowesttol) ? tol : lowesttol;
	float rkmaxx=0;
	float rkmaxy=0;
	std::vector<double> klocs_d1(numlocs);
	std::vector<double> klocs_d2(numlocs);
	std::vector<double> a1(numlocs);
	std::vector<double> a2(numlocs);
	std::vector<complex <double>> qc(numlocs);
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
		//qc[a]=complex<double> (q[a],0); //make q complex for finufft1d3
		qc[a]=q[a];
		klocs_d1[a]=klocs_d1_[a]; // to ensure that program does not alter input
		klocs_d2[a]=klocs_d2_[a];
		a1[a]=a1_[a];
		a2[a]=a2_[a];
	}
	if (ifl==1)
	{
		rkmaxx=pi*rkmaxx;
		rkmaxy=pi*rkmaxy;
		for(int a=0;a<numlocs;a++)
		{
			klocs_d1[a]=klocs_d1[a]*pi;
			klocs_d2[a]=klocs_d2[a]*pi;
			a1[a]=a1[a]*pi;
			a2[a]=a2[a]*pi;
		}
	}
	double rsamp;
	int nx,ny;

	if (quad==0)
	{
		rsamp=2;
		nx=ceil(rsamp*round(rkmaxx+3));
		ny=ceil(rsamp*round(rkmaxy+3));
		//create grid
		std::vector<double> allxx(nx*ny);
		std::vector<double> allyy(nx*ny);
		std::vector<double> allww(nx*ny);
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
	    nufft_opts opts; finufft_default_opts(&opts);
		if (newtol<1e-9)
		{
			opts.upsampfac=2;
		}
		//h_at_xx will be complex
		std::vector<complex <double>> h_at_xxyy(nx*ny);
		int ier1=finufft2d3(numlocs,klocs_d1.data(),klocs_d2.data(),qc.data(),-1,newtol,nx*ny,allxx.data(),allyy.data(),h_at_xxyy.data(),opts);
		if (ier1 != 0)
		{
			cout<<"Error: call 1 to finufft2d3 failed: "<<ier1<<"\n";
			return ier1;
		}
		std::vector<complex <double>> weighted(nx*ny);
		for(int a=0;a<(nx*ny);a++)
		{
			weighted[a]=0.25*h_at_xxyy[a]*allww[a];
		}

		std::vector<complex <double>> wtrans(numlocs);
		int ier2=finufft2d3(nx*ny,allxx.data(),allyy.data(),weighted.data(),1,newtol,numlocs,a1.data(),a2.data(),wtrans.data(),opts);
		if (ier2 != 0)
		{
			cout<<"Error: call 2 to finufft2d3 failed: "<<ier2<<"\n";
			return ier2;
		}

		//make wtrans real (double array)
		for(int a=0;a<numlocs;a++)
		{
			res[a]=wtrans[a];
		}
	}
	else // trapezoidal quadrature
	{
		rsamp=3;
		nx=ceil(rsamp*round(rkmaxx+3));
		ny=ceil(rsamp*round(rkmaxy+3));
		int a=-1;
		int b=1;
		int e=25; // increase (up to 60) to impose higher accuracy; will increase runtime
		std::vector<double> constants(e);
		for (int i=0;i<e;i++)
			constants[i]=allconstants[e-1][i];
		if ((nx % 2) != 0)
			nx=nx+1; // ensure even so that 0 is a quadrature point
		int n=nx;
		double h=(double (b-a))/n;
		int aind=e;
		int bind=aind+n;
		int numunif=n+(2*e)+1;
		// Get constants, double vector of length e 
		std::vector<double> xx(numunif);
		std::vector<double> ww_x(numunif);
		int msx=numunif;
		for (int i=0;i<numunif;i++) // Initialize all to 0
			ww_x[i]=0;
		double curr=a-(e*h);
		for (int i=0;i<numunif;i++)
		{
			xx[i]=curr;
			curr=curr+h;
		}
		ww_x[aind]=0.5;
		ww_x[bind]=0.5;
		for (int i=(aind+1);i<bind;i++)
		{
			ww_x[i]=1;
		}
		for (int k=1;k<=e;k++) 
		{
			ww_x[aind-k] = ww_x[aind-k] - constants[k-1];
		    ww_x[aind+k] = ww_x[aind+k] + constants[k-1];
		    ww_x[bind-k] = ww_x[bind-k] + constants[k-1];
		    ww_x[bind+k] = ww_x[bind+k] - constants[k-1];
		}
		for (int i=0;i<numunif;i++)
			ww_x[i]=ww_x[i]*h;

		if ((ny % 2) != 0)
			ny=ny+1; // ensure even so that 0 is a quadrature point
		n=ny;
		h=(double (b-a))/n;
		aind=e; 
		bind=aind+n;
		numunif=n+(2*e)+1;
		int msy=numunif;
		std::vector<double> yy(numunif);
		std::vector<double> ww_y(numunif);
		for (int i=0;i<numunif;i++) // Initialize all to 0
			ww_y[i]=0;
		curr=a-(e*h);
		for (int i=0;i<numunif;i++)
		{
			yy[i]=curr;
			curr=curr+h;
		}
		ww_y[aind]=0.5;
		ww_y[bind]=0.5;
		for (int i=(aind+1);i<bind;i++)
		{
			ww_y[i]=1;
		}
		for (int k=1;k<=e;k++) 
		{
			ww_y[aind-k] = ww_y[aind-k] - constants[k-1];
		    ww_y[aind+k] = ww_y[aind+k] + constants[k-1];
		    ww_y[bind-k] = ww_y[bind-k] + constants[k-1];
		    ww_y[bind+k] = ww_y[bind+k] - constants[k-1];
		}
		for (int i=0;i<numunif;i++)
			ww_y[i]=ww_y[i]*h;


		double actual_unif_spacex=xx[1]-xx[0];
		double actual_unif_spacey=yy[1]-yy[0];
		double Lx=xx[0];
		double Ly=yy[0];

		double DU1x,DU1y;
		if ((msx%2)==0)
			DU1x=-msx/2;
		else
			DU1x=(-msx+1)/2;
		double translationx=(DU1x-(Lx/actual_unif_spacex));
		if ((msy%2)==0)
			DU1y=-msy/2;
		else
			DU1y=(-msy+1)/2;
		double translationy=(DU1y-(Ly/actual_unif_spacey));

		// Make call to finufft1d1
		nufft_opts opts; finufft_default_opts(&opts);

		std::vector<complex<double>> h_at_xxyy(msx*msy); // length?
		
		complex<double> icomp=-1;
		icomp=sqrt(icomp);
		std::vector<double> newklocsx(numlocs);
		std::vector<double> newklocsy(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newklocsx[i]=klocs_d1[i]*actual_unif_spacex;
		}
		for (int i=0;i<numlocs;i++)
		{
			newklocsy[i]=klocs_d2[i]*actual_unif_spacey;
		}
		std::vector<complex<double>> strengths(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			strengths[i]=q[i]*exp(icomp*((actual_unif_spacex*klocs_d1[i]*translationx)+(actual_unif_spacey*klocs_d2[i]*translationy))); 
		}
		int ier=finufft2d1(numlocs,newklocsx.data(),newklocsy.data(),strengths.data(),-1,newtol,msx,msy,h_at_xxyy.data(),opts);

		translationx=(-1*DU1x)*actual_unif_spacex+Lx;
		translationy=(-1*DU1y)*actual_unif_spacey+Ly;

		// need to make allww!
		int counter=0;
		std::vector<complex <double>> newh_at_xxyy(msx*msy);
		for (int i=0;i<msy;i++)
		{
			for (int j=0;j<msx;j++)
			{
				newh_at_xxyy[counter]=h_at_xxyy[counter]*ww_x[j]*ww_y[i];
				counter=counter+1;
			}
		}
		std::vector<double> newa1(numlocs);
		std::vector<double> newa2(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newa1[i]=a1[i]*actual_unif_spacex;
			newa2[i]=a2[i]*actual_unif_spacey;
		}
		std::vector<complex<double>> temp(numlocs);
		ier=finufft2d2(numlocs,newa1.data(),newa2.data(),temp.data(),1,newtol,msx,msy,newh_at_xxyy.data(),opts); 
		//do exp thing
		for (int i=0;i<numlocs;i++)
		{
			res[i]=0.25*exp(i*(translationx*a1[i]+translationy*a2[i]))*temp[i];
		}
	}

	return 0;
}


int sincsq2d(int ifl,int numlocs,double *a1_,double *a2_,double *klocs_d1_,double *klocs_d2_,complex<double> *q,double tol, complex<double> *res,int quad)
{
	/*
	Computes res[j] = sum sinc^2(a1[k]-klocs_d1_[j]) * sinc^2(a2[k]-klocs_d2_[j]) * q[j]
	             	   k

	Inputs:
		ifl = sinc convention
			0: sinc(x) = sin(x)/x
			1: sinc(x)=sin(pi*x)/(pi*x)
		a1_ = (real) evaluation locations
		a2_ = (real) evaluation locations
		klocs_ = (real) sample locations
		q = sample strengths
		tol = requested precision
		quad = quadrature mode; 0 for legendre, 1 for corrected trapezoidal

	Returns:
		0: if success
		Other: error code, as returned by finufft (see finufft documentation)
	*/

	double pi=4*atan(1);
	double lowesttol=1e-15; // best precision that finufft can take
	double newtol=(tol > lowesttol) ? tol : lowesttol;
	float rkmaxx=0;
	float rkmaxy=0;
	std::vector<double> klocs_d1(numlocs);
	std::vector<double> klocs_d2(numlocs);
	std::vector<double> a1(numlocs);
	std::vector<double> a2(numlocs);
	std::vector<complex <double>> qc(numlocs);
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
		qc[a]=q[a]; //make q complex for finufft1d3
		klocs_d1[a]=klocs_d1_[a]; // to ensure that program does not alter input
		klocs_d2[a]=klocs_d2_[a];
		a1[a]=a1_[a];
		a2[a]=a2_[a];
	}
	if (ifl==1)
	{
		rkmaxx=pi*rkmaxx;
		rkmaxy=pi*rkmaxy;
		for(int a=0;a<numlocs;a++)
		{
			klocs_d1[a]=klocs_d1[a]*pi;
			klocs_d2[a]=klocs_d2[a]*pi;
			a1[a]=a1[a]*pi;
			a2[a]=a2[a]*pi;
		}
	}
	double rsamp;
	int nx,ny;

	if (quad==0)
	{
		rsamp=2;
		nx=ceil(rsamp*round(rkmaxx+3));
		ny=ceil(rsamp*round(rkmaxy+3));
		std::vector<double> xx(nx*2);
		std::vector<double> wwx(nx*2);
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
		std::vector<double> yy(ny*2);
		std::vector<double> wwy(ny*2);
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
		std::vector<double> allxx(4*nx*ny);
		std::vector<double> allyy(4*nx*ny);
		std::vector<double> allww(4*nx*ny);
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
		nufft_opts opts; finufft_default_opts(&opts);
		if (newtol<1e-9)
		{
			opts.upsampfac=2;
		}
		tempx=3;
		std::vector<complex <double>> h_at_xxyy(4*nx*ny);
		for (int a=0;a<(4*nx*ny);a++)
			h_at_xxyy[a]=0;
		int ier1=finufft2d3(numlocs,klocs_d1.data(),klocs_d2.data(),qc.data(),-1,newtol,4*nx*ny,allxx.data(),allyy.data(),h_at_xxyy.data(),opts);
		if (ier1 != 0)
		{
			cout<<"Error: call 1 to finufft2d3 failed: "<<ier1<<"\n";
			return ier1;
		}
		std::vector<complex <double>> weighted(4*nx*ny);
		for(int a=0;a<(4*nx*ny);a++)
		{
			weighted[a]=(1.0/16)*h_at_xxyy[a]*allww[a];
		}
		std::vector<complex <double>> wtrans(numlocs);
		int ier2=finufft2d3(4*nx*ny,allxx.data(),allyy.data(),weighted.data(),1,newtol,numlocs,a1.data(),a2.data(),wtrans.data(),opts);
		if (ier2 != 0)
		{
			cout<<"Error: call 2 to finufft2d3 failed: "<<ier2<<"\n";
			return ier2;
		}

		//make wtrans real (double array)
		for(int a=0;a<numlocs;a++)
		{
			res[a]=wtrans[a];
		}
	}
	else
	{
		rsamp=3;
		nx=ceil(rsamp*round(rkmaxx+3));
		ny=ceil(rsamp*round(rkmaxy+3));
		int a=-2;
		int b=0;
		int e=21; // increase (up to 60) to impose higher accuracy; will increase runtime
		std::vector<double> constants(e);
		for (int i=0;i<e;i++)
			constants[i]=allconstants[e-1][i];
		if ((nx % 2) != 0)
			nx=nx+1; // ensure even so that 0 is a quadrature point
		int n=nx;
		double h=(double (b-a))/n;
		int aind=e;
		int zind=aind+n;
		int bind=zind+n;
		int numunif=(2*n)+(2*e)+1;
		// Get constants, double vector of length e 
		std::vector<double> xx(numunif);
		std::vector<double> ww_x(numunif);
		int msx=numunif;
		for (int i=0;i<numunif;i++) // Initialize all to 0
			ww_x[i]=0;
		double curr=a-(e*h);
		for (int i=0;i<numunif;i++)
		{
			xx[i]=curr;
			curr=curr+h;
		}
		std::vector<double> leftvec(numunif);
		std::vector<double> rightvec(numunif);
		std::vector<double> trianglevec(numunif);
		std::vector<double> ww_trap(numunif);
		std::vector<double> ww_left(numunif);
		std::vector<double> ww_right(numunif);
		std::vector<double> ww(numunif);
		for (int i=0;i<numunif;i++)
		{
			leftvec[i]=0;
			rightvec[i]=0;
			trianglevec[i]=0;
			ww_trap[i]=0;
			ww_left[i]=0;
			ww_right[i]=0;
		}
		for (int i=0;i<numunif;i++)
		{
			leftvec[i]=2+xx[i];
			rightvec[i]=2-xx[i];
			trianglevec[i]=2-abs(xx[i]);
		}
		ww_trap[aind]=0.5;
		ww_trap[bind]=0.5;
		for (int i=aind+1;i<bind;i++)
			ww_trap[i]=1;
		for (int i=aind;i<=bind;i++)
			ww_trap[i]=ww_trap[i]*trianglevec[i]; // incorporate real values
		for (int k=1;k<=e;k++)
		{
			ww_left[aind-k] = ww_left[aind-k] - leftvec[aind-k]*constants[k-1];
		    ww_left[aind+k] = ww_left[aind+k] + leftvec[aind+k]*constants[k-1];
		    ww_left[zind-k] = ww_left[zind-k] + leftvec[zind-k]*constants[k-1];
		    ww_left[zind+k] = ww_left[zind+k] - leftvec[zind+k]*constants[k-1];
		    ww_right[zind-k] =  ww_right[zind-k]- rightvec[zind-k]*constants[k-1];
		    ww_right[zind+k] =  ww_right[zind+k]+ rightvec[zind+k]*constants[k-1];
		    ww_right[bind-k] =  ww_right[bind-k]+ rightvec[bind-k]*constants[k-1];
		    ww_right[bind+k] =  ww_right[bind+k]- rightvec[bind+k]*constants[k-1];
		}
		for (int i=0;i<numunif;i++)
			ww_x[i]=h*(ww_trap[i]+ww_left[i]+ww_right[i]);

		if ((ny % 2) != 0)
			ny=ny+1; // ensure even so that 0 is a quadrature point
		n=ny;
		h=(double (b-a))/n;
		aind=e;
		zind=aind+n;
		bind=zind+n;
		numunif=(2*n)+(2*e)+1;
		// Get constants, double vector of length e 
		std::vector<double> yy(numunif);
		std::vector<double> ww_y(numunif);
		int msy=numunif;
		for (int i=0;i<numunif;i++) // Initialize all to 0
			ww_y[i]=0;
		curr=a-(e*h);
		for (int i=0;i<numunif;i++)
		{
			yy[i]=curr;
			curr=curr+h;
		}

		leftvec.resize(numunif);
		rightvec.resize(numunif);
		trianglevec.resize(numunif);
		ww_trap.resize(numunif);
		ww_left.resize(numunif);
		ww_right.resize(numunif);
		ww.resize(numunif);
		for (int i=0;i<numunif;i++)
		{
			leftvec[i]=0;
			rightvec[i]=0;
			trianglevec[i]=0;
			ww_trap[i]=0;
			ww_left[i]=0;
			ww_right[i]=0;
		}
		for (int i=0;i<numunif;i++)
		{
			leftvec[i]=2+yy[i];
			rightvec[i]=2-yy[i];
			trianglevec[i]=2-abs(yy[i]);
		}
		ww_trap[aind]=0.5;
		ww_trap[bind]=0.5;
		for (int i=aind+1;i<bind;i++)
			ww_trap[i]=1;
		for (int i=aind;i<=bind;i++)
			ww_trap[i]=ww_trap[i]*trianglevec[i]; // incorporate real values
		for (int k=1;k<=e;k++)
		{
			ww_left[aind-k] = ww_left[aind-k] - leftvec[aind-k]*constants[k-1];
		    ww_left[aind+k] = ww_left[aind+k] + leftvec[aind+k]*constants[k-1];
		    ww_left[zind-k] = ww_left[zind-k] + leftvec[zind-k]*constants[k-1];
		    ww_left[zind+k] = ww_left[zind+k] - leftvec[zind+k]*constants[k-1];
		    ww_right[zind-k] =  ww_right[zind-k]- rightvec[zind-k]*constants[k-1];
		    ww_right[zind+k] =  ww_right[zind+k]+ rightvec[zind+k]*constants[k-1];
		    ww_right[bind-k] =  ww_right[bind-k]+ rightvec[bind-k]*constants[k-1];
		    ww_right[bind+k] =  ww_right[bind+k]- rightvec[bind+k]*constants[k-1];
		}
		for (int i=0;i<numunif;i++)
			ww_y[i]=h*(ww_trap[i]+ww_left[i]+ww_right[i]);

		double actual_unif_spacex=xx[1]-xx[0];
		double actual_unif_spacey=yy[1]-yy[0];
		double Lx=xx[0];
		double Ly=yy[0];

		double DU1x,DU1y;
		if ((msx%2)==0)
			DU1x=-msx/2;
		else
			DU1x=(-msx+1)/2;
		double translationx=(DU1x-(Lx/actual_unif_spacex));
		if ((msy%2)==0)
			DU1y=-msy/2;
		else
			DU1y=(-msy+1)/2;
		double translationy=(DU1y-(Ly/actual_unif_spacey));

		// Make call to finufft1d1
		nufft_opts opts; finufft_default_opts(&opts);

		std::vector<complex<double>> h_at_xxyy(msx*msy); 
		complex<double> icomp=-1;
		icomp=sqrt(icomp);
		std::vector<double> newklocsx(numlocs);
		std::vector<double> newklocsy(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newklocsx[i]=klocs_d1[i]*actual_unif_spacex;
		}
		for (int i=0;i<numlocs;i++)
		{
			newklocsy[i]=klocs_d2[i]*actual_unif_spacey;
		}
		std::vector<complex<double>> strengths(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			strengths[i]=q[i]*exp(icomp*((actual_unif_spacex*klocs_d1[i]*translationx)+(actual_unif_spacey*klocs_d2[i]*translationy))); 
		}
		int ier=finufft2d1(numlocs,newklocsx.data(),newklocsy.data(),strengths.data(),-1,newtol,msx,msy,h_at_xxyy.data(),opts);

		translationx=(-1*DU1x)*actual_unif_spacex+Lx;
		translationy=(-1*DU1y)*actual_unif_spacey+Ly;

		// need to make allww!

		int counter=0;
		std::vector<complex <double>> newh_at_xxyy(msx*msy);
		for (int i=0;i<msy;i++)
		{
			for (int j=0;j<msx;j++)
			{
				newh_at_xxyy[counter]=h_at_xxyy[counter]*ww_x[j]*ww_y[i];
				counter=counter+1;
			}
		}
		std::vector<double> newa1(numlocs);
		std::vector<double> newa2(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newa1[i]=a1[i]*actual_unif_spacex;
			newa2[i]=a2[i]*actual_unif_spacey;
		}
		std::vector<complex<double>> temp(numlocs);
		ier=finufft2d2(numlocs,newa1.data(),newa2.data(),temp.data(),1,newtol,msx,msy,newh_at_xxyy.data(),opts); 
		//do exp thing
		for (int i=0;i<numlocs;i++)
		{
			res[i]=(0.0625)*exp(i*(translationx*a1[i]+translationy*a2[i]))*temp[i];
		}
	}
	return 0;
}



