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

// Contains sinc1d and sincsq1d

int sinc1d(int ifl,int numlocs,double *a1_,double *klocs_, complex<double> *q,double tol, complex<double> *res,int quad)
{
	/*  
	Computes res[j] = sum sinc(a1[k]-klocs_[j]) * q[j]
	             	   k

	Inputs:
		ifl = sinc convention
			0: sinc(x) = sin(x)/x
			1: sinc(x)=sin(pi*x)/(pi*x)
		a1_ = (real) evaluation locations
		klocs_ = (real) sample locations
		q = sample strengths
		tol = requested precision
		quad = quadrature mode; 0 for legendre, 1 for corrected trapezoidal

	Returns:
		0: if success
		Other: error code, as returned by finufft (see finufft documentation)
	*/

	double pi=4*atan(1); 
	double lowesttol=1e-15; //best precision that finufft can take
	double newtol=(tol > lowesttol) ? tol : lowesttol;

	std::vector<double> klocs(numlocs);
	std::vector<double> a1(numlocs);

	// Find rkmax
	double rkmax=0;
	double temp,temp2;
	std::vector<complex<double>> qc(numlocs);

	for(int a=0;a<numlocs;a++)
	{
		temp=abs(klocs_[a]);
		temp2=abs(a1_[a]);
		if (temp>rkmax)
			rkmax=temp;
		if (temp2>rkmax)
			rkmax=temp2;
		qc[a]=q[a];
		klocs[a]=klocs_[a]; // to ensure that program does not alter input
		a1[a]=a1_[a];
	}

	// adjust based on ifl convention
	if (ifl==1)
	{
		rkmax=pi*rkmax;
		for(int a=0;a<numlocs;a++)
		{
			klocs[a]=klocs[a]*pi;
			a1[a]=a1[a]*pi;
		}
	}

	// set nx
	double rsamp;
	int nx;

	if (quad==0)
	{
		rsamp=2;
		nx=ceil(rsamp*round(rkmax+3));
		// get lgwt weights
		std::vector<double> xx(nx);
		std::vector<double> ww(nx);
		for(int a=1;a<=nx;a++)
		{
			fastgl::QuadPair p = fastgl::GLPair((size_t)nx,(size_t)a);
			xx[a-1]=p.x();
			ww[a-1]=p.weight;
		}

		// call 1 to finufft; h_at_xx will be complex
                nufft_opts* opts = new nufft_opts;
                finufft_default_opts(opts);
		if (newtol<1e-9)
		{
			opts->upsampfac=2;
		}
		std::vector<complex<double>> h_at_xx(nx);

		int ier1=finufft1d3(numlocs,klocs.data(),qc.data(),-1,newtol,nx,xx.data(),h_at_xx.data(),opts);
		if (ier1 >1)
		{
			cout<<"Error: call 1 to finufft1d3 failed: "<<ier1<<"\n";
			return ier1;
		}

		// call 2 to finufft
		// source strengths will need to be complex!
		std::vector<complex<double>> weighted(nx);
		for(int a=0;a<nx;a++)
		{
			weighted[a]=0.5*h_at_xx[a]*ww[a];
		}
		std::vector<complex<double>> wtrans(numlocs);
		int ier2=finufft1d3(nx,xx.data(),weighted.data(),1,newtol,numlocs,a1.data(),wtrans.data(),opts);
		if (ier2 >1)
		{
			cout<<"Error: call 2 to finufft1d3 failed: "<<ier2<<"\n";
			return ier2;
		}

		//make wtrans real (double array) ONLY if q is real
		for(int a=0;a<numlocs;a++)
		{
			res[a]=wtrans[a];
		}
	}
	else
	{
		rsamp=3;
		nx=ceil(rsamp*round(rkmax+3));
		int a=-1;
		int b=1;
		int e=25; // increase (up to 60) to impose higher accuracy; will increase runtime
		if ((nx % 2) != 0)
			nx=nx+1; // ensure even so that 0 is a quadrature point
		int n=nx;
		double h=(double (b-a))/n;
		int aind=e;
		int bind=aind+n;
		int numunif=n+(2*e)+1;
		// Get constants, double vector of length e 
		std::vector<double> xx(numunif);
		std::vector<double> ww(numunif);
		for (int i=0;i<numunif;i++) // Initialize all to 0
			ww[i]=0;
		double curr=a-(e*h);
		for (int i=0;i<numunif;i++)
		{
			xx[i]=curr;
			curr=curr+h;
		}
		ww[aind]=0.5;
		ww[bind]=0.5;
		for (int i=(aind+1);i<bind;i++)
		{
			ww[i]=1;
		}
		std::vector<double> constants(e);
		for (int i=0;i<e;i++)
			constants[i]=allconstants[e-1][i];

		for (int k=1;k<=e;k++) 
		{
			ww[aind-k] = ww[aind-k] - constants[k-1];
		    ww[aind+k] = ww[aind+k] + constants[k-1];
		    ww[bind-k] = ww[bind-k] + constants[k-1];
		    ww[bind+k] = ww[bind+k] - constants[k-1];
		}

		for (int i=0;i<numunif;i++)
			ww[i]=ww[i]*h;
		// Transform sources and modes to proper form for finufft1d1 and finufft1d2
		int ms=numunif;
		double actual_unif_space=xx[1]-xx[0];
		double L=xx[0];
		double DU1;
		if ((ms%2)==0)
			DU1=-ms/2;
		else
			DU1=(-ms+1)/2;
		double translation=(DU1-(L/actual_unif_space));

		// Make call to finufft1d1
                nufft_opts* opts; opts = new nufft_opts;
                finufft_default_opts(opts);
		std::vector<complex<double>> h_at_xx(ms);
		complex<double> icomp=-1;
		icomp=sqrt(icomp);
		std::vector<double> newklocs(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newklocs[i]=klocs[i]*actual_unif_space;
		}
		std::vector<complex<double>> strengths(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			strengths[i]=q[i]*exp(icomp*actual_unif_space*klocs[i]*translation); 
		}
		
		finufft1d1(numlocs,newklocs.data(),strengths.data(),-1,newtol,ms,h_at_xx.data(),opts);
		translation=((-1*DU1)*actual_unif_space)+L;
		std::vector<complex<double>> temp(numlocs);
		// Make call to finufft1d2, write to temp
		std::vector<double> newa1(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newa1[i]=a1[i]*actual_unif_space;
		}
		std::vector<complex<double>> newh_at_xx(ms);
		for (int i=0;i<ms;i++)
		{
			newh_at_xx[i]=h_at_xx[i]*ww[i];
		}

		finufft1d2(numlocs,newa1.data(),temp.data(),1,newtol,ms,newh_at_xx.data(),opts);
		for (int i=0;i<numlocs;i++)
		{
			res[i]=0.5*exp(i*translation*a1[i])*temp[i];
		}
	}
	return 0;
}

int sincsq1d(int ifl,int numlocs,double *a1_,double *klocs_, complex<double> *q,double tol, complex<double> *res,int quad)
{
	/*  
	Computes res[j] = sum sinc^2(klocs_[k]-klocs_[j]) * q[j]
	             	   k

	Inputs:
		ifl = sinc convention
			0: sinc(x) = sin(x)/x
			1: sinc(x)=sin(pi*x)/(pi*x)
		a1_ = (real) evaluation locations
		klocs_ = (real) sample locations
		q = sample strengths
		tol = requested precision
		quad = quadrature mode; 0 for legendre, 1 for corrected trapezoidal

	Returns:
		0: if success
		Other: error code, as returned by finufft (see finufft documentation)
	*/
	double pi=4*atan(1);
	double lowesttol=1e-15; //hopefully 1e-16, which is best precision that finufft can take?
	double newtol=(tol > lowesttol) ? tol : lowesttol;

	std::vector<double> klocs(numlocs);
	std::vector<double> a1(numlocs);
	
	// find rkmax
	double rkmax=0;
	double temp;
	std::vector<complex<double>> qc(numlocs);
	for(int a=0;a<numlocs;a++)
	{
		temp=abs(klocs_[a]);
		if (temp>rkmax)
		{
			rkmax=temp;
		}
		qc[a]=q[a];
		klocs[a]=klocs_[a]; // to ensure that program does not alter input
		a1[a]=a1_[a];
	}

	// adjust based on ifl convention
	if (ifl==1)
	{
		rkmax=pi*rkmax;
		for(int a=0;a<numlocs;a++)
		{
			klocs[a]=klocs[a]*pi;
			a1[a]=a1[a]*pi;
		}
	}

	// set nx
	double rsamp;
	int nx;
	
	if (quad==0) // Legendre Quadrature
	{
		rsamp=2;
		nx=ceil(rsamp*round(rkmax+3));
		// get lgwt weights
		std::vector<double> xx(2*nx);
		std::vector<double> ww(2*nx);
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
                nufft_opts* opts = new nufft_opts;
                finufft_default_opts(opts);	
		if (newtol<1e-9)
		{
			opts->upsampfac=2;
		}	
		
		std::vector<complex<double>> h_at_xx(2*nx);
		int ier1=finufft1d3(numlocs,klocs.data(),qc.data(),-1,newtol,2*nx,xx.data(),h_at_xx.data(),opts);
		if (ier1 >1)
		{
			cout<<"Issue: call 1 to finufft1d3 failed\n";
			return ier1;
		}

		// call 2 to finufft
		// source strengths will need to be complex!
		std::vector<complex<double>> weighted(2*nx);
		for(int a=0;a<2*nx;a++)
		{
			weighted[a]=0.25*h_at_xx[a]*ww[a]*(2-abs(xx[a]));
		}
		std::vector<complex<double>> wtrans(numlocs);
		int ier2=finufft1d3(2*nx,xx.data(),weighted.data(),1,newtol,numlocs,a1.data(),wtrans.data(),opts);
		if (ier2 >1)
		{
			cout<<"Issue: call 2 to finufft1d3 failed\n";
			return ier2;
		}
		//make wtrans real (double array) ONLY if q is real
		for(int a=0;a<numlocs;a++)
		{
			res[a]=wtrans[a];
		}
	}
	else // trapezoidal rule
	{

		rsamp=3;
		nx=ceil(rsamp*round(rkmax+3));
		int a=-2;
		int b=0;
		int e=21; // increase (up to 60) to impose higher accuracy; will increase runtime
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

		for (int i=(aind+1);i<bind;i++) // remove this?
		{
			ww[i]=1;
		}
		std::vector<double> constants(e);
		for (int i=0;i<e;i++)
			constants[i]=allconstants[e-1][i];

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
			ww[i]=h*(ww_trap[i]+ww_left[i]+ww_right[i]);

		// Transform sources and modes to proper form for finufft1d1 and finufft1d2
		int ms=numunif;
		double actual_unif_space=xx[1]-xx[0];
		double L=xx[0];
		double DU1;
		if ((ms%2)==0)
			DU1=-ms/2;
		else
			DU1=(-ms+1)/2;
		double translation=(DU1-(L/actual_unif_space));

		// Make call to finufft1d1
		nufft_opts* opts; opts = new nufft_opts;
                finufft_default_opts(opts);
		std::vector<complex<double>> h_at_xx(ms);
		complex<double> icomp=-1;
		icomp=sqrt(icomp);
		std::vector<double> newklocs(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newklocs[i]=klocs[i]*actual_unif_space;
		}
		std::vector<complex<double>> strengths(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			strengths[i]=q[i]*exp(icomp*actual_unif_space*klocs[i]*translation); 
		}
		
		finufft1d1(numlocs,newklocs.data(),strengths.data(),-1,newtol,ms,h_at_xx.data(),opts);
		translation=((-1*DU1)*actual_unif_space)+L;
		std::vector<complex<double>> temp(numlocs);
		// Make call to finufft1d2, write to temp
		std::vector<double> newa1(numlocs);
		for (int i=0;i<numlocs;i++)
		{
			newa1[i]=a1[i]*actual_unif_space;
		}
		std::vector<complex<double>> newh_at_xx(ms);
		for (int i=0;i<ms;i++)
		{
			newh_at_xx[i]=h_at_xx[i]*ww[i];
		}

		finufft1d2(numlocs,newa1.data(),temp.data(),1,newtol,ms,newh_at_xx.data(),opts);
		for (int i=0;i<numlocs;i++)
		{
			res[i]=0.25*exp(i*translation*a1[i])*temp[i];
		}
	}
	
	return 0;
}

