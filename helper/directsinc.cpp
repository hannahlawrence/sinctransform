#include <iostream>
#include <iomanip>
#include "sincutil.hpp"
#include "directsinc.hpp"

// Direct calculations of sinc transformations in 1 dimension

// macro to handle the origin
#define SINC(x) ((abs(x)<1e-8) ? 1.0 : sin(x)/x)

void directsinc1d(int ifl,int numlocs,int numeval,double *a1,double *klocs, complex<double>* q,complex<double> *ans)
// alex clarified version
{
	double pi=4*atan(1);
	double diff;
	// init. to zeros
	for(int a=0;a<numeval;a++)
	  ans[a]=0.0;
	for(int b=0;b<numeval;b++) {
	  for(int a=0;a<numlocs;a++) {
	    diff=a1[b]-klocs[a];
	    if(ifl==1)
	      diff=diff*pi;
	    ans[b] += q[a] * SINC(diff);   // note macro above, no pr needed	
	  }
	}
}

void directsincsq1d(int ifl, int numlocs,int numeval,double *a1,double *klocs, complex<double>* q,complex<double> *ans)
{
	double pi=4*atan(1);
	double diff,x;
	// init. to zeros
	for(int a=0;a<numeval;a++)
		ans[a]=0;
	for(int b=0;b<numeval;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			diff=a1[b]-klocs[a];
			if(ifl==1)
				diff=diff*pi;
			x=SINC(diff);
			ans[b] += q[a] * (x*x);
		}
	}
}

// Direct calculations of sinc transformations in 2 dimensions

void directsinc2d(int ifl,int numlocs,int numeval,double *a1,double *a2,double *klocs_d1,double *klocs_d2, complex<double>* q, complex<double> *ans)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2;
	// init. to zeros
	for(int a=0;a<numeval;a++)
		ans[a]=0;
	for(int b=0;b<numeval;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			diff_d1=a1[b]-klocs_d1[a];
			diff_d2=a2[b]-klocs_d2[a];
			if(ifl==1)
			{
				diff_d1=diff_d1*pi;
				diff_d2=diff_d2*pi;
			}
			ans[b] += q[a]*SINC(diff_d1)*SINC(diff_d2);
		}
	}
}
void directsincsq2d(int ifl, int numlocs,int numeval,double *a1,double *a2,double *klocs_d1,double *klocs_d2, complex<double>* q, complex<double> *ans)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2,x1,x2;
	// init. to zeros
	for(int a=0;a<numeval;a++)
		ans[a]=0;
	for(int b=0;b<numeval;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			diff_d1=a1[b]-klocs_d1[a];
			diff_d2=a2[b]-klocs_d2[a];
			if(ifl==1)
			{
				diff_d1=diff_d1*pi;
				diff_d2=diff_d2*pi;
			}
			x1=SINC(diff_d1);
			x2=SINC(diff_d2);
			ans[b] += q[a]*x1*x1*x2*x2;
		}
	}
}

// Direct calculations of sinc transformations in 3 dimensions

void directsinc3d(int ifl,int numlocs,int numeval,double *a1,double *a2,double *a3,double *klocs_d1,double *klocs_d2,double *klocs_d3, complex<double>* q,complex<double> *ans)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2,diff_d3;
	// init. to zeros
	for(int a=0;a<numeval;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numeval;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			diff_d1=a1[b]-klocs_d1[a];
			diff_d2=a2[b]-klocs_d2[a];
			diff_d3=a3[b]-klocs_d3[a];
			if(ifl==1)
			{
				diff_d1=diff_d1*pi;
				diff_d2=diff_d2*pi;
				diff_d3=diff_d3*pi;
			}
			ans[b] += q[a]*SINC(diff_d1)*SINC(diff_d2)*SINC(diff_d3);
		}
	}
}

void directsincsq3d(int ifl, int numlocs,int numeval,double *a1,double *a2,double *a3,double *klocs_d1,double *klocs_d2,double *klocs_d3, complex<double>* q,complex<double> *ans)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2,diff_d3,x1,x2,x3;
	// init. to zeros
	for(int a=0;a<numeval;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numeval;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			diff_d1=a1[b]-klocs_d1[a];
			diff_d2=a2[b]-klocs_d2[a];
			diff_d3=a3[b]-klocs_d3[a];
			if(ifl==1)
			{
				diff_d1=diff_d1*pi;
				diff_d2=diff_d2*pi;
				diff_d3=diff_d3*pi;
			}
			x1=SINC(diff_d1);
			x2=SINC(diff_d2);
			x3=SINC(diff_d3);
			ans[b] += q[a]*x1*x1*x2*x2*x3*x3;
		}
	}
}
