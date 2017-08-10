#include <iostream>
#include <iomanip>
#include "sincutil.hpp"
#include "directsinc.hpp"

// Direct calculations of sinc transformations in 1 dimension

void directsinc1d(int ifl,int numlocs,double *klocs, double* q,double *ans, double pr)
{
	double pi=4*atan(1);
	double diff;
	// init. to zeros
	for(int a=0;a<numlocs;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numlocs;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			if(equals(a,b,1e-5)| equals(klocs[a],klocs[b],pr))
			{
				ans[b]=ans[b]+q[a];
			}
			else
			{
				diff=klocs[b]-klocs[a];
				if(ifl==1)
				{
					diff=diff*pi;
				}
				ans[b]=ans[b]+(q[a]*sin(diff)/diff);
			}
		}
	}
}
void directsincsq1d(int ifl, int numlocs,double *klocs, double* q,double *ans, double pr)
{
	double pi=4*atan(1);
	double diff;
	// init. to zeros
	for(int a=0;a<numlocs;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numlocs;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			if(equals(a,b,1e-5) || equals(klocs[a],klocs[b],pr))
			{
				ans[b]=ans[b]+q[a];
			}
			else
			{
				diff=klocs[b]-klocs[a];
				if(ifl==1)
				{
					diff=diff*pi;
				}
				ans[b]=ans[b]+(q[a]*pow(sin(diff)/diff,2));
			}
		}
	}
}

// Direct calculations of sinc transformations in 2 dimensions

void directsinc2d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2, double* q, double *ans, double pr)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2,a1,a2;
	// init. to zeros
	for(int a=0;a<numlocs;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numlocs;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			if(equals(a,b,1e-5) || equals(klocs_d1[a],klocs_d1[b],pr) || equals(klocs_d2[a],klocs_d2[b],pr))
			{
				ans[b]=ans[b]+q[a];
			}
			else
			{
				diff_d1=klocs_d1[b]-klocs_d1[a];
				diff_d2=klocs_d2[b]-klocs_d2[a];
				if(ifl==1)
				{
					diff_d1=diff_d1*pi;
					diff_d2=diff_d2*pi;
				}
				if(diff_d1==0)
				{
					a1=1;
				}
				else
				{
					a1=sin(diff_d1)/diff_d1;
				}
				if(diff_d2==0)
				{
					a2=1;
				}
				else
				{
					a2=sin(diff_d2)/diff_d2;
				}
				ans[b]=ans[b]+(q[a]*a1*a2);
			}
		}
	}
}
void directsincsq2d(int ifl, int numlocs,double *klocs_d1,double *klocs_d2, double* q, double *ans, double pr)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2,a1,a2;
	// init. to zeros
	for(int a=0;a<numlocs;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numlocs;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			if(equals(a,b,1e-5) || equals(klocs_d1[a],klocs_d1[b],pr) || equals(klocs_d2[a],klocs_d2[b],pr))
			{
				ans[b]=ans[b]+q[a];
			}
			else
			{
				diff_d1=klocs_d1[b]-klocs_d1[a];
				diff_d2=klocs_d2[b]-klocs_d2[a];
				if(ifl==1)
				{
					diff_d1=diff_d1*pi;
					diff_d2=diff_d2*pi;
				}
				if(diff_d1==0)
				{
					a1=1;
				}
				else
				{
					a1=pow(sin(diff_d1)/diff_d1,2);
				}
				if(diff_d2==0)
				{
					a2=1;
				}
				else
				{
					a2=pow(sin(diff_d2)/diff_d2,2);
				}
				ans[b]=ans[b]+(q[a]*a1*a2);
			}
		}
	}
}

// Direct calculations of sinc transformations in 3 dimensions

void directsinc3d(int ifl,int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3, double* q,double *ans, double pr)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2,diff_d3,a1,a2,a3;
	// init. to zeros
	for(int a=0;a<numlocs;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numlocs;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			if(equals(a,b,1e-5) || equals(klocs_d1[a],klocs_d1[b],pr) || equals(klocs_d2[a],klocs_d2[b],pr) || equals(klocs_d3[a],klocs_d3[b],pr))
			{
				ans[b]=ans[b]+q[a];
			}
			else
			{
				diff_d1=klocs_d1[b]-klocs_d1[a];
				diff_d2=klocs_d2[b]-klocs_d2[a];
				diff_d3=klocs_d3[b]-klocs_d3[a];
				if(ifl==1)
				{
					diff_d1=diff_d1*pi;
					diff_d2=diff_d2*pi;
					diff_d3=diff_d3*pi;
				}
				if(diff_d1==0)
				{
					a1=1;
				}
				else
				{
					a1=sin(diff_d1)/diff_d1;
				}
				if(diff_d2==0)
				{
					a2=1;
				}
				else
				{
					a2=sin(diff_d2)/diff_d2;
				}
				if(diff_d3==0)
				{
					a3=1;
				}
				else
				{
					a3=sin(diff_d3)/diff_d3;
				}
				ans[b]=ans[b]+(q[a]*a1*a2*a3);
			}
		}
	}
}
void directsincsq3d(int ifl, int numlocs,double *klocs_d1,double *klocs_d2,double *klocs_d3, double* q,double *ans, double pr)
{
	double pi=4*atan(1);
	double diff_d1,diff_d2,diff_d3,a1,a2,a3;
	// init. to zeros
	for(int a=0;a<numlocs;a++)
	{
		ans[a]=0;
	}
	for(int b=0;b<numlocs;b++)
	{
		for(int a=0;a<numlocs;a++)
		{
			if(equals(a,b,1e-5) || equals(klocs_d1[a],klocs_d1[b],pr) || equals(klocs_d2[a],klocs_d2[b],pr) || equals(klocs_d3[a],klocs_d3[b],pr))
			{
				ans[b]=ans[b]+q[a];
			}
			else
			{
				diff_d1=klocs_d1[b]-klocs_d1[a];
				diff_d2=klocs_d2[b]-klocs_d2[a];
				diff_d3=klocs_d3[b]-klocs_d3[a];
				if(ifl==1)
				{
					diff_d1=diff_d1*pi;
					diff_d2=diff_d2*pi;
					diff_d3=diff_d3*pi;
				}
				if(diff_d1==0)
				{
					a1=1;
				}
				else
				{
					a1=pow(sin(diff_d1)/diff_d1,2);
				}
				if(diff_d2==0)
				{
					a2=1;
				}
				else
				{
					a2=pow(sin(diff_d2)/diff_d2,2);
				}
				if(diff_d3==0)
				{
					a3=1;
				}
				else
				{
					a3=pow(sin(diff_d3)/diff_d3,2);
				}
				ans[b]=ans[b]+(q[a]*a1*a2*a3);
			}
		}
	}
}
