#include <iostream>
#include <iomanip>
#include "sincutil.hpp"
#include "directsinc.hpp"
#include "sinctransform.hpp"

int main()
{
	cout<<"Testing: 3D\n";
	double precisions[]={1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15};
	double pr,err,start;
	double klb=-1;
	double kub=1;
	double qlb=-1;
	double qub=1;
	int numlocs=10000;
	int ifl=1;
	int s_err;
	cout<<"Sinc:\n";
	for(int a=0;a<14;a++)
	{
		pr=precisions[a];

		double* klocs_d1=(double*)malloc(sizeof(double)*numlocs);
		double* klocs_d2=(double*)malloc(sizeof(double)*numlocs);
		double* klocs_d3=(double*)malloc(sizeof(double)*numlocs);	
		double* q=(double*)malloc(sizeof(double)*numlocs);
		randarr(klb,kub,numlocs,klocs_d1);
		randarr(klb,kub,numlocs,klocs_d2);
		randarr(klb,kub,numlocs,klocs_d3);
		randarr(qlb,qub,numlocs,q);

		double* myout=(double*)malloc(sizeof(double)*numlocs);
		double* corr=(double*)malloc(sizeof(double)*numlocs);
		
		start=clock();
		directsinc3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,corr,1e-14); 
		cout<<"Testprog time: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		start=clock();
		s_err=sinc3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,pr,myout); 
		cout<<"Runtime: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		err=geterr(myout,corr,numlocs);
		cout<<"Rq: "<<pr<<" "; // Requested precision
		cout<<"Err: "<<err<<"\n"; // Error compared to direct calculation
		free(myout);
		free(corr);
		free(klocs_d1);
		free(klocs_d2);
		free(klocs_d3);
		free(q);		
	}

	cout<<"Sincsq:\n";
	for(int a=0;a<14;a++)
	{
		pr=precisions[a];

		double* klocs_d1=(double*)malloc(sizeof(double)*numlocs);
		double* klocs_d2=(double*)malloc(sizeof(double)*numlocs);
		double* klocs_d3=(double*)malloc(sizeof(double)*numlocs);	
		double* q=(double*)malloc(sizeof(double)*numlocs);
		randarr(klb,kub,numlocs,klocs_d1);
		randarr(klb,kub,numlocs,klocs_d2);
		randarr(klb,kub,numlocs,klocs_d3);
		randarr(qlb,qub,numlocs,q);

		double* myout=(double*)malloc(sizeof(double)*numlocs);
		double* corr=(double*)malloc(sizeof(double)*numlocs);

		start=clock();
		directsincsq3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,corr,1e-14); 
		cout<<"Testprog time: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		start=clock();
		s_err=sincsq3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,pr,myout); 
		cout<<"Runtime: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		err=geterr(myout,corr,numlocs);
		cout<<"Rq: "<<pr<<" "; // Requested precision
		cout<<"Err: "<<err<<"\n"; // Error compared to direct calculation
		free(myout);
		free(corr);
		free(klocs_d1);
		free(klocs_d2);
		free(klocs_d3);
		free(q);		
	}
}