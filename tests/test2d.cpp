#include <iostream>
#include <iomanip>
#include "sincutil.hpp"
#include "directsinc.hpp"
#include "sinctransform.hpp"

int main()
{
	cout<<"---Testing: 2D---\n\n";
	double precisions[]={1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15};
	double pr,err,start;
	double klb=-10;
	double kub=10;
	double qlb=-10;
	double qub=10;
	int numlocs=10000;
	int ifl=1;
	int s_err;

	cout<<"Sinc with "<<numlocs<<" samples:\n";

	double* klocs_d1=(double*)malloc(sizeof(double)*numlocs);
	double* klocs_d2=(double*)malloc(sizeof(double)*numlocs);
	complex<double>* q=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	randarr(klb,kub,numlocs,klocs_d1);
	randarr(klb,kub,numlocs,klocs_d2);
	randcarr(qlb,qub,numlocs,q);

	complex<double>* corr=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	start=clock();
	directsinc2d(ifl,numlocs,klocs_d1,klocs_d2,q,corr,1e-14); 
	cout<<"Direct calculation: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. \n";


	for(int a=0;a<14;a++)
	{
		pr=precisions[a];

		complex<double>* myout=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);

		start=clock();
		s_err=sinc2d(ifl,numlocs,klocs_d1,klocs_d2,q,pr,myout); 
		cout<<"Runtime: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		err=getcerr(myout,corr,numlocs);
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation
		free(myout);	
	}
	free(corr);

	cout<<"\nSincsq with "<<numlocs<<" samples:\n";

	corr=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	start=clock();
	directsincsq2d(ifl,numlocs,klocs_d1,klocs_d2,q,corr,1e-14); 
	cout<<"Direct calculation: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. \n";

	for(int a=0;a<14;a++)
	{
		pr=precisions[a];

		complex<double>* myout=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);

		start=clock();
		s_err=sincsq2d(ifl,numlocs,klocs_d1,klocs_d2,q,pr,myout); 
		cout<<"Runtime: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		err=getcerr(myout,corr,numlocs);
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation
		free(myout);		
	}
	free(corr);
	free(klocs_d1);
	free(klocs_d2);
	free(q);
}