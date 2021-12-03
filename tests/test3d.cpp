#include <iostream>
#include <iomanip>
#include <chrono>
#include "sincutil.hpp"
#include "directsinc.hpp"
#include "sinctransform.hpp"

int main()
{
	cout<<"---Testing: 3D---\n\n";
	double precisions[]={1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12};
        int nprec = 11;
	double pr,err;
	double klb=-1;
	double kub=1;
	double qlb=-1;
	double qub=1;
	int numlocs=5000;
	int ifl=1;
	int s_err;

	cout<<"Sinc with "<<numlocs<<" samples:\n";

	double* klocs_d1=(double*)malloc(sizeof(double)*numlocs);
	double* klocs_d2=(double*)malloc(sizeof(double)*numlocs);
	double* klocs_d3=(double*)malloc(sizeof(double)*numlocs);	
	complex<double>* q=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
	randarr(klb,kub,numlocs,klocs_d1);
	randarr(klb,kub,numlocs,klocs_d2);
	randarr(klb,kub,numlocs,klocs_d3);
	randcarr(qlb,qub,numlocs,q);

	complex<double>* corr=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
        auto start = chrono::system_clock::now();
	directsinc3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,corr,1e-14); 
        chrono::duration<double> dur = chrono::system_clock::now() - start;
        cout<<setprecision(3);
	cout<<"Direct calculation: "<<dur.count()<<" sec. \n";

	for(int a=0;a<nprec;a++)
	{
		pr=precisions[a];		

		complex<double>* myout=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
                start = chrono::system_clock::now();
		s_err=sinc3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,pr,myout); 
                dur = chrono::system_clock::now() - start;
                cout<<setprecision(3);
		cout<<"Runtime: "<<dur.count()<<" sec. ";

		err=getcerr(myout,corr,numlocs); //PUT BACK
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation 
		free(myout);
			
	}

	free(corr);	

	cout<<"\nSincsq with "<<numlocs<<" samples:\n";
	corr=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
        start = chrono::system_clock::now();
	directsincsq3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,corr,1e-14);
        dur = chrono::system_clock::now() - start;
        cout<<setprecision(3);
	cout<<"Direct calculation: "<<dur.count()<<" sec. \n";

	for(int a=0;a<nprec;a++)
	{
		pr=precisions[a];

		complex<double>* myout=(complex<double>*)malloc(sizeof(complex<double>)*numlocs);
                start = chrono::system_clock::now();
		s_err=sincsq3d(ifl,numlocs,klocs_d1,klocs_d2,klocs_d3,q,pr,myout); 
                dur = chrono::system_clock::now() - start;
                cout<<setprecision(3);
		cout<<"Runtime: "<<dur.count()<<" sec. ";

		err=getcerr(myout,corr,numlocs);
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation
		free(myout);
	}
	free(klocs_d1);
	free(klocs_d2);
	free(klocs_d3);
	free(q);
	free(corr);		

	return s_err;
}
