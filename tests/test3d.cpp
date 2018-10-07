#include <iostream>
#include <iomanip>
#include "sincutil.hpp"
#include "directsinc.hpp"
#include "sinctransform.hpp"
#include <vector>

int main()
{
	cout<<"---Testing: 3D---\n\n";
	double precisions[]={1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15};
	double pr,err,start;
	double klb=-1;
	double kub=1;
	double qlb=-1;
	double qub=1;
	int numlocs=50000;
	int numeval=500;
	int ifl=1;
	int quad=1;
	int s_err;

	cout<<"Sinc with "<<numlocs<<" samples:\n";

	std::vector<double> klocs_d1(numlocs);
	std::vector<double> klocs_d2(numlocs);
	std::vector<double> klocs_d3(numlocs);
	std::vector<double> a1(numlocs);
	std::vector<double> a2(numlocs);
	std::vector<double> a3(numlocs);
	std::vector<complex <double>> q(numlocs);
	randarr(klb,kub,numlocs,klocs_d1.data());
	randarr(klb,kub,numlocs,klocs_d2.data());
	randarr(klb,kub,numlocs,klocs_d3.data());
	randarr(klb,kub,numlocs,a1.data());
	randarr(klb,kub,numlocs,a2.data());
	randarr(klb,kub,numlocs,a3.data());
	randcarr(qlb,qub,numlocs,q.data());

	std::vector<complex <double>> corr(numlocs);
	start=clock();
	directsinc3d(ifl,numlocs,numeval,a1.data(),a2.data(),a3.data(),klocs_d1.data(),klocs_d2.data(),klocs_d3.data(),q.data(),corr.data()); 
	cout<<"Direct calculation: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. \n";

	for(int a=0;a<14;a++)
	{
		pr=precisions[a];		

		std::vector<complex <double>> myout(numlocs);
		start=clock();
		s_err=sinc3d(ifl,numlocs,a1.data(),a2.data(),a3.data(),klocs_d1.data(),klocs_d2.data(),klocs_d3.data(),q.data(),pr,myout.data(),quad); 
		cout<<"Runtime: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		err=getcerr(myout.data(),corr.data(),numeval); //PUT BACK
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation 
	}

	cout<<"\nSincsq with "<<numlocs<<" samples:\n";
	start=clock();
	directsincsq3d(ifl,numlocs,numeval,a1.data(),a2.data(),a3.data(),klocs_d1.data(),klocs_d2.data(),klocs_d3.data(),q.data(),corr.data()); 
	cout<<"Direct calculation: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. \n";

	for(int a=0;a<14;a++)
	{
		pr=precisions[a];

		std::vector<complex <double>> myout(numlocs);
		start=clock();
		s_err=sincsq3d(ifl,numlocs,a1.data(),a2.data(),a3.data(),klocs_d1.data(),klocs_d2.data(),klocs_d3.data(),q.data(),pr,myout.data(),quad); 
		cout<<"Runtime: "<<setprecision(6)<<(clock()-start)/(double) CLOCKS_PER_SEC<<" sec. ";

		err=getcerr(myout.data(),corr.data(),numeval);
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation
	}	

	return s_err;
}
