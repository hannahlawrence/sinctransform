#include <iostream>
#include <iomanip>
#include <chrono>
#include "sincutil.hpp"
#include "directsinc.hpp"
#include "sinctransform.hpp"
#include <vector>

int main()
{
        cout<<"---Testing 1D---\n\n";
	double precisions[]={1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12};
        int nprec = 11;
	double pr,err;
	double klb=-10;
	double kub=10;
	double qlb=-10;
	double qub=10;
	int numlocs=3000;
	int numeval=numlocs;
	int ifl=1;
	int quad=1;
	int s_err;
	cout<<"Sinc with "<<numlocs<<" samples:\n";

	std::vector<double> klocs(numlocs);
	std::vector<double> a1(numlocs);
	std::vector<complex <double>> q(numlocs);
	randarr(klb,kub,numlocs,klocs.data());
	randarr(klb,kub,numlocs,a1.data());
	randcarr(qlb,qub,numlocs,q.data());

	std::vector<complex<double>> corr(numlocs);
        auto start = chrono::system_clock::now();
	directsinc1d(ifl,numlocs,numeval,a1.data(),klocs.data(),q.data(),corr.data()); 
        chrono::duration<double> dur = chrono::system_clock::now() - start;
        cout<<setprecision(3);
	cout<<"Direct calculation: "<<dur.count()<<" sec. \n";
	for(int a=0;a<nprec;a++)
	{
		pr=precisions[a];
		std::vector<complex<double>> myout(numlocs);

                start = chrono::system_clock::now();
		s_err=sinc1d(ifl,numlocs,a1.data(),klocs.data(),q.data(),pr,myout.data(),quad);
                dur = chrono::system_clock::now()-start;
                cout<<setprecision(3);            // why needed?
		cout<<"Runtime: "<<dur.count()<<" sec. ";

		err=getcerr(myout.data(),corr.data(),numeval);
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation
	}

	cout<<"\nSincsq with "<<numlocs<<" samples:\n";

        start = chrono::system_clock::now();
	directsincsq1d(ifl,numlocs,numeval,a1.data(),klocs.data(),q.data(),corr.data());
        dur = chrono::system_clock::now() - start;
	cout<<"Direct calculation: "<<dur.count()<<" sec. \n";

	for(int a=0;a<nprec;a++)
	{
		pr=precisions[a];

		std::vector<complex<double>> myout(numlocs);

		start = chrono::system_clock::now();
		s_err=sincsq1d(ifl,numlocs,a1.data(),klocs.data(),q.data(),pr,myout.data(),quad); 
                dur = chrono::system_clock::now() - start;
                cout<<setprecision(3);            // why needed?
		cout<<"Runtime: "<<dur.count()<<" sec. ";

		err=getcerr(myout.data(),corr.data(),numeval);
		cout<<"Requested precision: "<<pr<<" "; // Requested precision
		cout<<"Error: "<<err<<"\n"; // Error compared to direct calculation	
	}
	return s_err;
}
