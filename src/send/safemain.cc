#include "Parton.h"
#include "Quench.h"
#include "Random.h"
#include "Wake.h"
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iostream>
#include <assert.h>
#include <vector>

#include "global.h"

//#define DO_PRINT

using std::vector;
using namespace std;

void TreeDoer(int, int, vector<Parton> &);
void read_nuclear(void);
void set_b_range(int);
void read_hydro(int);
void gxy(double &, double &, double &, numrand &);
void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha);
void do_wake(vector<Quench> quenched, vector<Wake> &wake, numrand &nr);

int main(int argc, char** argv)
{
	assert(argc==6);
	int nhat=atoi(argv[1]);
	int N=atoi(argv[2]);
	int cent=atoi(argv[3]);
	double kappa=atof(argv[4]);
	double alpha=atof(argv[5]);

	cout << " nhat= " << nhat << " N= " << N << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << endl;

	//Read nuclear thickness function
	read_nuclear();

	//Read hydro file, event averaged
	read_hydro(cent);

	//Initialize Random Seed
	numrand nr(1346);
	//cout << " rando= " << nr.rando() << endl;

	//Set impact parameter range and normalisation for overlap function of thickness functions
	set_b_range(cent);
	
	//Event Loop
	int count=0;
	do {
		//Declare partons vector
		vector<Parton> partons;
	
		//Generate the parton tree from pythia
		TreeDoer(nhat, count, partons);
		for (unsigned int i = 0; i < partons.size(); i++)
                {
                	//cout << " Parton " << i << endl;
			//partons[i].display();
			//cout << endl;
		}
		
		//Generate b,x,y
		//cout << " GetIr= " << nr.GetIr() << endl;
		gxy(x, y, b, nr);
		//cout << " GetIr= " << nr.GetIr() << endl;
		cout << " xcre= " << x << " ycre= " << y << " b= " << b << endl;
		//Do energy loss, return a vector of quenched partons
		//Create vector of quenched partons initially equal to vacuum partons
		
		//vector<Quench> quenched {Quench(partons[0])};	//Want to use this initialization form, but how?
		vector<Quench> quenched;	
		for (unsigned int i = 0; i<partons.size(); i++)
		{
			quenched.push_back ( Quench ( partons[i] ) );
		}
		
		for (unsigned int i = 0; i < quenched.size(); i++)
                {
                        //cout << " Pre-quenched Parton " << i << endl;
                        //quenched[i].display();
			//cout << endl;
                }
		do_eloss(partons, quenched, x, y, nr, kappa, alpha);
		#ifdef DO_PRINT
		for (unsigned int i = 0; i < quenched.size(); i++)
                {
			cout << " Parton " << i << endl;
                        partons[i].display();
                        cout << endl;
                        cout << " Quenched Parton " << i << endl;
                        quenched[i].display();
			cout << endl;
			if (quenched[i].GetD1()!=-1)
			{
				int d1=quenched[i].GetD1();
				int d2=quenched[i].GetD2();
				cout << " D1 Parton " << endl;
				quenched[d1].display();
				cout << endl;
				cout << " D2 Parton " << endl;
                                quenched[d2].display();
                                cout << endl;
				cout << " SumPx= " << quenched[d1].GetInhP()[0]+quenched[d2].GetInhP()[0] << endl;
				cout << " SumPy= " << quenched[d1].GetInhP()[1]+quenched[d2].GetInhP()[1] << endl;
				cout << " SumPz= " << quenched[d1].GetInhP()[2]+quenched[d2].GetInhP()[2] << endl;
				cout << " SumEn= " << quenched[d1].GetInhP()[3]+quenched[d2].GetInhP()[3] << endl;
			}
                }
		#endif
		//Do back-reaction, return a vector of wake hadrons
		vector<Wake> wake;
		do_wake(quenched,wake,nr);
		//Hadronize in pythia, return a vector of pythia hadrons

		//Print all hadrons into a file
	
		partons.clear();
		quenched.clear();	
		count+=1;
		cout << " \n NEXT EVENT \n";
	} while (count<N);

	return 0;
}
