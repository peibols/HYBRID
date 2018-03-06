#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
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
#include <algorithm>

#include "vector_operators.h"

//#define DO_PRINT
//#define DO_WAKE

using std::vector;
using namespace std;
using namespace Pythia8;

bool safe_jet(vector<Parton>, double);
void read_nuclear(void);
void set_b_range(int);
void read_hydro(int);
void gxy(double &, double &, double &, numrand &);
void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha, int tmethod);
void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);

int main(int argc, char** argv)
{
	assert(argc==7);
	int nhat=atoi(argv[1]);
	int N=atoi(argv[2]);
	int cent=atoi(argv[3]);
	double kappa=atof(argv[4]);
	double alpha=atof(argv[5]);
	int tmethod=atoi(argv[6]);

	cout << " nhat= " << nhat << " N= " << N << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << " tmethod= " << tmethod << endl;

	//Partons Output File
	char outPart[100];
	sprintf(outPart,"./PART.dat");
	ofstream part_file;
	part_file.open (outPart);

	//Wake Output File
        char outWake[100];
        sprintf(outWake,"./WAKE.dat");
        ofstream wake_file;
        wake_file.open (outWake);
	
	//Read nuclear thickness function
	read_nuclear();

	//Read hydro file, event averaged
	read_hydro(cent);

	//Initialize Random Seed
	numrand nr(1346);
	//cout << " rando= " << nr.rando() << endl;

	//Set impact parameter range and normalisation for overlap function of thickness functions
	set_b_range(cent);
	
	//Generate the parton tree from pythia
        Pythia pythia;

	//Read cmnd file
        pythia.readFile("setup_pythia.cmnd");

	//Set Random Seed
        ostringstream seedstring;
        seedstring << "Random:seed = " << 33;
        pythia.readString(seedstring.str().c_str());

	//Set pTHatMin
        double nhatmin[]={0.0,3.5,7.0,15.0,30.0,50.0,80.0,120.0,170.0,220.0,280.0};
        ostringstream hatmin;
        hatmin << "PhaseSpace:pTHatMin = " << nhatmin[nhat];
        pythia.readString(hatmin.str().c_str());

	//JetSafePt Values
	double safept[]={0.,10.,18.,34.,60.,100.,153.,227.,305.,400.,480.};

	//Print parameters used
	part_file << "pTHatMin = " << nhatmin[nhat] << ", #Events = " << N << ", Centrality = " << cent << ", Alpha = " << alpha << ", Kappa = " << kappa << " Tc= " << tmethod << endl;
	wake_file << "pTHatMin = " << nhatmin[nhat] << ", #Events = " << N << ", Centrality = " << cent << ", Alpha = " << alpha << ", Kappa = " << kappa << " Tc= " << tmethod << endl;

	//Initialize PYTHIA
        pythia.init();
	
	//Event Loop
	int count=0;
	do {
		if (!pythia.next()) continue;

		//Declare partons vector
		vector<Parton> partons;
	
		//Find Final Particles
        	for (int i = 0; i < pythia.event.size(); i++)
        	{
        	        if (pythia.event[i].isFinal())
        	        {
        	                vector<double> p;
				for (unsigned int j=1; j<4; j++) p.push_back(pythia.event[i].p()[j]);
				p.push_back(pythia.event[i].p()[0]);
				//Simply store remnants (virt 0, mother 0, daugthers -1)
        	                if (pythia.event[i].status() == 63)
        	                {
					partons.push_back ( Parton ( p, 0., pythia.event[i].m(), 0, -1, -1, pythia.event[i].id(), "rem", pythia.event[i].col(), pythia.event[i].acol(), true ) );
        	                        continue;
        	                }
	
	                        //Find first non-trivial mother, skipping carbon copies
	                        int use = i;
	                        int m1 = 0;
	                        int m2 = 0;
	                        do
	                        {
	                                m1 = pythia.event[use].mother1();
	                                m2 = pythia.event[use].mother2();
	                                if (m1==m2) use = m1;
	                        } while (m1==m2);
	
        	                //Compute virtuality
        	                double virt=sqrt(abs(pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2()));
				//cout << " virt= " << virt << endl;	
				//cout << " virt2= " << pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2() << endl;

	                        //Add it to partons array
	                        partons.push_back ( Parton ( p, virt, pythia.event[i].m(), m1, -1, -1, pythia.event[i].id(), "ps", pythia.event[i].col(), pythia.event[i].acol(), false ) );

	                        //If mother is ISR initiator, say mother is -1
	                        if (pythia.event[m1].status()==-41)
	                        {
	                                partons[partons.size()-1].SetMom(-1);
	                                partons[partons.size()-1].SetOrig("isr");
	                                partons[partons.size()-1].SetIsDone(true);
	                        }
	
	                        //If mother is hard scattering initiator, say mother is -1
	                        if (pythia.event[m1].status()==-21)
	                        {
	                                //cout << " ahora " << endl;
	                                partons[partons.size()-1].SetMom(-1);
	                                partons[partons.size()-1].SetOrig("hs");
	                                partons[partons.size()-1].SetIsDone(true);
	                        }
	                }
	        }

		//Check if there is a jet above the safe pt cut using only final partons
		bool is_safe=safe_jet(partons, safept[nhat]);	
		if (is_safe==false) continue;
		cout << " PtHat= " << pythia.info.pTHat() << endl;  

		//Reconstruct tree, using momentum of daughters to get the one for mothers
        	int changes=0;
        	do
        	{
        	        changes=1;
        	        unsigned int ps = partons.size();
        	        for (unsigned int i = 0; i < ps; i++)
        	        {
        	                if (partons[i].GetIsDone()==false)
        	                {
        	                        for (unsigned int j = 0; j < ps; j++)
        	                        {
        	                                if (partons[i].GetMom()==partons[j].GetMom() && i!=j && partons[j].GetIsDone()==false)
        	                                {
        	                                        int mom=partons[i].GetMom();
	
	                                                //Set mother momentum and virtuality
							vector<double> p=partons[i].vGetP()+partons[j].vGetP();
							double virt=sqrt(abs(pow(p[3],2.)-pow(p[0],2.)-pow(p[1],2.)-pow(p[2],2.)-pow(pythia.event[mom].m(),2.)));
	
	                                                //Find first non-trivial mother of the mother
	                                                int use = mom;
	                                                int m1 = 0;
	                                                int m2 = 0;
	                                                do {
	                                                        m1 = pythia.event[use].mother1();
	                                                        m2 = pythia.event[use].mother2();
	                                                        if (m1==m2) use = m1;
	                                                } while (m1==m2);
	
	                                                //Fill it in partons array
	                                                partons.push_back ( Parton ( p, virt, pythia.event[mom].m(), m1, i, j, pythia.event[mom].id(), "ps", pythia.event[mom].col(), pythia.event[mom].acol(), false ) );
	                                                //Update mother of daughters to point to position in partons array instead of pythia list, and declare as done
	                                                partons[i].SetMom(partons.size()-1);
	                                                partons[j].SetMom(partons.size()-1);
	                                                partons[i].SetIsDone(true);
	                                                partons[j].SetIsDone(true);
	
	                                                //If mother is ISR initiator, say mother of mother is -1
	                                                if (pythia.event[m1].status()==-41)
	                                                {
	                                                        partons[partons.size()-1].SetMom(-1);
	                                                        partons[partons.size()-1].SetOrig("isr");
	                                                        partons[partons.size()-1].SetIsDone(true);
	                                                }
	
	                                                //If mother is hard scattering initiator, say mother of mother is -1
	                                                if (pythia.event[m1].status()==-21)
	                                                {
	                                                        partons[partons.size()-1].SetMom(-1);
	                                                        partons[partons.size()-1].SetOrig("hs");
	                                                        partons[partons.size()-1].SetIsDone(true);
	                                                }
	                                                changes=0;
	                                                break;
	                                        }
	                                }
	                        }
	                }
	        } while (changes==0);

		//Generate b,x,y
		double x,y,b;
		gxy(x, y, b, nr);
		cout << " xcre= " << x << " ycre= " << y << " b= " << b << endl;
		
		//Create vector of quenched partons initially equal to vacuum partons
		//vector<Quench> quenched {Quench(partons[0])};	//Want to use this initialization form, but how?
		vector<Quench> quenched;	
		for (unsigned int i = 0; i<partons.size(); i++)
		{
			quenched.push_back ( Quench ( partons[i] ) );
		}
		
		do_eloss(partons, quenched, x, y, nr, kappa, alpha, tmethod);
		
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
	#ifdef DO_WAKE
		do_wake(quenched,partons,wake,nr);
	#endif
		//Hadronize in pythia, return a vector of pythia hadrons
		//How to do without initializing pythia again?

		//Print partons into a file
		part_file << "EVENT# = " << count << ", PARTICLES# = " << partons.size() << ", X = " << x << ", Y = " << y << ", B = " << b << endl;
		for (unsigned int i=0; i<quenched.size(); i++)
		{
			part_file << i << " " << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].vGetP()[3] << " ";
			part_file << quenched[i].vGetP()[0] << " " << quenched[i].vGetP()[1] << " " << quenched[i].vGetP()[2] << " " << quenched[i].vGetP()[3] << " ";
			part_file << partons[i].GetQ() << " " << partons[i].GetMom() << " " << partons[i].GetId() << " ";
			part_file << partons[i].GetOrig() << " " << partons[i].GetCol() << " " << partons[i].GetAcol() << endl; 
		}

	#ifdef DO_WAKE
		//Print wake hadrons into a file
		wake_file << "EVENT# = " << count << ", PARTICLES# = " << wake.size() << ", X = " << x << ", Y = " << y << ", B = " << b << endl;
		for (unsigned int i=0; i<wake.size(); i++)
		{
			wake_file << i << " " << wake[i].vGetP()[0] << " " << wake[i].vGetP()[1] << " " << wake[i].vGetP()[2] << " " << wake[i].vGetP()[3] << " ";
			wake_file << wake[i].GetMass() << " " << wake[i].GetCharge() << " " << wake[i].GetId() << " " << wake[i].GetStatus() << endl; 
		}
	#endif
		
		
		partons.clear();
		quenched.clear();	
		wake.clear();
		count+=1;
		cout << " \n NEXT EVENT \n";
	} while (count<N);

	part_file.close();
	wake_file.close();

	return 0;
}
