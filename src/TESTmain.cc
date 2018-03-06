//./main 5 10000 00_05 0. 0.35 0

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

#include "global.h"
#include "vector_operators.h"

//#define DO_PRINT
//#define DO_WAKE

char Sfile[100];	//Source File name: one file per event

using std::vector;
using namespace std;
using namespace Pythia8;

bool safe_jet(vector<Parton>, double);
void read_nuclear(int, std::string);
void read_hydro(int, std::string);
void gxy(double &, double &, numrand &);
void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha, int tmethod);
void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);
void jet_obs (vector<Parton> partons, vector<Quench> quenched, std::string cent, int tmethod, double alpha, double R);

int main(int argc, char** argv)
{
	//Cent: 00_05, 05_10

	assert(argc==7);
	int nhat=atoi(argv[1]);
	int Nev=atoi(argv[2]);
	std::string cent=argv[3];
	double kappa=atof(argv[4]);
	double alpha=atof(argv[5]);
	int tmethod=atoi(argv[6]);

	cout << " nhat= " << nhat << " N= " << Nev << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << " tmethod= " << tmethod << endl;

	//Partons Output File
	char outPart[100];
	sprintf(outPart,"../results/%s/test_T%iA%iPARTONS.dat",cent.c_str(),tmethod,int(alpha*100.));
	ofstream part_file;
	part_file.open (outPart);

#ifdef DO_WAKE
	//Wake Output File
        char outWake[100];
        sprintf(outWake,"../results/%s/test_WAKE.dat",cent.c_str());
        ofstream wake_file;
        wake_file.open (outWake);
#endif	
	//Initialize Random Seed
	numrand nr(1350);
	//cout << " rando= " << nr.rando() << endl;

	//Read Ncoll File
	ifstream ncoll_file;
	char collFile[100];
	sprintf(collFile,"../ncoll_weight/PbPb/ncoll%s.dat",cent.c_str());
	ncoll_file.open(collFile);
	double ncoll_ev[50];
	for (unsigned int i=0; i<50; i++)
	{
		int bla;
		ncoll_file >> bla >> ncoll_ev[i];
	}	

	//Generate the parton tree from pythia
        Pythia pythia;

	//Read cmnd file
        pythia.readFile("setup_pythia.cmnd");

	//Set Random Seed
        ostringstream seedstring;
        seedstring << "Random:seed = " << 37;
        pythia.readString(seedstring.str().c_str());

	//Set pTHatMin
        double nhatmin[]={0.0,3.5,7.0,15.0,30.0,50.0,80.0,120.0,170.0,220.0,280.0};
        ostringstream hatmin;
        hatmin << "PhaseSpace:pTHatMin = " << nhatmin[nhat];
        pythia.readString(hatmin.str().c_str());

	//JetSafePt Values
	//double safept[]={0.,10.,18.,34.,60.,100.,153.,227.,305.,400.,480.};
	double safept[]={0.,10.,18.,34.,70.,100.,153.,227.,305.,400.,480.};

	//Print parameters used
	part_file << "pTHatMin = " << nhatmin[nhat] << ", #Events = " << Nev << ", Centrality = " << cent << ", Alpha = " << alpha << ", Kappa = " << kappa << " Tc= " << tmethod << endl;
#ifdef DO_WAKE	
	wake_file << "pTHatMin = " << nhatmin[nhat] << ", #Events = " << Nev << ", Centrality = " << cent << ", Alpha = " << alpha << ", Kappa = " << kappa << " Tc= " << tmethod << endl;
#endif
	//Initialize PYTHIA
        pythia.init();

	//Total Event Loop
	int N=50;
	int totcount=0;
	do {
		//cout << " Totcount+1= " << totcount+1 << endl;	
		if (totcount+1==38 && cent=="05_10") { totcount+=1; continue; }
		//Read initial energy density
		read_nuclear(totcount+1, cent);
		
		//Read hydro file, event averaged
		read_hydro(totcount+1, cent);

		//Determine #events for this hydro
		int Nhyd=int(ncoll_ev[totcount]*double(Nev));
		if (int((ncoll_ev[totcount]*double(Nev)-double(Nhyd))*10.)>=5) Nhyd+=1;	
		//Force to do 50 jobs per hydro - for mayank sources
		Nhyd=50;
		cout << " Doing " << Nhyd << " events for hydro event " << totcount+1 << endl;
		//Hydro Event Loop
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
			//cout << " PtHat= " << pythia.info.pTHat() << endl;  
	
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
	
			//Generate x,y
			double x,y;
			gxy(x, y, nr);
			cout << " xcre= " << x << " ycre= " << y << endl;
			//Print creation point, and four momentum of hs partons in Source file
			ofstream source_file;
	        	sprintf(Sfile,"../results/%s/sources/hyd%i/SOURCE%i.dat",cent.c_str(),totcount+1,count);
	        	source_file.open(Sfile);
			source_file << "X_cre= " << x << " Y_cre= " << y << endl;
			for (unsigned int i = 0; i<partons.size(); i++)
			{
				if (partons[i].GetOrig()=="hs") source_file << "p_x = " << partons[i].vGetP()[0] << " p_y= " << partons[i].vGetP()[1] << " p_z = " << partons[i].vGetP()[2] << " p_e = " << partons[i].vGetP()[3] << endl; 
			}
			
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
	
		#ifdef DO_WAKE
			//Do back-reaction, return a vector of wake hadrons
			vector<Wake> wake;
			do_wake(quenched,partons,wake,nr);
		#endif
			//Hadronize in pythia, return a vector of pythia hadrons
			//How to do without initializing pythia again?
	
			//Print partons into a file
			part_file << "EVENT# = " << count << " PARTICLES# = " << partons.size() << " NHYD# = " << totcount+1 << " X= " << x << " Y= " << y << endl;
			for (unsigned int i=0; i<quenched.size(); i++)
			{
				//Eloss sanity check
				/*
				cout << " lambPx= " << quenched[i].vGetP()[0]/partons[i].vGetP()[0] << " ";
				cout << " lambPy= " << quenched[i].vGetP()[1]/partons[i].vGetP()[1] << " ";
				cout << " lambPz= " << quenched[i].vGetP()[2]/partons[i].vGetP()[2] << " ";
				cout << " lambEn= " << quenched[i].vGetP()[3]/partons[i].vGetP()[3] << endl;
				*/
				part_file << i << " " << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].vGetP()[3] << " ";
				part_file << quenched[i].vGetP()[0] << " " << quenched[i].vGetP()[1] << " " << quenched[i].vGetP()[2] << " " << quenched[i].vGetP()[3] << " ";
				part_file << partons[i].GetQ() << " " << partons[i].GetMom() << " " << partons[i].GetD1() << " " << partons[i].GetD2() << " ";
				part_file << partons[i].GetId() << " " << partons[i].GetOrig() << " " << partons[i].GetCol() << " " << partons[i].GetAcol() << endl; 
			}
	
		#ifdef DO_WAKE
			//Print wake hadrons into a file
			wake_file << "EVENT# = " << count << ", PARTICLES# = " << wake.size() << ", X = " << x << ", Y = " << y << endl;
			for (unsigned int i=0; i<wake.size(); i++)
			{
				wake_file << i << " " << wake[i].vGetP()[0] << " " << wake[i].vGetP()[1] << " " << wake[i].vGetP()[2] << " " << wake[i].vGetP()[3] << " ";
				wake_file << wake[i].GetMass() << " " << wake[i].GetCharge() << " " << wake[i].GetId() << " " << wake[i].GetStatus() << endl; 
			}
		#endif

			//Compute jet observables
			//jet_obs(partons,quenched,cent,tmethod+1000,alpha,0.3);			
			
			partons.clear();
			quenched.clear();	
		#ifdef DO_WAKE
			wake.clear();
		#endif
			source_file.close();	
			count+=1;
			//cout << " \n NEXT EVENT \n";
		} while (count<Nhyd);
	
		totcount+=1;
	} while (totcount<N);

	part_file.close();
#ifdef DO_WAKE
        wake_file.close();
#endif

	return 0;
}
