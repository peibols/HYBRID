#include "Parton.h"
#include "Quench.h"
#include "Hadron.h"
#include "Random.h"
#include "Wake.h"
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <sstream>

#include "global.h"
#include "vector_operators.h"

//#define DO_PRINT
#define DO_WAKE
#define JET_TOOLS

char Sfile[100];	//Source File name: one file per event

using std::vector;
using namespace std;

void read_nuclear(int, std::string);
void read_hydro(int, std::string);

void init_tree(int njob);
void do_tree(vector<Parton> &partons, double &weight, double &cross);

void gxy(double &, double &, numrand &);

void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha, int tmethod);

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);

void init_lund();
void do_lund(vector<Parton> partons, vector<Quench> quenched, vector<Hadron> &vhadrons, vector<Hadron> &qhadrons);

void jet_obs (vector<Parton> partons, vector<Quench> quenched, vector<Hadron> vhadrons, vector<Hadron> qhadrons, std::string part_or_had, std::string cent, int tmethod, double alpha, double R, std::string system, double weight, double cross, std::string pdf);

void shower_analysis(vector<Parton> partons);

int main(int argc, char** argv)
{
	//Cent: 00_05, 05_10 in PbPb
	//Cent: 0-10 in pPb

	assert(argc==9);
	int njob=atoi(argv[1]);
	int Nev=atoi(argv[2]);
	std::string cent=argv[3];
	double kappa=atof(argv[4]);
	double alpha=atof(argv[5]);
	int tmethod=atoi(argv[6]);
	std::string system=argv[7];
	std::stringstream ss(argv[8]);	
	bool do_quench;
	if (!(ss >> std::boolalpha >> do_quench)) { cout << " Wrong bool variable for do_quench "; exit(0); }
	
	std::string pdf="obsolete";

	cout << " njob= " << njob << " N= " << Nev << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << " tmethod= " << tmethod << " system= " << system << " pdf= " << pdf << endl;

	//Partons Output File
	char outPart[100];
	sprintf(outPart,"./PARTONS_spec%i.dat",njob);
	ofstream part_file;
	part_file.open (outPart);

	//Hadrons Output File
	char outHad[100];
        sprintf(outHad,"./HADRONS_spec%i.dat",njob);
	ofstream had_file;
        had_file.open (outHad);

#ifdef JET_TOOLS
	char JToutHad[100];
	if (do_quench==1) sprintf(JToutHad,"./HYBRID_Hadrons_kappa_%.3f_K_%.3f_run%i.out",alpha,kappa,njob);
	else sprintf(JToutHad,"./HYBRID_Hadrons_pp_run%i.out",njob);
	ofstream hjt_file;
	hjt_file.open (JToutHad);
	
	char JToutPart[100];
	if (do_quench==1) sprintf(JToutPart,"./HYBRID_Partons_kappa_%.3f_K_%.3f_run%i.out",alpha,kappa,njob);
	else sprintf(JToutPart,"./HYBRID_Partons_pp_run%i.out",njob);
	ofstream pjt_file;
	pjt_file.open (JToutPart);
#endif

#ifdef DO_WAKE
	//Wake Output File
        char outWake[100];
        sprintf(outWake,"./WAKE_spec%i.dat",njob);
	ofstream wake_file;
        wake_file.open (outWake);
#endif	
	//Initialize Random Seed
	numrand nr(1346+njob);
	//cout << " rando= " << nr.rando() << endl;

	//Read Ncoll File
	ifstream ncoll_file;
	char collFile[100];
	//sprintf(collFile,"/gs/project/cqn-654-ad/peibols/hybrid_ebe/ncoll_weight/%s/ncoll%s.dat",system.c_str(),cent.c_str());
	sprintf(collFile,"/home/peibols/projects/rrg-jeon-ac/peibols/hybrid/ncoll_weight/%s/ncoll%s.dat",system.c_str(),cent.c_str());
	ncoll_file.open(collFile);
	double ncoll_ev[100];
	int N=1;
	if (cent=="0-100") N=100;
	for (int i=0; i<N; i++)
	{
		int bla;
		ncoll_file >> bla >> ncoll_ev[i];
	}	

	//Print parameters used
	had_file << "#Events = " << Nev << ", Centrality = " << cent << ", Alpha = " << alpha << ", Kappa = " << kappa << " Tc= " << tmethod << endl;
	part_file << "#Events = " << Nev << ", Centrality = " << cent << ", Alpha = " << alpha << ", Kappa = " << kappa << " Tc= " << tmethod << endl;
#ifdef DO_WAKE	
	wake_file << "#Events = " << Nev << ", Centrality = " << cent << ", Alpha = " << alpha << ", Kappa = " << kappa << " Tc= " << tmethod << endl;
#endif

	//Total Event Loop
	int totcount=0;
	do {

	int Nhyd;
	if (do_quench) {	
		//Read initial energy density
		read_nuclear(totcount+1, cent);
	
		//cout << " Totcount+1= " << totcount+1 << endl;	
		//Read hydro file, event averaged
		if (totcount+1==70 && cent=="0-100" && system=="pPb") { totcount+=1; continue; } //Skip non-existent hydro file in minbias pPb
		if (totcount+1==38 && cent=="05_10" && system=="PbPb") { totcount+=1; continue; } //Skip non-existent hydro file in PbPb
		if (totcount+1==5 && cent=="0-10" && system=="pPb") { totcount+=1; continue; } //Skip non-existent hydro file in pPb
		read_hydro(totcount+1, cent);

		//Determine #events for this hydro
		Nhyd=int(ncoll_ev[totcount]*double(Nev));
		if (int((ncoll_ev[totcount]*double(Nev)-double(Nhyd))*10.)>=5) Nhyd+=1;	
		cout << " Doing " << Nhyd << " events for hydro event " << totcount+1 << endl;
		Nhyd=Nev;	//One smooth hydro setup
	}
	else {
		Nhyd=Nev;
	}
		//Hydro Event Loop
		int count=0;
		do {
		        //Generate PYTHIA tree
		        //Declare partons vector
		        vector<Parton> partons;
			if (totcount==0 && count==0) init_tree(njob);
			double weight, cross;
			do_tree(partons,weight,cross);
			//cout << " partons.size= " << partons.size() << endl;
	
			//Create vector of quenched partons initially equal to vacuum partons
			//vector<Quench> quenched {Quench(partons[0])};	//Want to use this initialization form, but how?
			vector<Quench> quenched;	
			for (unsigned int i = 0; i<partons.size(); i++)
			{
				quenched.push_back ( Quench ( partons[i] ) );
			}

			// Shower Analysis
			//shower_analysis(partons);

		double x,y;
		ofstream source_file;	
		if (do_quench) {	
			//Generate x,y
			gxy(x, y, nr);
			cout << " xcre= " << x << " ycre= " << y << endl;

			//Print creation point, and four momentum of hs partons in Source file
	        	sprintf(Sfile,"./results/sources/SOURCE_Hyd%i_Ev%i.dat",totcount,count);
	        	source_file.open(Sfile);
			source_file << "X_cre= " << x << " Y_cre= " << y << endl;
			for (unsigned int i = 0; i<partons.size(); i++)
			{
				if (partons[i].GetOrig()=="hs") source_file << "p_x = " << partons[i].vGetP()[0] << " p_y= " << partons[i].vGetP()[1] << " p_z = " << partons[i].vGetP()[2] << " p_e = " << partons[i].vGetP()[3] << endl; 
			}
			
			do_eloss(partons, quenched, x, y, nr, kappa, alpha, tmethod);
		}
		else {
			x=0., y=0.;	
		}		

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
			vector<Hadron> vhadrons, qhadrons;
			if (totcount==0 && count==0) init_lund();
			do_lund(partons,quenched,vhadrons,qhadrons);
			//cout << " Vac Hadron size= " << vhadrons.size() << " Med Hadron size= " << qhadrons.size() << endl;
		#ifndef JET_TOOLS
			//Print partons into a file
			part_file << "EVENT# = " << count << ", PARTICLES# = " << partons.size() << ", X = " << x << ", Y = " << y << ", Weight = " << weight << ", Cross = " << cross << endl;
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
		#endif
		#ifdef JET_TOOLS
		//Print partonic
		if (do_quench==1) {
		  pjt_file << "# event " << count << endl; 
		  pjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << endl; 
		  for (unsigned int i=0; i<quenched.size(); i++) {
		    if (quenched[i].GetD1()==-1 && quenched[i].vGetP()[3]!=0.) pjt_file << quenched[i].vGetP()[0] << " " << quenched[i].vGetP()[1] << " " << quenched[i].vGetP()[2] << " " << quenched[i].GetMass() << " " << quenched[i].GetId() << " " << 0 << endl;
		  }
		}
		else {
		  pjt_file << "# event " << count << endl; 
		  pjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << endl; 
		  for (unsigned int i=0; i<partons.size(); i++) {
		    if (partons[i].GetD1()==-1) pjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << 0 << endl;
		  }
		}
	        pjt_file << "end" << endl;
		#endif
		#ifndef JET_TOOLS
			//Print shower hadrons into a file
			had_file << "VAC EVENT# = " << count << ", PARTICLES# = " << vhadrons.size() << ", X = " << x << ", Y = " << y << ", Weight = " << weight << ", Cross = " << cross << endl;
			for (unsigned int i=0; i<vhadrons.size(); i++)
                        {
				had_file << i << " " << vhadrons[i].vGetP()[0] << " " << vhadrons[i].vGetP()[1] << " " << vhadrons[i].vGetP()[2] << " " << vhadrons[i].vGetP()[3] << " ";
                                had_file << vhadrons[i].GetId() << " " << vhadrons[i].GetOrig() << " " << vhadrons[i].GetCharge() << " " << vhadrons[i].GetMass() << endl;
                        }
	
		if (do_quench) {
			had_file << "MED EVENT# = " << count << ", PARTICLES# = " << qhadrons.size() << ", X = " << x << ", Y = " << y << ", Weight = " << weight << ", Cross = " << cross << endl;
                        for (unsigned int i=0; i<qhadrons.size(); i++)
                        {
                                had_file << i << " " << qhadrons[i].vGetP()[0] << " " << qhadrons[i].vGetP()[1] << " " << qhadrons[i].vGetP()[2] << " " << qhadrons[i].vGetP()[3] << " ";
                                had_file << qhadrons[i].GetId() << " " << qhadrons[i].GetOrig() << " " << qhadrons[i].GetCharge() << " " << qhadrons[i].GetMass() << endl;
                        }
		}
		#endif
	
		#ifdef DO_WAKE
		#ifndef JET_TOOLS
			//Print wake hadrons into a file
			wake_file << "EVENT# = " << count << ", PARTICLES# = " << wake.size() << ", X = " << x << ", Y = " << y << ", Weight = " << weight << ", Cross = " << cross << endl;
			for (unsigned int i=0; i<wake.size(); i++)
			{
				wake_file << i << " " << wake[i].vGetP()[0] << " " << wake[i].vGetP()[1] << " " << wake[i].vGetP()[2] << " " << wake[i].vGetP()[3] << " ";
				wake_file << wake[i].GetMass() << " " << wake[i].GetCharge() << " " << wake[i].GetId() << " " << wake[i].GetStatus() << endl; 
			}
		#endif
		#ifdef JET_TOOLS
		//Print hadronic
		if (do_quench==1) {
		  hjt_file << "# event " << count << endl; 
		  hjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << endl; 
		  for (unsigned int i=0; i<qhadrons.size(); i++) {
		    hjt_file << qhadrons[i].vGetP()[0] << " " << qhadrons[i].vGetP()[1] << " " << qhadrons[i].vGetP()[2] << " " << qhadrons[i].GetMass() << " " << qhadrons[i].GetId() << " " << 0 << endl;
		  }
		  for (unsigned int i=0; i<wake.size(); i++)
		  {
		    int ide_jt;
		    if (int(wake[i].GetStatus())==1) ide_jt=1;
		    else ide_jt=2;
			int wake_id;
			double wake_ch=wake[i].GetCharge();
			if (wake[i].GetMass()<0.5) {
			 if (wake_ch==0.) wake_id=111;
			 else if (wake_ch==1.) wake_id=211;
			 else wake_id=-211;
			}
			else {
			 if (wake_ch==1.) wake_id=2212;
			 else wake_id=-2212;
			}
		    hjt_file << wake[i].vGetP()[0] << " " << wake[i].vGetP()[1] << " " << wake[i].vGetP()[2] << " " << wake[i].vGetP()[3] << " ";
		    hjt_file << wake[i].GetMass() << " " << wake_id << " " << ide_jt << endl; 
		  }
		}
		else {
		  hjt_file << "# event " << count << endl; 
		  hjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << endl; 
		  for (unsigned int i=0; i<vhadrons.size(); i++) {
		    hjt_file << vhadrons[i].vGetP()[0] << " " << vhadrons[i].vGetP()[1] << " " << vhadrons[i].vGetP()[2] << " " << vhadrons[i].GetMass() << " " << vhadrons[i].GetId() << " " << 0 << endl;
		  }
		}
	        hjt_file << "end" << endl;
		#endif
		#endif

			//Compute jet observables, choose part or had
			std::string part_or_had = "part";
			//jet_obs(partons,quenched,vhadrons,qhadrons,part_or_had,cent,tmethod,alpha,0.3,system,weight,cross,pdf);

			vhadrons.clear();
			qhadrons.clear();
			partons.clear();
			quenched.clear();	
		#ifdef DO_WAKE
			wake.clear();
		#endif
		if (do_quench) {
			source_file.close();
		}	
			count+=1;
			//cout << setprecision(9) << " cross= " << info.sigmaGen() << endl;
			//cout << " \n NEXT EVENT \n";
		} while (count<Nhyd);
	
		totcount+=1;
	} while (totcount<N);

	part_file.close();
#ifdef DO_WAKE
        wake_file.close();
#endif

#ifdef JET_TOOLS
	pjt_file.close();
	hjt_file.close();
#endif
	return 0;
}
