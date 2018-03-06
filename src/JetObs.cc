#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <vector>
#include <fstream>

#include "Parton.h"
#include "Quench.h"
#include "Hadron.h"

using std::vector;
using namespace fastjet;
using namespace std;

void jet_obs (vector<Parton> partons, vector<Quench> quenched, vector<Hadron> vhadrons, vector<Hadron> qhadrons, std::string part_or_had, std::string cent, int tmethod, double alpha, double R, std::string system, double weight, double cross, std::string pdf);
void jet_spec(vector<PseudoJet> jets, int v_or_m, double weight);
void had_spec(vector<Hadron> hadrons, int v_or_m, double weight, std::string system);
void part_spec(vector<PseudoJet> particles, int v_or_m, double weight); 
void jet_aj_aco(vector<PseudoJet> jets, int v_or_m, int L, int S, double weight);
void select_dijet(vector<PseudoJet> jets, int&, int&);

double my_pi=3.14159265359;

double ptbin=5.;
double jetacut=2.1;
double jet_spec_hist[4][2][100];
double jet_spec_histtwo[4][2][100];

double leadcut[5]={100.,126.,158.,200.,100000.};
double subcut=25.;
double dphicut=7.*my_pi/8.;
double jet_naj_hist[4][2][17][4];
double jet_xaj_hist[4][2][17][4];
double jet_aco_hist[4][2][30];
double jet_naj_histtwo[4][2][17][4];
double jet_xaj_histtwo[4][2][17][4];
double jet_aco_histtwo[4][2][30];

double hetacut=1.;
double bincut[33]={1.,1.1,1.2,1.4,1.6,1.8,2.,2.2,2.4,3.2,4.,4.8,5.6,6.4,7.2,9.6,12.,14.4,19.2,24.,28.8,35.2,41.6,48.,60.8,73.6,86.4,103.6,180.,280.,400.,600.,1000.};
double hadspec[31][2], hadspectwo[31][2];

int nev=0;
double weightsum=0.;
double crosssum=0.;

int nr;

void jet_obs (vector<Parton> partons, vector<Quench> quenched, vector<Hadron> vhadrons, vector<Hadron> qhadrons, std::string part_or_had, std::string cent, int tmethod, double alpha, double R, std::string system, double weight, double cross, std::string pdf)
{
	//R bin
	nr=int(10.*R)-2;

	//Initialize histos if first event
	if (nev==0)
	{
		//Change eta cuts for RHIC (rest to be implemented)
		if (system == "AuAu") {
			hetacut=0.35;
			jetacut=0.35;
		}
		for (unsigned int k=0; k<4; k++)
		{
			for (unsigned int i=0; i<2; i++)
			{
        			for (unsigned int j=0; j<100; j++)
        			{
                			jet_spec_hist[k][i][j]=0.;
        				jet_spec_histtwo[k][i][j]=0.;
				}
				for (unsigned int j=0; j<17; j++)
				{
					for (unsigned int b=0; b<4; b++)
					{
						jet_naj_hist[k][i][j][b]=0.;
						jet_xaj_hist[k][i][j][b]=0.;
						jet_naj_histtwo[k][i][j][b]=0.;
                                                jet_xaj_histtwo[k][i][j][b]=0.;
					}
				}
				for (unsigned int j=0; j<30; j++)
				{
					jet_aco_hist[k][i][j]=0.;
					jet_aco_histtwo[k][i][j]=0.;
				}
			}
		}
		for (unsigned int i=0; i<2; i++) {
			for (unsigned int a=0; a<31; a++) {
    				hadspec[a][i]=0.;
    				hadspectwo[a][i]=0.;
  			}
		}
	}

	//Accumulate weight and cross
	weightsum+=weight;
	crosssum+=cross;
	
	vector<PseudoJet> particles;
	vector<PseudoJet> particlesq;
	if (part_or_had == "part") {
		for (unsigned int i=0; i<partons.size(); i++)
		{
			//Include only final particles - also remnants
			if (partons[i].GetD1()==-1)
			{ 
				particles.push_back( PseudoJet ( partons[i].vGetP()[0], partons[i].vGetP()[1], partons[i].vGetP()[2], partons[i].vGetP()[3]) );
				particlesq.push_back( PseudoJet ( quenched[i].vGetP()[0], quenched[i].vGetP()[1], quenched[i].vGetP()[2], quenched[i].vGetP()[3]) );
			}
		}
	}
	else if (part_or_had == "had") {
		for (unsigned int i=0; i<vhadrons.size(); i++)
                {
			//Include only final vac hadrons
			if (vhadrons[i].GetD1()==-1)
                        {
                                particles.push_back( PseudoJet ( vhadrons[i].vGetP()[0], vhadrons[i].vGetP()[1], vhadrons[i].vGetP()[2], vhadrons[i].vGetP()[3]) );
			}
		}
		for (unsigned int i=0; i<qhadrons.size(); i++)
                {
			//Include only final med hadrons
			if (qhadrons[i].GetD1()==-1)
                        {
                                particlesq.push_back( PseudoJet ( qhadrons[i].vGetP()[0], qhadrons[i].vGetP()[1], qhadrons[i].vGetP()[2], qhadrons[i].vGetP()[3]) );
                        }
                }
	}
	else {
		cout << " EXITING, non-standard string for particle type= " << part_or_had << endl;
		exit(0);
	}	


	// choose a jet definition
	JetDefinition jet_def(antikt_algorithm, R);

  	// run the clustering, extract the jets
  	ClusterSequence cs(particles, jet_def);
  	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(10.));

	ClusterSequence csq(particlesq, jet_def);
        vector<PseudoJet> jetsq = sorted_by_pt(csq.inclusive_jets(10.));

	/*
	for (unsigned int i=0; i<jets.size(); i++)
	{
		cout << "V jet#" << i << " has pt= " << jets[i].pt() <<  " and eta= " << jets[i].eta() << endl;
	}
	cout << "\n";
	for (unsigned int i=0; i<jetsq.size(); i++)
        {
                cout << "M jet#" << i << " has pt= " << jetsq[i].pt() <<  " and eta= " << jetsq[i].eta() << endl;
        }
	*/

	//Select Leading and Subleading Jets
	int LV=-1, SV=-1;
	select_dijet(jets, LV, SV);
	int LM=-1, SM=-1;
	select_dijet(jetsq, LM, SM);

	//Jet RAA
	jet_spec(jets, 0, weight);
	jet_spec(jetsq, 1, weight);

	//Hadron or Parton RAA
	if (part_or_had=="part") {
		part_spec(particles, 0, weight);
		part_spec(particlesq, 1, weight);
	}
	else if (part_or_had=="had") {
		had_spec(vhadrons, 0, weight, system);
		had_spec(qhadrons, 1, weight, system);
	}

	//AJ and Aco
	if (LV!=-1 && SV!=-1) jet_aj_aco(jets, 0, LV, SV, weight);
	if (LM!=-1 && SM!=-1) jet_aj_aco(jetsq, 1, LM, SM, weight);
	//cout << " LV= " << LV << " SV= " << SV << endl;
	//cout << " LM= " << LM << " SM= " << SM << endl;

	
	//HISTOS
	//Print on file every time
	ostringstream odir;
	odir << "../results_" << system.c_str() << "/" << cent.c_str() << "/" << part_or_had.c_str() << "_" << pdf.c_str() << "_tmethod_" << tmethod << "_alpha_" << int(100.*alpha);
	ostringstream pRaaF, RaaF, AcoF;
	
	//Hadron or Parton RAA
	ofstream pRAA;
        pRaaF << odir.str().c_str() << "particle_RAA.dat";
        pRAA.open (pRaaF.str().c_str());
	for (unsigned a=0; a<31; a++) {
		double praaerr=0.;
		if (hadspec[a][0]!=0.) praaerr=sqrt(pow(sqrt(hadspectwo[a][1])/hadspec[a][0],2.)+pow(hadspec[a][1]/pow(hadspec[a][0],2.)*sqrt(hadspectwo[a][0]),2.));
    		pRAA << (bincut[a]+bincut[a+1])/2. << " " << hadspec[a][0] << " " << hadspectwo[a][0] << " " << hadspec[a][1] << " " << hadspectwo[a][1] << " ";
		if (hadspec[a][0]!=0.) pRAA << hadspec[a][1]/hadspec[a][0] << " " << praaerr << " "; 
		else pRAA << 0. << " " << 0. << " ";
		pRAA << (bincut[a+1]-bincut[a]) << endl;
  	}
	pRAA << "END " << nev+1 << " " << crosssum << " " << weightsum << endl;
	pRAA.close();

	//Jet RAA
	ofstream RAA;
	RaaF << odir.str().c_str() << "RAA_R0" << int(R*10.) << ".dat";
	RAA.open (RaaF.str().c_str());
	for (unsigned int i=0; i<100; i++)
	{
		double binsize=ptbin;
		double raa_err=0.;
		double numerr=sqrt(jet_spec_histtwo[nr][1][i]);
		double denerr=sqrt(jet_spec_histtwo[nr][0][i]);
		if (jet_spec_hist[nr][0][i]!=0.) raa_err=sqrt(pow(numerr/jet_spec_hist[nr][0][i],2.)+pow(jet_spec_hist[nr][1][i]/jet_spec_hist[nr][0][i]/jet_spec_hist[nr][0][i]*denerr,2.));
		RAA << ptbin*double(i)+ptbin/2. << " " << jet_spec_hist[nr][0][i] << " " << jet_spec_histtwo[nr][0][i] << " "
		<< jet_spec_hist[nr][1][i] << " " << jet_spec_histtwo[nr][1][i] << " ";
		if (jet_spec_hist[nr][0][i]!=0.) RAA << jet_spec_hist[nr][1][i]/jet_spec_hist[nr][0][i] << " " << raa_err << " " << binsize << endl;
		else RAA << 0. << " " << 0. << " " << binsize << endl;
	}
	RAA << "END " << nev+1 << " " << crosssum << " " << weightsum << endl;
	RAA.close();

	//ASYM
	for (unsigned int j=0; j<4; j++)
	{
		ofstream ASYM;
		ostringstream AsyF;
		AsyF << odir.str().c_str() << "ASYM_R0" << int(R*10.) << "BIN" << j << ".dat";
        	ASYM.open (AsyF.str().c_str());
		double v_totasym=0., m_totasym=0.;
		for (unsigned int i=0; i<17; i++)
		{
			v_totasym+=jet_naj_hist[nr][0][i][j];
			m_totasym+=jet_naj_hist[nr][1][i][j];
		}
		for (unsigned int i=0; i<17; i++)
		{
			double binsize=1./17.;
			double vajerr=sqrt(jet_naj_histtwo[nr][0][i][j]);
			double vxjerr=sqrt(jet_xaj_histtwo[nr][0][i][j]);
			double majerr=sqrt(jet_naj_histtwo[nr][1][i][j]);
                        double mxjerr=sqrt(jet_xaj_histtwo[nr][1][i][j]);
			ASYM << double(i)*binsize+binsize/2. << " " << jet_naj_hist[nr][0][i][j] << " " << vajerr << " " << v_totasym << " ";
			ASYM << jet_naj_hist[nr][1][i][j] << " " << majerr << " " << m_totasym << " "; 
			ASYM << jet_xaj_hist[nr][0][i][j] << " " << vxjerr << " " << v_totasym << " ";
			ASYM << jet_xaj_hist[nr][1][i][j] << " " << mxjerr << " " << m_totasym << " " << binsize << endl;
		}
		ASYM << "END " << nev+1 << " " << crosssum << " " << weightsum << endl;
		ASYM.close();
	}

	//ACO
	ofstream ACO;
	AcoF << odir.str().c_str() << "ACO_R0" << int(R*10.) << ".dat";
        ACO.open (AcoF.str().c_str());
	double v_totaco=0., m_totaco=0.;
	for (unsigned int i=0; i<30; i++)
	{
		v_totaco+=jet_aco_hist[nr][0][i];
		m_totaco+=jet_aco_hist[nr][1][i];
	}
	for (unsigned int i=0; i<30; i++)
	{
		double binsize=pi/30.;
		double vacerr=sqrt(jet_aco_histtwo[nr][0][i]);
		double macerr=sqrt(jet_aco_histtwo[nr][1][i]);
		ACO << double(i)*binsize+binsize/2. << " " << jet_aco_hist[nr][0][i] << " " << vacerr << " " << v_totaco << " ";
		ACO << jet_aco_hist[nr][1][i] << " " << macerr << " " << m_totaco << " " << binsize << endl;
	}
	ACO << "END " << nev+1 << " " << crosssum << " " << weightsum << endl;
	ACO.close();

	nev+=1;
	//cout << "NEXT EVENT" << endl;	

	particles.clear();
	particlesq.clear();
	jets.clear();
	jetsq.clear();
}

//Jet Spectrum
void jet_spec(vector<PseudoJet> jets, int v_or_m, double weight)
{
	for (unsigned int i=0; i<jets.size(); i++)
        {
		if (abs(jets[i].eta())<jetacut)
                {
			int nj=int(jets[i].pt()/ptbin);		
                	jet_spec_hist[nr][v_or_m][nj]+=weight;
			jet_spec_histtwo[nr][v_or_m][nj]+=pow(weight,2.);
		}
	}
}

//Hadron Spectrum
void had_spec(vector<Hadron> hadrons, int v_or_m, double weight, std::string system)
{
	for (unsigned int i=0; i<hadrons.size(); i++) {
		if (fabs(hadrons[i].GetEta())<hetacut) {
			double pt=hadrons[i].GetPt();
			int npt=-1000;
			for (unsigned a=0; a<32; a++) {
          			if (pt>=bincut[a] && pt<bincut[a+1]) npt=int(a);
        		}
			if (npt!=-1000) {
				if (system=="AuAu" && hadrons[i].GetId()==111) hadspec[npt][v_or_m]+=weight, hadspectwo[npt][v_or_m]+=pow(weight,2.);
				if ((system=="pPb" || system=="PbPb") && hadrons[i].GetCharge()!=0.) hadspec[npt][v_or_m]+=weight, hadspectwo[npt][v_or_m]+=pow(weight,2.);
			}
		}
	}
}

//Parton Spectrum
void part_spec(vector<PseudoJet> particles, int v_or_m, double weight)
{
        for (unsigned int i=0; i<particles.size(); i++) {
                if (fabs(particles[i].eta())<hetacut) {
                        double pt=particles[i].pt();
                        int npt=-1000;
                        for (unsigned a=0; a<32; a++) {
                                if (pt>=bincut[a] && pt<bincut[a+1]) npt=int(a);
                        }
                        if (npt!=-1000) {
                                hadspec[npt][v_or_m]+=weight, hadspectwo[npt][v_or_m]+=pow(weight,2.);
                        }
                }
        }
}



//Jet Asymmetry and Acoplanarity
void jet_aj_aco(vector<PseudoJet> jets, int v_or_m, int L, int S, double weight)
{
	if (jets[L].pt()>=leadcut[0] && jets[S].pt()>=subcut)
	{
		double acodelta=acos((jets[L].px()*jets[S].px()+jets[L].py()*jets[S].py())/jets[L].pt()/jets[S].pt());
		if (acodelta<0.) acodelta+=2.*my_pi;
		jet_aco_hist[nr][v_or_m][int(acodelta/my_pi*30.)]+=weight;
		jet_aco_histtwo[nr][v_or_m][int(acodelta/my_pi*30.)]+=pow(weight,2.);
		//cout << " acodelta= " << acodelta << endl;
		if (acodelta>dphicut)
		{
			int nb=-1;
			for (unsigned int a=0; a<4; a++)
			{
				if (jets[L].pt()>=leadcut[a] && jets[L].pt()<leadcut[a+1]) nb=a;
			}
			if (nb<0) cout << " NEGA nb WTFFFFFF !!!! " << endl;
			int n_aj=int(17.*(jets[L].pt()-jets[S].pt())/(jets[L].pt()+jets[S].pt()));
			int x_aj=int(17.*jets[S].pt()/jets[L].pt());
			jet_naj_hist[nr][v_or_m][n_aj][nb]+=weight;
			jet_xaj_hist[nr][v_or_m][x_aj][nb]+=weight;	
			jet_naj_histtwo[nr][v_or_m][n_aj][nb]+=pow(weight,2.);
                        jet_xaj_histtwo[nr][v_or_m][x_aj][nb]+=pow(weight,2.);
		}
	}	
}

//Select Dijets within EtaCut
void select_dijet(vector<PseudoJet> jets, int &L, int &S)
{
	double leadpt=0., subpt=0.;
	for (unsigned int i=0; i<jets.size(); i++)
        {
                if (abs(jets[i].eta())<=jetacut)
                {
                        if (jets[i].pt()>leadpt) subpt=leadpt, leadpt=jets[i].pt(), S=L, L=int(i);
                        else if (jets[i].pt()>subpt) subpt=jets[i].pt(), S=int(i);
                }
        }
}
