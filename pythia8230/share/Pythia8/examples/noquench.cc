#include "fastjet/ClusterSequence.hh"
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"

#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <assert.h>

using namespace Pythia8;
using namespace fastjet;

double my_pi=3.14159265359;

void select_dijet(vector<PseudoJet> jets, int&, int&);
double dphicut, jetacut, leadpt, subpt;

//Args: N, seed, sqrts (in GeV), ISR(off=0, on=1), MPI(off=0, on=1)
int main(int argc, char** argv) {

  assert(argc==6);

  int nEvent  = atoi(argv[1]);
  int isisr = atoi(argv[4]);
  int ismpi = atoi(argv[5]);

  //Jet Spec
  double jbincut[51]={0.};
  double ptbin[51]={0.};
  jbincut[0]=0.;
  ptbin[0]=5.;
  for (unsigned a=1; a<=20; a++) {
    ptbin[a]=5.;
    jbincut[a]=jbincut[a-1]+ptbin[a];
  }
  for (unsigned a=21; a<=30; a++) {
    ptbin[a]=10.;
    jbincut[a]=jbincut[a-1]+ptbin[a];
  }
  for (unsigned a=31; a<=36; a++) {
    ptbin[a]=50.;
    jbincut[a]=jbincut[a-1]+ptbin[a];
  }
  for (unsigned a=37; a<=41; a++) {
    ptbin[a]=100.;
    jbincut[a]=jbincut[a-1]+ptbin[a];
  }
  for (unsigned a=42; a<=50; a++) {
    ptbin[a]=500.;
    jbincut[a]=jbincut[a-1]+ptbin[a];
  }
  double ptmin=0.;
  double jetspec[500][4], jetspectwo[500][4];
  for (unsigned a=0; a<500; a++) {
    for (unsigned c=0; c<4; c++) {
        jetspec[a][c]=0.;
        jetspectwo[a][c]=0.;
    }
  }

  //Had Spec
  double hadspec[31], hadspectwo[31];
  for (unsigned a=0; a<31; a++) {
    hadspec[a]=0.;
    hadspectwo[a]=0.;
  }

  //FF
  double fragcut=1.;
  double ff[20][500][4], fftwo[20][500][4];
  double jetff[500][4], jetfftwo[500][4];
  for (unsigned c=0; c<4; c++) {
    for (unsigned a=0; a<500; a++) {
      jetff[a][c]=0., jetfftwo[a][c]=0.;
      for (unsigned b=0; b<20; b++) {
        ff[b][a][c]=0., fftwo[b][a][c]=0.;
      }
    }
  }

  //Asymmetry
  double jet_naj_hist[17][4];
  double jet_xaj_hist[17][4];
  double jet_naj_histtwo[17][4];
  double jet_xaj_histtwo[17][4];
  for (unsigned a=0; a<17; a++) {
    for (unsigned b=0; b<4; b++) {
      jet_naj_hist[a][b]=0.;
      jet_xaj_hist[a][b]=0.;
      jet_naj_histtwo[a][b]=0.;
      jet_xaj_histtwo[a][b]=0.;
    }
  }

  double bincut[32];
  bincut[0]=1.;
  bincut[1]=1.1;
  bincut[2]=1.2;
  bincut[3]=1.4;
  bincut[4]=1.6;
  bincut[5]=1.8;
  bincut[6]=2.;
  bincut[7]=2.2;
  bincut[8]=2.4;
  bincut[9]=3.2;
  bincut[10]=4.;
  bincut[11]=4.8;
  bincut[12]=5.6;
  bincut[13]=6.4;
  bincut[14]=7.2;
  bincut[15]=9.6;
  bincut[16]=12.;
  bincut[17]=14.4;
  bincut[18]=19.2;
  bincut[19]=24.;
  bincut[20]=28.8;
  bincut[21]=35.2;
  bincut[22]=41.6;
  bincut[23]=48.;
  bincut[24]=60.8;
  bincut[25]=73.6;
  bincut[26]=86.4;
  bincut[27]=103.6;
  bincut[28]=180.;
  bincut[29]=280.;
  bincut[30]=400.;
  bincut[31]=600.;

  Pythia pythia;
  Info& info = pythia.info;
  Event& event = pythia.event;

  string pdfSet = "LHAPDF5:cteq6ll.LHpdf";
  pythia.readString("PDF:pSet = " + pdfSet);
  pythia.readString("PDF:extrapolate = on");
  
  pythia.readString("Beams:idB = 2212");
  pythia.readString("SoftQCD:all = off");
  pythia.readString("HardQCD:all = on");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("PartonLevel:FSR = off");
  if (isisr==1) pythia.readString("PartonLevel:ISR = on");
  else pythia.readString("PartonLevel:ISR = off");
  if (ismpi==1) pythia.readString("PartonLevel:MPI = on");
  else pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("TimeShower:pTmin = 1.");
  pythia.readString("Random:setSeed= on");
  pythia.readString("111:mayDecay = off");  

  //Sqrts 
  ostringstream sqrtstring;
  double sqrts=atof(argv[3]);
  sqrtstring << "Beams:eCM = " << sqrts;
  pythia.readString(sqrtstring.str().c_str());
  double jet_etacut=2., had_etacut=1.;
  dphicut=2.*my_pi/3.;
  jetacut=2.;
  leadpt=120., subpt=30.;
  double fptbin=5.;
  if (sqrts==200.) fptbin=2.;
  if (sqrts==200.) jet_etacut=0.35, had_etacut=0.35, dphicut=my_pi-0.4, jetacut=0.35, leadpt=20., subpt=10.;
  cout << " Jet Etacut= " << jet_etacut << " Had Etacut= " << had_etacut << endl;
  cout << " DPhicut= " << dphicut << " Asym etacut= " << jetacut << " LeadPtcut= " << leadpt << " SubPtcut= " << subpt << endl;  

  //Seed
  ostringstream seedstring; 
  int see=33+atoi(argv[2]);
  int core=atoi(argv[2]);
  seedstring << "Random:seed = " << see;
  pythia.readString(seedstring.str().c_str());
 
  ostringstream hatmin, hatmax;
  double nhatmin = 0.;
  double nhatmax = 0.;
  hatmin << "PhaseSpace:pTHatMin = " << nhatmin;
  hatmax << "PhaseSpace:pTHatMax = " << nhatmax;
  pythia.readString(hatmin.str().c_str());
  pythia.readString(hatmax.str().c_str());

  pythia.readString("PhaseSpace:bias2Selection = on");
  if (sqrts==200.) pythia.readString("PhaseSpace:bias2SelectionPow = 6.");
  else pythia.readString("PhaseSpace:bias2SelectionPow = 4.");
  pythia.readString("PhaseSpace:bias2SelectionRef = 10.");
  if (sqrts==200.) pythia.readString("PhaseSpace:bias2SelectionRef = 1.");

  pythia.init();

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    double weight = info.weight();

    int qcharge[100000], qide[100000];
    int k=0;
    vector<PseudoJet> particles;
    for (int i=0; i<pythia.event.size(); i++) {
      if (pythia.event[i].isFinal()) {
	double Px=pythia.event[i].px();
	double Py=pythia.event[i].py();
	double Pz=pythia.event[i].pz();
	double E=pythia.event[i].e();
	int charge=pythia.event[i].charge();
	int ide=pythia.event[i].id();
	qcharge[k]=charge;	
	qide[k]=ide;

	particles.push_back( PseudoJet(   Px,  Py,  Pz, E) );
	double pt=sqrt(pow(Px,2.)+pow(Py,2.));
	double raper=1./2.*log((sqrt(pow(pt,2.)+pow(Pz,2.))+Pz)/(sqrt(pow(pt,2.)+pow(Pz,2.))-Pz));
	int npt=-1000;
	for (unsigned a=0; a<31; a++) {
	  if (pt>=bincut[a] && pt<bincut[a+1]) npt=int(a);
	}
	if (npt!=-1000 && fabs(raper)<=had_etacut) {
	  if (sqrts==200. && int(ide)==111) hadspec[npt]+=weight, hadspectwo[npt]+=pow(weight,2.);
	  if (sqrts!=200. && charge!=0.) hadspec[npt]+=weight, hadspectwo[npt]+=pow(weight,2.);
	}

	k+=1;
      }
    }

    for (int rad=2; rad<=5; rad++) {
      double R=double(rad)/10.;
      JetDefinition jet_def(antikt_algorithm, R);
      ClusterSequence cs(particles, jet_def);
      vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
      //Jet Spec
      for (unsigned int nj=0; nj<jets.size(); nj++) {
        int npt=-1000;
        for (unsigned a=0; a<50; a++) {
          if (jets[nj].pt()>=jbincut[a] && jets[nj].pt()<jbincut[a+1]) npt=int(a);
        }
        if (npt!=-1000 && fabs(jets[nj].eta())<=jet_etacut) {
	  jetspec[npt][rad-2]+=weight;
	  jetspectwo[npt][rad-2]+=pow(weight,2.);
	}
      }
      //FF
      for (unsigned int nj=0; nj<jets.size(); nj++) {
	int npt=int(jets[nj].pt()/fptbin);
	if (jets[nj].pt()>ptmin && fabs(jets[nj].eta())<jet_etacut && npt<500) {
	  jetff[npt][rad-2]+=weight, jetfftwo[npt][rad-2]+=pow(weight,2.);
	  for (unsigned int py=0; py<particles.size(); py++) {
	    double deltaR=particles[py].delta_R(jets[nj]);
	    if (deltaR<R && particles[py].pt()>fragcut) {
              double zet=(particles[py].px()*jets[nj].px()+particles[py].py()*jets[nj].py()+particles[py].pz()*jets[nj].pz())/((pow(jets[nj].px(),2.)+pow(jets[nj].py(),2.)+pow(jets[nj].pz(),2.)));
	      int na=int(4.*log(1./zet));
	      if (na>=0 && na<20) {
                if (sqrts==200. && qide[py]==111) ff[na][npt][rad-2]+=weight, fftwo[na][npt][rad-2]+=pow(weight,2.); 
	        if (sqrts!=200. && qcharge[py]!=0) ff[na][npt][rad-2]+=weight, fftwo[na][npt][rad-2]+=pow(weight,2.);
	      }	
	    }
          }
        }
      }
      //Asym
      int L=-1, S=-1;
      select_dijet(jets, L, S);
      if (L!=-1 && S!=-1) {
        if (jets[L].pt()>=leadpt && jets[S].pt()>=subpt) {
          double acodelta=acos((jets[L].px()*jets[S].px()+jets[L].py()*jets[S].py())/jets[L].pt()/jets[S].pt());
          if (acodelta<0.) cout << " acodelta= " << acodelta << endl;
          if (acodelta>dphicut) {
	    int n_aj=int(17.*(jets[L].pt()-jets[S].pt())/(jets[L].pt()+jets[S].pt()));
            int x_aj=int(17.*jets[S].pt()/jets[L].pt());
	    jet_naj_hist[n_aj][rad-2]+=weight;
            jet_xaj_hist[x_aj][rad-2]+=weight;
	    jet_naj_histtwo[n_aj][rad-2]+=pow(weight,2.);
            jet_xaj_histtwo[x_aj][rad-2]+=pow(weight,2.);
          }
        }
      } 

      jets.clear();
    }
    particles.clear();      
 
  } //End event loop

  //Histos
  double sigmaNorm = (info.sigmaGen() / info.weightSum());
  for (int rad=2; rad<=5; rad++) {
    //Jet Spec
    std::ostringstream fg;
    fg << "../results/sqrts_" << sqrts << "_ISR_" << isisr << "_MPI_" << ismpi << "_JETR0" << rad << "spec" << core << ".dat";
    std::ofstream outfile5(fg.str().c_str(),std::ios_base::binary);
    for (unsigned a=0; a<50; a++) {
      double binwidth=ptbin[a+1];
      outfile5 << (jbincut[a]+jbincut[a+1])/2. << " " << jetspec[a][rad-2] << " " << jetspectwo[a][rad-2] << " " << binwidth << " " << sigmaNorm << endl;
    }
    outfile5 << "END" << " " << std::setprecision(9) << info.sigmaGen() << " " << info.weightSum() << endl;
    outfile5.close();

    //Asym
    std::ostringstream fas;
    fas << "../results/sqrts_" << sqrts << "_ISR_" << isisr << "_MPI_" << ismpi << "_ASYMR0" << rad << "spec" << core << ".dat";
    std::ofstream outfileAS(fas.str().c_str(),std::ios_base::binary);
    double v_totasym=0.;
    for (unsigned int i=0; i<17; i++)
    {
      v_totasym+=jet_naj_hist[i][rad-2];
    }
    for (unsigned int i=0; i<17; i++)
    {
      double binsize=1./17.;
      outfileAS << double(i)*binsize+binsize/2. << " " << jet_naj_hist[i][rad-2] << " " << jet_naj_histtwo[i][rad-2] << " " << jet_xaj_hist[i][rad-2] << " " << jet_xaj_histtwo[i][rad-2] << " " << v_totasym << " " << binsize << " " << sigmaNorm << endl;
    }
    outfileAS << "END" << " " << std::setprecision(9) << info.sigmaGen() << " " << info.weightSum() << endl;
    outfileAS.close();

    //FF
    std::ostringstream fff;
    fff << "../results/sqrts_" << sqrts << "_ISR_" << isisr << "_MPI_" << ismpi << "_FFR0" << rad << "spec" << core << ".dat";
    std::ofstream outfile7(fff.str().c_str(),std::ios_base::binary);
    double binsize=0.25;
    for (unsigned int jbin=0; jbin<500; jbin++) {
      if (jetff[jbin][rad-2]!=0.) {
        outfile7 << jbin << " " << jetff[jbin][rad-2] << " " << jetfftwo[jbin][rad-2] << " ";
        for (unsigned e=0; e<20; e++) {
          outfile7 << ff[e][jbin][rad-2] << " " << fftwo[e][jbin][rad-2] << " ";
        }
        outfile7 << " " << binsize << " " << sigmaNorm << endl;
      }
    }
    outfile7 << "END" << " " << std::setprecision(9) << info.sigmaGen() << " " << info.weightSum() << " jet_ptbin= " << fptbin << endl;
    outfile7.close();
  }
  //Hadron Spec
  std::ostringstream fg;
  fg << "../results/sqrts_" << sqrts << "_ISR_" << isisr << "_MPI_" << ismpi << "_HADRONspec" << core << ".dat";
  std::ofstream outfile5(fg.str().c_str(),std::ios_base::binary);
  for (unsigned a=0; a<31; a++) {
    outfile5 << (bincut[a]+bincut[a+1])/2. << " " << hadspec[a] << " " << hadspectwo[a] << " " << (bincut[a+1]-bincut[a]) << " " << sigmaNorm << endl;
  }
  outfile5 << "END" << " " << std::setprecision(9) << info.sigmaGen() << " " << info.weightSum() << endl;
  outfile5.close();

  pythia.stat();

  return 0;
}

//Select Dijets within EtaCut
void select_dijet(vector<PseudoJet> jets, int &L, int &S)
{
        double ileadpt=0., isubpt=0.;
        for (unsigned int i=0; i<jets.size(); i++)
        {
                if (abs(jets[i].eta())<=jetacut)
                {
                        if (jets[i].pt()>ileadpt) isubpt=ileadpt, ileadpt=jets[i].pt(), S=L, L=int(i);
                        else if (jets[i].pt()>isubpt) isubpt=jets[i].pt(), S=int(i);
                }
        }
} 
