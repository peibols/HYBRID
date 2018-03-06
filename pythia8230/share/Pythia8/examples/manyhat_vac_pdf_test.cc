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

//Args: N, seed, sqrts 
int main(int argc, char** argv) {

  assert(argc==4);

  std::string sqrts=argv[3];
  std::string tpdf="VAC";

  int nEvent  = atoi(argv[1]);
  int core = atoi(argv[2]); //Not used

  double ptmin=0.;
  double ptbin=2.;
  double jetspec[500][4], p_jetspec[500][4];
  for (unsigned a=0; a<500; a++) {
    for (unsigned c=0; c<4; c++) {
      jetspec[a][c]=0.;
      p_jetspec[a][c]=0.;
    }
  }
  double hadspec[31], p_hadspec[31];
  for (unsigned a=0; a<31; a++) {
    hadspec[a]=0.;
    p_hadspec[a]=0.;
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

  string pdfSet = "LHAPDF5:cteq6ll.LHpdf";

  Pythia pythia;
  Info& info = pythia.info;
  Event& event = pythia.event;

  pythia.readString("PDF:pSet = " + pdfSet);
  pythia.readString("PDF:extrapolate = on");
  
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("SoftQCD:all = off");
  pythia.readString("HardQCD:all = on");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("TimeShower:pTmin = 1.");
  pythia.readString("Random:setSeed= on");
  pythia.readString("111:mayDecay = off");  
  
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 6.");
  pythia.readString("PhaseSpace:bias2SelectionRef = 1.");

  int nBin=6;
  for (int iBin = 0; iBin < nBin; ++iBin) {

    ostringstream hatmin, hatmax;
    double nhatmin[7]={1.,8.,16.,24.,32.,40.,0.};
    hatmin << "PhaseSpace:pTHatMin = " << nhatmin[iBin];
    hatmax << "PhaseSpace:pTHatMax = " << nhatmin[iBin+1];
    pythia.readString(hatmin.str().c_str());
    pythia.readString(hatmax.str().c_str());

    pythia.init();

    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      double weight = info.weight();

      vector<PseudoJet> particles;
      for (int i=0; i<pythia.event.size(); i++) {
        if (pythia.event[i].isFinal()) {
	  double Px=pythia.event[i].px();
	  double Py=pythia.event[i].py();
	  double Pz=pythia.event[i].pz();
	  double E=pythia.event[i].e();
	  int charge=pythia.event[i].charge();
	  int ide=pythia.event[i].id();

	  particles.push_back( PseudoJet(   Px,  Py,  Pz, E) );
	  double pt=sqrt(pow(Px,2.)+pow(Py,2.));
	  double raper=1./2.*log((sqrt(pow(pt,2.)+pow(Pz,2.))+Pz)/(sqrt(pow(pt,2.)+pow(Pz,2.))-Pz));
	  int npt=-1000;
	  for (unsigned a=0; a<31; a++) {
	    if (pt>=bincut[a] && pt<bincut[a+1]) npt=int(a);
	  }
	  if (npt!=-1000 && fabs(raper)<=0.35) {
	    if (sqrts=="020" && int(ide)==111) p_hadspec[npt]+=weight;
	    if (sqrts!="020" && charge!=0.) p_hadspec[npt]+=weight;
	  }
        }
      }

      for (int rad=2; rad<=5; rad++) {
        double R=double(rad)/10.;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
        for (unsigned int nj=0; nj<jets.size(); nj++) {
	  int npt=int(jets[nj].pt()/ptbin);
	  if (npt<500 && fabs(jets[nj].eta())<=0.35) {
	    p_jetspec[npt][rad-2]+=weight;
	  }
        }
        jets.clear();
      }
      particles.clear();      

      //Fill total histograms
      double sigmaNorm = (info.sigmaGen() / info.weightSum());
      for (unsigned a=0; a<500; a++) {
        for (unsigned b=0; b<4; b++) {
          jetspec[a][b]+=sigmaNorm*p_jetspec[a][b];
	  p_jetspec[a][b]=0.;
        }
      }
      for (unsigned a=0; a<31; a++) {
        hadspec[a]+=sigmaNorm*p_hadspec[a];
	p_hadspec[a]=0.;
      }
 
    } //End iEvent loop

  } //End iBin loop

  //Histos
  double sigmaNorm = (info.sigmaGen() / info.weightSum());
  for (int rad=2; rad<=5; rad++) {
    std::ostringstream fg;
    fg << "../../results/ManyHat_" << tpdf << "_N" << nEvent << "_JETR0" << rad << "spec" << core << ".dat";
    std::ofstream outfile5(fg.str().c_str(),std::ios_base::binary);
    for (unsigned a=0; a<500; a++) {
      outfile5 << double(a)*ptbin+ptbin/2. << " " << jetspec[a][rad-2]*sigmaNorm << " " << ptbin << endl;
    }
    outfile5.close();
  }
  std::ostringstream fg;
  fg << "../../results/ManyHat_" << tpdf << "_N" << nEvent << "_HADRONspec" << core << ".dat";
  std::ofstream outfile5(fg.str().c_str(),std::ios_base::binary);
  for (unsigned a=0; a<31; a++) {
    outfile5 << (bincut[a]+bincut[a+1])/2. << " " << hadspec[a]*sigmaNorm << " " << (bincut[a+1]-bincut[a]) << endl;
  }
  outfile5.close(); 

  pythia.stat();

  return 0;
} 
