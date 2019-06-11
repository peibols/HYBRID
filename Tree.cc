#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include <iostream>
#include <vector>
#include <fstream>

#include "Parton.h"

#include "vector_operators.h"

using std::vector;
using namespace std;
using namespace Pythia8;

bool safe_jet(vector<Parton>, double);
void init_tree(int);
void do_tree(vector<Parton>&, double&, double&);

Pythia pythia;

void init_tree(int njob)
{
  //Read cmnd file
  ostringstream pythiaset;
  //pythiaset << "setup_pythia_" << pdf << "_" << system << ".cmnd";
  pythiaset << "setup_pythia.cmnd";
  pythia.readFile(pythiaset.str());

  //Set PDF (not used now)
  //std::string pdfSet = "LHAPDF5:cteq6ll.LHpdf";
  //pythia.readString("PDF:pSet = " + pdfSet);
  //pythia.readString("PDF:extrapolate = on");
	
  //Set Random Seed
  int seed=33+njob;
  ostringstream seedstring;
  seedstring << "Random:seed = " << seed;
  pythia.readString(seedstring.str().c_str());

  //Set pTHatMin (not used now)
  //double nhatmin[]={0.0,3.5,7.0,15.0,30.0,50.0,80.0,120.0,170.0,220.0,280.0};
  //ostringstream hatmin;
	
  //JetSafePt Values (not used now)
  //double safept[]={0.,10.,18.,34.,60.,100.,153.,227.,305.,400.,480.};

  pythia.init();
  cout << " TREE INITIALISED \n";
}

void do_tree(vector<Parton> &partons, double &weight, double &cross)
{
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;
  Info& info = pythia.info;
  
  if (!pythia.next()) {
    return;
  }

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
  //bool is_safe=safe_jet(partons, safept[nhat]);	
  //if (is_safe==false) continue;
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
  //cout << " Partons size= " << partons.size() << endl;

  //Extract weight and cross section of the event
  weight=info.weight();
  cross=info.sigmaGen();
	
}
