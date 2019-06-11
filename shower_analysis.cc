#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "Parton.h"

#include "vector_operators.h"

using std::vector;
using namespace fastjet;
using namespace std;

double ****numparts_hist;
double totaljets[22][2]={{0.}};

class Times: public Parton
{
  private:
    double _ti;
    double _deltat;
    int _HJ;

  public :
    Times();
    Times(Parton partons) : Parton(partons), _ti(0.), _deltat(0.), _HJ(-1) {};
    ~Times(){};

    void SetTi(double ti) {_ti=ti;}
    double GetTi() const {return _ti;}

    void SetDeltaT(double deltat) {_deltat=deltat;}
    double GetDeltaT() const {return _deltat;}

    void SetHJ(int HJ) {_HJ=HJ;}
    int GetHJ() const {return _HJ;}
    
};

void shower_analysis(vector<Parton> partons);

bool first_event=true;
int totalmine=0;

void shower_analysis(vector<Parton> partons)
{

  // Histos
  double timestep=1.;
  int maxbins=100;
  int maxjetbins=22;
  double jetbincuts[23]={0.,10.,20.,30.,40.,50.,60.,70.,80.,100.,150.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,1200.,1500.,3000.};
  if (first_event==true)
  {
    cout << " RESETTING \n \n";
    first_event=false;
    // Allocate memory for multidim array
    numparts_hist = (double ****)malloc(maxbins * sizeof(double ***)); 
    for (int i=0; i<maxbins; i++) numparts_hist[i]=(double ***)malloc(2 * sizeof(double **));
    for (int i=0; i<maxbins; i++)
    {
      for (int j=0; j<2; j++)
      {
        numparts_hist[i][j]=(double **)malloc(maxjetbins * sizeof(double *));
      }
    }
    for (int i=0; i<maxbins; i++)
    {
      for (int j=0; j<2; j++)
      {
        for (int k=0; k<maxjetbins; k++)
        {
          numparts_hist[i][j][k]=(double *)malloc(maxjetbins * sizeof(double));
        }
      }
    }
    // Set multidim array to zero
    for (int i=0; i<maxbins; i++)
    {
      for (int j=0; j<2; j++)
      {
        for (int k=0; k<maxjetbins; k++)
        {
          for (int l=0; l<2; l++)
          {
            numparts_hist[i][j][k][l]=0.;
          }
        }
      }
    }
  } 

    
  // Minimum virtuality for formation time of final particles
  double minvirt=0.001;

  vector<Times> tparton;
  int HJ_is[2]={-1};
  int HardSpecies[2]={-1};
  for (unsigned int i=0; i<partons.size(); i++)
  {
    tparton.push_back(Times(partons[i]));
    if (partons[i].GetOrig()=="hs" && HJ_is[0]!=-1)
    {
      HJ_is[1]=tparton.size()-1;
      HardSpecies[1]=partons[i].GetId();
    }
    else if (partons[i].GetOrig()=="hs")
    {
      HJ_is[0]=tparton.size()-1;
      HardSpecies[0]=partons[i].GetId();
    }
  }

  //cout << " AKI HJ[0]= " << HJ_is[0] << " HJ[1]= " << HJ_is[1] << "\n \n \n";

  for (unsigned int i=0; i<tparton.size(); i++)
  {
    if (tparton[i].GetD1()==-1)
    {
      if (partons[i].GetOrig()=="rem") continue;
      //Stablish particular chain
      vector<int> chain;
      chain.push_back(i);
      int j=i;
      while (true)
      {
        if (tparton[j].GetOrig()=="hs" || tparton[j].GetDeltaT()!=0.) break;
        else
        {
          int mom=tparton[j].GetMom();
          chain.push_back(mom);
          j=mom;
        }  
      }

      //Add up times
      std::reverse(chain.begin(),chain.end());
      for (unsigned int k=0; k<chain.size(); k++)
      {
        int ind=chain[k];

        double form_time= 0.2 * 2. * tparton[ind].vGetP()[3] / pow(max(tparton[ind].GetQ(),minvirt),2.);
        tparton[ind].SetDeltaT(form_time);

        if (tparton[ind].GetOrig()=="hs")
        {
	  tparton[ind].SetTi(0.);
	  tparton[ind].SetHJ(ind);
	}
	else if (k==0);
	else
        {
          tparton[ind].SetTi(tparton[chain[k-1]].GetTi()+tparton[chain[k-1]].GetDeltaT());
          tparton[ind].SetHJ(tparton[chain[k-1]].GetHJ());
        }
      }
    }

  }

  for (unsigned int k=0; k<2; k++)
  {
    double prov_numparts_hist[100]={0.};
    for (unsigned int i=0; i<tparton.size(); i++)
    {
      //cout << " myHJ now for parton " << i << " is= " << tparton[i].GetHJ() << endl;
      if (tparton[i].GetHJ()!=HJ_is[k]) continue;
      double tini=tparton[i].GetTi();
      double tfin=tini+tparton[i].GetDeltaT();
      int first_bin=int(tini/timestep);
      int last_bin=min(maxbins-1,int(tfin/timestep));
      if (last_bin<0) last_bin=maxbins-1;
      if (last_bin!=maxbins-1 && first_bin!=last_bin) last_bin-=1;
      //cout << " my hard is " << tparton[i].GetHJ() << " " << " my daughter is " << tparton[i].GetD1() << " ";
      //cout << " First bin= " << first_bin << " tini= " << tini << " Last_bin= " << last_bin << " tfin= " << tfin << endl;
      for (int j=first_bin; j<=last_bin; j++)
      {
        prov_numparts_hist[j]+=1.;
      }     
    }
    // Fill Final Histo
    int is_qg=-1;
    if (HardSpecies[k]==21)
    {
      is_qg=0;
    }
    else if (abs(HardSpecies[k])<=6)
    {
      is_qg=1;
    }
    else
    {
      cout << " Initiated not by quark nor gluon! \n \n";
      exit(0);
    }
    double hardpt=tparton[HJ_is[k]].GetPt();
    int njet=-1000;
    for (unsigned int a=0; a<22; a++)
    {
      if (hardpt>=jetbincuts[a] && hardpt<jetbincuts[a+1]) njet=a;  
    }
    //njet =3;
    if (njet==-1000) continue;
    for (int j=0; j<maxbins; j++)
    {
      //cout << " njet= " << njet << endl;
      //cout << prov_numparts_hist[j] << endl;
      numparts_hist[j][0][njet][is_qg]+=prov_numparts_hist[j];
      //cout << numparts_hist[j][0][njet]/(totaljets[njet]+1.) << endl;
      numparts_hist[j][1][njet][is_qg]+=pow(prov_numparts_hist[j],2.);
    }
    totalmine+=1;
    totaljets[njet][is_qg]+=1.;
    //cout << " totalmine= " << totalmine << " totaljets [njet]= " << totaljets[njet] << " njet= " << njet << endl;
  }

  // Print Histo

  // <#Part> vs t, for different jet pt
  for (int jb=0; jb<22; jb++)
  {
    char outPart[100];
    sprintf(outPart,"./results/numpart_vs_time/for_jetbin_%i.dat",jb);
    ofstream part_file;
    part_file.open (outPart);
    
    for (int i=0; i<maxbins; i++)
    {
      part_file << double(i)*timestep+timestep/2. << " ";
      for (int qg=0; qg<2; qg++)
      {
        //cout << " <num>= " << numparts_hist[i][0][jb]/totaljets[jb] << " totaljets= " << totaljets[jb] << endl;
        double num_err=sqrt((numparts_hist[i][1][jb][qg]-totaljets[jb][qg]*pow(numparts_hist[i][0][jb][qg]/totaljets[jb][qg],2.))/totaljets[jb][qg]/(totaljets[jb][qg]-1.));
        part_file << numparts_hist[i][0][jb][qg]/totaljets[jb][qg] << " " << num_err << " " << totaljets[jb][qg] << endl;
      }
    }
  }

  // <#Part> at certain times vs jet pt
  for (int qg=0; qg<2; qg++)
  {
    char outFix[100];
    sprintf(outFix,"./results/numpart_vs_time/fixed_vs_jetpt_%i.dat",qg);
    ofstream fix_file;
    fix_file.open (outFix);
    // Fix times(fm/c): 1, 3, 5, 10, 15, 20, 40;
    int fix_times[7]={1,3,5,10,15,20,40};
    for (unsigned int jj=0; jj<22; jj++)
    {
      fix_file << (jetbincuts[jj]+jetbincuts[jj+1])/2. << " ";
      for (unsigned int tt=0; tt<7; tt++)
      {
        int tm=fix_times[tt];
        double t_err=sqrt((numparts_hist[tm][1][jj][qg]-totaljets[jj][qg]*pow(numparts_hist[tm][0][jj][qg]/totaljets[jj][qg],2.))/totaljets[jj][qg]/(totaljets[jj][qg]-1.));
        if (totaljets[jj][qg]!=0.) fix_file << numparts_hist[tm][0][jj][qg]/totaljets[jj][qg] << " " << t_err << " ";
        else fix_file << 0. << " " << 0. << " ";
      }
      fix_file << endl; 
    }
  }
  

  cout << " NEXT \n";



// End Program
}
























