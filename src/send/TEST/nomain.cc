#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Parton.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <stdlib.h>

using std::vector;
using namespace std;
using namespace Pythia8;

void TreeDoer (vector<Parton>&);

int main(int argc, char** argv)
{
	assert(argc==4);
	int nhat=atoi(argv[1]);
	int N=atoi(argv[2]);
	int seed=atoi(argv[3]);

	//Initialize Pythia
	Pythia pythia;

	//Read cmnd file
	pythia.readFile("setup_pythia.cmnd");

	//Set Random Seed
	ostringstream seedstring;
        seedstring << "Random:seed = " << seed;
        pythia.readString(seedstring.str().c_str());
	
	//Set pTHatMin
        double nhatmin[]={0.0,3.5,7.0,15.0,30.0,50.0,80.0,120.0,170.0,220.0,280.0};
        ostringstream hatmin;
        hatmin << "PhaseSpace:pTHatMin = " << nhatmin[nhat];
        pythia.readString(hatmin.str().c_str());

	vector<Parton> partons;
	cout << "parton size= " << partons.size()<< endl;
	TreeDoer(partons);
	cout << "parton size= " << partons.size()<< endl;
	for (unsigned int i = 0; i < partons.size(); i++)
        {
        	cout << " Parton " << i << " has mom= " << partons[i].GetMom() << " is done or not " << partons[i].GetIsDone() << " orig= " << partons[i].GetOrig();
        	cout << " and has daughters = " << partons[i].GetD1() << " and  " << partons[i].GetD2() << endl;
        }
	
	return 0;
}
