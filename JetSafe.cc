#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <vector>

#include "Parton.h"

using std::vector;
using namespace fastjet;
using namespace std;

bool safe_jet (vector<Parton> partons, double safept);

int plushund=0;

bool safe_jet (vector<Parton> partons, double safept) {
	
	bool safe=false;

	vector<PseudoJet> particles;
	for (unsigned int i=0; i<partons.size(); i++)
	{
		particles.push_back( PseudoJet ( partons[i].vGetP()[0], partons[i].vGetP()[1], partons[i].vGetP()[2], partons[i].vGetP()[3]) );
	}
	
	// choose a jet definition
	double R = 0.3;
	JetDefinition jet_def(antikt_algorithm, R);

  	// run the clustering, extract the jets
  	ClusterSequence cs(particles, jet_def);
  	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  	for (unsigned int i=0; i<jets.size(); i++)
	{
		//cout << " JetPt= " << jets[i].pt() << endl;
		if (abs(jets[i].eta())<2. && jets[i].pt()>=safept) {
			safe=true;
			if (jets[i].pt()>100.) plushund+=1;
			cout << " Plushund= " << plushund << endl;
			break; 
		}
	}
	return safe;
}
