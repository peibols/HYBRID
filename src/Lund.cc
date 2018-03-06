#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include <iostream>
#include <vector>
#include <fstream>

#include "Parton.h"
#include "Quench.h"
#include "Hadron.h"

using std::vector;
using namespace std;
using namespace Pythia8;

Pythia pythia;
Event& event      = pythia.event;
ParticleData& pdt = pythia.particleData;

void do_lund(vector<Parton> partons, vector<Quench> quenched, vector<Hadron> &vhadrons, vector<Hadron> &qhadrons);
void init_lund();

void init_lund()
{
	pythia.readString("Random:setSeed = on");
        pythia.readString("Random:seed = 0");

        pythia.readString("Print:quiet = on");
        pythia.readString("111:mayDecay = off");
        pythia.readString("ProcessLevel:all = off");

	pythia.init();
}

void do_lund(vector<Parton> partons, vector<Quench> quenched, vector<Hadron> &vhadrons, vector<Hadron> &qhadrons)
{

	//Hadronize vacuum
	event.reset();
	int colsum=0;
	for (unsigned int i = 0; i < partons.size(); i++) {
        	//If final, introduce in PYTHIA
		if (partons[i].GetD1()==-1) {
			int ide=partons[i].GetId();
			vector<double> p=partons[i].vGetP();
			double px=p[0];
			double py=p[1];
			double pz=p[2];
			double mm=pdt.m0(ide);
			int col=partons[i].GetCol();
			int acol=partons[i].GetAcol();
			double ee=sqrtpos(px*px+py*py+pz*pz+mm*mm);
			event.append(int(ide),23,int(col),int(acol),px,py,pz,ee,mm);
			colsum+=col-acol;
		}
	}
	if (colsum!=0.) cout << " Sumadecolor= " << colsum << endl;
	if (!pythia.next()) {
        	cout << " Event generation aborted prematurely, owing to error" << endl;
        }

	for (int i=0; i<pythia.event.size(); ++i) {
		if (pythia.event[i].isFinal()) {
			//Fill vhadrons vector
			vector<double> p;
                        for (unsigned int j=1; j<4; j++) p.push_back(pythia.event[i].p()[j]);
                        p.push_back(pythia.event[i].p()[0]);
			vhadrons.push_back( Hadron ( Parton ( p, 0., pythia.event[i].m(), 0, -1, -1, pythia.event[i].id(), "lund", 0, 0, true), pythia.event[i].charge(), -1 ) );
		}
	}

	//Hadronize medium
	event.reset();
        colsum=0;
        for (unsigned int i = 0; i < quenched.size(); i++) {
		//If final, introduce in PYTHIA
		if (quenched[i].GetD1()==-1) {
                        int ide=quenched[i].GetId();
                        vector<double> p=quenched[i].vGetP();
                        double px=p[0];
                        double py=p[1];
                        double pz=p[2];
                        double mm=pdt.m0(ide);
                        int col=quenched[i].GetCol();
                        int acol=quenched[i].GetAcol();
                        double ee=sqrtpos(px*px+py*py+pz*pz+mm*mm);
                        event.append(int(ide),23,int(col),int(acol),px,py,pz,ee,mm);
                        colsum+=col-acol;
                }
        }
        if (colsum!=0.) cout << " Sumadecolor= " << colsum << endl;
        if (!pythia.next()) {
                cout << " Event generation aborted prematurely, owing to error" << endl;
        }

	for (int i=0; i<pythia.event.size(); ++i) {
                if (pythia.event[i].isFinal()) {
			//Fill qhadrons vector
			vector<double> p;
                        for (unsigned int j=1; j<4; j++) p.push_back(pythia.event[i].p()[j]);
                        p.push_back(pythia.event[i].p()[0]);
                        qhadrons.push_back( Hadron ( Parton ( p, 0., pythia.event[i].m(), 0, -1, -1, pythia.event[i].id(), "lund", 0, 0, true), pythia.event[i].charge(), -1 ) );
                }
        }

//End program
}
