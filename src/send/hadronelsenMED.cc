#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

using namespace Pythia8;
ofstream geneal;
int main(int argc, char** argv) {
	assert(argc==3);
	char InpF[100], OutF[100];
	sprintf(InpF,"../../%s.out",argv[1]);
	sprintf(OutF,"../../%shad.out",argv[1]);	

	ifstream infile;
	geneal.open (OutF);
        infile.open (InpF);
	double px, py, pz, E, qpx, qpy, qpz, qe, qg, ini, col, acol, ide;

	int nEvent = atoi(argv[2]);

	Pythia pythia;
	Event& event      = pythia.event;
	ParticleData& pdt = pythia.particleData;

	pythia.readString("111:mayDecay = off");
	pythia.readString("ProcessLevel:all = off");

	pythia.init();

	for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
		event.reset();
		double colsum=0.;
		double mm, ee;
		std::string sQx,sQy,sQz,sQe;
		while (infile >> px >> py >> pz >> E >> sQx >> sQy >> sQz >> sQe >> qg >> ini >> col >> acol >> ide) {
			qpx = strtod(sQx.c_str(), NULL);
                        qpy = strtod(sQy.c_str(), NULL);
                        qpz = strtod(sQz.c_str(), NULL);
			qe = strtod(sQe.c_str(), NULL);
                        if (qpx!=qpx || qpy!=qpy || qpz!=qpz || qe!=qe ) {
                                cout << " BUG STILL THEEEEERE\n";
                        }
                        if (qpx!=qpx) qpx=0.;
                        if (qpy!=qpy) qpy=0.;
                        if (qpz!=qpz) qpz=0.;
                        if (qe!=qe) qe=0.;
			if (px==0.000123 && py==0.) {
				break;
			}
			else {
				if (qe==0.) {
					qpx=px/100000.;
					qpy=py/100000.;
					qpz=pz/100000.;
				}
				mm=pdt.m0(int(ide));
				ee=sqrtpos(qpx*qpx+qpy*qpy+qpz*qpz+mm*mm);
				event.append(int(ide),23,int(col),int(acol),qpx,qpy,qpz,ee,mm);
				colsum+=col-acol;
				//cout << " Col= " << col << " Acol= " << acol << "\n";
			}
		}
		if (colsum!=0.) cout << " Sumadecolor= " << colsum << " at count= " << iEvent << " ";
		if (!pythia.next()) {
			cout << " Event generation aborted prematurely, owing to error at count= " << iEvent << "\n";
			//break;
		}

		for (int i=0; i<pythia.event.size(); ++i) {
			if (pythia.event[i].isFinal()) {
				geneal << pythia.event[i].id() << " " << pythia.event[i].px() << " " << pythia.event[i].py() << " " << pythia.event[i].pz() << " " << pythia.event[i].e() << " " << pythia.event[i].charge() << "\n";
			}
		}
		if (pythia.event.size()!=0) geneal << "0.000123" << " 0. " << " 0. " << " 0. " << " 0. "  << iEvent << "\n" ;
	//End Event Loop
	}
	geneal.close();
	return 0;
//End Program
}
