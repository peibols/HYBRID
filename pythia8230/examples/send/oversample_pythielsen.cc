#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <assert.h>

using namespace std;
using namespace Pythia8;

//Args: N, seed, sqrts (in GeV)
int main(int argc, char** argv) {

assert(argc==4);

char outFile[100];
sprintf(outFile,"../../tree_%s.txt",argv[3]);
ofstream geneal;
geneal.open (outFile);

Pythia pythia;
Info& info = pythia.info;
Event& event = pythia.event;

pythia.readString("Beams:idB = 2212");
pythia.readString("SoftQCD:all = off");
pythia.readString("HardQCD:all = on");
pythia.readString("HadronLevel:all = off");
pythia.readString("PartonLevel:ISR = on");
pythia.readString("PartonLevel:MPI = off");
pythia.readString("TimeShower:pTmin = 1.");
pythia.readString("Random:setSeed= on");

//#Events
int N=atoi(argv[1]);

//Sqrts
ostringstream sqrtstring;
double sqrts=atof(argv[3]);
sqrtstring << "Beams:eCM = " << sqrts;
pythia.readString(sqrtstring.str().c_str());

//Seed
ostringstream seedstring;
double see=33+atoi(argv[2]);
seedstring << "Random:seed = " << see; 
pythia.readString(seedstring.str().c_str());

//PDF
string pdfSet = "LHAPDF5:cteq6ll.LHpdf";
pythia.readString("PDF:pSet = " + pdfSet);
pythia.readString("PDF:extrapolate = on");

ostringstream hatmin, hatmax;
double nhatmin = 1.;
double nhatmax = 0.;
hatmin << "PhaseSpace:pTHatMin = " << nhatmin;
hatmax << "PhaseSpace:pTHatMax = " << nhatmax;
pythia.readString(hatmin.str().c_str());
pythia.readString(hatmax.str().c_str());

int finals, counter, m1, m2, use, qg, mam, ini;
double ultra[1000][7];

pythia.readString("PhaseSpace:bias2Selection = on");
pythia.readString("PhaseSpace:bias2SelectionPow = 6.");
pythia.readString("PhaseSpace:bias2SelectionRef = 1.");

//Initialize
pythia.init();

int count=0;
do {
	if (pythia.event.size()>950) {
		cout << " size= " << pythia.event.size();
	}
	finals = 0;
	counter = 0;
	for (unsigned e=0;e<1000; e++) {
		for (unsigned f=0;f<7; f++) {
			ultra[e][f]=0.;
		}
	}

	if (!pythia.next()) continue;
	for (int i = 0; i < pythia.event.size(); ++i) {
		if (pythia.event[i].isFinal()) {
			if (pythia.event[i].mother2() == 0) continue;
			if (pythia.event.size()>1000) continue;
			finals+=1;
			use = i;
			for (unsigned k = 0; k<30; k++){
				m1 = pythia.event[use].mother1();
				m2 = pythia.event[use].mother2();
				if (m1!=m2) break;
				use = m1;
			}
			ultra[i][0]=pythia.event[i].px();
			ultra[i][1]=pythia.event[i].py();
			ultra[i][2]=pythia.event[i].pz();
			ultra[i][3]=pythia.event[i].e();
			ultra[i][4]=m1;
		}
	}
	donantli:
	for (unsigned l = 0; l<1000; l++) {
		for (unsigned m=0; m<1000; m++) {
			if (ultra[m][4]==ultra[l][4] && m!=l && ultra[m][5]==0 && ultra[m][4]!=0) {
				counter +=1;
				use = int(ultra[m][4]);
				ultra[use][0]=ultra[m][0]+ultra[l][0];
				ultra[use][1]=ultra[m][1]+ultra[l][1];
				ultra[use][2]=ultra[m][2]+ultra[l][2];
				ultra[use][3]=ultra[m][3]+ultra[l][3];
				int early = use;
				for ( unsigned y = 0; y<30; y++) {
					m1 = pythia.event[use].mother1();
					m2 = pythia.event[use].mother2();
					if (m1!=m2) break;
					use = m1;
				}
				if (pythia.event[m1].status()!=-41) {
					ultra[early][4]=m1;
				}
				else ultra[early][4]=3;
				ultra[m][5]=1;
				ultra[l][5]=1;
				if (counter < finals) goto donantli;
			}
		}
	}
	ultra[3][0]=0.;
	ultra[3][1]=0.;
	ultra[3][2]=0.;
	ultra[3][3]=0.;
	int col, acol, ide;
	for (unsigned c=3; c<1000; c++) {
		if (pythia.event[ultra[c][4]].status()==-41) ultra[c][4]=3;
		if ((ultra[c][4]!=0 || c==3) && pythia.event[c].status()!=-41) {
			ultra[c][6]=abs(sqrt(abs(pow(ultra[c][3],2.)-pow(ultra[c][0],2.)-pow(ultra[c][1],2.)-pow(ultra[c][2],2.)-pow(pythia.event[c].m(),2.))));
			ide=abs(pythia.event[c].id());
			if (ide >=1 && ide <=6) qg=1;
			if (ide==21) qg=2;
			if (ide==22) {
				qg=3;
			}
			ini=0;
			use=c;
			for (unsigned w=0; w<50; w++) {
				mam=pythia.event[use].mother1();
				if (pythia.event[mam].status()==-23) {
					ini=mam;
					ide=abs(pythia.event[mam].id());
					if (ide >=1 && ide <=6) ini*=1.;
					if (ide==21) ini*=-1.;
					break;
				}
				use=mam;
			}
			ide=pythia.event[c].id();
			col=pythia.event[c].col();
			acol=pythia.event[c].acol();
			geneal << c << " " << ultra[c][0] << " " << ultra[c][1] << " " << ultra[c][2] << " " << ultra[c][3] << " " << ultra[c][6] << " " << int(ultra[c][4]) << " " << qg << " " << ini << " " << col << " " << acol << " " << ide << "\n";
		}
	}
	for (int c=0; c<pythia.event.size(); c++) {
		if (pythia.event[c].status()==63) {
			ide=pythia.event[c].id();
			col=pythia.event[c].col();
			acol=pythia.event[c].acol();
			double px=pythia.event[c].px();
			double py=pythia.event[c].py();
			double pz=pythia.event[c].pz();
			double pe=pythia.event[c].e();
			double pm=pythia.event[c].m();
			double pq=sqrtpos(pe*pe-px*px-py*py-pz*pz-pm*pm);
			geneal << c << " " << px << " " << py << " " << pz << " " << pe << " " << pq << " " << 0 << " " << 0 << " " << 0 << " "  << col << " "  << acol << " "  << ide << "\n";
		}
	}
	double weight = info.weight();
	geneal << "NEXT" << " 0. " << " 0. " << " 0. " << "0.00" << count << " 0. " << " " << weight << " 0. " << " 0. " << " 0. " << " 0. " << " 0. " << "\n" ;
	count+=1;
} while (count<N);

geneal << "END" << " " << std::setprecision(9) << info.sigmaGen() << " " << info.weightSum() << endl;
geneal.close();
pythia.stat();

return 0;

} //End Program





