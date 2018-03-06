#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Parton.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>

using std::vector;
using namespace std;
using namespace Pythia8;

void TreeDoer(vector<Parton>& partons)
{
	//Pythia::pythia.init();

	//Find Final Particles
	for (int i = 0; i < pythia.event.size(); i++)
	{
        	if (pythia.event[i].isFinal())
                {        
			//Simply store remnants
			if (pythia.event[i].status() == 63)
			{
				partons.push_back ( Parton ( pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e(), 0., pythia.event[i].m(), 0, -1, -1, pythia.event[i].id(), "rem", pythia.event[i].col(), pythia.event[i].acol(), true ) );
				continue;
			}
			
			//Find first non-trivial mother
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
			double virt=abs(sqrt(pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2()));
			
			//Add it to partons array
			partons.push_back ( Parton ( pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e(), virt, pythia.event[i].m(), m1, -1, -1, pythia.event[i].id(), "", pythia.event[i].col(), pythia.event[i].acol(), false ) );	

			//If mother is ISR initiator, say mother is 3
			if (pythia.event[m1].status()==-41)
                        {
                  	    	partons[partons.size()-1].SetMom(3);
				partons[partons.size()-1].SetOrig("isr");
                                partons[partons.size()-1].SetIsDone(true);
                        }

			//If mother is hard scattering initiator, say mother is 3
			if (pythia.event[m1].status()==-21)
                        {
				//cout << " ahora " << endl;
                                partons[partons.size()-1].SetMom(3);
				partons[partons.size()-1].SetOrig("hs");
                                partons[partons.size()-1].SetIsDone(true);
                        }

		}
	}
		
	/*
	cout << "#Partons in array = " << partons.size() << endl;
	for (unsigned int i = 0; i < partons.size(); i++)
	{
		cout << " Parton " << i << " has mom= " << partons[i].GetMom() << " is done or not= " << partons[i].GetIsDone() << endl;
	}
	*/

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

						//Set mother momentum
						double px=partons[i].GetPx()+partons[j].GetPx();
						double py=partons[i].GetPy()+partons[j].GetPy();
						double pz=partons[i].GetPz()+partons[j].GetPz();
						double en=partons[i].GetEn()+partons[j].GetEn();
						double virt=abs(sqrt(pow(en,2.)-pow(px,2.)-pow(py,2.)-pow(pz,2.)-pow(pythia.event[mom].m(),2.)));

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
						partons.push_back ( Parton ( px, py, pz, en, virt, pythia.event[mom].m(), m1, i, j, pythia.event[mom].id(), "", pythia.event[mom].col(), pythia.event[mom].acol(), false ) );
						//Update mother of daughters to point to position in partons array instead of pythia list, and declare as done
						partons[i].SetMom(partons.size()-1);
						partons[j].SetMom(partons.size()-1);
						partons[i].SetIsDone(true);
						partons[j].SetIsDone(true);

						//If mother is ISR initiator, say mother of mother is 3
						if (pythia.event[m1].status()==-41)
						{
							partons[partons.size()-1].SetMom(3);
							partons[partons.size()-1].SetOrig("isr");
							partons[partons.size()-1].SetIsDone(true);
						}

						//If mother is hard scattering initiator, say mother of mother is 3
						if (pythia.event[m1].status()==-21)
						{
							partons[partons.size()-1].SetMom(3);
							partons[partons.size()-1].SetOrig("hs");
							partons[partons.size()-1].SetIsDone(true);
						}
							//cout << " Hola " << partons.size() << " daugh1= " << i << " daugh2= " << j << " mom= " << mom << " array mom= " << partons[i].GetMom() << endl;
						changes=0;
						break;
					}
				}
			}
		}
	} while (changes==0);
	/*
	cout << " I aki partons size = " << partons.size() << endl;
	for (unsigned int i = 0; i < partons.size(); i++)
	{
		cout << " Parton " << i << " has mom= " << partons[i].GetMom() << " is done or not " << partons[i].GetIsDone() << " orig= " << partons[i].GetOrig();
		cout << " and has daughters = " << partons[i].GetD1() << " and  " << partons[i].GetD2() << endl;
	}
	*/
	/*
	//Write to Output File
	for (unsigned int i = 0; i < partons.size(); i++)
	{
		geneal << partons[i].GetPx() << " " << partons[i].GetPy() << " " << partons[i].GetPz() << " " << partons[i].GetEn() << " " << partons[i].GetQ() << " ";
		geneal << partons[i].GetMom() << " " << partons[i].GetPx() << " "   	
	}
	*/
}
