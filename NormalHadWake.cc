#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Random.h"
#include "Quench.h"
#include "Wake.h"

#include "vector_operators.h"

using std::vector;
using namespace std;

double transcut=0.;		//Lower threshold for transverse mass for back-reaction
double basesig=0.65;		//Starting gaussian width for pass function in Metropolis
int Nrun=800000;	        //Maximum iterations of Metropolis
double maxptsq=pow(3.5,2.);	//Maximum squared pt for MC
double maxrap=2.5;		//Maximum absolute rapidity for MC
double tole=0.25;		//Non-conservation tolerance in GeV
double masspi=0.1396;		//Pion mass
double masspro=0.938;		//Proton mass
double masstra[2]={masspi,masspro};

double normcoop[2]={30.,30.};		//Norm for OneBody MC, adaptable 

double maxcooper[2]={0.,0.};		//Record highest value of OneBody MC

int toomuch=0;			//How many times we gave up in Metropolis

int totalsave=0;

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);
double rapid(double pt, double pz);
int set_charge(int spe, numrand &nr);
double thermal(int spe, double ptrand);
void one_body(vector<Wake> &wake, vector<double> delta, vector<double> momback, double ptlost, double mtlost, double raplost, numrand &nr, int spe, int mode);
vector<double> vec_abs(vector<double> p);

double passdist[4][20][1];
double phirap[201][201][2];
double phisto[201][2];
double rbin=2.*maxrap/200.;
double risto[201][2];
double ptisto[200][2];
double nmed=0.;
double ntwomed=0.;
double nisto[100];
double Yield[2];

#define DO_PRINT

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr)
{
    if (totalsave==0)
    {
	for (unsigned a=0; a<4; a++) {
		for (unsigned b=0; b<20; b++) {
			for (unsigned c=0; c<1; c++) {
				passdist[a][b][c]=0.;
			}
		}
	}

	//Back Histos Ini
	for (unsigned a=0; a<201; a++) {
		for (unsigned b=0; b<201; b++) {
			for (unsigned c=0; c<2; c++) {
				phirap[a][b][c]=0.;
			}
		}
	}
	for (unsigned a=0; a<201; a++) {
		for (unsigned b=0; b<2; b++) {
			phisto[a][b]=0.;
		}
	}
	for (unsigned a=0; a<201; a++) {
		for (unsigned b=0; b<2; b++) {
			risto[a][b]=0.;
		}
	}
	for (unsigned a=0; a<200; a++) {
		for (unsigned b=0; b<2; b++) {
			ptisto[a][b]=0.;
		}
	}
	for (unsigned a=0; a<100; a++) {
		nisto[a]=0.;
	}
	Yield[0]=0.;
	Yield[1]=0.;
    }

        double emismatch=0.;
	for (unsigned int i=0; i<quenched.size(); i++)
	{
	        if (partons[i].GetD1()!=-1) continue;
		vector<double> delta = partons[i].vGetP()-quenched[i].vGetP();
		if (delta[3]<=0.) continue;				//Swiftly skip guys who didn't lose energy
		
		double ptlost=sqrt(pow(delta[0],2.)+pow(delta[1],2.));	//Lost Pt
		double raplost=rapid(ptlost,delta[2]);	//Want this or pseudorapidity?
		if (raplost!=raplost) {cout << " RapLost NaN: Dx= " << delta[0] << " Dy= " << delta[1] << " Dz= " << delta[2] << " De= " << delta[3] << endl; continue;}
		double mtlost=delta[3]/cosh(raplost);			//Lost transverse mass
		mtlost=delta[3];					//Boost to wake frame
		
		//Reasonable conditions to do back-reaction
		//if (delta[3]>=0. && delta[3]<3000. && mtlost>transcut && fabs(raplost)<3.5)
		if (delta[3]>=0. && delta[3]<10000000. && mtlost>transcut)
		{
			//Declare wake vector for this particle
			vector<Wake> pwake;
			//Declare variables here because of the goto:
			vector<double> momback (4,0.);
			vector<double> pmomback (4,0.);
			vector<double> dif (4,0.);
			vector<double> old_dif(4,0.);
			double msigma, pass, newpass, passrand;
			int spe, mode;
			int runi, encallao;
			int numenc=0;
			double clocklim=0.1;
                        int tooclock=0;
		
			//GoTo flag
			thisis:
			clock_t startClock = clock();
		
			pwake.clear();
			for (unsigned int j=0; j<4; j++) momback[j]=0.;
			runi=0, encallao=0;
		
			//Initial List of Hadrons, approximately satisfy energy conservation
			do
			{
				//Proton or Pion
        			if (nr.rando()<=0.05) spe=1;    //is proton
        			else spe=0;                     //is pion
				//OneBody Dist: last arg=-1 because it is creation of initial list
				one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, -1);
				momback+=pwake[pwake.size()-1].vGetP()*pwake[pwake.size()-1].GetStatus();
			} while(momback[3]<delta[3]);

			//Absolute difference wrt Wake Momentum
			dif=vec_abs(delta-momback);
			
			//Select random particle, flip it and see whether it improves conservation
			msigma=sqrt(pow(dif[0],2.)+pow(dif[1],2.)+pow(dif[2],2.)+pow(dif[3],2.))/sqrt(log(2.));
			
			do
			{
				//Select random particle
				mode=int(double(pwake.size())*nr.rando());
			
				//Get the species it was
				spe=pwake[mode].GetId();
			
				//Generate new particle, same species and status than selected
				one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, mode);
			
				//Define pass function	
				pass=exp((-pow(dif[0],2.)-pow(dif[1],2.)-pow(dif[2],2.)-pow(dif[3],2.))/pow(msigma,2.));
				
				//Update sum of 4momentum: remove previous particle, add new particle
				pmomback=momback+(pwake[pwake.size()-1].vGetP()-pwake[mode].vGetP())*pwake[mode].GetStatus();
				
				//Update difference wrt Lost Momentum
				dif=vec_abs(delta-pmomback);
				
				//Compute new pass function
				newpass=exp((-pow(dif[0],2.)-pow(dif[1],2.)-pow(dif[2],2.)-pow(dif[3],2.))/pow(msigma,2.));
		
				//If conservation improved, accept substitution
				if (newpass>pass)
				{
					pwake.erase(pwake.begin()+mode);
					runi+=1;
					encallao=0;
					momback=pmomback;
					//cout << " changed! " << endl;
				}
				else
				{
					/*
					//If conservation not improved, accept with probability newpass/pass
					passrand=nr.rando();
					if (newpass/pass > passrand)
					//if (pass/newpass > passrand)
					{
						pwake.erase(pwake.begin()+mode);
						runi+=1;
						encallao=0;
						momback=pmomback;
					}
					else
					*/
					{
						pwake.erase(pwake.begin()+pwake.size()-1);
						dif=vec_abs(delta-momback);
						encallao+=1;
					}
				}
			
				//If truly stuck, start over. If 5 times truly stuck, give up
				if (encallao>50000)
				{
					cout << " ENCALLAO !!!! \n \n";
					numenc+=1;
					if (numenc>5) {toomuch+=1; break;}
					goto thisis;
				}
		
				//If time spent surpasses clocklim, start over, and wait longer for next attempt.
				clock_t endClock = clock();
				if (double((endClock - startClock)) / CLOCKS_PER_SEC>clocklim)
				{
					//cout << " CLOCK !!! \n";
					clocklim+=0.02;
					tooclock+=1;
		  			
					double difx=-momback[0]+delta[0];
		  			double dify=-momback[1]+delta[1];
		  			double difz=-momback[2]+delta[2];
		  			double dife=-momback[3]+delta[3];
        				double remass2=dife*dife-difx*difx-dify*dify-difz*difz;
					
					//cout << " Got to= " << " DifX= " << difx << " DifY= " << dify << " DifZ= " << difz << " DifE= " << dife << endl;
					if (fabs(difx)<3.*tole && fabs(dify)<3.*tole && fabs(difz)<8.*tole && remass2>0.) {
					  //Add extra particle with remnant momentum
					  //cout << " OUTSIDE TOLERANCE \n";
					  //cout << " remnant mass= " << remass2 << " " << sqrt(remass2) << endl;

			                  if (fabs(remass2-masstra[0]*masstra[0])<fabs(remass2-masstra[1]*masstra[1])) spe=0;
					  else spe=1;
			  		  //cout << " spe= " << spe << endl;
					  
					  //if (nr.rando()<=0.05) spe=1;    //is proton
        			          //else spe=0;                     //is pion
					 
					  int charge=set_charge(spe, nr);
 
					  double stat;
					  if (dife<0.) stat=-1.;
					  else stat=1.;
					
					  vector<double> p;
					  double tdife=sqrt(difx*difx+dify*dify+difz*difz+masstra[spe]*masstra[spe]);
					  p.push_back(stat*difx), p.push_back(stat*dify), p.push_back(stat*difz), p.push_back(tdife);
					  pwake.push_back( Wake ( p, masstra[spe], charge, spe, stat ) );
					  momback+=pwake[pwake.size()-1].vGetP()*stat;
					  dif=vec_abs(delta-momback);
					  break; 
					}
				
					if (tooclock>100) {toomuch+=1; break;}
                                        goto thisis;
				}

			} while(runi<Nrun && (dif[0]>tole || dif[1]>tole || dif[2]>tole || dif[3]>tole));
			
			if (runi>=Nrun) toomuch+=1;
			
			//If outside tolerance
			if (dif[0]!=0.) {

		  	  double difx=-momback[0]+delta[0];
		  	  double dify=-momback[1]+delta[1];
		  	  double difz=-momback[2]+delta[2];
		  	  double dife=-momback[3]+delta[3];
        		  double remass2=dife*dife-difx*difx-dify*dify-difz*difz;
        		
			  //cout << " WITHIN TOLERANCE \n";
			  //cout << " remnant mass= " << remass2 << " " << sqrt(remass2) << endl;
			  if (fabs(remass2-masstra[0]*masstra[0])<fabs(remass2-masstra[1]*masstra[1])) spe=0;
			  else spe=1;
			  //cout << " spe= " << spe << endl;
			  
			  //if (nr.rando()<=0.05) spe=1;    //is proton
        	          //else spe=0;                     //is pion
					 
			  int charge=set_charge(spe, nr);
 
			  double stat;
			  if (dife<0.) stat=-1.;
			  else stat=1.;
			  
			  vector<double> p;
			  double tdife=sqrt(difx*difx+dify*dify+difz*difz+masstra[spe]*masstra[spe]);
			  p.push_back(stat*difx), p.push_back(stat*dify), p.push_back(stat*difz), p.push_back(tdife);
			  pwake.push_back( Wake ( p, masstra[spe], charge, spe, stat ) );
			  momback+=pwake[pwake.size()-1].vGetP()*stat;
			  dif=vec_abs(delta-momback);

			}
		#ifdef DO_PRINT
			{
		//	  cout << " runi= " << runi << endl;
		//	  cout << " DeltaX= " << delta[0] << " DeltaY= " << delta[1] << " DeltaZ= " << delta[2] << " DeltaE= " << delta[3] << endl;
		//	  cout << " DifX= " << dif[0] << " DifY= " << dif[1] << " DifZ= " << dif[2] << " DifE= " << dif[3] << endl;
		//	  cout << " rel e= " << delta[3]-momback[3] << endl;
		//	  emismatch+=delta[3]-momback[3];
		//	  cout << " pwake.size()= " << pwake.size() << endl << endl;
			}
		#endif
			//Fill total wake with pwake
			int Npart=0;
			for (unsigned int k=0; k<pwake.size(); k++)
			{
                        	double ptrand=sqrt(pwake[k].vGetP()[0]*pwake[k].vGetP()[0]+pwake[k].vGetP()[1]*pwake[k].vGetP()[1]);
                                double mtrand=sqrt(pow(ptrand,2.)+pow(pwake[k].GetMass(),2.));
				double phirand=acos((pwake[k].vGetP()[0]*delta[0]+pwake[k].vGetP()[1]*delta[1])/ptrand/ptlost);
                                //raprand=asinh(bomz[c]/mtrand);
				double raprand=asinh(pwake[k].vGetP()[2]/ptrand);
                                int nphi=int(phirand/(2.*3.141592)*200.);
                                int nrap=0;
                                raprand-=raplost;
                                if (raprand>0.) nrap=100+int((raprand+rbin/2.)/rbin);
                                if (raprand<0.) nrap=100+int((raprand-rbin/2.)/rbin);
                                int npt=int(200.*(ptrand/sqrt(maxptsq)));
				spe=0;
				if (pwake[k].GetMass()>0.5) spe=1;
                                if (npt<200) ptisto[npt][spe]+=1.*pwake[k].GetStatus();
                                if (nrap<201 && nrap>-1) risto[nrap][spe]+=1.*pwake[k].GetStatus();
                                if (nphi<201) phisto[nphi][spe]+=1.*pwake[k].GetStatus();
                                if (nrap<201 && nphi<201 && nrap>-1) {
					phirap[nphi][nrap][spe]+=1.*pwake[k].GetStatus();
				}
				Yield[spe]+=int(1.*pwake[k].GetStatus());
				Npart+=int(1.*pwake[k].GetStatus());
				wake.push_back ( pwake[k] );
			}
			nmed+=Npart;
			pwake.clear();
		}
	}
	cout << " emismatch= " << emismatch << endl;
	cout << " wake.size= " << wake.size() << " toomuch= " << toomuch << endl << endl;
	totalsave+=1;
#ifdef DO_PRINT
	//cout << " maxcooper pion= " << maxcooper[0] << " maxcooper prot= " << maxcooper[1] << endl;
	//cout << " normcoop pion= " << normcoop[0] << " normcoop prot= " << normcoop[1] << endl;
	/*
	for (unsigned int l=0; l<wake.size(); l++)
	{
		cout << " Wake " << l << endl;
		wake[l].display();
		cout << endl;
	}
	*/
	//BackFiles
	ofstream phifile;
        phifile.open ("dists/Phidebug.txt");
        ofstream rapfile;
        rapfile.open ("dists/Rapdebug.txt");
        ofstream ptfile;
        ptfile.open ("dists/Ptdebug.txt");
        ofstream nfile;
        nfile.open ("dists/Ndebug.txt");
	ofstream raphi;
	raphi.open ("dists/Raphi.txt");
	ofstream cons;
        cons.open ("dists/Cons.txt");
	for (unsigned a=0; a<200; a++) {
		for (unsigned b=0; b<200; b++) {
			raphi << 2.*3.141592/400.+double(a)/200.*2.*3.141592 << " " << -maxrap+double(b)*rbin << " " << phirap[a][b][0] << " " << phirap[a][b][1] <<  "\n";
		}
	}
	for (unsigned a=0; a<200; a++) {
                phifile << 2.*3.141592/400.+double(a)/200.*2.*3.141592 << " " << phisto[a][0]/double(totalsave) << " " << sqrt(phisto[a][0])/double(totalsave) << " " << phisto[a][1]/double(totalsave) << " " << sqrt(phisto[a][1])/double(totalsave) << "\n";
        }
        for (unsigned a=0; a<200; a++) {
                rapfile << -maxrap+double(a)*rbin << " " << risto[a][0]/double(totalsave) << " " << sqrt(risto[a][0])/double(totalsave) << " " << risto[a][1]/double(totalsave) << sqrt(risto[a][1])/double(totalsave) << "\n";
        }
        for (unsigned a=0; a<200; a++) {
                ptfile << double(a)/200.*sqrt(maxptsq)+sqrt(maxptsq)/400. << " " << ptisto[a][0]/double(totalsave) << " " << sqrt(ptisto[a][0])/double(totalsave) << " " << ptisto[a][1]/double(totalsave) << " " << sqrt(ptisto[a][1])/double(totalsave) << "\n";
        }
        for (unsigned a=0; a<100; a++) {
                nfile << -49+double(a) << " " << nisto[a]/double(totalsave) << " " << sqrt(nisto[a])/totalsave << "\n";
        }
        double nerr=sqrt(ntwomed+double(totalsave)*pow(nmed/double(totalsave),2.)-2.*pow(nmed,2.)/double(totalsave))/sqrt(double(totalsave)*(double(totalsave)-1.));
        cout << " <N>= " << nmed/double(totalsave) << " " << nerr << "\n";
	cout << " <Pion Yield>= " << Yield[0]/double(totalsave) << " <Proton Yield>= " << Yield[1]/double(totalsave) << " <Proton/Pion Yield> " << Yield[1]/Yield[0] << "\n";

	raphi.close();
        phifile.close();
        rapfile.close();
        ptfile.close();
        nfile.close();

#endif
}

void one_body(vector<Wake> &wake, vector<double> delta, vector<double> momback, double ptlost, double mtlost, double raplost, numrand &nr, int spe, int mode)
{
	double mc=0.;
	double cooper=0.;
	double randian=1.;
	double pxrand, pyrand, raprand, ptrand;
	double phirand, mtrand, phidif, rapdif;
	double mphidif;
	do {
	        {	 
		  ptrand=max(sqrt(maxptsq*nr.rando()),0.000001);	
		  phirand=2.*3.141592*nr.rando();
		  raprand=maxrap*(-1.+2.*nr.rando());
		  pxrand=ptrand*cos(phirand);
              	  pyrand=ptrand*sin(phirand);
                  mtrand=sqrt(pow(ptrand,2.)+pow(masstra[spe],2.));
		  phidif=acos((delta[0]*pxrand+delta[1]*pyrand)/(ptlost*ptrand));
		  //rapdif=raplost-raprand;
		  rapdif=raprand;
	        }
		
		double Temp=thermal(spe,ptrand);

                mtrand=sqrt(pow(ptrand,2.)+pow(masstra[spe],2.));

		//Usual expression
		cooper=exp(-mtrand/Temp*cosh(rapdif))*mtrand/pow(Temp,5.)*cosh(rapdif)*
		(ptrand*3.*ptlost/mtlost*cos(phidif)+mtrand*cosh(rapdif))/normcoop[spe];

		//Second order expression
                //cooper=exp(-mtrand/Temp*cosh(rapdif))*mtrand/pow(Temp,5.)*cosh(rapdif)*
                //(ptrand*3.*ptlost/mtlost*cos(phidif)+mtrand*cosh(rapdif)+pow(ptrand*ptlost/mtlost*cos(phidif),2.)*9./2./mtrand/cosh(rapdif))/normcoop[spe];

		//Resummed expression
		//cooper=exp(-mtrand/Temp*cosh(rapdif))*pow(mtrand,2.)/pow(Temp,5.)*pow(cosh(rapdif),2.)*
                //exp(ptrand*3.*ptlost/mtlost*cos(phidif)/mtrand/cosh(rapdif))/normcoop[spe];

		if (cooper!=cooper) { 
			cout << " Cooper NaN \n", cooper=0.; cout << " mtrand= " << mtrand << " Temp= " << Temp << " rapdif= " << rapdif << " phidif= " << phidif << " ptrand= " << ptrand << endl; 
			cout << " pxrand= " << pxrand << " pyrand= " << pyrand << endl;
			//exit(0);
		}
		//Check maximum value in situ
		if (fabs(cooper)>maxcooper[spe]) maxcooper[spe]=fabs(cooper);
		//Adapt normalization so that maximum cannot be greater than 1
		if (fabs(cooper)>1.) 
		{
                	normcoop[spe]*=(fabs(cooper)+0.0001);
                        cooper=0.;
                }
		if (mode==-1) mc=fabs(cooper);
		else mc=cooper*wake[mode].GetStatus();
		
		randian=nr.rando();
	
	} while(mc<randian);
	
	//Set status
	double status;
	if (cooper>0.) status=1.;
	else status=-1.;
	
	//Set charge
	int charge=set_charge(spe, nr);
	
	//Fill wake vector with this hadron
	vector<double> p;
	p.push_back(pxrand), p.push_back(pyrand), p.push_back(ptrand*sinh(raprand+raplost)), p.push_back(sqrt(pow(ptrand*cosh(raprand+raplost),2.)+pow(masstra[spe],2.)));
	wake.push_back( Wake ( p, masstra[spe], charge, spe, status ) );	
}

double thermal(int spe, double ptrand)
{
	double temp;
	if (spe==0)
	{
		temp=0.211501*pow(ptrand,0.275362);
		if (temp>0.4) temp=0.4;
		if (temp<0.19) temp=0.19;
	}
	else {
		temp=0.33*pow(ptrand,0.3);
                if (temp>0.4) temp=0.4;
                if (temp<0.15) temp=0.15;
	}
	return temp;
}

int set_charge(int spe, numrand &nr)
{
	double ranchar=nr.rando();
        int charge;
        if (spe==0)
        {
                if (ranchar>2./3.) charge=0;
                else if (ranchar>1./3.) charge=1;
                else charge=-1;
        }
        else
        {
                if (ranchar>1./2.) charge=1;
                else charge=-1;
        }
	return charge;
}

double rapid(double pt, double pz)
{
        return 1./2.*log((sqrt(pow(pt,2.)+pow(pz,2.))+pz)/(sqrt(pow(pt,2.)+pow(pz,2.))-pz));
}
