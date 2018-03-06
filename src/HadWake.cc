#include <vector>
#include <cmath>
#include <iostream>

#include "Random.h"
#include "Quench.h"
#include "Wake.h"

#include "vector_operators.h"

using std::vector;
using namespace std;

double transcut=0.5;		//Lower threshold for transverse mass for back-reaction
double basesig=0.45;		//Starting gaussian width for pass function in Metropolis
int Nrun=4000;			//Maximum iterations of Metropolis
double maxptsq=pow(5.,2.);	//Maximum squared pt for MC
double maxrap=5.;		//Maximum absolute rapidity for MC
double tole=0.4;		//Non-conservation tolerance in GeV
double masspi=0.1396;		//Pion mass
double masspro=0.938;		//Proton mass
double masstra[2]={masspi,masspro};

double normcoop[2]={30.,30.};		//Norm for OneBody MC, adaptable 

double maxcooper[2]={0.,0.};		//Record highest value of OneBody MC

int toomuch=0;			//How many times we gave up in Metropolis

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);
double rapid(double pt, double pz);
int set_charge(int spe, numrand &nr);
double thermal(int spe, double ptrand);
void one_body(vector<Wake> &wake, vector<double> delta, double ptlost, double mtlost, double raplost, numrand &nr, int spe, int mode);
vector<double> vec_abs(vector<double> p);

//#define DO_PRINT

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr)
{
	for (unsigned int i=0; i<quenched.size(); i++)
	{
		vector<double> delta = partons[i].vGetP()-quenched[i].vGetP();
		if (delta[3]<=0.) continue;				//Swiftly skip guys who didn't lose energy
		double ptlost=sqrt(pow(delta[0],2.)+pow(delta[1],2.));	//Lost Pt
		double raplost=rapid(ptlost,delta[2]);	//Want this or pseudorapidity?
		if (raplost!=raplost) {cout << " RapLost NaN: Dx= " << delta[0] << " Dy= " << delta[1] << " Dz= " << delta[2] << " De= " << delta[3] << endl; continue;}
		double mtlost=delta[3]/cosh(raplost);			//Lost transverse mass
		//Reasonable conditions to do back-reaction
		if (delta[3]>=0. && delta[3]<1000. && mtlost>transcut && fabs(raplost)<3.5)
		{
			//Declare wake vector for this particle
			vector<Wake> pwake;
			//Declare variables here because of the goto:
			vector<double> momback (4,0.);
			vector<double> pmomback (4,0.);
			vector<double> dif (4,0.);
			double msigma, pass, newpass, passrand;
			int spe, mode;
			int runi, encallao;
			int numenc=0;
			double clocklim=0.05;
                        int tooclock=0;
			//GoTo flag
			thisis:
			clock_t startClock = clock();
			for (unsigned int j=0; j<4; j++) momback[j]=0.;
			runi=0, encallao=0;
			//Initial List of Hadrons, approximately satisfy energy conservation
			do
			{
				//Proton or Pion
        			if (nr.rando()<=0.05) spe=1;    //is proton
        			else spe=0;                     //is pion
				//OneBody Dist: last arg=-1 because it is creation of initial list
				one_body(pwake, delta, ptlost, mtlost, raplost, nr, spe, -1);
				momback+=pwake[pwake.size()-1].vGetP()*pwake[pwake.size()-1].GetStatus();
			} while(momback[3]<delta[3]);		
			//Absolute difference wrt Wake Momentum
			dif=vec_abs(delta-momback);
			//Select random particle, flip it and see whether it improves conservation
			msigma=basesig;
			//cout << "BEF pwake.size()= " << pwake.size() << endl;
			do
			{
				//Select random particle
				mode=int(double(pwake.size())*nr.rando());
				//Get the species it was
				spe=pwake[mode].GetId();
				//Generate new particle, same species and status than selected
				one_body(pwake, delta, ptlost, mtlost, raplost, nr, spe, mode);
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
				}
				else
				{
					//If conservation not improved, accept with probability newpass/pass
					passrand=nr.rando();
					if (newpass/pass > passrand)
					{
						pwake.erase(pwake.begin()+mode);
						runi+=1;
						encallao=0;
						momback=pmomback;
					}
					else {
						pwake.erase(pwake.begin()+pwake.size()-1);
						dif=vec_abs(delta-momback);
						encallao+=1;
					}
				}
				//Tricks to improve convergence
				//Adapt width of gaussian in pass function. Suddenly open, and rapidly close again
				if (encallao>10000) msigma=1./double(encallao-9999);
				if (encallao>30000) msigma=2./double(encallao-29999);
				if (encallao>40000) msigma=2.5/double(encallao-39999);	
				if (msigma<basesig) msigma=basesig;	//Never below minimum;
				//If conservation reasonably good, don't open width
				if (sqrt(pow(dif[0],2.)+pow(dif[1],2.)+pow(dif[2],2.)+pow(dif[3],2.))<15.) msigma=basesig;
				//If truly stuck, start over. If 5 times truly stuck, give up
				if (encallao>50000)
				{
					numenc+=1;
					if (numenc>5) {toomuch+=1; break;}
					goto thisis;
				}
				//If time spent surpasses clocklim, start over, and wait longer for next attempt.
				clock_t endClock = clock();
				if (double((endClock - startClock)) / CLOCKS_PER_SEC>clocklim)
				{
					//cout << " CLOCK !!! \n";
					clocklim+=0.05;
					tooclock+=1;
					//If 20 times start over for clock, give up
					if (tooclock>20) {toomuch+=1; break;}
                                        goto thisis;
				}
			} while(runi<Nrun && (dif[0]>tole || dif[1]>tole || dif[2]>tole || dif[3]>tole));
			if (runi>=Nrun) toomuch+=1;
		#ifdef DO_PRINT
			cout << " DeltaX= " << delta[0] << " DeltaY= " << delta[1] << " DeltaZ= " << delta[2] << " DeltaE= " << delta[3] << endl;
			cout << " DifX= " << dif[0] << " DifY= " << dif[1] << " DifZ= " << dif[2] << " DifE= " << dif[3] << endl;
			cout << " pwake.size()= " << pwake.size() << endl;
		#endif
		//cout << "AFT pwake.size()= " << pwake.size() << endl;
			//Fill total wake with pwake
			for (unsigned int k=0; k<pwake.size(); k++)
			{
				wake.push_back ( pwake[k] );
			}
			pwake.clear();
		}
	}
#ifdef DO_PRINT
	cout << " wake.size= " << wake.size() << " toomuch= " << toomuch << endl;
	cout << " maxcooper pion= " << maxcooper[0] << " maxcooper prot= " << maxcooper[1] << endl;
	cout << " normcoop pion= " << normcoop[0] << " normcoop prot= " << normcoop[1] << endl;
	/*
	for (unsigned int l=0; l<wake.size(); l++)
	{
		cout << " Wake " << l << endl;
		wake[l].display();
		cout << endl;
	}
	*/
#endif
}

void one_body(vector<Wake> &wake, vector<double> delta, double ptlost, double mtlost, double raplost, numrand &nr, int spe, int mode)
{
	double mc=0.;
	double cooper=0.;
	double randian=1.;
	double pxrand, pyrand, raprand, ptrand;
	do {
		ptrand=sqrt(maxptsq*nr.rando());	
		double phirand=2.*3.141592*nr.rando();
		raprand=maxrap*(-1.+2.*nr.rando());
		pxrand=ptrand*cos(phirand);
              	pyrand=ptrand*sin(phirand);
                double mtrand=sqrt(pow(ptrand,2.)+pow(masstra[spe],2.));
		double phidif=acos((delta[0]*pxrand+delta[1]*pyrand)/(ptlost*ptrand));
		double rapdif=raplost-raprand;
	
		double Temp=thermal(spe,ptrand);

		//CHECK!
		cooper=exp(-mtrand/Temp*cosh(rapdif))*mtrand/pow(Temp,5.)*cosh(rapdif)*
		(ptrand*3.*ptlost/mtlost*cos(phidif)+mtrand*cosh(rapdif))/normcoop[spe];

		if (cooper!=cooper) cout << " Cooper NaN \n", cooper=0.;
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
	p.push_back(pxrand), p.push_back(pyrand), p.push_back(ptrand*sinh(raprand)), p.push_back(sqrt(pow(ptrand*cosh(raprand),2.)+pow(masstra[spe],2.)));
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
        //CHECK!
        return 1./2.*log((sqrt(pow(pt,2.)+pow(pz,2.))+pz)/(sqrt(pow(pt,2.)+pow(pz,2.))-pz));
}
