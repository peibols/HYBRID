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

double transcut=0.5;		//Lower threshold for transverse mass for back-reaction
double basesig=0.65;		//Starting gaussian width for pass function in Metropolis
int Nrun=800000;			//Maximum iterations of Metropolis
double maxptsq=pow(3.5,2.);	//Maximum squared pt for MC
double maxrap=2.5;		//Maximum absolute rapidity for MC
double tole=0.4;		//Non-conservation tolerance in GeV
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

        int mom1=-1000, mom2=-1000;
	//Group all partons in terms of their originator
	for (unsigned int i=0; i<partons.size(); i++)
	{
          if (partons[i].GetOrig()=="hs" && mom1==-1000) mom1=int(i);
	  else if (partons[i].GetOrig()=="hs") { mom2=int(i); break; }

	}
	vector<double> delta1 (4,0.);
	vector<double> delta2 (4,0.);
        for (unsigned int i=0; i<partons.size(); i++)
	{
	  if (partons[i].GetD1()!=-1) continue;
	  if (partons[i].GetOrig()=="rem") continue;
	  int fp=int(i);
	  if (fp==mom1) delta1+=(partons[i].vGetP()-quenched[i].vGetP());
	  if (fp==mom2) delta2+=(partons[i].vGetP()-quenched[i].vGetP());
	  while (true) {
	    int tmom=partons[fp].GetMom();
            if (tmom==mom1) { delta1+=(partons[i].vGetP()-quenched[i].vGetP()); break; }
            if (tmom==mom2) { delta2+=(partons[i].vGetP()-quenched[i].vGetP()); break; }
	    if (tmom==-1) break;
	    fp=tmom;
	  }
	}

	for (unsigned int i=0; i<2; i++)
	{
		//vector<double> delta=partons[i].vGetP()*0.2;
		vector<double> delta;
		if (i==0) delta=delta1;
		else delta=delta2;
		if (delta[3]<=0.) continue;				//Swiftly skip guys who didn't lose energy
		cout << " DOING DeltaX= " << delta[0] << " DeltaY= " << delta[1] << " DeltaZ= " << delta[2] << " DeltaE= " << delta[3] << endl;
		double ptlost=sqrt(pow(delta[0],2.)+pow(delta[1],2.));	//Lost Pt
		double raplost=rapid(ptlost,delta[2]);	//Want this or pseudorapidity?
		if (raplost!=raplost) {cout << " RapLost NaN: Dx= " << delta[0] << " Dy= " << delta[1] << " Dz= " << delta[2] << " De= " << delta[3] << endl; continue;}
		double mtlost=delta[3]/cosh(raplost);			//Lost transverse mass
		mtlost=delta[3];					//Boost to wake frame
		//Reasonable conditions to do back-reaction
		if (delta[3]>=0. && delta[3]<3000. && mtlost>transcut && fabs(raplost)<2.5)
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
			//cout << " Difs= " << delta[0]-momback[0] << " " << delta[1]-momback[1] << " " << delta[2]-momback[2] << " " << delta[3]-momback[3] << endl;
			//Select random particle, flip it and see whether it improves conservation
			msigma=basesig;
			msigma=sqrt(pow(dif[0],2.)+pow(dif[1],2.)+pow(dif[2],2.)+pow(dif[3],2.))/sqrt(log(2.));
			//cout << "BEF pwake.size()= " << pwake.size() << endl;
			do
			{
				//Select random particle
				mode=int(double(pwake.size())*nr.rando());
				//Get the species it was
				spe=pwake[mode].GetId();
				//Generate new particle, same species and status than selected
				one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, mode);
				//Define pass function	
				//msigma=sqrt(pow(dif[0],2.)+pow(dif[1],2.)+pow(dif[2],2.)+pow(dif[3],2.))/sqrt(log(2.));
				//msigma=sqrt(pow(dif[0],2.)+pow(dif[1],2.))/sqrt(log(2.));
				pass=exp((-pow(dif[0],2.)-pow(dif[1],2.)-pow(dif[2],2.)-pow(dif[3],2.))/pow(msigma,2.));
				//pass=pow(pow(dif[0],2.)+pow(dif[1],2.)+pow(dif[2],2.)+pow(dif[3],2.),3.);
				//pass=exp((-pow(dif[0],2.)-pow(dif[1],2.))/pow(msigma,2.));
				//Update sum of 4momentum: remove previous particle, add new particle
				pmomback=momback+(pwake[pwake.size()-1].vGetP()-pwake[mode].vGetP())*pwake[mode].GetStatus();
				//Update difference wrt Lost Momentum
				dif=vec_abs(delta-pmomback);
				//Compute new pass function
				newpass=exp((-pow(dif[0],2.)-pow(dif[1],2.)-pow(dif[2],2.)-pow(dif[3],2.))/pow(msigma,2.));
				//newpass=pow(pow(dif[0],2.)+pow(dif[1],2.)+pow(dif[2],2.)+pow(dif[3],2.),3.);
				//newpass=exp((-pow(dif[0],2.)-pow(dif[1],2.))/pow(msigma,2.));
				//If conservation improved, accept substitution
				if (newpass>pass)
				//if (newpass<pass)
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
				//Tricks to improve convergence
				
				//Adapt width of gaussian in pass function. Suddenly open, and rapidly close again
				if (encallao>10000) msigma=1./double(encallao-9999);
				if (encallao>30000) msigma=2./double(encallao-29999);
				if (encallao>40000) msigma=2.5/double(encallao-39999);	
				if (msigma<basesig) msigma=basesig;	//Never below minimum;
				//If conservation reasonably good, don't open width
				if (sqrt(pow(dif[0],2.)+pow(dif[1],2.)+pow(dif[2],2.)+pow(dif[3],2.))<15.) msigma=basesig;
			
				//cout << " Difs= " << delta[0]-momback[0] << " " << delta[1]-momback[1] << " " << delta[2]-momback[2] << " " << delta[3]-momback[3] << endl;	
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
					clocklim+=0.01;
					tooclock+=1;
					//If 20 times start over for clock, give up
					if (tooclock>60) {toomuch+=1; break;}
                                        goto thisis;
				}
			} while(runi<Nrun && (dif[0]>tole || dif[1]>tole || dif[2]>tole || dif[3]>tole));
			if (runi>=Nrun) toomuch+=1;
		#ifdef DO_PRINT
			if (dif[0]>tole || dif[1]>tole || dif[2]>tole || dif[3]>tole)
			{
			  cout << " runi= " << runi << endl;
			  cout << " DeltaX= " << delta[0] << " DeltaY= " << delta[1] << " DeltaZ= " << delta[2] << " DeltaE= " << delta[3] << endl;
			  cout << " DifX= " << dif[0] << " DifY= " << dif[1] << " DifZ= " << dif[2] << " DifE= " << dif[3] << endl;
			  cout << " pwake.size()= " << pwake.size() << endl;
			}
		#endif
		//cout << "AFT pwake.size()= " << pwake.size() << endl;
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
			pwake.clear();
		}
	}
	totalsave+=1;
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
	double mphidif;
	do {
		ptrand=sqrt(maxptsq*nr.rando());	
		double phirand=2.*3.141592*nr.rando();
		if (mode!=-1)
		{
		  //mphidif=acos((delta[0]*momback[0]+delta[1]*momback[1])/(ptlost*sqrt(pow(momback[0],2.)+pow(momback[1],2.))));	
		  //phirand=atan2(delta[1],delta[0])-mphidif+3.141592;
		  //phirand=atan2(delta[1]-momback[1],delta[0]-momback[0]);
		  //mphidif=phirand-atan2(delta[1],delta[0]);
		}
		raprand=maxrap*(-1.+2.*nr.rando());
		pxrand=ptrand*cos(phirand);
              	pyrand=ptrand*sin(phirand);
		//cout << " MomPx= " << momback[0]-delta[0] << " thisPx= " << pxrand << " MomPy= " << momback[1]-delta[1] << " thisPy= " << pyrand << endl;
                double mtrand=sqrt(pow(ptrand,2.)+pow(masstra[spe],2.));
		double phidif=acos((delta[0]*pxrand+delta[1]*pyrand)/(ptlost*ptrand));
		//if (mode!=-1) cout << " mphidif= " << mphidif << " phidif= " << phidif << endl;
		double rapdif=raplost-raprand;
		rapdif=raprand;
	
		double Temp=thermal(spe,ptrand);

		//Usual expression
		cooper=exp(-mtrand/Temp*cosh(rapdif))*mtrand/pow(Temp,5.)*cosh(rapdif)*
		(ptrand*3.*ptlost/mtlost*cos(phidif)+mtrand*cosh(rapdif))/normcoop[spe];

		//Second order expression
                //cooper=exp(-mtrand/Temp*cosh(rapdif))*mtrand/pow(Temp,5.)*cosh(rapdif)*
                //(ptrand*3.*ptlost/mtlost*cos(phidif)+mtrand*cosh(rapdif)+pow(ptrand*ptlost/mtlost*cos(phidif),2.)*9./2./mtrand/cosh(rapdif))/normcoop[spe];

		//Resummed expression
		//cooper=exp(-mtrand/Temp*cosh(rapdif))*pow(mtrand,2.)/pow(Temp,5.)*pow(cosh(rapdif),2.)*
                //exp(ptrand*3.*ptlost/mtlost*cos(phidif)/mtrand/cosh(rapdif))/normcoop[spe];

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
        //CHECK!
        return 1./2.*log((sqrt(pow(pt,2.)+pow(pz,2.))+pz)/(sqrt(pow(pt,2.)+pow(pz,2.))-pz));
}
