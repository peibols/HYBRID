#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "Random.h"

#include "global.h"
#include "vector_operators.h"

using std::vector;
using namespace std;

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod);
double normalise(vector<double> &p);
vector<double> vec_prod(vector<double> a, vector<double> b);
double gQ(double Del, numrand &nr);
void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr);
double gVx(double tau, double x, double y, double eta);
double gVy(double tau, double x, double y, double eta);
double gVz(double tau, double x, double y, double eta);
double gT(double tau, double x, double y, double eta);

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod)
{
	//Source File, append
	ofstream source_file;	
	source_file.open(Sfile, ios::app);

	double Tc;
	if (tmethod==0) Tc=0.170;
	else Tc=0.145;

	double tot=pos[3]+tof;		//Final time

	double ei=p[3];			//Initial energy

	double f_dist=0.;		//Traversed distance in Fluid Frame
	double l_dist=0.;		//Traversed distance in Lab Frame

	double CF;
	if (id==21) CF=pow(9./4.,1./3.);	//If gluon, color charge dependence is ratio of casimirs to power 1/3
	else CF=1.;

	int marker=0;		//If one, exit loop
	double step=0.02;	//Time step in LAB frame

	vector<double> w = p/p[3];	//4-velocity
	//cout << " x= " << pos[0] << " y= " << pos[1] << " z= " << pos[2] << " t= " << pos[3] << endl;
	//cout << " px= " << p[0] << " py= " << p[1] << " pz= " << p[2] << " en= " << p[3] << endl;
	//cout << " wx= " << w[0] << " wy= " << w[1] << " wz= " << w[2] << " we= " << w[3] << endl;
	do {
		//Keep 4momentum before applying quenching this step
		vector<double> p_prev = p;

		if (pos[3]==tot) marker=1;
		if (pos[3]>tot) cout << " Warning: Went beyond tot= " << tot << " t= " << pos[3] << endl;	
		//Proper time
		double tau=sqrt(pos[3]*pos[3]-pos[2]*pos[2]);
		if (tau!=tau) tau=0., cout << " TAU Not a number z= " << pos[2] << " t= " << pos[3] << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << "\n";

		//Rapidity
		double eta = 1./2.*log((pos[3]+pos[2])/(pos[3]-pos[2]));
                //if (abs(eta)>=10.) cout << " eta= " << eta << endl;
		if (eta!=eta && tau>0.) {
                	cout << " Eta is NaN= " << eta << " t= " << pos[3] << " z= " << pos[2] << endl;
                        eta=0.;
                }

		int will_hot=0;	//Advance variable (to reach hot zones)
		if (tau>=0.6)	//Smooth profile starting time
		{
			vector<double> v;
			double vx=gVx(tau,pos[0],pos[1],eta);
			double vy=gVy(tau,pos[0],pos[1],eta);
			//double vz=gVz(tau,pos[0],pos[1],eta);
			//**BOOST INVARIANT**
			double vz=pos[2]/pos[3];
			double frap=atanh(vz);
			vx/=cosh(frap);
			vy/=cosh(frap);
			//*******************
			v.push_back(vx), v.push_back(vy), v.push_back(vz), v.push_back(1.);
			
			double v2=pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.);
			double w2=pow(w[0],2.)+pow(w[1],2.)+pow(w[2],2.);
			double vscalw=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
			if (v2>=1.) v2=0.999999999, cout << " V2 >= 1 \n";
			double lore=1./sqrt(1.-v2);

			l_dist+=step;
			f_dist+=step*sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));

			double temp=gT(tau,pos[0],pos[1],eta);
			
			//Safe way to exit the plasma: check whether temperature will be above Tc in the next 1000 steps
			if (temp<Tc)
			{
				//Check whether temperature increases in its way
				for (unsigned int j=1; j<1000; j++)
				{
					vector<double> tpos=pos+w*step*double(j);
					if (tpos[3]>tot) break;
					tau=sqrt(tpos[3]*tpos[3]-tpos[2]*tpos[2]);
					eta = 1./2.*log((tpos[3]+tpos[2])/(tpos[3]-tpos[2]));
					if (gT(tau,tpos[0],tpos[1],eta)>Tc)
					{
						//cout << " Will Hot! from time= " << pos[3] << " at j= " << j << endl;
						will_hot=int(j);
						break;
					}
				}
				if (will_hot==0)
				{
					pos+=w*(tot-pos[3]);
					marker=1;
				}
			}

			//Now broad&quench
			if (p[3]>0. && temp>=Tc)
			{
				//Broadening
				if (kappa!=0.) trans_kick(w,w2,v,p,temp,vscalw,lore,step,kappa,nr);
				//Quenching
				if (alpha!=0.)
				{
					double Efs=ei*lore*(1.-vscalw);
					double tstop=0.2*pow(Efs,1./3.)/(2.*pow(temp,4./3.)*alpha)/CF;
					double beta=tstop/f_dist;
					if (beta>1.)
					{
						double intpiece=Efs*step*4./(3.141592)*(1./(beta*tstop*sqrt(beta*beta-1.)));
						//Update 4momentum
						double quench=(p[3]-intpiece)/p[3];
						p*=quench;
					}
					else
					{
						p[3]=0.;
					}
				}
			} 
		}

		//Source File: Write spacetime coordinates if there was quenching
		if (p[3]!=p_prev[3]) source_file << tau << " " << pos[0] << " " << pos[1] << " " << eta << " ";

		//If parton gets totally quenched, exit 
		if (p[3]<=0.)
		{
			marker=1;
			for (unsigned int i=0; i<4; i++) p[i]=0.;
		}
		else {
			//Manually protect very soft particles from getting kicks that yield velocities greater than 1
			for (unsigned int i=0; i<3; i++) { if (p[i]>p[3]) p[i]=0.99999*p[3]; }
			//Update kinematical quantities, with the possibility of advancing to hot regions
			w=p/p[3];
			double tstep=max(double(will_hot),1.)*step;
			if (pos[3]+tstep>tot)
                	{
                	        tstep=tot-pos[3];
                	}
			if (marker!=1) pos+=w*tstep;
		}
		
		//Source File: write Jmu if there was quenching
		if (p[3]!=p_prev[3]) source_file << -p[3]+p_prev[3] << " " << -p[0]+p_prev[0] << " " << -p[1]+p_prev[1] << " " << -p[2]+p_prev[2] << endl;
	} while (marker==0);
}

void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr)
{
	vector<double> e1 = vec_prod(w,v);
	double Ne1=normalise(e1);

	double Nw=sqrt(w2);
	vector<double> l = vec_prod(w,e1)/Nw;

	double uscalW=lore*(1.-vscalw);
	double uscall=lore*(-v[0]*l[0]-v[1]*l[1]-v[2]*l[2]);
	double W2=1.-w2;

	vector<double> Wp = w-v*lore*W2/uscalW;

	double Nalpha=-uscall*uscalW/(pow(uscalW,2.)-W2);
	double NN=1.+W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2);
	//In some rare situations, this norm squared can be negative. Only do kick otherwise
	if (sqrt(NN)!=sqrt(NN)) cout << " negative NN " << endl;
	else {
		vector<double> e2 = (l+Wp*Nalpha)/sqrt(NN);

		double Ef=p[3]*lore*(1.-vscalw);
		double wf2=1.-W2/pow(lore*(1.-vscalw),2.);
		double DelQ2=kappa*pow(temp,3.)*lore*(1.-vscalw)*step*5.;	//Should this step be in f_frame?

		double qfac=0.;
		//Only do kick if energy is greater than temperature
		if (Ef>temp)
		{
			do
			{
				qfac=gQ(DelQ2,nr);
			} while(qfac>Ef*sqrt(wf2));
		}

		double qbeta=sqrt(1.-qfac*qfac/Ef/Ef/wf2)-1.;

		double qphi=2.*3.141592654*nr.rando(); 

		vector<double> e = e1*cos(qphi)+e2*sin(qphi); 

		vector<double> Wt = (w-v*uscalW*lore)/lore/(1.-vscalw);

		//Update 4momentum
		p+=Wt*qbeta*Ef+e*qfac;
	}
}

double gQ(double Del, numrand &nr) //Need to double check how optimized this is
{

        double gaussq;
        double qfac;
        double gaussmax=sqrt(2./Del/exp(1.));

        qhatelsen:
        qfac=5.*sqrt(Del/2.)*nr.rando();
        gaussq=2.*qfac/Del*exp(-qfac*qfac/Del)/gaussmax;
        double nrand=nr.rando();
        if (nrand>gaussq) goto qhatelsen;

	return qfac;
}
