#include <vector>
#include <iostream>
#include <cmath>
#include "Random.h"

#include "vector_operators.h"

using std::vector;
using namespace std;

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod);
double normalise(vector<double> &p);
vector<double> vec_prod(vector<double> a, vector<double> b);
double gQ(double Del, numrand &nr);
void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr);
double gVx(double tau, double x, double y);
double gVy(double tau, double x, double y);
double gT(double tau, double x, double y);

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod)
{
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
	double step=0.01;	//Time step in LAB frame

	vector<double> w = p/p[3];	//4-velocity
	//cout << " px= " << p[0] << " py= " << p[1] << " pz= " << p[2] << " en= " << p[3] << endl;
	//cout << " wx= " << w[0] << " wy= " << w[1] << " wz= " << w[2] << " we= " << w[3] << endl;
	do {
		if (pos[3]==tot) marker=1;
		
		double tau=sqrt(pos[3]*pos[3]-pos[2]*pos[2]);
		if (tau!=tau) tau=0., cout << " TAU Not a number z= " << pos[2] << " t= " << pos[3] << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << "\n";

		if (tau>=0.6)
		{
			//Boost invariant fluid velocities
			//double vz=z/t;
			//double frap=atanh(vz);
			//double vx=gVx(tau,x,y)/cosh(frap);
			//double vy=gVy(tau,x,y)/cosh(frap);

			vector<double> v;
			v.push_back(gVx(tau,pos[0],pos[1])/cosh(atanh(pos[2]/pos[3]))), v.push_back(gVy(tau,pos[0],pos[1])/cosh(atanh(pos[2]/pos[3]))), v.push_back(pos[2]/pos[3]), v.push_back(1.);
			double v2=pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.);
			double w2=pow(w[0],2.)+pow(w[1],2.)+pow(w[2],2.);
			double vscalw=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
			if (v2>=1.) v2=0.999999999, cout << " V2 >= 1 \n";
			double lore=1./sqrt(1.-v2);

			l_dist+=step;
			f_dist+=step*sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));

			double temp=gT(tau,pos[0],pos[1]);
			
			if (temp<Tc)
			{
				//Check for partons that are created cold, but point towards hot. If point outwards, set temp=0 and let if fly until tot
				double cospart=pos[0]*w[0]+pos[1]*w[1];
				if (cospart>0.)
				{
					temp=0.;
					pos+=w*(tot-pos[3]);
				}
			}

			//Check temperature exit condition
			if (temp==0.) marker=1;
		
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
		//If parton gets totally quenched, exit 
		if (p[3]<=0.)
		{
			marker=1;
			for (unsigned int i=0; i<4; i++) p[i]=0.;
		}
		else {
			//Update kinematical quantities
			//Manually protect very soft particles from getting kicks that yield velocities greater than 1
			for (unsigned int i=0; i<3; i++) { if (p[i]>p[3]) p[i]=0.99999*p[3]; }
			w=p/p[3];
			if (pos[3]+step>tot)
                	{
                	        step=tot-pos[3];
                	}
			pos+=w*step;
		}
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
