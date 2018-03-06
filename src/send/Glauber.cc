#include <math.h>
#include <fstream>
#include <iostream>

#include "Random.h"

using namespace std;

double TA[4000];		//Nuclear Thickness, as a function of distance squared from nucleus center
double d2step;			//Step in distance squared
double norm;			//Normalization 

double bmin, bmax;		//Min and Max Impact Parameter

void gxy(double &x, double &y, double &b, numrand &nr);
void set_b_range(int cent);
void read_nuclear();
double gTAA(double x, double y, double b);

void read_nuclear()
{
	//Read Thickness function
        std::ifstream glauber("../hydro502/TAb2LL.dat");
        double b2;
        for (unsigned int a=0; a<4000; a++)
        {
                glauber >> b2 >> TA[a];
                if (a==1) d2step=b2;
        }
}

void set_b_range(int cent)
{
	double imp_param[]={0., 3.5, 4.94, 6.98, 8.55, 9.88, 11.04, 12.09, 13.05};
        bmin=imp_param[cent];
        bmax=imp_param[cent+1];
	//Set normalization for Overlap function gTAA
	norm=1.;
        norm=gTAA(0.,0.,bmin);
	
}

void gxy(double &x, double &y, double &b, numrand &nr) {

        double rho,phi;
        double P;

        naiguels:
        b=sqrt((bmax*bmax-bmin*bmin)*nr.rando()+bmin*bmin);
        rho=sqrt(150.*nr.rando());
        phi=2.*3.141592654*nr.rando();
        x=rho*cos(phi);
        y=rho*sin(phi);
        P=nr.rando();
        if(P>gTAA(x,y,b)) goto naiguels;
}

double gTAA(double x, double y, double b)
{
        int il, irr;
        double rho2;
	double use;

        rho2=pow(x+b/2.,2.)+y*y;
        il=int(rho2/d2step);
        rho2=pow(x-b/2.,2.)+y*y;
        irr=int(rho2/d2step);
        use=0.;
        if(il<4000 && irr<4000) {
        	use=TA[il]*TA[irr]/norm;
        }
        return use;
}
