#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <assert.h>

#include "Random.h"

#include "global.h"

using std::vector;
using namespace std;

double glaub[200][200];
int g_maxx=200;
int g_maxy=200;
double g_deltax=0.098650;
double g_deltay=0.098650;

void read_nuclear(int nhyd, std::string cent);
void gxy(double &x, double &y, numrand &nr);
double gGlaub(double x, double y);

void read_nuclear(int nhyd, std::string cent)
{
	char glauFile[200];
        sprintf(glauFile,"/gs/project/cqn-654-ad/mayank/events_for_jets/PbPb_%s_2p76/job-%i/results/u_field_1.dat",cent.c_str(),nhyd);
        ifstream initial (glauFile);

	assert(!initial.fail());

	cout << " Reading Initial Energy Density..." << endl;
	string s, sub;
	//First line: some crap
	getline(initial,s);
        //Rest of lines: crap, x, y, edensity(MeV/fm^3), 7 x crap
	double crap, x, y, edens;
	double maxedens=0.;
	int ix=0, iy=0;
	do {
		getline(initial,s);
		if (s.empty()) continue;
		istringstream iss(s);
		iss >> crap >> x >> y >> edens;
		if (edens>maxedens) maxedens=edens;
		//cout << " x= " << x << " y= " << y << " edens= " << edens << endl;
		glaub[ix][iy]=edens;
		iy+=1;
		if (iy==g_maxy) ix+=1, iy=0;	
	} while(ix<g_maxx);
	cout << " Finish Reading Initial Energy Density." << endl;
	cout << " Max Energy Density = " << maxedens << endl;
	//Set Maximum to 1
	for (int i=0; i<g_maxx; i++)
	{
		for (int j=0; j<g_maxy; j++)
		{
			//double xt = -9.855135e+00+g_deltax*double(i);
			//double yt = -9.855135e+00+g_deltay*double(j);
			//cout << xt << " " << yt << " " << gGlaub(xt,yt) << endl;
			glaub[i][j]/=maxedens;
		}
	}
}

void gxy(double &x, double &y, numrand &nr) {

        double rho,phi;
	double P;

        naiguels:
        rho=sqrt(150.*nr.rando());
        phi=2.*3.141592654*nr.rando();
        x=rho*cos(phi);
        y=rho*sin(phi);
	P=nr.rando();
        if(P>gGlaub(x,y)) goto naiguels;
}

double gGlaub(double x, double y)
{
	double gdens=0.;
	
	int ix, dx, iy, dy;
	
	if (x>=0.) {
		ix = int(x/g_deltax)+(g_maxx)/2;
		dx = (x - double(ix-(g_maxx)/2)*g_deltax)/g_deltax;
	}
	else {
		ix = int(x/g_deltax)+(g_maxx)/2-1;
                dx = (x - double(ix-(g_maxx)/2)*g_deltax)/g_deltax;
	}
	
	if (y>=0.) {
                iy = int(y/g_deltay)+(g_maxy)/2;
                dy = (y - double(iy-(g_maxy)/2)*g_deltay)/g_deltay;
        }
        else {
                iy = int(y/g_deltay)+(g_maxy)/2-1;
                dy = (y - double(iy-(g_maxy)/2)*g_deltay)/g_deltay;
        }

	if (ix<0 || ix>=g_maxx-1 || iy<0 || iy>=g_maxy-1) return gdens; 
	gdens=glaub[ix][iy]*(1.-dx)*(1.-dy);
	gdens+=glaub[ix][iy+1]*(1.-dx)*dy;
	gdens+=glaub[ix+1][iy]*dx*(1.-dy);
	gdens+=glaub[ix+1][iy+1]*dx*dy;

	if (gdens>1.) cout << " gGlaub not properly normalised: gdens = " << gdens << endl;
	return gdens;
}
