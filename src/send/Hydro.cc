#include <fstream>
#include <iostream>

double hydrot[200][101][101];   //Hydro Temperature
double hydroe[200][101][101];   //Hydro Energy Density 
double hydrox[200][101][101];   //Hydro Vx 
double hydroy[200][101][101];   //Hydro Vy

//Hydro Grid
int maxx=100;
int maxy=100;
double deltat=0.1;
double deltax=0.3;
double deltay=0.3;

double dt, dx, dy;
int it, ix, iy;

//Hydro Ini time
double tau0=0.6;

void read_hydro(int cent);
void getGrid(double tau, double x, double y, double* dt, double* dx, double* dy, int* it, int* ix, int* iy);
double gVx(double tau, double x, double y);
double gVy(double tau, double x, double y);
double gT(double tau, double x, double y);

using namespace std;

void read_hydro(int cent)
{
        //Choose Event Averaged Hydro File
        std::string hyf[]={"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70"};
        char inpFile[100];
        sprintf(inpFile,"../hydro502/hydroinfoPlaintxtHuichaoFormat_C%s.dat",hyf[cent].c_str());
        std::ifstream hydfile(inpFile);
        if (!hydfile.is_open()) cout << " Hydro file not found! ";
        //Read Hydro
        double enedat, tdat, vxdat, vydat;
        double tou, hor, ver;
	//int it, ix, iy;
	//double dt, dx, dy;
        while (hydfile >> hor >> ver >> tou >> enedat >> tdat >> vxdat >> vydat)
        {
                it = int((tou+deltat/2.-tau0)/deltat);
                ix = int((hor+deltax*maxx/2.+deltax/2.)/deltax);
                iy = int((ver+deltay*maxy/2.+deltay/2.)/deltay);

                hydrot[it][ix][iy]=tdat;
                hydroe[it][ix][iy]=enedat;
                hydrox[it][ix][iy]=vxdat;
                hydroy[it][ix][iy]=vydat;
        }
}

void getGrid(double tau, double x, double y, double* udt, double* udx, double* udy, int* uit, int* uix, int* uiy)
{
        *uit=int((tau-tau0)/deltat);
        *udt=(tau-tau0-double(*uit)*deltat)/deltat;

        if (y>0.) {
                *uiy = int(y/deltay)+maxy/2;
                *udy = (y - double(*uiy-maxy/2.)*deltay)/deltay;
        }
        else {
                *uiy = int(y/deltay)+maxy/2-1;
                *udy = (y - double(*uiy-maxy/2.)*deltay)/deltay;
        }

        if (x>0.) {
                *uix = int(x/deltax)+maxx/2;
                *udx = (x - double(*uix-maxx/2.)*deltax)/deltax;
        }
        else {
                *uix = int(x/deltax)+maxx/2-1;
                *udx = (x - double(*uix-maxx/2.)*deltax)/deltax;
        }
}

double gT(double tau, double x, double y)
{
        double gete=0.;
        double tau1=18.5;

        if (tau>tau1) {
                return gete;
        }
	//cout << " BEF " << dt << " " << dx << " " << dy << " " << it << " " << ix << " " << iy << "\n";	
        getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " AFT " << dt << " " << dx << " " << dy << " " << it << " " << ix << " " << iy << "\n";
        if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gete;
        gete=hydrot[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
        gete+=hydrot[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
        gete+=hydrot[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
        gete+=hydrot[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
        gete+=hydrot[it+1][ix+1][iy]*dx*(1.-dy)*dt;
        gete+=hydrot[it][ix+1][iy+1]*(1.-dt)*dy*dx;
        gete+=hydrot[it+1][ix][iy+1]*(1.-dx)*dt*dy;
        gete+=hydrot[it+1][ix+1][iy+1]*dx*dt*dy;

	//Return temp in GeV	
        return gete*0.2;
}

double gVx(double tau, double x, double y)
{
        double gvelx=0.;
        double tau1=18.5;

        if (tau>tau1) {
                return gvelx;
        }

        getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
        if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gvelx;
        gvelx=hydrox[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
        gvelx+=hydrox[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
        gvelx+=hydrox[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
        gvelx+=hydrox[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
        gvelx+=hydrox[it+1][ix+1][iy]*dx*(1.-dy)*dt;
        gvelx+=hydrox[it][ix+1][iy+1]*(1.-dt)*dy*dx;
        gvelx+=hydrox[it+1][ix][iy+1]*(1.-dx)*dt*dy;
        gvelx+=hydrox[it+1][ix+1][iy+1]*dx*dt*dy;

        return gvelx;
}

double gVy(double tau, double x, double y)
{
        double gvely=0.;
        double tau1=18.5;

        if (tau>tau1) {
                return gvely;
        }

        getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
        if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gvely;
        gvely=hydroy[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
        gvely+=hydroy[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
        gvely+=hydroy[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
        gvely+=hydroy[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
        gvely+=hydroy[it+1][ix+1][iy]*dx*(1.-dy)*dt;
        gvely+=hydroy[it][ix+1][iy+1]*(1.-dt)*dy*dx;
        gvely+=hydroy[it+1][ix][iy+1]*(1.-dx)*dt*dy;
        gvely+=hydroy[it+1][ix+1][iy+1]*dx*dt*dy;

        return gvely;
}
