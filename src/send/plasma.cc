// HYDRO QUENCHING

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

using namespace std;
int ir;
double TA[4000], step, bmin, bmax;
double norm=1.;
double tau0, deltat, deltax, deltay;
int maxt, maxx, maxy, it, ix, iy;
double hydrot[200][101][101], hydroe[200][101][101], hydrox[200][101][101], hydroy[200][101][101];
double dt, dx, dy;
double PrNn=0.;
double rando(int*);
double gTAA(double, double, double);
double gWN(double, double, double);
double gT(double, double, double);
double gVy(double, double, double);
double gVx(double, double, double);
double gE(double, double, double);
double getGrid(double, double, double, double *, double *, double *, int *, int *, int *);
double gQ(double);
double fuckedN=0.;
double totalN=0.;

double alpha;
double kappa;
int qmethod; //1 Coll, 2 Radiative, 3 Strong
int tmethod;

double tquench=0.6;
//ARGUMENTS: File, Alpha, Kappa, N, Model, Tmethod, Cent, Rand
int main (int argc, char** argv) {
	assert(argc==9);
	char inpFile[100], outFile[100];
	sprintf(outFile,"./%s.out",argv[1]);
        sprintf(inpFile,"./SAFEnoupTREEhatBINdogaN%s.txt",argv[1]);
	alpha=atof(argv[2]);
	kappa=atof(argv[3]);
	int N=atoi(argv[4]);
	qmethod=atoi(argv[5]);
	tmethod=atoi(argv[6]);
	int Cent=atoi(argv[7]);
	ir=atoi(argv[8]);

	//Input Files
        std::ifstream infile1(inpFile);
	std::ifstream infile2("./TAb2LL.dat");

	ofstream check;
	check.open (outFile);

	double P[1000][12], num, px, py, pz, E, q, mom, qge, ini, col, acol, ide;
	double x,y,z,b;
	double pexe, peye, peze;
	double gxy(double*, double*, double*);
	double gdE(double, double, double, double, double, double *, double *, double *, double, double, double *, double *, double *);
	int count=0;
	double ene=5.;

	ifstream infile3;
        if (Cent==0) {          //0-5%
                bmin=0.;
                bmax=3.5;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C0-5.dat");
        }
        if (Cent==1) {          //5-10%
                bmin=3.5;
                bmax=4.94;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C5-10.dat");
        }
        if (Cent==2) {          //10-20%
                bmin=4.94;
                bmax=6.98;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C10-20.dat");
        }
        if (Cent==3) {          //20-30%
                bmin=6.98;
                bmax=8.55;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C20-30.dat");
        }
	if (Cent==4) {          //30-40%
                bmin=8.55;
                bmax=9.88;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C30-40.dat");
        }
        if (Cent==5) {          //40-50%
                bmin=9.88;
                bmax=11.04;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C40-50.dat");
        }
        if (Cent==6) {          //50-60%
                bmin=11.04;
                bmax=12.09;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C50-60.dat");
        }
        if (Cent==7) {          //60-70%
                bmin=12.09;
                bmax=13.05;
                infile3.open ("./hydroinfoPlaintxtHuichaoFormat_C60-70.dat");
        }

	//Reading T(t)
	double b2;
	for (unsigned a=0; a<4000; a++) {
		infile2 >> b2 >> TA[a];
		if (a == 1) step = b2;
	}

	//Reading HydroTemp
	double enedat, tdat, vxdat, vydat;
	double tou, hor, ver;
	maxt=179;
	maxx=100;
	maxy=100;
	tau0=0.6;
	deltat=0.1;
	deltax=0.3;
	deltay=0.3;

	tou=0.6;
	while (infile3 >> hor >> ver >> tou >> enedat >> tdat >> vxdat >> vydat)
	{
		it = int((tou+deltat/2.-tau0)/deltat);
          	ix = int((hor+deltax*maxx/2.+deltax/2.)/deltax);
           	iy = int((ver+deltay*maxy/2.+deltay/2.)/deltay);

		hydrot[it][ix][iy]=tdat;
		hydroe[it][ix][iy]=enedat;
		hydrox[it][ix][iy]=vxdat;
		hydroy[it][ix][iy]=vydat;
	}

	int event=0;
	clock_t startClock = clock();
	//Event Loop
	do {
		//Getting initial x, y and b
		gxy(&x,&y,&b);
		//cout << " X= " << x << " Y= " << y << " B= " << b << "\n";
		double xcre=x;
		double ycre=y;

		//Getting tree info
		int FinId[500];
		for (unsigned a=0; a<1000; a++) {
			P[a][0]=0.;
			P[a][1]=0.;
			P[a][2]=0.;
			P[a][3]=0.;
			P[a][4]=0.;
			P[a][5]=0.;
			P[a][6]=0.;
			P[a][7]=0.;
			P[a][8]=0.;
			P[a][11]=0.;
		}
		//Store information from tree and write remnants
		while (infile1 >> num >> px >> py >> pz >> E >> q >> mom >> qge >> ini >> col >> acol >> ide)
		{
			if (num!=0.000123) {
				P[int(num)][0]=px;
				P[int(num)][1]=py;
				P[int(num)][2]=pz;
				P[int(num)][3]=E;
				P[int(num)][4]=q;	//Virtuality
				P[int(num)][5]=mom;
				P[int(num)][6]=qge;	//Quark or Gluon
				P[int(num)][7]=ini;
				P[int(num)][9]=col;
				P[int(num)][10]=acol;
				P[int(num)][11]=ide;
				//Remnants
				if (int(q+0.9)==0 && E!=0. && mom==0.) {
					check << px << " " << py << " " << pz << " " << E  << " " << px << " " << py << " " << pz << " " << E << " " << " 0. " << " " 
					<< int(qge) << " " << int(ini) << " " << int(col) << " " << int(acol) << " " << int(ide) << "\n";
				}
			}
			else {
				event+=1;
				break;
			}
		}
		int k=0;
		//Assign final status
		for (unsigned a=4; a<1000; a++) {
			if (P[int(a)][5]!=0.) {
				int canfin=1;
				for (unsigned b=4; b<1000; b++) {
					if (int(P[int(b)][5])==int(a)) {
						canfin=0;
						break;
					}
				}
				if (canfin==1) FinId[k]=int(a), k+=1;
			}
		}

		int fam[30];
		double  Q[1000][13];
		//Energy Loss Loop
		for (int cr=0; cr<k; cr++) {
			int ind=FinId[cr];
			int fr=0;
			double frac;
			double acum;
			int mymomisko=0;
			fam[fr]=ind;
			for (unsigned a=0; a<50; a++) {
				int mom=int(P[ind][5]);
				if (P[mom][8]==1. && mom!=3) {
					frac=Q[mom][3]/P[mom][3];
					Q[ind][0]=Q[mom][9];		//Px
					Q[ind][1]=Q[mom][10];		//Py
					Q[ind][2]=Q[mom][11];		//Pz
					Q[ind][3]=P[ind][3]*frac;	//E
					if (Q[mom][3]==0.) mymomisko=1;
					x=Q[mom][4];		// X
					y=Q[mom][5];		// Y
					z=Q[mom][6];		// Z
					acum=Q[mom][7];		// t
					break;
				}
				if (mom==3) {
					Q[ind][0]=P[ind][0];
					Q[ind][1]=P[ind][1];
					Q[ind][2]=P[ind][2];
					Q[ind][3]=P[ind][3];
					x=xcre;
					y=ycre;
					z=0.;
					acum=0.;
					break;
				}
				if (P[mom][8]==0.) {
					fr+=1;
					fam[fr]=mom;
					ind=mom;
				}
			}
			for (int w=fr; w>-1; w--) {
				if (mymomisko==1) {
					for (int ku=w; ku>-1; ku--) {
                                                Q[fam[ku]][3]=0.;
                                                P[fam[ku]][8]=1.;
                                        }
                                        break;
				}
				pexe=Q[fam[w]][0];
				peye=Q[fam[w]][1];
				peze=Q[fam[w]][2];
				double we=Q[fam[w]][3];
				double previous=acum;
				double tauf=0.2*2.*P[fam[w]][3]/pow(P[fam[w]][4],2.);
				if (fam[w]==FinId[cr]) tauf=10000.;
				acum=previous+tauf;
				double phi=atan2(peye,pexe);	//Obsolete
				double dE;
				if (abs(P[fam[w]][11])<=6. || P[fam[w]][11]==21.) {
					dE=gdE(we,P[fam[w]][6],P[fam[w]][4],acum-previous,previous,&x,&y,&z,phi,b,&pexe,&peye,&peze);
					Q[fam[w]][0]=pexe;
                                        Q[fam[w]][1]=peye;
                                        Q[fam[w]][2]=peze;
					if (we-dE<=0.)  {
                                        	for (int ko=w; ko>-1; ko--) {
                                        	        Q[fam[ko]][3]=0.;
                                        	        P[fam[ko]][8]=1.;
                                        	}
                                        	break;
                                	}
				}
				else {
					dE=0.;
					x+=Q[fam[w]][0]/we*(acum-previous);
					y+=Q[fam[w]][1]/we*(acum-previous);
					z+=Q[fam[w]][2]/we*(acum-previous);
				}

				Q[fam[w]][3]=we-dE;
				Q[fam[w]][4]=x;
				Q[fam[w]][5]=y;
				Q[fam[w]][6]=z;
				Q[fam[w]][7]=acum;
				P[fam[w]][8]=1.;		//Parton done
				//Propagate Quench&Broad
				if (fam[w]!=FinId[cr]) {
					//cout << " NOTFINAL-> EnerPYTHIA= " << P[fam[w]][3] << " EnerQUEN= " << Q[fam[w]][3] << "\n";
					double lambda=Q[fam[w]][3]/P[fam[w]][3];
					double qux=Q[fam[w]][0];
					double quy=Q[fam[w]][1];
					double quz=Q[fam[w]][2];
					double qmod=sqrt(pow(qux,2.)+pow(quy,2.)+pow(quz,2.));
					if (qmod!=qmod) cout << " QMOD= " << qmod << "\n ";
					double ux=P[fam[w]][0];
					double uy=P[fam[w]][1];
					double uz=P[fam[w]][2];
					double modmom=sqrt(pow(ux,2.)+pow(uy,2.)+pow(uz,2.));
					if (modmom!=modmom) cout << " MODMOM= " << modmom << "\n ";
					qux/=qmod;
					quy/=qmod;
					quz/=qmod;
					ux/=modmom;
					uy/=modmom;
					uz/=modmom;
					if (qmod!=qmod || qmod==0.) {
						qux=0.;
						quy=0.;
						quz=0.;
						qmod=0.;
					}
					double xaxis=uy*quz-uz*quy;
					double yaxis=uz*qux-ux*quz;
					double zaxis=ux*quy-uy*qux;
					double modaxis=sqrt(pow(xaxis,2.)+pow(yaxis,2.)+pow(zaxis,2.));
					if (modaxis!=0.) {
						xaxis/=modaxis;
						yaxis/=modaxis;
						zaxis/=modaxis;
					}
					double angle=0.;
					if (ux*qux+uy*quy+uz*quz>=1.) angle=0.;
					else angle=acos(ux*qux+uy*quy+uz*quz);
					if (angle!=angle) cout << " Cos= " << ux*qux+uy*quy+uz*quz <<  " Acos(cos)= " << acos(ux*qux+uy*quy+uz*quz) << "\n";
					//Who's the other son?
					int oson=1002;
					for (unsigned a=0; a<1000; a++) {
						if (P[a][5]==fam[w] && int(a)!=fam[w-1]) oson=int(a);
					}
					if (oson==1002) cout << " DIDN'T FIND OTHER SON!!! \n";
					double modoson=sqrt(pow(P[oson][0],2.)+pow(P[oson][1],2.)+pow(P[oson][2],2.));
					double modthis=sqrt(pow(P[fam[w-1]][0],2.)+pow(P[fam[w-1]][1],2.)+pow(P[fam[w-1]][2],2.));
					double thisx=P[fam[w-1]][0]/modthis;
					double thisy=P[fam[w-1]][1]/modthis;
					double thisz=P[fam[w-1]][2]/modthis;
					double otx=P[oson][0]/modoson;
					double oty=P[oson][1]/modoson;
					double otz=P[oson][2]/modoson;
					double thisscal=xaxis*thisx+yaxis*thisy+zaxis*thisz;
					double usex=cos(angle)*thisx+sin(angle)*(yaxis*thisz-zaxis*thisy)+(1.-cos(angle))*xaxis*thisscal;
					double usey=cos(angle)*thisy+sin(angle)*(zaxis*thisx-xaxis*thisz)+(1.-cos(angle))*yaxis*thisscal;
					double usez=cos(angle)*thisz+sin(angle)*(xaxis*thisy-yaxis*thisx)+(1.-cos(angle))*zaxis*thisscal;
					double oscal=xaxis*otx+yaxis*oty+zaxis*otz;
					double xoson=cos(angle)*otx+sin(angle)*(yaxis*otz-zaxis*oty)+(1.-cos(angle))*xaxis*oscal;
					double yoson=cos(angle)*oty+sin(angle)*(zaxis*otx-xaxis*otz)+(1.-cos(angle))*yaxis*oscal;
					double zoson=cos(angle)*otz+sin(angle)*(xaxis*oty-yaxis*otx)+(1.-cos(angle))*zaxis*oscal;
					double lamp=qmod/modmom;
					//cout << " ModThis Bef= " << modthis << " ModThis Aft= " << sqrt(pow(usex,2.)+pow(usey,2.)+pow(usez,2.)) << "\n";
					//cout << " ModOson Bef= " << modoson << " ModOson Aft= " << sqrt(pow(xoson,2.)+pow(yoson,2.)+pow(zoson,2.)) << "\n";
					usex*=modthis*lamp;
					usey*=modthis*lamp;
					usez*=modthis*lamp;
					xoson*=modoson*lamp;
					yoson*=modoson*lamp;
					zoson*=modoson*lamp;
					//cout << " Angle= " << angle << " ";
					//cout << " MomCons Px= " << Q[fam[w]][0]-usex-xoson << " Py= " << Q[fam[w]][1]-usey-yoson << " Pz= " << Q[fam[w]][2]-usez-zoson << " E= " << Q[fam[w]][3] -P[fam[w-1]][3]*lambda-P[oson][3]*lambda << "\n";
					//cout << " This Son Q2 Bef= " << pow(P[fam[w-1]][3],2.)-pow(modthis,2.) << "\n";
					//cout << " Other Son Q2 Bef= " << pow(P[oson][3],2.)-pow(modoson,2.) << "\n";
					Q[fam[w-1]][0]=usex;
                                        Q[fam[w-1]][1]=usey;
                                        Q[fam[w-1]][2]=usez;
					Q[fam[w-1]][3]=P[fam[w-1]][3]*lambda;
					//cout << " This Son Q2 Aft= " << pow(Q[fam[w-1]][3],2.)-pow(usex,2.)-pow(usey,2.)-pow(usez,2.) << "\n";
                                        //cout << " Other Son Q2 Aft= " << pow(P[oson][3]*lambda,2.)-pow(xoson,2.)-pow(yoson,2.)-pow(zoson,2.) << "\n";
					Q[fam[w-1]][4]=x;
					Q[fam[w-1]][5]=y;
					Q[fam[w-1]][6]=z;
					Q[fam[w]][9]=xoson;
					Q[fam[w]][10]=yoson;
					Q[fam[w]][11]=zoson;
				}
			}
			check << P[fam[0]][0] << " " << P[fam[0]][1] << " " << P[fam[0]][2] << " " << P[fam[0]][3]  << " " 
			<< Q[fam[0]][0] << " " << Q[fam[0]][1] << " " << Q[fam[0]][2] << " " << Q[fam[0]][3] << " " 
			<< 0. << " " << int(P[fam[0]][6]) << " " << int(P[fam[0]][7]) << " " << int(P[fam[0]][9]) 
			<< " " << int(P[fam[0]][10]) << " " << int(P[fam[0]][11]) << "\n";
			//cout << " EnerPYTHIA= " << P[fam[0]][3] << " EnerQUEN= " << Q[fam[0]][3] << "\n";
		}

		count +=1;
		//cout << " NEXT EVENT  \n";
		int mult=int(double(count)/(ene*100.));
		if (mult>0) {
			cout << " " << count << "\n";
			ene+=5.;
		}

		if (k!=0) check << "0.000123" << " " << " 0. " << " " << " 0. " << " " << xcre << " " << ycre << " " << b << " " << " 0. "  << " " << " 0. " << " " << count  << " " << " 0. " << " " << " 0. " << " " << " 0. " << " " << " 0. " << " " <<  " 0. " << "\n";
	//End Event Loop
	} while (count < N);
	clock_t endClock = clock();
	cout << " Time= ";
	printf("%ld", (endClock - startClock) / CLOCKS_PER_SEC);
	cout  << " #FuckedN= " << fuckedN/totalN << " ";
	cout << " #Events= " << event << " ";

	check.close();
//End Program
}
//------------------gdE-----------------------------//
double gdE(double E, double M, double gamma, double t, double t0, double *x, double *y, double *z, double phi, double p, double *pexe, double *peye, double *peze) {

	double ti, tih, tot;
	double proper, properh;
	double xp, yp, zp;
	double xph, yph, zph;
	double Temp, Temph;
	double Ev, CF;
	double Tc;
	if (tmethod==0) Tc=170./1000.;
	if (tmethod==1) Tc=145./1000.;
	double stepe=0.01;
	double wx, wy, wz;
	double vx, vy, vz;
	double lore, vscalw, v2;
	int marker=0;
	double i=0.;
	double quench=0.;
	double gede=0.;
	double dist=0.;
	double distlab=0.;
	if (int(M)==1) CF=1.;
	if (int(M)==2) {
		if (qmethod==3) CF=pow(9./4.,1./3.);
		else CF=pow(9./4.,1.);
	}
	Ev=E;
	tot=t0+t;
	int noproblem=1;
	int saidbef=0, saidaft=0;
	//cout << "\n New Rung--- Ei= " << E << "\n";
	//cout << " gamma= " << gamma << "\n";
	//cout << " xini= " << *x << " yini= "<< *y << " zini= " << *z << " tauini= " << sqrt(t0*t0-pow(*z,2.))<< "\n";
	//cout << " Initial Q2= " << E*E-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.) << "\n";
	do {
		ti=t0+i*stepe;
		tih=ti+stepe;
		if (tih>tot) {
			marker=1;
			tih=tot;
			stepe=tih-ti;
		}
		wx=*pexe/Ev;
		wy=*peye/Ev;
		wz=*peze/Ev;

		xp=*x;
		yp=*y;
		zp=*z;
		proper=sqrt(ti*ti-zp*zp);

		xph=*x+wx*stepe;
                yph=*y+wy*stepe;
                zph=*z+wz*stepe;
                properh=sqrt(tih*tih-zph*zph);

                *x=xph;
                *y=yph;
                *z=zph;

                if (proper!=proper) proper=0.;
                if (properh!=properh) properh=0.;
		if (proper>=tquench && properh>=tquench) {	//Need both proper and properh > tquench !!!
			//extract fluid velocity with LAB coordinates
			vx=gVx(proper,xp,yp);
			vy=gVy(proper,xp,yp);
			double frap=0.;
			vz=*z/tih;
			frap=atanh(vz);
			vx/=cosh(frap);
			vy/=cosh(frap);
			v2=pow(vx,2.)+pow(vy,2.)+pow(vz,2.);
			vscalw=vx*wx+vy*wy+vz*wz;
			if (v2>=1.) {
				v2=0.9999999;
				cout << " V2 > 1 WTFFF !!! \n";
			}
			lore=1./sqrt(1.-v2);
			double infi=stepe*sqrt(wx*wx+wy*wy+wz*wz+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
			if (infi!=infi) infi=0.;
			dist+=infi;
			distlab+=stepe;
			Temp=gT(proper,xp,yp);
			Temph=gT(properh,xph,yph);
			if (Temph<Tc) {
				double cospart=xph*wx+yph*wy;
				if (cospart>0.) {
					Temp=0.;
					Temph=0.;
					*x=xph+wx*(tot-tih);
					*y=yph+wy*(tot-tih);
					*z=zph+wz*(tot-tih);
				}
			}
			if (Temph==0.) marker=1;
			if (Ev<=0.) marker=1;
		}
		if (Ev>0.) {
			if (proper>=tquench && properh>=tquench && (Temp+Temph)/2.>Tc+0.00000001) {
				//Broadening
				double e1x=wy*vz-wz*vy;
				double e1y=wz*vx-wx*vz;
				double e1z=wx*vy-wy*vx;
				double Ne1=sqrt(pow(e1x,2.)+pow(e1y,2.)+pow(e1z,2.));
				if (Ne1!=Ne1) cout << " PROBLEMMMM WITH NE1 ";
				e1x=e1x/Ne1;
				e1y=e1y/Ne1;
				e1z=e1z/Ne1;

				double Nw=sqrt(pow(wx,2.)+pow(wy,2.)+pow(wz,2.));
				if (Nw!=Nw) cout << " PROBLEMMMM WITH NW ";
				double lx=(wy*e1z-wz*e1y)/Nw;
				double ly=(wz*e1x-wx*e1z)/Nw;
				double lz=(wx*e1y-wy*e1x)/Nw;

				//cout << " Modul L= " << lx*lx+ly*ly+lz*lz << " \n";
				double uscalW=lore*(1.-vx*wx-vy*wy-vz*wz);
				double uscall=lore*(-vx*lx-vy*ly-vz*lz);
				double W2=1.-wx*wx-wy*wy-wz*wz;
				//if (wx>1. || wy>1. || wz>1.) cout << " at i= " << i << " W2 loker= " << W2 << " with Ev= " << Ev << " wx= " << wx << " wy= " << wy << " wz= " << wz << "\n";
				double Wpt=1.-lore*W2/uscalW;
				double Wpx=wx-vx*lore*W2/uscalW;
				double Wpy=wy-vy*lore*W2/uscalW;
				double Wpz=wz-vz*lore*W2/uscalW;

				double NN=1.+W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2);
				double Nalpha=-uscall*uscalW/(pow(uscalW,2.)-W2);
				if (v2>1.-W2) {
					//cout << " vfluid= " << v2 << " vparton= " << 1.-W2 << " ";
					//cout << " nume/deno= " << W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2) << " ";
					//cout << " W2= " << W2 << "\n";
				}
				if (sqrt(NN)!=sqrt(NN)) {
					noproblem=0;
					cout << " Cos= " << (vx*wx+vy*wy+wz*vz)/sqrt(v2)/sqrt(1.-W2) << " ";
					//cout <<  " deno N= " << -pow(uscalW,2.)+W2 << " Ev= " << Ev << "\n";
					//cout << " nume/deno= " << W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2) << " ";
					//cout << " vfluid= " << v2 << " vparton= " << 1.-W2 << " W2= " << W2 << " lore= " << lore << "\n";
				}
				else {
					double e2t=Nalpha*Wpt/sqrt(NN);
					double e2x=(lx+Nalpha*Wpx)/sqrt(NN);
					double e2y=(ly+Nalpha*Wpy)/sqrt(NN);
					double e2z=(lz+Nalpha*Wpz)/sqrt(NN);
					if (sqrt(NN)!=sqrt(NN)) {
						cout << " e2t= " << e2t << " e2x= " << e2x << " e2y= " << e2y << " e2z= " << e2z << endl;
					}

					double Ef=Ev*lore*(1.-vscalw);
					double wf2=1.-W2/pow(lore*(1.-vscalw),2.);
					double DelQ2=kappa*pow((Temp+Temph)/2.,3.)*lore*(1.-vscalw)*stepe*5.;
					double qfac=0.;
					int bucle=0;
					if (Ef>(Temp+Temph)/2.) {
						do {
							bucle+=1;
							if (kappa!=0.) qfac=gQ(DelQ2);
							if (bucle>1000) cout << " Ef= " << Ef << " ";
						} while (qfac>Ef*sqrt(wf2) && bucle<1000);
					}
					double qbeta=sqrt(1.-qfac*qfac/Ef/Ef/wf2)-1.;
					if (qbeta!=qbeta) {
						qbeta=-1.;
						qfac=Ef*sqrt(wf2);
						if (qfac!=qfac) cout << " QFAC LOKKKKKKERRRR ";
					}
					double qphi=2.*3.141592654*rando(&ir);

					double et=sin(qphi)*e2t;
					double ex=cos(qphi)*e1x+sin(qphi)*e2x;
					double ey=cos(qphi)*e1y+sin(qphi)*e2y;
					double ez=cos(qphi)*e1z+sin(qphi)*e2z;

					double Wtt=(1.-(uscalW)*lore)/lore/(1.-vscalw);
					double Wtx=(wx-(uscalW)*vx*lore)/lore/(1.-vscalw);
					double Wty=(wy-(uscalW)*vy*lore)/lore/(1.-vscalw);
					double Wtz=(wz-(uscalW)*vz*lore)/lore/(1.-vscalw);

					if (W2<0. && saidbef==0) {
						//cout << " Q² bef= " << Ev*Ev-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.);
						saidbef=1;
					}
					//cout << " Echange= " << qbeta*Ef*Wpt+qfac*et << " ";
					if (qbeta*Ef*Wtt+qfac*et!=qbeta*Ef*Wtt+qfac*et) cout << " Ener WRONG ";
					if (qbeta*Ef*Wtx+qfac*ex!=qbeta*Ef*Wtx+qfac*ex) cout << " Px WRONG";
					if (qbeta*Ef*Wty+qfac*ey!=qbeta*Ef*Wty+qfac*ey) cout << " Py WRONG";
					if (qbeta*Ef*Wtz+qfac*ez!=qbeta*Ef*Wtz+qfac*ez) cout << " Pz WRONG";
					Ev+=qbeta*Ef*Wtt+qfac*et;
					if (Ev<=0.) {
						//cout << " at i= " << i << " Ev= " << Ev << " ";
						Ev=0.;
						marker=1;
					}
					*pexe+=qbeta*Ef*Wtx+qfac*ex;
					*peye+=qbeta*Ef*Wty+qfac*ey;
					*peze+=qbeta*Ef*Wtz+qfac*ez;
					if (W2<0. && saidaft==0) {
						//cout << " Q² aft= " << Ev*Ev-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.) << "\n";
						saidaft=1;
					}
				}
				//Coll
				if (qmethod==1 && Ev>0.) {
					double intpiece=alpha*CF*stepe*5.*pow((Temp+Temph)/2.,2.);
					quench=(Ev-intpiece)/Ev;
					Ev-=intpiece;
					if (Ev<0.) Ev=0.;
				}
				//Radiative
				if (qmethod==2 && Ev>0.) {
					double intpiece=5.*alpha*CF*stepe*dist*5.*pow((Temp+Temph)/2.,3.);
					quench=(Ev-intpiece)/Ev;
					Ev-=intpiece;
					if (Ev<0.) Ev=0.;
				}
				//Strong
				if (qmethod==3 && Ev>0.) {
					double Efs=E*lore*(1.-vscalw);
					double tstop=0.2*pow(Efs,1./3.)/(2.*pow((Temp+Temph)/2.,4./3.)*alpha)/CF;
					double beta=tstop/dist;
					if (beta>1.) {
						double intpiece=Efs*stepe*4./(3.141592)*(1./(beta*tstop*sqrt(beta*beta-1.)));
						quench=(Ev-intpiece)/Ev;
						Ev-=intpiece;
						if (intpiece<0.) cout << " WTF INTPIECE= " << intpiece << " ";
						if (quench>1.) cout << " QuenchFactor= " << quench << "  and IntPiece= " << intpiece << " and Ev= " << Ev << " ";
					}
					else {
						Ev=0.;
						quench=0.;
						marker=1;
					}
					if (Ev<0.) {
						Ev=0.;
						quench=0.;
						marker=1;
					}
				}
				if (quench>1.) cout << " QuenchFactor= " << quench << " ";
				if (Ev<=0.) quench=0.;
				*pexe*=quench;
				*peye*=quench;
				*peze*=quench;
			}
		}
		else {
			gede=E;
			Ev=0.;
			marker=1;
		}
		i+=1.;
	} while (marker!=1);
	//cout << " Final Q2= " << Ev*Ev-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.) << "\n";
	//cout << " Final 4mom= " << Ev << " " << *pexe << " " << *peye << " " << *peze << " quench= " << quench << "\n";
	gede=E-Ev;
	if (noproblem!=1) fuckedN+=1.;
	totalN+=1.;
	//cout << " xfin= " << *x << " yfin= " << *y << " zfin= " << *z <<" Temp fin= " << (Temp+Temph)/2. << "\n";
	//cout << " DistLab= " << distlab << " DistFluid= " << dist << "\n";
	//cout << " dE= " << gede << " ";
	return gede;
}
//------------------getGrid-------------------------//
double getGrid(double tau, double x, double y, double *dt, double *dx, double *dy, int *it, int *ix, int *iy) {

	*it=int((tau-tau0)/deltat);
	*dt=(tau-tau0-double(*it)*deltat)/deltat;

	if (y>0.) {
		*iy = int(y/deltay)+maxy/2;
		*dy = (y - double(*iy-maxy/2.)*deltay)/deltay;
	}
	else {
		*iy = int(y/deltay)+maxy/2-1;
		*dy = (y - double(*iy-maxy/2.)*deltay)/deltay;
	}

	if (x>0.) {
		*ix = int(x/deltax)+maxx/2;
		*dx = (x - double(*ix-maxx/2.)*deltax)/deltax;
	}
	else {
		*ix = int(x/deltax)+maxx/2-1;
		*dx = (x - double(*ix-maxx/2.)*deltax)/deltax;
	}

	return 0;
}

//------------------gT------------------------------//
double gT(double tau, double x, double y) {

	double gete=0.;
	double tau1=18.5;

	if (tau>tau1) {
		//cout << " No data for tau > " << tau1 << " fm ";
		return gete;
	}

	getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
	if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gete;
	gete=hydrot[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
	gete+=hydrot[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
	gete+=hydrot[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
	gete+=hydrot[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
	gete+=hydrot[it+1][ix+1][iy]*dx*(1.-dy)*dt;
	gete+=hydrot[it][ix+1][iy+1]*(1.-dt)*dy*dx;
	gete+=hydrot[it+1][ix][iy+1]*(1.-dx)*dt*dy;
	gete+=hydrot[it+1][ix+1][iy+1]*dx*dt*dy;

	return gete*0.2;
}

//------------------gE------------------------------//
double gE(double tau, double x, double y) {

	double gener=0.;
	double tau1=18.5;

	if (tau>tau1) {
		//cout << " No data for tau > " << tau1 << " fm ";
		return gener;
	}

	getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
	if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gener;
	gener=hydroe[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
	gener+=hydroe[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
	gener+=hydroe[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
	gener+=hydroe[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
	gener+=hydroe[it+1][ix+1][iy]*dx*(1.-dy)*dt;
	gener+=hydroe[it][ix+1][iy+1]*(1.-dt)*dy*dx;
	gener+=hydroe[it+1][ix][iy+1]*(1.-dx)*dt*dy;
	gener+=hydroe[it+1][ix+1][iy+1]*dx*dt*dy;

	return gener;
}

//------------------gVx------------------------------//
double gVx(double tau, double x, double y) {

        double gvelx=0.;
        double tau1=18.5;

        if (tau>tau1) {
                //cout << " No data for tau > " << tau1 << " fm ";
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

//------------------gVy------------------------------//
double gVy(double tau, double x, double y) {

        double gvely=0.;
        double tau1=18.5;

        if (tau>tau1) {
                //cout << " No data for tau > " << tau1 << " fm ";
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

//------------------gQ------------------------------//
double gQ(double Del) {

	double gaussq;
	double qfac;
	double gaussmax=sqrt(2./Del/exp(1.));

	qhatelsen:
	qfac=5.*sqrt(Del/2.)*rando(&ir);
	gaussq=2.*qfac/Del*exp(-qfac*qfac/Del)/gaussmax;
	double nrand=rando(&ir);
	if (nrand>gaussq) goto qhatelsen;

return qfac;
}

//------------------gxy-----------------------------//
double gxy(double *x, double *y, double *b) {

	double rho,phi;
	double P;

        naiguels:
	*b=sqrt((bmax*bmax-bmin*bmin)*rando(&ir)+bmin*bmin);
        norm=1.;
        norm=gTAA(0.,0.,bmin);
	rho=sqrt(150.*rando(&ir));
	phi=2.*3.141592654*rando(&ir);
	*x=rho*cos(phi);
	*y=rho*sin(phi);
	P=rando(&ir);
	if(P>gTAA(*x,*y,*b)) goto naiguels;
	return 0;
}

//-----------------gTAA------------------------------//
double gTAA(double x, double y, double b)
{
	int il, irr;
	double rho2, use;

	rho2=pow(x+b/2.,2.)+y*y;
	il=int(rho2/step);
	rho2=pow(x-b/2.,2.)+y*y;
	irr=int(rho2/step);
	use=0.;
	if(il<4000 && irr<4000) {
		use=TA[il]*TA[irr]/norm;
	}
	return use;
}

//--------------Random Number-----------------------//
double rando (int *ir)
{
	double da=16807.;
	double db=2147483647.;
	double dc=2147483648.;
	double usran;

	*ir=int(abs(fmod(da*(*ir),db)+0.5));
	usran=double(*ir)/dc;
	return usran;
}
