/************************************************************************
GreensTD - time-dependent Greens function method
TWS June 2012
Based on GreensV3
Variable array sizes and bounds, using Numerical Recipes utilities.
No nondimensionalization.  Lengths, diameters in microns, times in s.
Flows in nanoliters/min. Diffusivities in cm2/s, converted to micron^2/s

Main variables:
mat --- matrix for vessel strengths
gamma1 --- intravascular resistance, varies with vessel diameter
alphaconv, betaconv,gammaconv,omegaconv,zetaconv,xiconv
	--- matrices giving dependence of vessel levels on b.c.s, previous values and source strengths
lseg --- segment length
ds --- subsegment length
qv --- subsegment source strength
cv --- concentration in the subsegment
qt --- tissue source strength
ct --- concentration in the tissue
q --- flow rate, qq = abs(q)
Special parameters for oxygen
p50, fn --- parameters in the Hill equation
cs --- red blood cell oxygen binding capacity in cm^3 O2/cm^3
Updated October 2016 for impermeable solutes (alternative method for G0)
**************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void putrank(void);
void initgreensTD();
void convectmatrix();
void diffusionmatrix();
void blood(float c, float hem, float *p, float *pp);
float bloodconc(float p,float h);
void tissrate(int nsp, float *p, float *mtiss, float *mptiss);
float erfcc(float x);
double bicgstab(double **a, double *b, double *x, int n, double eps, int itmax);
void contour(int itime, float time, int method);

double bicgstabBLASD(double **a, double *b, double *x, int n, double eps, int itmax);
float bicgstabBLASSt(float *b, float *x, int n, float eps, int itmax,float diff);
void tissueGPU1c(int isp);
void tissueGPU2c();
void tissueGPU3c();
void tissueGPUcopy();

void greensTD(int irun)
{
	extern int nmax,ntpts,ninterval;
	extern int mxx,myy,mzz,nnt,nnv,nseg,nsp,nnod,nnodbc;
	extern int contourmethod;
	extern int *mainseg,*permsolute,*segtyp;
	extern int *segname,*nspoint,*istart,*nodout,*bcnod,*nodtyp,*bcnodname,*nodname;
	extern int *nresis,*oxygen,*diffsolute;
	extern int *intervalnst;
	extern int **nodseg;
	extern int ***nbou,**tisspoints;
	extern int nntGPU,nnvGPU,useGPU;//needed for GPU version

	extern float p50,cs,req,fac,flowfac,errfac,totalq,N_D,N_P;
	extern float v,vol,vdom,tlength,tarea,pi1,alx,aly,alz;
	extern float *axt,*ayt,*azt,*ds,*diff,*alphab,*alphat,*cmin,*cmax,*cmean;
	extern float *crefb,*creft,*cinitb,*cinitt,*lambda;
	extern float *diam,*rseg,*volseg,*q,*qq,*hd,*bchd,*qvfac;
	extern float *x,*lseg,*mtiss,*mptiss;
	extern float *epsvesselq,*epstissueq,*epsvesselc,*epstissuec,*p,*pv,*dcdp;
	extern float *g0,*g0prev,*g0previt,*ctct,*qtsum,*qtpsum,*qvsum;
	extern float *intervaldt,*tpts,deltat,*cumulerror,*cevlhs,*cevrhs;
	extern float **start,**end,**scos,**ax,**bcp,**bcpbase;
	extern float **tissparam,**resisdiam,**resis;
	extern float **cv,**ct,**cvprev,**ctprev,**cvprevit,**ctprevit;
	extern float **qv,**qt,**qtp,**qvprevit,**qtprevit,**qvseg,**cvseg,**gamma1;
	extern float **gtt,**gbartt,**tconcs,**gbarvtsum,**gbarttsum,**gttsum,*gttsum1,*gbarttsum1;
	extern float **alphaconv,**betaconv,**gammaconv,**omegaconv,**zetaconv,**xiconv;
	extern float ***gbarvv,***gvt,***gbarvt;
	extern double **mat,*rhs,*matx;
	extern float *h_ct,*h_ctprev,*h_qt,*h_qtp,*h_pv,*h_qv,*h_gtt,*h_gbartt;	//needed for GPU version

	int i,j,k,ix,iy,iz,jx,jy,jz,iseg,inod,nt,ineg,ihigh,isp,imaxerr,inodbc,jxyz,iinterval,isubint;
	int isn,itime,kmain,itp,jtp,convflag;
	int greensverbose = 1;	//set to 0 for terse output, 1 for verbose
	int bicgstabit = 2000; //parameter for biconjugate gradient method.  March 2010

	float duration,time,pathlength,meanflow,rcbar,lambda1,den;
	float dif,err,matfac,rhs2,perm,cvmin,cvmax;
	float qtqvsum,sumt,sumtprev,sumv,flowsumin,flowsumout,concin,concout;
	float gbarttratesum,g0fac,ctsum,ctprevsum;	//added October 2016 for non-permeable solutes
	double bicgstaberr = 1.e-12; //parameter for biconjugate gradient method

	FILE *ofp,*ofp1,*ofp2,*ofp3;
	clock_t tstart,tfinish,tstart1,tfinish1;
	char varname[3],fname1[80],fname2[80];

//use sequentially numbered files
	sprintf(fname1,"GreensWashout%03i.txt",irun);
	sprintf(fname2,"GreensLog%03i.txt",irun);


	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5)	//setup mainseg (must be done after setuparrays2)	
		for(i=0; i<nspoint[iseg]; i++) mainseg[istart[iseg] + i] = iseg;
	for(i=1; i<=nnv; i++){	//identify vessel points
		iseg = mainseg[i];
		isn = i - istart[iseg];
		for(j=1; j<=3; j++)	ax[j][i] = start[j][iseg] + scos[j][iseg]*ds[iseg]*(isn+0.5);
		volseg[i] = pi1*SQR(rseg[iseg])*ds[iseg];
	}
	for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){	//index tissue points
		nt = nbou[i][j][k];
		if(nt > 0){
			tisspoints[1][nt] = i;
			tisspoints[2][nt] = j;
			tisspoints[3][nt] = k;
		}
	}
	vdom = nnt*vol;
	den = sqrt(vdom/tlength);
	if(greensverbose){
		printf("Tissue node points used in simulation = %i\n", nnt);
		printf("V of simulation region/V of cubic region = %f\n", nnt/(1.0*mxx*myy*mzz));
		printf("Total vessel length = %g micron\n",tlength);
		printf("Total vessel area = %g micron^2\n",tarea);
		rcbar = tarea/tlength/2./pi1;
		printf("Equivalent vessel radius = %f micron\n", rcbar);
		printf("Tissue volume = %g micron^3\n",vdom);
		printf("Sqrt(Tissue Volume/vessel length) = %f\n", den);
		printf("Total inflow to network = %f nl/min\n", totalq);
		printf("Capillary density = %f /mm2\n", tlength/vdom*1.e6);
		printf("Perfusion = %f cm3/cm3/min (uncorrected for path length effect)\n", totalq/vdom*1.e6);
		pathlength = 0.;
		meanflow = 0.;
		for(iseg=1; iseg<=nseg; iseg++){
			pathlength += fabs(q[iseg])*lseg[iseg]/totalq;
			meanflow += fabs(q[iseg])*lseg[iseg]/tlength;
		}
		printf("Mean flow-weighted path length = %f micron\n", pathlength);
		printf("Mean length-weighted flow rate = %f nl/min\n", meanflow);
		printf("Note: Perfusion = (mean flow)*(capillary density)/(mean path length)\n");
		for(isp=1; isp<=nsp; isp++){
			if(diffsolute[isp]){
			N_D = 2.*pi1*pathlength*diff[isp]/meanflow/flowfac;
				printf("Solute %i: dimensionless diffusivity N_D = %f\n",isp,N_D);
			}
			if(permsolute[isp]){
			N_P = 2.*pi1*rcbar*pathlength/resis[1][isp]/meanflow/flowfac;
				printf("Solute %i: dimensionless permeability N_P = %f\n",isp,N_P);
			}			
		}
	}
//Calculate intravascular or wall transport resistance.  Zero unless specified in IntravascRes.dat.
//If not oxygen, assume value from data is 1/(wall permeability in um/s)
	for(isp=1; isp<=nsp; isp++)	for(iseg=1; iseg<=nseg; iseg++)	gamma1[iseg][isp] = 0.;
	for(isp=1; isp<=nsp; isp++)	if(nresis[isp] != 0) for(iseg=1; iseg<=nseg; iseg++){
		gamma1[iseg][isp] = resis[1][isp];
		for(j=2; j<=nresis[isp]; j++) if(diam[iseg] <= resisdiam[j][isp] && diam[iseg] > resisdiam[j-1][isp])
			gamma1[iseg][isp] = resis[j-1][isp] + (resis[j][isp]-resis[j-1][isp])
			*(diam[iseg]-resisdiam[j-1][isp])/(resisdiam[j][isp]-resisdiam[j-1][isp]);
		if(diam[iseg] > resisdiam[nresis[isp]][isp]) gamma1[iseg][isp] = resis[nresis[isp]][isp];
		if(oxygen[isp] != 1) gamma1[iseg][isp] = gamma1[iseg][isp]/pi1/diam[iseg];
	}

	initgreensTD();
	putrank();

	ofp1 = fopen(fname2, "w");	//create log file
	fprintf(ofp1,"GreensLog.txt\n");
	perm = 0.;
	if(resis[1][1] != 0.) perm = 1./resis[1][1];
	fprintf(ofp1,"%f %f %f delta_t diff perm\n",intervaldt[1],diff[1],perm);
	fclose(ofp1);
	ofp2 = fopen("GreensConc.txt", "w");	//create vessel conc. file
	fprintf(ofp2,"GreensConc.txt\n");
	fclose(ofp2);
	ofp3 = fopen(fname1, "w");	//create washout concentration file
	fprintf(ofp3,"itime   time      ");
	for(isp=1; isp<=nsp; isp++){
		fprintf(ofp3," cumulerror%i ",isp);
		if(permsolute[isp]) fprintf(ofp3,"   concin%i      concout%i  ",isp,isp);
		if(diffsolute[isp]) fprintf(ofp3,"   ambient%i  ",isp);
	}
	fprintf(ofp3,"\n");
	fprintf(ofp3,"  0     0.0     ");
	for(isp=1; isp<=nsp; isp++){
		fprintf(ofp3,"     0.0     ");
		if(permsolute[isp]) fprintf(ofp3,"     0.0          0.0     ");
		if(diffsolute[isp]) fprintf(ofp3,"     0.0     ");
	}
	fprintf(ofp3,"\n");
	fclose(ofp3);	

	for(isp=1; isp<=nsp; isp++){	//data for plot with initial values
		for(iseg=1; iseg<=nseg; iseg++) cvseg[iseg][isp] = 0.;
		for(i=1; i<=nnv; i++){
			iseg = mainseg[i];
			if(oxygen[isp]){
				blood(cv[i][isp],hd[iseg],&pv[i],&dcdp[i]);
				cvseg[iseg][isp] += pv[i]/nspoint[iseg];
			}
			else cvseg[iseg][isp] += cv[i][isp]/nspoint[iseg];
		}
	}
	time = 0.;
	itime = 0;
	contourmethod = 1;
	contour(itime,time,contourmethod);

	for(isp=1; isp<=nsp; isp++)	cumulerror[isp] = 0.;
	flowsumin = 0.;
	flowsumout = 0.;
	for(inodbc=1; inodbc<=nnodbc; inodbc++){
		inod = bcnod[inodbc];
		iseg = nodseg[1][inod];
		if(nodout[inod] == 1) flowsumin += qq[iseg];
		else flowsumout += qq[iseg];
	}
	if(fabs(flowsumin - flowsumout) > 0.001) printf("***Error: input and output flow not equal\n");

	tstart = clock();
	deltat = 0.;
//********************** start of time loop *****************************
	for(iinterval=1; iinterval<=ninterval; iinterval++) for(isubint=1; isubint<=intervalnst[iinterval]; isubint++){
		if((iinterval == 1 && isubint == 1) || fabs(intervaldt[iinterval] - deltat) > 0.001*deltat){	//reset deltat, recompute coefficients
			tstart1 = clock();
			deltat = intervaldt[iinterval];
			diffusionmatrix();
			convectmatrix();
			tfinish1 = clock();
			duration = (float)(tfinish1 - tstart1)/CLOCKS_PER_SEC;
			if(greensverbose) printf("%2.3f seconds for matrix setup\n", duration);
			if(useGPU) tissueGPUcopy();	//this has to be done after gtt and gttbar are set up
		}
		itime++;
		time += deltat;
		if(fabs(time - tpts[itime]) > 1.e-6) printf("*** Error in time points list\n");
		tstart1 = clock();
		if(greensverbose) printf("\n----- itime = %i time = %f -----\n",itime,time);
		else printf(" %i",itime);
		for(isp=1; isp<=nsp; isp++){	//'prev' variables are values from the previous time step 
			for(itp=1; itp<=nnt; itp++)	ctprev[itp][isp] = ct[itp][isp];
			for(i=1; i<=nnv; i++) cvprev[i][isp] = cv[i][isp];
			for(inodbc=1; inodbc<=nnodbc; inodbc++)	bcp[inodbc][isp] = bcpbase[inodbc][isp]*tconcs[itime][isp];
			g0prev[isp] = g0[isp];
		}
//********************** start of main iteration loop *****************************
		for(kmain=1; kmain<nmax; kmain++){
			if(greensverbose) printf("\n----- kmain = %i -----\n",kmain);
			if(kmain == 1){			//initialize intravascular solute levels
				for(isp=1; isp<=nsp; isp++)	if(permsolute[isp]){
					for(i=1; i<=nnv; i++){
						iseg = mainseg[i];
						cv[i][isp] = 0.;
						for(inodbc=1; inodbc<=nnodbc; inodbc++) cv[i][isp] += gammaconv[i][inodbc]*bcp[inodbc][isp];
						for(j=1; j<=nnv; j++) cv[i][isp] += betaconv[i][j]*cvprev[j][isp] - alphaconv[i][j]*qv[j][isp]*deltat/volseg[j];
					}
				}
			}
			for(isp=1; isp<=nsp; isp++){	//'previt' variables are values from the previous iteration 
				for(itp=1; itp<=nnt; itp++){
					ctprevit[itp][isp] = ct[itp][isp];
					qtprevit[itp][isp] = qt[itp][isp];
				}
				if(permsolute[isp]) for(i=1; i<=nnv; i++){
					cvprevit[i][isp] = cv[i][isp];
					qvprevit[i][isp] = qv[i][isp];
				}
				g0previt[isp] = g0[isp];
			}			
//generate linear system to be solved for changes in qv and g0
			for(isp=1; isp<=nsp; isp++)	if(permsolute[isp]){
				for(i=1; i<=nnv; i++){
					iseg = mainseg[i];
					if(oxygen[isp]) blood(cv[i][isp],hd[iseg],&pv[i],&dcdp[i]);
					else{
						pv[i] = cv[i][isp];
						dcdp[i] = 1.;
					}
//right hand side of linear system
					cevlhs[i] = g0[isp];
					for(itp=1; itp<=nnt; itp++) cevlhs[i] += vol*gvt[i][itp][isp]*(ctprev[itp][isp] - g0prev[isp])
						+ gbarvt[i][itp][isp]*(qt[itp][isp] - vol/deltat*(g0[isp] - g0prev[isp]));
					for(j=1; j<=nnv; j++) cevlhs[i] += gbarvv[i][j][isp]*qv[j][isp];
					if(oxygen[isp]) cevrhs[i] = alphat[isp]*pv[i];
					else cevrhs[i] = alphat[isp]*cv[i][isp];		//corrected to include alphat. April 2016
					cevrhs[i] -= alphat[isp]*qv[i][isp]*gamma1[iseg][isp]/ds[iseg];
					rhs[i] = cevrhs[i] - cevlhs[i];
//matrix of linear system
					for(j=1; j<=nnv; j++) mat[i][j] = gbarvv[i][j][isp] + alphat[isp]/dcdp[i]*alphaconv[i][j]*deltat/volseg[j];
					mat[i][i] += alphat[isp]*gamma1[iseg][isp]/ds[iseg];
					mat[i][nnv+1] = 1. - vol/deltat*gbarvtsum[i][isp];
					mat[nnv+1][i] = 1. - vol/deltat*gbarvtsum[i][isp];	//note gbartv = transpose of gbarvt
					matx[i] = 0.;
				}
				mat[nnv+1][nnv+1] = -vol/deltat*(nnt - vol/deltat*gbarttsum1[isp]);
				qtqvsum = 0.;
				sumtprev = 0.;
				sumt = nnt*g0[isp];
				for(i=1; i<=nnv; i++){
					sumt += gbarvtsum[i][isp]*qv[i][isp];
					qtqvsum += qv[i][isp];
				}
				for(itp=1; itp<=nnt; itp++){
					qtqvsum += qt[itp][isp];
					sumtprev += ctprev[itp][isp];
					sumt += vol*gttsum[itp][isp]*(ctprev[itp][isp] - g0prev[isp])
						+ gbarttsum[itp][isp]*(qt[itp][isp] - vol/deltat*(g0[isp] - g0prev[isp]));
				}
				rhs[nnv+1] = -qtqvsum + vol/deltat*(sumt - sumtprev);
				matx[nnv+1] = 0.;
				matfac = 1./mat[1][1];	//scale matrix to avoid very small values
				rhs2 = 0.;
				for(i=1; i<=nnv+1; i++){
					for(j=1; j<=nnv+1; j++) mat[i][j] *= matfac;
					rhs[i] *= matfac;
					rhs2 += SQR(rhs[i]);
				}
				if(rhs2 > 1.e-12){	//if RHS is not zero, solve linear system: sum mat[i][j]*qv[j] = rhs[i]
					if(useGPU) bicgstabBLASD(mat, rhs, matx, nnv+1, bicgstaberr, bicgstabit);
					else bicgstab(mat, rhs, matx, nnv+1, bicgstaberr, bicgstabit);
				}
				qvsum[isp] = 0.;
				if(kmain <= 5) lambda1 = 1.;
				else if(kmain <= 10) lambda1 = 0.5;	//if convergence is slow, underrelax
				else if(kmain <= 15) lambda1 = 0.2;	
				else if(kmain <= 20) lambda1 = 0.1;
				else if(kmain <= 25) lambda1 = 0.05;
				else lambda1 = 0.02;
				for(i=1; i<=nnv; i++){
					qv[i][isp] += lambda1*matx[i];
					qvsum[isp] += qv[i][isp];
				}
				g0[isp] += lambda1*matx[nnv+1];
			}
//update blood solute levels
			for(isp=1; isp<=nsp; isp++)	if(permsolute[isp]){
				ineg = 0;
				cvmin = 0.;
				cvmax = 0.;
				ihigh = 0;
				for(i=1; i<=nnv; i++){
					iseg = mainseg[i];
					cv[i][isp] = 0.;
					for(inodbc=1; inodbc<=nnodbc; inodbc++) cv[i][isp] += gammaconv[i][inodbc]*bcp[inodbc][isp];
					for(j=1; j<=nnv; j++) cv[i][isp] += betaconv[i][j]*cvprev[j][isp] - alphaconv[i][j]*qv[j][isp]*deltat/volseg[j];
					if(cv[i][isp] < 0.){	//check for solute concentrations outside physical range and give warnings
						ineg++;
						cvmin = FMIN(cvmin,cv[i][isp]);
					}
					if(oxygen[isp]) if(cv[i][isp] > cs*hd[iseg]){
						ihigh++;
						cvmax = FMAX(cvmax,cv[i][isp]);
					}
					if(oxygen[isp]) blood(cv[i][isp],hd[iseg],&pv[i],&dcdp[i]);
					else dcdp[i] = 1.;
				}
				if(ineg > 0 && greensverbose) printf("*** Warning: Solute %i cblood is negative in %i segments, minimum is %f\n",isp,ineg,cvmin);
				if(ihigh > 0 && greensverbose) printf("*** Warning: Solute %i cblood is high in %i segments, maximum is %f\n",isp,ihigh,cvmax);	
			}
//update tissue solute levels
			for(isp=1; isp<=nsp; isp++){
				if(diffsolute[isp]){
					if(useGPU){
						for(itp=1; itp<=nnt; itp++){
							ct[itp][isp] = g0[isp];
							if(permsolute[isp]) for(i=1; i<=nnv; i++) ct[itp][isp] += gbarvt[i][itp][isp]*qv[i][isp];
							h_qt[itp-1] = qt[itp][isp] - vol/deltat*(g0[isp] - g0prev[isp]);
							h_ctprev[itp-1] = vol*(ctprev[itp][isp] - g0prev[isp]);
						}
						tissueGPU1c(isp);	//Compute contribution from tissue source strengths and previous solute levels
						for(itp=1; itp<=nnt; itp++)	ct[itp][isp] += h_ct[itp-1];
					}
					else{
						for(itp=1; itp<=nnt; itp++){
							ct[itp][isp] = g0[isp];
							if(permsolute[isp]) for(i=1; i<=nnv; i++) ct[itp][isp] += gbarvt[i][itp][isp]*qv[i][isp];
							ix = tisspoints[1][itp];
							iy = tisspoints[2][itp];
							iz = tisspoints[3][itp];
							for(jtp=1; jtp<=nnt; jtp++){
								jx = tisspoints[1][jtp];
								jy = tisspoints[2][jtp];
								jz = tisspoints[3][jtp];
								jxyz = 1 + abs(ix - jx) + abs(iy - jy)*mxx + abs(iz - jz)*mxx*myy;
								ct[itp][isp] += vol*gtt[jxyz][isp]*(ctprev[jtp][isp] - g0prev[isp])
									+ gbartt[jxyz][isp]*(qt[jtp][isp] - vol/deltat*(g0[isp] - g0prev[isp]));
							}
						}
					}
				}
				else for(itp=1; itp<=nnt; itp++) ct[itp][isp] = ctprev[itp][isp] + qt[itp][isp]/vol*deltat;	//added April 2015
			}
//update tissue reaction rates
			for(isp=1; isp<=nsp; isp++){
				qtsum[isp] = 0.;
				qtpsum[isp] = 0.;
			}
			for(itp=1; itp<=nnt; itp++){
				for(isp=1; isp<=nsp; isp++) ctct[isp] = 0.5*(ct[itp][isp] + ctprev[itp][isp]);	//average over time step
				tissrate(nsp,ctct,mtiss,mptiss);
				for(isp=1; isp<=nsp; isp++){
					qt[itp][isp] = lambda[isp]*mtiss[isp]*vol + (1. - lambda[isp])*qt[itp][isp];	//underrelax here for nonlinear problems
					qtp[itp][isp] = mptiss[isp]*vol;	//added October 2016
					qtsum[isp] += qt[itp][isp];
					qtpsum[isp] += qtp[itp][isp];
				}
			}
			for(isp=1; isp<=nsp; isp++){
				if(permsolute[isp] == 0){	//update G0 for non-permeable solutes. Added October 2016.
					gbarttratesum = 0.;
					ctsum = 0.;
					ctprevsum = 0.;
					for(itp=1; itp<=nnt; itp++){
						ctsum += ct[itp][isp];
						ctprevsum += ctprev[itp][isp];
						ix = tisspoints[1][itp];
						iy = tisspoints[2][itp];
						iz = tisspoints[3][itp];
						for(jtp=1; jtp<=nnt; jtp++){
							jx = tisspoints[1][jtp];
							jy = tisspoints[2][jtp];
							jz = tisspoints[3][jtp];
							jxyz = 1 + abs(ix - jx) + abs(iy - jy)*mxx + abs(iz - jz)*mxx*myy;
							gbarttratesum += gbartt[jxyz][isp]*qtp[itp][isp];
						}
					}
					g0fac = (1. - gbarttratesum/nnt/2.)/(1. - gbarttsum1[isp]*vol/deltat/nnt);
					g0[isp] -= (qtsum[isp] - vol/deltat*(ctsum - ctprevsum))/(qtpsum[isp]/2. - nnt*vol/deltat)*g0fac;
				}
			}
			if(greensverbose) for(isp=1; isp<=nsp; isp++){
				printf("Solute %i: qtsum = %f",isp,qtsum[isp]);
				if(permsolute[isp]) printf(" qvsum = %f",qvsum[isp]);
				if(diffsolute[isp]) printf(" g0 = %f",g0[isp]);
				printf("\n");
			}
			//Convergence based on changes in cv, qv, ct, qt and g0
			if(kmain > 1){	//force at least two iterations to allow convergence testing
				convflag = 1;
				for(isp=1; isp<=nsp; isp++){
					err = 0.;
					if(permsolute[isp]) for(i=1; i<=nnv; i++){
						dif = fabs(cv[i][isp] - cvprevit[i][isp])/epsvesselc[isp];
						if(dif > err){
							imaxerr = i;
							err = dif;
							strcpy(varname, "cv");
						}
						dif = fabs(qv[i][isp] - qvprevit[i][isp])/epsvesselq[isp];
						if(dif > err){
							imaxerr = i;
							err = dif;
							strcpy(varname, "qv");
						}
					}
					for(itp=1; itp<=nnt; itp++){
						dif = fabs(ct[itp][isp] - ctprevit[itp][isp])/epstissuec[isp];
						if(dif > err){
							imaxerr = itp;
							err = dif;
							strcpy(varname, "ct");
						}
						dif = fabs(qt[itp][isp] - qtprevit[itp][isp])/epstissueq[isp];
						if(dif > err){
							imaxerr = itp;
							err = dif;
							strcpy(varname, "qt");
						}
					}
					dif = fabs(g0[isp] - g0previt[isp])/epstissuec[isp];
					if(dif > err){
						imaxerr = 0;
						err = dif;
						strcpy(varname, "g0");
					}
					if(greensverbose && err > 0.) printf("Solute %i: maximum relative error = %f in %s[%i]\n",isp,err,varname,imaxerr);
					if(err > 1.) convflag = 0;
				}
				if(convflag) goto mainconv;
			}
		}
		if(greensverbose) printf("*** Did not converge\n");
		ofp1 = fopen(fname2, "a");
		fprintf(ofp1,"*** Did not converge at time step %i\n", itime);
		fclose(ofp1);
		mainconv:
//********************** end of main loop *****************************
		tfinish1 = clock();
		duration = (float)(tfinish1 - tstart1)/CLOCKS_PER_SEC;
		if(greensverbose) printf("itime = %i, %2.3f seconds for step\n", itime,duration);
		for(isp=1; isp<=nsp; isp++){	//data for plot with updated values
			for(iseg=1; iseg<=nseg; iseg++) cvseg[iseg][isp] = 0.;
			if(permsolute[isp]) for(i=1; i<=nnv; i++){
				iseg = mainseg[i];
				if(oxygen[isp]) cvseg[iseg][isp] += pv[i]/nspoint[iseg];
				else cvseg[iseg][isp] += cv[i][isp]/nspoint[iseg];
			}
		}
		contour(itime,time,contourmethod);

		//calculate outflow concentration, for washout simulations
		ofp3 = fopen(fname1, "a");
		fprintf(ofp3,"%3i %12f ", itime, time);
		for(isp=1; isp<=nsp; isp++){
			concin = 0.;
			concout = 0.;
			if(permsolute[isp]) for(inodbc=1; inodbc<=nnodbc; inodbc++){
				inod = bcnod[inodbc];
				iseg = nodseg[1][inod];
				if(nodout[inod] == 1) concin += bcp[inodbc][isp]*qq[iseg]/flowsumin;	//inflow boundary nodes
				else{																	//outflow boundary nodes
					bcp[inodbc][isp] = 0.;
					for(j=1; j<=nnodbc; j++) bcp[inodbc][isp] += xiconv[inodbc][j]*bcp[j][isp];
					for(j=1; j<=nnv; j++) bcp[inodbc][isp] += omegaconv[inodbc][j]*cvprev[j][isp]
						- zetaconv[inodbc][j]*qv[j][isp]*deltat/volseg[j];
					concout += bcp[inodbc][isp]*qq[iseg]/flowsumout;
				}
			}
//Check conservation of solute
			sumt = 0.;														//tissue sum = previous + inputs - current
			for(itp=1; itp<=nnt; itp++) sumt += ctprev[itp][isp] - ct[itp][isp] + qt[itp][isp]*deltat/vol;
			if(permsolute[isp]) for(i=1; i<=nnv; i++) sumt += qv[i][isp]*deltat/vol;
			sumt = sumt/creft[isp]/nnt;
			if(fabs(sumt) > 0.0001)			//check solute conservation in tissue
				printf("*** Tissue conservation test: solute %i sum = %f\n", isp, sumt);
			sumv = 0.;
			if(permsolute[isp]){
				sumv = (concin*flowsumin - concout*flowsumout)*flowfac*deltat/vol;				//vessel sum = previous + inputs - current - outputs
				for(i=1; i<=nnv; i++) sumv += (cvprev[i][isp] - cv[i][isp])*volseg[i]/vol - qv[i][isp]*deltat/vol;
				sumv = sumv/creft[isp]/nnt;
				if(fabs(sumv) > 0.0001)			//check solute conservation in vessel
					printf("*** Vessel conservation test: solute %i sum = %f\n", isp, sumv);
			}
			cumulerror[isp] += sumt + sumv;
			if(fabs(cumulerror[isp]) > 0.001)	//check cumulative solute conservation (larger error tolerance)
				printf("*** Cumulative conservation test: solute %i sum = %f\n", isp, cumulerror[isp]);	
			fprintf(ofp3,"%12f ", cumulerror[isp]);
			if(permsolute[isp]) fprintf(ofp3,"%12f %12f ", concin, concout);
			if(diffsolute[isp]) fprintf(ofp3,"%12f ", g0[isp]);
		}
		fprintf(ofp3,"\n");
		fclose(ofp3);		
	}
//********************** end of time loop *****************************
	tfinish = clock();
	duration = (float)(tfinish - tstart)/CLOCKS_PER_SEC;
	printf("\n%i time steps, %3.1f seconds for main loop\n", itime,duration);

	ofp1 = fopen(fname2, "a");
	fprintf(ofp1,"%i time steps, %3.1f seconds for main loop\n", itime,duration);
	fclose(ofp1);

//output files describe state of system at final time point - useful for steady state
//general output file
	ofp = fopen("GreensRes.out", "w");
	for(isp=1; isp<=nsp; isp++) fprintf(ofp,"g0[%i] = %f\n", isp,g0[isp]);
	fprintf(ofp,"\n");
	for(isp=1; isp<=nsp; isp++) if(permsolute[isp]){
		fprintf(ofp,"Solute %i: qtsum = %f, qvsum = %f\n",isp,qtsum[isp],qvsum[isp]);
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
			qvseg[iseg][isp] = 0.;
			cvseg[iseg][isp] = 0.;
		}
		fprintf(ofp,"subseg# segment#       cv       qv\n");
		for(i=1; i<=nnv; i++){
			iseg = mainseg[i];
			if(permsolute[isp]) fprintf(ofp,"%8i %8i %10.4f %10.4f\n", i,iseg,cv[i][isp],qv[i][isp]);
			qvseg[iseg][isp] += qv[i][isp];
			cvseg[iseg][isp] += cv[i][isp]/nspoint[iseg];
		}
		fprintf(ofp,"segname length      cvseg    qvseg     gamma\n");
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5)
			fprintf(ofp,"%4i %10.4f %10.4f %10.4f %10.4f\n",
			segname[iseg],lseg[iseg],cvseg[iseg][isp],qvseg[iseg][isp],gamma1[iseg][isp]);
	}
	fclose(ofp);

	ofp = fopen("TissueSources.out", "w");
	ofp1 = fopen("TissueLevels.out", "w");
	for(isp=1; isp<=nsp; isp++){
		fprintf(ofp,"\nTissue source strengths for solute %i",isp);
		fprintf(ofp1,"\nTissue levels for solute %i",isp);
		cmax[isp] = -1.e8;
		cmean[isp] = 0.;
		cmin[isp] = 1.e8;
		for(itp=1; itp<=nnt; itp++){
			if(oxygen[isp]) ct[itp][isp] /= alphat[isp];	//convert to PO2 values for oxygen
			if(itp%10 == 1) fprintf(ofp,"\n");
			if(itp%10 == 1) fprintf(ofp1,"\n");
			fprintf(ofp," %12f", qt[itp][isp]);
			fprintf(ofp1," %12f", ct[itp][isp]);
			cmean[isp] += ct[itp][isp];
			cmax[isp] = FMAX(ct[itp][isp],cmax[isp]);
			cmin[isp] = FMIN(ct[itp][isp],cmin[isp]);
		}
		cmean[isp] = cmean[isp]/nnt;
		fprintf(ofp1,"\n%12f %12f %12f Solute %i: cmean, cmin, cmax\n", cmean[isp],cmin[isp],cmax[isp],isp);
	}
	fclose(ofp);
	fclose(ofp1);

	ofp = fopen("VesselSources.out", "w");
	ofp1 = fopen("VesselLevels.out", "w");
	for(isp=1; isp<=nsp; isp++) if(permsolute[isp]){
		fprintf(ofp,"\nVessel source strengths for solute %i",isp);
		fprintf(ofp1,"\nVessel levels for solute %i",isp);
		for(i=1; i<=nnv; i++){
			if(i%10 == 1) fprintf(ofp,"\n");
			if(i%10 == 1) fprintf(ofp1,"\n");
			fprintf(ofp," %12f", qv[i][isp]);
			fprintf(ofp1," %12f", cv[i][isp]);
		}
	}
	fclose(ofp);
	fclose(ofp1);
}