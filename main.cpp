/************************************************************************
Main program to call greensTD
Time-dependent Greens function method, June 9, 2012. Updated August 2017.
See greensTD.cpp for description of changes.
Cuda 10.1 version, September 2019
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#if defined(__linux__)
// Requires c++17 support, should be included in all current linux releases
#include <experimental/filesystem> 
namespace fs = std::experimental::filesystem::v1;
#elif defined(__APPLE__)
// Requires removal of the -lstdc++fs flag from makefile
#include <filesystem>
namespace fs = std::filesystem;
#elif defined(_WIN32)    //Windows version
#include <Windows.h>
#endif

void input(int irun, int n_step, int n_diff, int n_perm);
void analyzenet();
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void greensTD(int irun);
void histogram();
void setuparrays0();
void setuparrays1(int nseg, int nnod);
void setuparrays2(int nnv, int nnt);
void cmgui(float *segvar);
void washoutrate(int irun);

void bicgstabBLASDinit(int nnvGPU);
void bicgstabBLASDend(int nnvGPU);
void tissueGPUinit(int nntGPU, int nnvGPU);
void tissueGPUend(int nntGPU, int nnvGPU);

int max=100,nmax,initgreens,ninterval;
int mxx,myy,mzz,nnt,nseg,nnod,nnodfl,nnv,nsp,nnodbc,nodsegm;
int slsegdiv,nsl1,nsl2,ntpts,contourmethod;
int *mainseg,*permsolute,*nodrank,*nodtyp,*nodout,*bcnodname,*bcnod,*bctyp;
int *nodname,*segname,*segtyp,*nspoint,*istart,*nl,*nk,*ista,*iend,*nresis;
int *intervalnst,*oxygen,*diffsolute;
int **segnodname,**nodseg,**tisspoints,**nodnod,**intervalin;
int ***nbou;

float fn,c,alphaO2,p50,cs,cext,hext,req,q0fac,totalq,flowfac=1.e6/60.;
float plow,phigh,clowfac,chighfac,pphighfac;
float pi1 = atan(1.)*4.,fac = 1./4./pi1,fact,fact1,fact2;
float lb,maxl,v,vol,vdom,errfac,tlength,tarea,alx,aly,alz;
float deltat,xmax,ymax,scalefac,N_D,N_P;
float *axt,*ayt,*azt,*ds,*diff,*alphab,*alphat,*cmin,*cmax,*cmean,*crefb,*creft,*cinitb,*cinitt,*lambda;
float *diam,*rseg,*q,*qq,*hd,*segc,*bcprfl,*bchd,*nodvar,*segvar,*qvfac;
float *tpts,*intervaldt,*transit,*cumulerror;
float *x,*y,*lseg,*volseg,*ss,*mtiss,*mptiss;
float *g0,*g0prev,*g0previt,*ctct,*qtsum,*qtpsum,*qvsum,*gttsum1,*gbarttsum1;
float *xsl0,*xsl1,*xsl2,*clmin,*clint,*cl,*cevlhs,*cevrhs;
float *epsvesselq,*epstissueq,*epsvesselc,*epstissuec,*p,*pv,*dcdp;
float **start,**end,**scos,**ax,**cnode,**resisdiam,**resis,**bcp,**bcpbase;
float **cv,**ct,**cvprev,**ctprev,**cvprevit,**ctprevit;
float **qv,**qt,**qtp,**qvprevit,**qtprevit,**qvseg,**cvseg;
float **tconcs,**gbarvtsum,**gbarttsum,**gttsum;
float **alphaconv,**betaconv,**gammaconv,**omegaconv,**zetaconv,**xiconv;
float **tissparam,**gamma1,**gtt,**gbartt,**zv;
float ***gbarvv,***gvt,***gbarvt,***gbarvc,***gtc,***gbartc,***psl;
double **mat,*rhs,*matx;

//Needed for GPU version
int useGPU,nnvGPU,nntGPU;
double *h_x, *h_b, *h_a, *h_rs;
float *h_rst;
double *d_a, *d_x, *d_b, *d_res, *d_r, *d_rs, *d_v, *d_s, *d_t, *d_p, *d_er;
float *d_xt, *d_bt, *d_rest, *d_rt, *d_rst, *d_vt, *d_st, *d_tt, *d_pt, *d_ert;
int *h_tisspoints,*d_tisspoints;
float *h_ct,*h_ctprev,*h_qt,*h_qtp,*h_pv,*h_qv,*h_tissxyz,*h_vessxyz;
float *d_qt,*d_qtp,*d_ct,*d_ctprev,*d_qv,*d_pv;
float *h_gtt,*d_gtt,*h_gbartt,*d_gbartt;
float *d_tissxyz,*d_vessxyz;

int main(int argc, char *argv[])
{
	int iseg,inod,irun,nrun,n_step,n_diff,n_perm;
	FILE *ifp;

#if defined(__unix__) //Create a Current subdirectory if it does not already exist.
	if (!fs::exists("Current")) fs::create_directory("Current");
#elif defined(_WIN32)
	DWORD ftyp = GetFileAttributesA("Current\\");
	if (ftyp != FILE_ATTRIBUTE_DIRECTORY) system("mkdir Current");
#endif

//This file defines job. Use int n_step = int n_diff = int n_perm = 1 for single run.
	ifp = fopen("DefineRuns.dat", "r");	//"SoluteParams.dat"
	fscanf(ifp,"%i%*[^\n]", &n_diff);
	fscanf(ifp,"%i%*[^\n]", &n_perm);
	fscanf(ifp,"%i%*[^\n]", &n_step);
	nrun = n_diff*n_perm*n_step;
	fclose(ifp);

	for(irun=0;irun<nrun;irun++){

		input(irun,n_step,n_diff,n_perm);

		if(ntpts > 0){	//if no time points in TimeDep.dat, skip that case
			if(irun == 0){	//no need to set up arrays every time

				setuparrays0();

				setuparrays1(nseg,nnod);

				analyzenet();

				setuparrays2(nnv,nnt);

				if(useGPU){
					nntGPU = mxx*myy*mzz;	//this is the maximum possible number of tissue points
					nnvGPU = nnv + 1;	//linear system dimension is nnv + 1
					bicgstabBLASDinit(nnvGPU);
					tissueGPUinit(nntGPU, nnvGPU);
				}
			}
			for(iseg=1; iseg<=nseg; iseg++) segvar[iseg] = segname[iseg];
			for(inod=1; inod<=nnod; inod++) nodvar[inod] = nodname[inod];
			picturenetwork(nodvar,segvar,"NetNodesSegs.ps");

			greensTD(irun);

			for(iseg=1; iseg<=nseg; iseg++) segvar[iseg] = cvseg[iseg][1];
			for(inod=1; inod<=nnod; inod++) nodvar[inod] = nodname[inod];
			picturenetwork(nodvar,segvar,"NetNodesSolute1.ps");

			for(iseg=1; iseg<=nseg; iseg++) segvar[iseg] = log10(qq[iseg]);
			cmgui(segvar);
			histogram();
			washoutrate(irun);
			if(irun == nrun){
				if(useGPU){
					tissueGPUend(nntGPU, nnvGPU);
					bicgstabBLASDend(nnvGPU);
				}
			}
		}
	}
	return 0;
}