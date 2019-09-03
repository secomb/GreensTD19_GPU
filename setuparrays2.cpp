/************************************************************************
setuparrays2 - for GreensTD.  TWS June 2012
Set up arrays with dimensions nnv and nnt
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays2(int nnv, int nnt)
{
	extern int nsp,nnodbc,nsl1,nsl2;
	extern int *mainseg,**tisspoints;

	extern float *pv,*dcdp,*volseg;
	extern float **cv,**ct,**cvprev,**ctprev,**cvprevit,**ctprevit;
	extern float **qv,**qt,**qtp,**qvprevit,**qtprevit;
	extern float ***gbarvv,***gvt,***gbarvt;
	extern float ***gbarvc,***gtc,***gbartc;
	extern float **ax;
	extern float *cevlhs,*cevrhs;
	extern double **mat,*rhs,*matx;
	extern float *transit;
	extern float **alphaconv,**betaconv,**gammaconv,**omegaconv,**zetaconv,**xiconv;
	extern float **gbarvtsum,**gbarttsum,**gttsum,*gttsum1,*gbarttsum1;

	mainseg = ivector(1,nnv);
	tisspoints = imatrix(1,3,1,nnt);
	dcdp = vector(1,nnv);
	pv = vector(1,nnv);
	volseg = vector(1,nnv);

	cv = matrix(1,nnv,1,nsp);
	ct = matrix(1,nnt,1,nsp);
	cvprev = matrix(1,nnv,1,nsp);
	ctprev = matrix(1,nnt,1,nsp);
	cvprevit = matrix(1,nnv,1,nsp);
	ctprevit = matrix(1,nnt,1,nsp);
	qv = matrix(1,nnv,1,nsp);
	qt = matrix(1,nnt,1,nsp);
	qtp = matrix(1,nnt,1,nsp);	//added October 2016
	qvprevit = matrix(1,nnv,1,nsp);
	qtprevit = matrix(1,nnt,1,nsp);
	ax = matrix(1,3,1,nnv);

	gbarvv = f3tensor(1,nnv,1,nnv,1,nsp);
	gvt = f3tensor(1,nnv,1,nnt,1,nsp);
	gbarvt = f3tensor(1,nnv,1,nnt,1,nsp);

	gbarvc = f3tensor(1,nsl1*nsl2,1,nnv,1,nsp);
	gtc = f3tensor(1,nsl1*nsl2,1,nnt,1,nsp);
	gbartc = f3tensor(1,nsl1*nsl2,1,nnt,1,nsp);

	alphaconv = matrix(1,nnv,1,nnv);
	betaconv = matrix(1,nnv,1,nnv);
	gammaconv = matrix(1,nnv,1,nnodbc);
	omegaconv = matrix(1,nnodbc,1,nnv);
	zetaconv = matrix(1,nnodbc,1,nnv);
	xiconv = matrix(1,nnodbc,1,nnodbc);
	transit = vector(1,nnv);
	cevlhs = vector(1,nnv);
	cevrhs = vector(1,nnv);
	gbarvtsum = matrix(1,nnv,1,nsp);
	gbarttsum = matrix(1,nnt,1,nsp);
	gttsum = matrix(1,nnt,1,nsp);
	gttsum1 = vector(1,nsp);
	gbarttsum1 = vector(1,nsp);

	mat = dmatrix(1,nnv+1,1,nnv+1);
	rhs = dvector(1,nnv+1);
	matx = dvector(1,nnv+1);
}