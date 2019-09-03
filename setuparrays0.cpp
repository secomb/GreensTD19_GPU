/************************************************************************
setuparrays0 - for GreensTD.  TWS June 2012
Set up arrays with fixed dimensions
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays0()
{
	extern int nsp,mxx,myy,mzz;
	extern int ***nbou;
	extern float *cmin,*cmax,*cmean,*mtiss,*mptiss,*ctct,*qtsum,*qtpsum,*qvsum;
	extern float *p,*epsvesselq,*epstissueq,*epsvesselc,*epstissuec,*g0,*g0prev,*g0previt;
	extern float *qvfac; //May 2010
	extern float *cumulerror;
	extern float *x,*y,*ss,*axt,*ayt,*azt;
	extern float **gtt,**gbartt;

	nbou = i3tensor(1,mxx,1,myy,1,mzz);

	cmin = vector(1,nsp);
	cmax = vector(1,nsp);
	cmean = vector(1,nsp);
	mtiss = vector(1,nsp);
	mptiss = vector(1,nsp);
	g0prev = vector(1,nsp);
	g0previt = vector(1,nsp);
	ctct = vector(1,nsp);
	qtsum = vector(1,nsp);
	qtpsum = vector(1,nsp);
	qvsum = vector(1,nsp);
	p = vector(1,nsp);
	epsvesselq = vector(1,nsp);
	epstissueq = vector(1,nsp);
	epsvesselc = vector(1,nsp);
	epstissuec = vector(1,nsp);
	qvfac = vector(1,nsp);
	cumulerror = vector(1,nsp);

	x = vector(1,3);
	y = vector(1,3);
	ss = vector(1,3);

	axt = vector(1,mxx);
	ayt = vector(1,myy);
	azt = vector(1,mzz);

	gtt = matrix(1,mxx*myy*mzz,1,nsp);
	gbartt = matrix(1,mxx*myy*mzz,1,nsp);
}