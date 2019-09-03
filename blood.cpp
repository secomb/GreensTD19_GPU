/*****************************************************
Oxygen concentration and partial pressure - used for oxygen only
Calculate p = p(c) where p is the inverse function of f.
c = f(p) is the hill equation  of the oxygen dissociation curve,
c (input) oxygen concentration in red blood cell
p (output) partial oxygen pressure
pp (output) pp = f'(p) at the computed p
Version 3.0, May 17, 2011.
Revised TWS January 2012 to use consistent f(p) for all cases
******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float rtflsp(float (*func)(float), float x1, float x2, float xacc);
float rtbis(float (*func)(float), float x1, float x2, float xacc);
float func(float x);
float bloodconc(float p,float h);

void blood(float c, float h, float *p, float *pp)
{
	extern float fn,alphaO2,p50,cs,cext,hext;
	extern float plow,phigh,clowfac,chighfac,pphighfac;//added January 2012
	float pl,ph;
	float clow,chigh,pphigh;

	if(h < 1.e-6){	//changed from 1.e-3, January 2012
		*p = c/alphaO2;
		*pp = alphaO2;
		return;
	}
//changed for better behavior in severe hypoxia. 2 Oct 08
	if(c < 0.){
		*p = c/alphaO2;
		*pp = alphaO2;
		return;
	}
	clow = clowfac*h + alphaO2*plow;
	if(c < clow){
		*p = c*plow/clow;
		*pp = clow/plow;
		return;
	}
	chigh = chighfac*h + alphaO2*phigh;
	if(c < chigh){
		if(c/h/cs < 1.){
			ph = pow((c/h/cs)/(1.0-c/h/cs),1./fn)*p50;
			pl = 0.;	//0. here to be sure to bracket the root.  June 2009
		}
		else{
			ph = phigh;
			pl = plow;//changed August 2010
		}
		hext = h;
		cext = c;
		*p = rtflsp(func, pl, ph, 0.001f);//lower tolerance, August 2010
//if false position algorithm fails, use bisection!  June 2009.
		if(*p < 0.) *p = rtbis(func, pl, ph, 0.001f);
		*pp = cs*h*fn/p50*pow(*p/p50,(fn-1))/SQR(1. + pow(*p/p50,fn)) + alphaO2;
		return;
	}
	pphigh = pphighfac*h + alphaO2;
	*p = phigh + (c - chigh)/pphigh;
	*pp = pphigh;
	return;
}
//Using the false position method, find the root of a function func known to lie between x1 and x2.
//The root, returned as rtflsp, is refined until its accuracy is +- xacc.

float rtflsp(float (*func)(float), float x1, float x2, float xacc)
{
	int j;
	int maxit=30;
	float fl,fh,xl,xh,swap,dx,del,f,rtf;
	fl = (*func)(x1);
	fh = (*func)(x2);
	if(fl*fh > 0.0){
		printf("\n*** Warning: Root must be bracketed in rtflsp %f %f\n", x1,x2);
		return -1.0;
	}
	if (fl < 0.0){
		xl = x1;
		xh = x2;
	}
	else {
		xl = x2;
		xh = x1;
		swap = fl;
		fl = fh;
		fh = swap;
	}
	dx = xh - xl;
	for(j=1; j<=maxit; j++){
		rtf = xl + dx*fl/(fl - fh);
		f = (*func)(rtf);
		if (f < 0){
			del = xl - rtf;
			xl = rtf;
			fl = f;
		}
		else {
			del = xh - rtf;
			xh = rtf;
			fh = f;
		}
		dx = xh - xl;
		if(fabs(del) < xacc || f == 0.0) return rtf;
	}
//	printf("\n*** Warning: Maximum number of iterations exceeded in rtflsp %f %f\n", x1,x2);
	return -1.0;
}

float rtbis(float (*func)(float), float x1, float x2, float xacc)
{
	int j;
	int maxit=30;
	float dx,f,fmid,xmid,rtb;
	f = (*func)(x1);
	fmid = (*func)(x2);
	if(f*fmid >= 0.0){
		printf("*** Error: Root must be bracketed in rtbis %f %f\n", x1,x2);
		return 0.0;
	}
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for(j=1; j<=maxit; j++){
		fmid = (*func)(xmid=rtb+(dx *= 0.5));
		if(fmid <= 0.) rtb = xmid;
		if(fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	printf("*** Error: Maximum number of iterations exceeded in rtbis %f %f\n", x1,x2);
	return 0.0;
}

float func(float p)
{
	extern float cext,hext;
	float y;
	y = bloodconc(p,hext) - cext;
	return y;
}

float bloodconc(float p, float h)
{
	extern float fn,p50,alphaO2,cs;
	extern float plow,phigh,clowfac,chighfac,pphighfac;//added January 2012
	float cv;
	if(p < 0.) cv = alphaO2*p;  //This gives more consistent behavior when p<0
	else if(p < plow) cv = clowfac*h*p/plow + alphaO2*p;
	else if(p < phigh) cv = cs*h*(1.0 - 1.0/(1.0 + pow((p/p50),fn))) + alphaO2*p;
	else cv = (chighfac + (p - phigh)*pphighfac)*h + alphaO2*p;
	return cv;
}
float bloodconcp(float p, float h)
{
	extern float fn,p50,alphaO2,cs;
	extern float plow,phigh,clowfac,chighfac,pphighfac;//added January 2012
	float cvp;
	if(p < 0.) cvp = alphaO2;  //This gives more consistent behavior when p<0
	else if(p < plow) cvp = clowfac*h/plow + alphaO2;
	else if(p < phigh) cvp = cs*h*fn/p50*pow(p/p50,(fn-1))/SQR(1. + pow(p/p50,fn)) + alphaO2;
	else cvp = pphighfac*h + alphaO2;
	return cvp;
}