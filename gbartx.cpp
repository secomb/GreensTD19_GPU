/*****************************************************
gbartx - Evaluate interaction functions needed for GreensTD.  TWS, August 2012.
Note use of Greens function integrated over a spherical region for
for evaluating all interactions.
dist1: distance between centers
req1: equivalent radius of (larger) interacting volume
Updated January 2015 to include explicit dependence on req1 
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

float erfcc(float x);

void gbartx(int isp, float dist1, float req1, float *gtx, float *gbartx)
{
	extern float fact,fact1,fact2,pi1,deltat,*diff;
	float a1fact,a2fact,x0,kkx0,erfcc1,erfcc2,exp1,exp2,pi12;

	a1fact = req1/sqrt(fact);
	a2fact = SQR(req1)/fact;
	x0 = dist1/req1;
	pi12 = sqrt(pi1);
	if(x0 > 1.) kkx0 = 0.;
	else kkx0 = 2.;
	if(fabs(x0) < 1.e-6){	//coincident points
		erfcc1 = erfcc(a1fact);
		exp1 = exp(-a2fact);
		*gtx = fact1*(0.75*pi12/a1fact*(1. - erfcc1) - 1.5*exp1)/a2fact;
		*gbartx = -0.75*fact2/req1*(2./pi12/a1fact*exp1 - 1./a2fact + (1./a2fact - 2.)*erfcc1);
	}
	else if(fabs(x0) < 10.){	//close points
		erfcc1 = erfcc(a1fact*(x0 - 1));
		erfcc2 = erfcc(a1fact*(x0 + 1));
		exp1 = exp(-a2fact*SQR(x0 - 1));
		exp2 = exp(-a2fact*SQR(x0 + 1));
		*gtx = fact1*3./(8.*a2fact*a2fact*x0)*(pi12*a1fact*x0*(erfcc1 - erfcc2) - exp1 + exp2);
		*gbartx = fact2/dist1/4.*(1./pi12/a1fact*
			(((x0 - 2.)*(x0 + 1) + 1./a2fact)*exp2 - ((x0 + 2.)*(x0 - 1) + 1./a2fact)*exp1)
			- ((x0 - 2.)*SQR(x0 + 1) + 1.5*x0/a2fact)*erfcc2 + ((x0 + 2.)*SQR(x0 - 1) + 1.5*x0/a2fact)*erfcc1
			- kkx0*(x0 + 2.)*SQR(x0 - 1));
	}
	else{	//distant points
		*gtx = fact1*exp(-SQR(dist1)/fact);
		*gbartx = fact2/dist1*erfcc(dist1/sqrt(fact));
	}
}
