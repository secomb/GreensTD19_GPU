/****************************************************************************
initgreens - initial concentrations and source strengths
TWS August 2012
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void tissrate(int nsp, float *p, float *mtiss, float *mptiss);
float bloodconc(float p,float h);

void initgreensTD()
{
	extern int nnt,nnv,nsp;
	extern int *mainseg,*permsolute,*oxygen,*diffsolute;

	extern float vol,errfac,tlength;
	extern float *mtiss,*mptiss,*epsvesselq,*epstissueq,*epsvesselc,*epstissuec,*p;
	extern float *g0,*ds,*qq,*hd,*alphap,*alphat,*crefb,*creft,*cinitb,*cinitt,*qtsum,*qvsum;
	extern float **qt,**qv,**ct,**cv,**tissparam;

	int isp,i,itp;

	tissrate(nsp,cinitt,mtiss,mptiss);
	for(isp=1; isp<=nsp; isp++){
		qtsum[isp] = 0.;
		for(itp=1; itp<=nnt; itp++){
			qt[itp][isp] = mtiss[isp]*vol;
			qtsum[isp] += qt[itp][isp];
			ct[itp][isp] = cinitt[isp];
		}
		qvsum[isp] = 0.;
		if(permsolute[isp]) for(i=1; i<=nnv; i++){
			qv[i][isp] = -mtiss[isp]*vol*nnt*ds[mainseg[i]]/tlength;
			qvsum[isp] += qv[i][isp];
			cv[i][isp] = cinitb[isp];
		}
		g0[isp] = cinitt[isp];
	}
//set error bounds on qt, qv, ct, cv proportional to errfac
	tissrate(nsp,creft,mtiss,mptiss);
	for(isp=1; isp<=nsp; isp++){
		epstissueq[isp] = FMAX(fabs(mtiss[isp])*vol*errfac,0.001);
		epsvesselq[isp] = nnt*epstissueq[isp]/nnv;
		epstissuec[isp] = creft[isp]*errfac;
		epsvesselc[isp] = crefb[isp]*errfac;
	}
}