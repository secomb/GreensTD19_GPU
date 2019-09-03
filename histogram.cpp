/*****************************************************
Evaluate histograms of solute levels.  TWS December 07.
Version 2.0, May 1, 2010.
Version 3.0, Ma7 17,2011.
******************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include "nrutil.h"
#include <stdio.h>

void histogram(void)
{
	extern int mxx,myy,mzz,nnt,nnv,nseg,nnod,nsp,*oxygen;
	extern float *cmin,*cmax,*cmean,*creft,*alphat;
	extern float **ct;
	int i,j,itp,isp,nctop,ncbot,nstat,nstatm=101;
	float step,dev;
	float *ctsamp,*stat,*cumu,*mstat;
	FILE *ofp;
	ctsamp = vector(1,nstatm);
	stat = vector(1,nstatm);
	cumu = vector(1,nstatm);
	mstat = vector(1,nstatm);
		
	ofp = fopen("Histogram.out", "w");
	for(isp=1; isp<=nsp; isp++){
		step = creft[isp]/100;
		if(oxygen[isp]) step /= alphat[isp];
		nctop = cmax[isp]/step + 1.;
		if(nctop > 100) nctop = 100;
		ncbot = cmin[isp]/step;
		if(ncbot <= 2) ncbot = 0;
		nstat = nctop - ncbot + 1;
		dev = 0.;
		if(nstat > nstatm) printf("*** nstatm too small in histogram\n");
		for(i=1; i<=nstat; i++){
			ctsamp[i] = step*(i - 1 + ncbot);
			mstat[i] = 0;
		}
		for(itp=1; itp<=nnt; itp++){
			dev = dev + SQR(cmean[isp] - ct[itp][isp]);
			for(j=1; j<=nstat; j++)	if(ct[itp][isp] <= ctsamp[j]){
				mstat[j] = mstat[j] + 1;
				goto binned;
			}
			binned:;
		}
		dev = sqrt(dev/nnt);
		for(i=1; i<=nstat; i++) stat[i] = mstat[i]*100./nnt;
		cumu[1] = stat[1];
		for(i=2; i<=nstat; i++) cumu[i] = cumu[i-1] + stat[i];

		fprintf(ofp,"Histogram data for solute %i\n", isp);
		fprintf(ofp,"value  %% cumul. %%\n");
		for(i=1; i<=nstat; i++) fprintf(ofp,"%g %7.2f %7.2f\n", ctsamp[i],stat[i],cumu[i]);
		fprintf(ofp,"Mean = %f deviation = %f min = %f max  = %f\n", cmean[isp],dev,cmin[isp],cmax[isp]);
	}
	fclose(ofp);
	free_vector(ctsamp,1,nstatm);
	free_vector(stat,1,nstatm);
	free_vector(cumu,1,nstatm);
	free_vector(mstat,1,nstatm);
}