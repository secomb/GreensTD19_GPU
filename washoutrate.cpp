/************************************************************************
washoutrate - for GreensTD.  TWS February 2015
computes rate constant for washout decay 
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <cstdio>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "nrutil.h"

void washoutrate(int irun)
{
	extern int max;
	extern float N_D,N_P;
	int itime,ntime,idum,isum;
	float perm,diff,delta_t,time,concin,ambient,concout,cumulerror;
	float xsum,x2sum,ysum,xysum,afit,bfit,afit1,bfit1,ysum1,xysum1;
	char fname1[80],fname2[80],bb[100];
	FILE *ifp,*ofp;

	//use sequentially numbered files
    sprintf(fname1,"GreensLog%03i.txt",irun);
    sprintf(fname2,"GreensWashout%03i.txt",irun);

//read GreensLogxxx.txt
	ifp = fopen(fname1,"r");
	fgets(bb,max,ifp);	
	fscanf(ifp,"%f %f %f%*[^\n]", &delta_t,&diff,&perm);
	fscanf(ifp,"%i%*[^\n]", &ntime);
	fclose(ifp);
//read GreensWashoutxxx.txt
	ifp = fopen(fname2,"r");
	fgets(bb,max,ifp);
	xsum = 0.;
	ysum = 0.;
	ysum1 = 0.;
	x2sum = 0.;
	xysum = 0.;
	xysum1 = 0.;
	isum = 0;
	for(itime=0;itime<=ntime;itime++){	//sums for linear regression
		fscanf(ifp,"%i %f %i %f %f %f %f%*[^\n]", &idum,&time,&idum,&cumulerror,&concin,&concout,&ambient);
		if(itime >= ntime/2 && concout > 0.){	//use only second half of points to calculate lambda
			xsum += time;
			x2sum += SQR(time);
			ysum += log(concout);
			ysum1 += log(ambient);
			xysum += time*log(concout);
			xysum1 += time*log(ambient);
			isum++;
		}
	}
	fclose(ifp);
//linear regression
	afit = 0.;
	bfit = 0.;
	afit1 = 0.;
	bfit1 = 0.;
	if(isum > 0){
		bfit = (isum*xysum - xsum*ysum)/(isum*x2sum - xsum*xsum);
		afit = (ysum - bfit*xsum)/isum;
		bfit1 = (isum*xysum1 - xsum*ysum1)/(isum*x2sum - xsum*xsum);
		afit1 = (ysum1 - bfit*xsum)/isum;
	}
	if(irun) ofp = fopen("WashoutRate.txt","a");
	else{
		ofp = fopen("WashoutRate.txt","w");
		fprintf(ofp,"WashoutRate.txt\ndelta_t      diff     perm     N_D      N_P       lambda   lambda (ambient)\n");
	}
	fprintf(ofp,"%f %f %f %f %f %12.8f %12.8f\n",delta_t,diff,perm,N_D,N_P,bfit,bfit1);
	fclose(ofp);
}