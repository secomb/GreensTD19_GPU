/*****************************************************
Tissue uptake rates of solutes as a function of levels
Updated for Greens12TD - August 2012.
Note: uses concentration, not PO2
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>
void tissrate(int nsp, float *c, float *mtiss, float *mptiss)
{	
	extern int *oxygen;
	extern float **tissparam;

	int isp;

	#include "tissrate.cpp.dat"
}



/*

Some sample code for tissrate.cpp.dat
	int isp;
	float ccr,m0,kP;
	float kinstabA,kfA,krA,kAH,kinstabH,kfH,krH,kHP,kinstabM,kfM,krM,kMP,po2;	//For PR104 model

	for(isp=1; isp<=nsp; isp++){
		if(isp == 1){
			kP = tissparam[1][isp];
			mtiss[isp] = -kP*c[isp];
			mptiss[isp] = -kP;
		}
*/
/*
			m0 = tissparam[1][isp];
			ccr = tissparam[2][isp];
			if(c[isp] >= 0.){
				mtiss[isp] = -m0*c[isp]/(c[isp] + ccr);
				mptiss[isp] = -m0*ccr/SQR(c[isp] + ccr);
			}
			else{
				mtiss[isp] = 0.;
				mptiss[isp] = 0.;
			}
		}
		if(isp == 2){		//prodrug
*/

/*

if(nsp >= 3){
//isp = 2: PR104A - free
//isp = 3: PR104A - bound
		phi = tissparam[3][2];
		kfA = tissparam[1][2];	//note: dependence on oxygen c[1] can be introduced here
		kinstabA = tissparam[2][2];
		krA = tissparam[1][3];
		kAH = tissparam[2][3]*(0.023+0.0977/(0.1+po2)); // oxygen dependence of intracellular metabolism
		mtiss[2] = -(kfA*phi + kinstabA)*c[2] + phi*krA*c[3];
		mptiss[2] = -(kfA*phi + kinstabA);
		mtiss[3] = -(krA + kAH)*c[3] + kfA*c[2];
		mptiss[3] = -(krA + kAH);
	}
	if(nsp >= 5){
//isp = 4: PR104H - bound
//isp = 5: PR104H - free
		krH = tissparam[1][4];
		kHP = tissparam[2][4]; //*(0.2+8/(10+po2));  //*(0.2+0.08/(0.1+po2)) oxygen dependence of h-> M could be added here
		kfH = tissparam[1][5];
		kinstabH = tissparam[2][5];
		mtiss[4] = -(krH + kHP)*c[4] + kfH*c[5] + kAH*c[3];
		mptiss[4] = -(krH + kHP);
		mtiss[5] = -(kfH*phi + kinstabH)*c[5] + phi*krH*c[4];
		mptiss[5] = -(kfH*phi + kinstabH);
	}
	if(nsp >= 7){
//isp = 4: PR104M - bound ---added Jan 2012 KH
//isp = 5: PR104M - free
		krM = tissparam[1][6];
		kMP = tissparam[2][6];
		kfM = tissparam[1][7];
		kinstabM = tissparam[2][7];
		mtiss[6] = -(krM + kMP)*c[6] + kfM*c[7] + kHP*c[4];
		mptiss[6] = -(krM + kMP);
		mtiss[7] = -(kfM*phi + kinstabM)*c[7] + phi*krM*c[6];
		mptiss[7] = -(kfM*phi + kinstabM);
	}
}


		float gf0,gf1,gf2,aterm;
		case 1: //tracer
			mtiss[isp] = 0.;
			mptiss[isp] = 0.;
			break;

		case 2: //VEGF: non-permeable diffusible solute, based on Mac Gabhann and Popel - 2010
			if(c[1] <= 1.) mtiss[2] = 6.*tissparam[1][2];
			else if(c[1] <= 20.) mtiss[2] = (1. + 5.*pow((20. - c[1])/19.,3.))*tissparam[1][2];
			else mtiss[2] = tissparam[1][2];
			mtiss[2] -= tissparam[2][2]*c[2];
			mptiss[2] = -tissparam[2][2];

		case 2: //non-permeable diffusible solute produced in hypoxic regions - old version 
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			if(c[1] >= 0.) mtiss[isp] = gf0*ccr/(c[1] + ccr) - gf1*c[isp];
			else mtiss[isp] = gf0 - gf1*c[isp];
			mptiss[isp] = -gf1;
			break;
		case 2: //permeable solute delivered in blood with linear consumption in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(c[1] >= 0.) mtiss[isp] = - gf1*c[isp]*ccr/(c[1] + ccr);
			else mtiss[isp] =  - gf1*c[isp];
			mptiss[isp] = -gf1*ccr/(c[1] + ccr);
			break;

		case 2: //permeable solute delivered in blood with linear consumption in hypoxic regions - KOH, updated TWS July 2011
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(c[1] >= 0.){
				aterm = gf1*c[isp]*ccr/(c[1] + ccr);
				mtiss[isp] = -aterm;
				mptiss[isp] = -gf1*ccr/(c[1] + ccr);
			}
			else{
				aterm = gf1*c[isp];
				mtiss[isp] = -aterm;
				mptiss[isp] = -gf1;
			}
			
			break;
		case 3: //permeable diffusiblbe solute produced in hypoxic regions by solute 2
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(c[1] >= 0.) mtiss[isp] =  aterm;
			else mtiss[isp] =  aterm;
			mptiss[isp] = 0;
			break;
/*
		case 3: //permeable solute delivered in blood with linear consumption in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(c[1] >= 0.) mtiss[isp] = - gf1*c[isp]*ccr/(c[1] + ccr);
			else mtiss[isp] =  - gf1*c[isp];
			mptiss[isp] = -gf1*ccr/(c[1] + ccr);
			break;
		case 3: //permeable solute delivered in blood with linear consumption
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			mtiss[isp] = -gf1*c[isp];
			mptiss[isp] = -gf1;
			break;

		case 4: //non-diffusible solute produced in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			if(c[1] >= 0.) mtiss[isp] = gf0*ccr/(c[1] + ccr) - gf1*c[isp];
			else mtiss[isp] = gf0 - gf1*c[isp];
			mptiss[isp] = -gf1;
			break;

		case 4: //non-diffusible non permeable solute produced from case 2 in hypoxic regions
			gf0 = tissparam[1][isp];
			gf1 = tissparam[2][isp];
			gf2 = tissparam[3][isp];
			if(c[1] >= 0.) mtiss[isp] =  aterm - gf0*c[isp];
			else mtiss[isp] =  aterm - gf0*c[isp];
			mptiss[isp] = -gf0;
			break;

*/