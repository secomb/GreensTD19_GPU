/************************************************************************
diffusionmatrix - Coefficients for time-dependent Greens function method
TWS July 2012, updated March 2017
**************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float expBessi0(float x);
float erfcc(float x);
void gbartx(int isp, float dist1, float req1, float *gtx, float *gbartx);

void diffusionmatrix(void)
{
	extern int mxx,myy,mzz,nnt,nnv,nseg,nsp,nnodbc,slsegdiv,nsl1,nsl2;
	extern int contourmethod;
	extern int *mainseg,**tisspoints,*diffsolute,*permsolute;
	extern float req,fac,deltat,pi1,fact,fact1,fact2;
	extern float *axt,*ayt,*azt,*diff,*ds,*rseg;
	extern float **gtt,**gbartt;
	extern float ***gbarvv,***gvt,***gbarvt;
	extern float ***gbarvc,***gtc,***gbartc;
	extern float **gbarvtsum,**gbarttsum,**gttsum,*gttsum1,*gbarttsum1;
	extern float *xsl0,*xsl1,*xsl2;
	extern float **start,**scos,**ax,*x,*y;

	int isl1,isl2,isl12,gvtsegdiv,gvtsegdiv1;
	int i,j,j1,k,k1,iseg,jseg,ix,iy,iz,jx,jy,jz,jxyz,isp,it,itp,jtp;
	int nit=10000;	//number of points in time integration of gbarvv - integral converges slowly with nit!
	float dist1,dist2,fact3,fact4,b2fact,c2fact,tint,integrand,dummy,gbarvc1,a1fact,a2fact,req1;
	float gvt0,gbarvt0,gbarvv0;
	float xmin,ymin,xmax,ymax,**cos;
	//***** May need to increase up to 10.*req if instability occurs
	float lsegmax = 2.*req;	//length criterion for subdividing segments in calculation of gvv, gvt

	cos = matrix(1,3,1,3);
	printf("Computing diffusion coefficient matrices...");
	for(isp=1; isp<=nsp; isp++) if(diffsolute[isp]){
		fact = 4.*diff[isp]*deltat;
		fact1 = 1./(sqrt(fact*pi1)*fact*pi1);
		fact2 = 1./(4.*pi1*diff[isp]);
// tissue ~ tissue matrices gtt and gbartt, based on relative positions of tissue points
		printf("gtt,gbartt...");
		req1 = req;	//based on tissue volume
		for(jx=1; jx<=mxx; jx++) for(jy=1; jy<=myy; jy++) for(jz=1; jz<=mzz; jz++){
			jxyz = jx + (jy - 1)*mxx + (jz - 1)*mxx*myy;
			if(jxyz == 1){	//self-interaction
				a1fact = req/sqrt(fact);
				a2fact = SQR(req)/fact;
				gtt[jxyz][isp] = fact1*3./(16.*a2fact*a2fact*a2fact)*(1. - 6.*a2fact + (2.*a2fact - 1.)*exp(-4.*a2fact)
					+ 4.*sqrt(pi1)*a1fact*a2fact*(1. - erfcc(2.*a1fact)));
				gbartt[jxyz][isp] = 1./(4.*pi1*diff[isp]*req)*3./40*(16. - (16. - 10/a2fact)*(1. - erfcc(2.*a1fact))
					- 1./sqrt(pi1*a2fact)*(-1./SQR(a2fact) + 10./a2fact + (1./SQR(a2fact) - 6./a2fact + 8.)*exp(-4.*a2fact)));
			}
			else{
				dist2 = SQR(axt[1]-axt[jx]) + SQR(ayt[1]-ayt[jy]) + SQR(azt[1]-azt[jz]);
				dist1 = sqrt(dist2);
				gbartx(isp,dist1,req1,&gtt[jxyz][isp],&gbartt[jxyz][isp]);
			}
		}
		gttsum1[isp] = 0.;
		gbarttsum1[isp] = 0.;
		for(itp=1; itp<=nnt; itp++){
			ix = tisspoints[1][itp];
			iy = tisspoints[2][itp];
			iz = tisspoints[3][itp];
			gbarttsum[itp][isp] = 0.;
			gttsum[itp][isp] = 0.;
			for(jtp=1; jtp<=nnt; jtp++){
				jx = tisspoints[1][jtp];
				jy = tisspoints[2][jtp];
				jz = tisspoints[3][jtp];
				jxyz = 1 + abs(ix - jx) + abs(iy - jy)*mxx + abs(iz - jz)*mxx*myy;
				gbarttsum[itp][isp] += gbartt[jxyz][isp];
				gttsum[itp][isp] += gtt[jxyz][isp];
			}
			gttsum1[isp] += gttsum[itp][isp];
			gbarttsum1[isp] += gbarttsum[itp][isp];
		}
// vessel ~ vessel matrix gbarvv
		if(permsolute[isp]){
			printf("gbarvv...");
			for(i=1; i<=nnv; i++){
				iseg = mainseg[i];
				gvtsegdiv = ds[iseg]/lsegmax + 1.;
				req1 = pow(SQR(rseg[iseg])*ds[iseg]/gvtsegdiv*0.75,0.333333);	//based on vessel segment volume - Jan. 2015
				for(j=1; j<=nnv; j++){
					jseg = mainseg[j];
					gvtsegdiv1 = ds[jseg]/lsegmax + 1.;
					req1 = FMAX(pow(SQR(rseg[jseg])*ds[jseg]/gvtsegdiv1*0.75,0.333333),req1);	//based on other vessel segment volume
					if(i == j){	//self-interaction
						gbarvv[i][j][isp] = 0.;
						for(it=0; it<=nit; it++){	//integration using trapezium rule
							if(it == 0) integrand = 0.5/(2.*pi1*SQR(rseg[iseg])*ds[iseg]);
							else{
								tint = it*deltat/nit;
								fact3 = 4.*diff[isp]*tint;
								fact4 = 1./(sqrt(fact3*pi1)*fact3*pi1);
								b2fact = SQR(rseg[iseg])/fact3;
								c2fact = SQR(ds[iseg])/fact3;
								integrand = fact4/2./b2fact*(1. - expBessi0(2.*b2fact))
									/c2fact*(exp(-c2fact) - 1. + sqrt(pi1*c2fact)*(1. - erfcc(sqrt(c2fact))));
								if(it == nit) integrand *= 0.5;
							}
							gbarvv[i][j][isp] += integrand*deltat/nit;
						}					
					}
					else{			//subdivide each vessel segment in gvtsegdiv pieces - April 2015
						gbarvv[i][j][isp] = 0.;
						for(k=1; k<=gvtsegdiv; k++){
							for(j1=1; j1<=3; j1++) x[j1] = ax[j1][i] + scos[j1][iseg]*ds[iseg]*((k - 0.5)/gvtsegdiv - 0.5);
							for(k1=1; k1<=gvtsegdiv1; k1++){
								for(j1=1; j1<=3; j1++) y[j1] = ax[j1][j] + scos[j1][jseg]*ds[jseg]*((k1 - 0.5)/gvtsegdiv1 - 0.5);
								dist2 = SQR(x[1] - y[1]) + SQR(x[2] - y[2])	+ SQR(x[3] - y[3]);
								dist1 = sqrt(dist2);
								gbartx(isp,dist1,req1,&dummy,&gbarvv0);	//use tissue interaction for different vessels
								gbarvv[i][j][isp] += gbarvv0/gvtsegdiv/gvtsegdiv1;
							}
						}
					}
				}
			}
// tissue ~ vessel matrices gvt and gbarvt: subdivide each vessel segment in gvtsegdiv pieces - April 2015
			printf("gvt,gbarvt...");
			for(i=1; i<=nnv; i++){
				iseg = mainseg[i];
				gvtsegdiv = ds[iseg]/lsegmax + 1.;
				req1 = pow(SQR(rseg[iseg])*ds[iseg]/gvtsegdiv*0.75,0.333333);
				req1 = FMAX(req,req1);	//based on larger of tissue volume and segment volume - Jan 2015
				gbarvtsum[i][isp] = 0.;
				for(itp=1; itp<=nnt; itp++){
					gvt[i][itp][isp] = 0.;
					gbarvt[i][itp][isp] = 0.;
					for(k=1; k<=gvtsegdiv; k++){
						for(j=1; j<=3; j++) x[j] = ax[j][i] + scos[j][iseg]*ds[iseg]*((k - 0.5)/gvtsegdiv - 0.5);
						dist2 = SQR(x[1] - axt[tisspoints[1][itp]])
							+ SQR(x[2] - ayt[tisspoints[2][itp]])
							+ SQR(x[3] - azt[tisspoints[3][itp]]);
						dist1 = sqrt(dist2);
						gbartx(isp,dist1,req1,&gvt0,&gbarvt0);
						gvt[i][itp][isp] += gvt0/gvtsegdiv;
						gbarvt[i][itp][isp] += gbarvt0/gvtsegdiv;
					}
					gbarvtsum[i][isp] += gbarvt[i][itp][isp];
				}
			}
		}
	}
// vessel ~ contour matrix gbarvc, tissue ~ contour matrices gtc and gbartc
// this is needed only if contour method 2 is used, otherwise disable to save time
	if(contourmethod == 2){
		xmin = 0.;
		xmax = sqrt(SQR(xsl1[1]-xsl0[1]) + SQR(xsl1[2]-xsl0[2]) + SQR(xsl1[3]-xsl0[3]));
		ymin = 0.;
		ymax = sqrt(SQR(xsl2[1]-xsl0[1]) + SQR(xsl2[2]-xsl0[2]) + SQR(xsl2[3]-xsl0[3]));
		for(i=1; i<=3; i++){	//set up matrix of direction cosines
			cos[1][i] = (xsl1[i]-xsl0[i])/xmax;
			cos[2][i] = (xsl2[i]-xsl0[i])/ymax;
		}
		cos[3][1] = cos[1][2]*cos[2][3] - cos[1][3]*cos[2][2];
		cos[3][2] = cos[1][3]*cos[2][1] - cos[1][1]*cos[2][3];
		cos[3][3] = cos[1][1]*cos[2][2] - cos[1][2]*cos[2][1];
		for(isp=1; isp<=nsp; isp++) if(diffsolute[isp]){
			printf("gbarvc,gtc,gbartc...");
			for(isl1=1; isl1<=nsl1; isl1++)	for(isl2=1; isl2<=nsl2; isl2++){
				for(i=1; i<=3; i++) x[i] = xsl0[i] + (isl1-1)*(xsl1[i]-xsl0[i])/(nsl1-1) + (isl2-1)*(xsl2[i]-xsl0[i])/(nsl2-1);
				isl12 = isl1 + (isl2 - 1)*nsl1;
				for(i=1; i<=nnv; i++){	//add contributions from vessel sources.  Subdivide subsegments.
					gbarvc[isl12][i][isp] = 0.;
					iseg = mainseg[i];
					for(k=1; k<=slsegdiv; k++){
						for(j=1; j<=3; j++)	y[j] = ax[j][i] + scos[j][iseg]*ds[iseg]*(-0.5 + (k-0.5)/slsegdiv);
						dist2 = SQR(x[1] - y[1]) + SQR(x[2] - y[2]) + SQR(x[3] - y[3]);
						dist1 = sqrt(dist2);
						gbartx(isp,dist1,req1,&dummy,&gbarvc1);
						gbarvc[isl12][i][isp] += gbarvc1/slsegdiv;
					}
				}
				for(itp=1; itp<=nnt; itp++){	//add contributions from tissue sources
					dist2 = SQR(x[1] - axt[tisspoints[1][itp]])
						  + SQR(x[2] - ayt[tisspoints[2][itp]])
						  + SQR(x[3] - azt[tisspoints[3][itp]]);
					dist1 = sqrt(dist2);
					gbartx(isp,dist1,req1,&gtc[isl12][itp][isp],&gbartc[isl12][itp][isp]);
				}
			}
		}
	}
	printf("done\n");
}
