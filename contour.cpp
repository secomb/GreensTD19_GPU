/**********************************************************
contour.cpp - generate data for contour plot.  TWS Dec. 07
Version 3.0, May 17, 2011.
Produces a single postscript file with a page for each solute
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void contr_lines(FILE *ofp, int m, int n, float scalefac, int nl,
		   float xmin, float xmax, float ymin, float ymax, float *cl, float **zv);
void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl,
		   float xmin, float xmax, float ymin, float ymax, float *cl, float **zv,
		   int showscale, int lowcolor, int hatch, int plotcontour);

void contour(int itime, float time, int method)
 {
	extern int max,nsp,nseg,nnt,nnv,mxx,myy,mzz;
	extern int *segtyp,*nl,*ista,*iend,*oxygen,***nbou;
	extern int slsegdiv,nsl1,nsl2;
	extern float pi1,req,scalefac,vol;
	extern float *x,*p,*diam,*g0,**cnode,**cvseg,*alphat;
	extern float *xsl0,*xsl1,*xsl2,*clmin,*clint,*cl,**zv,***psl,*axt,*ayt,*azt;
	extern float ***gbarvc,***gtc,***gbartc,*g0,*g0prev,**qv,**qt;
	extern float **ct,**ctprev;

	int i,j,k,iseg,itp,isp,isl1,isl2,isl12,ilevel,nlevel = 100,ii,jj,kk;
	float xmin,ymin,xmax,ymax,xs1,ys1,xs2,ys2,**cos;
	float red,green,blue,xz,xzmin,xzmax,lamx,lamy,lamz;
	float diamfac = 1.,zcoord,zbottom,ztop,zmin,zmax;
	char fname[80];
	FILE *ofp;
	//method = 1: use current tissue field and interpolate
	//method = 2: use previous tissue field and source strengths.
	//Method 2 requires computation of gvarvc, gtc and gbartc in diffusionmatrix.cpp.

//create file name - need 3-digit frame number
    sprintf(fname,"Current\\Contour%03i.ps",itime);

	printf("Generating data for contour plots...");

	xmin = 0.;
	xmax = sqrt(SQR(xsl1[1]-xsl0[1]) + SQR(xsl1[2]-xsl0[2]) + SQR(xsl1[3]-xsl0[3]));
	ymin = 0.;
	ymax = sqrt(SQR(xsl2[1]-xsl0[1]) + SQR(xsl2[2]-xsl0[2]) + SQR(xsl2[3]-xsl0[3]));
	cos = matrix(1,3,1,3);
	for(i=1; i<=3; i++){	//set up matrix of direction cosines
		cos[1][i] = (xsl1[i]-xsl0[i])/xmax;
		cos[2][i] = (xsl2[i]-xsl0[i])/ymax;
	}
	cos[3][1] = cos[1][2]*cos[2][3] - cos[1][3]*cos[2][2];
	cos[3][2] = cos[1][3]*cos[2][1] - cos[1][1]*cos[2][3];
	cos[3][3] = cos[1][1]*cos[2][2] - cos[1][2]*cos[2][1];

//Determine range of z values
	zmin = 1.e6;
	zmax = -1.e6;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		zcoord = 0.;
		for(i=1; i<=3; i++)	zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]])/2.*cos[3][i];
		zmin = FMIN(zmin,zcoord-1.);
		zmax = FMAX(zmax,zcoord+1.);
	}
//Calculate P on a planar slice through the region
	for(isl1=1; isl1<=nsl1; isl1++)	for(isl2=1; isl2<=nsl2; isl2++){
		isl12 = isl1 + (isl2 - 1)*nsl1;
		for(isp=1; isp<=nsp; isp++){
			if(method == 1){	//compute using interpolation from current time step
				for(i=1; i<=3; i++) x[i] = xsl0[i] + (isl1-1)*(xsl1[i]-xsl0[i])/(nsl1-1) + (isl2-1)*(xsl2[i]-xsl0[i])/(nsl2-1);
				i = 0;
				while(x[1] > axt[i+1] && i < mxx) i++;
				if(i == 0) lamx = 1.;
				else if(i == mxx) lamx = 0.;
				else lamx = (x[1] - axt[i])/(axt[i+1] - axt[i]);
				j = 0;
				while(x[2] > ayt[j+1] && j < myy) j++;
				if(j == 0) lamy = 1.;
				else if(j == myy) lamy = 0.;
				else lamy = (x[2] - ayt[j])/(ayt[j+1] - ayt[j]);
				k = 0;
				while(x[3] > azt[k+1] && k < mzz) k++;
				if(k == 0) lamz = 1.;
				else if(k == mzz) lamz = 0.;
				else lamz = (x[3] - azt[k])/(azt[k+1] - azt[k]);
				psl[isl1][isl2][isp] = g0[isp];
				for(ii=0; ii<=1; ii++) for(jj=0; jj<=1; jj++) for(kk=0; kk<=1; kk++){
					if(i+ii>=1 && i+ii<=mxx && j+jj>=1 && j+jj<=myy && k+kk>=1 && k+kk<=mzz){
						itp = nbou[i+ii][j+jj][k+kk];
						if(itp != 0) psl[isl1][isl2][isp] += (1 - ii + (2*ii - 1)*lamx)*(1 - jj + (2*jj - 1)*lamy)
							*(1 - kk + (2*kk - 1)*lamz)*(ct[itp][isp] - g0[isp]);
					}
				}
			}
			else{		//compute using greens functions from previous time step
				psl[isl1][isl2][isp] = g0[isp];	//initialize to g0
				for(i=1; i<=nnv; i++) psl[isl1][isl2][isp] += gbarvc[isl12][i][isp]*qv[i][isp];
				for(itp=1; itp<=nnt; itp++) psl[isl1][isl2][isp] += vol*gtc[isl12][itp][isp]*(ctprev[itp][isp] - g0prev[isp])
					+ gbartc[isl12][itp][isp]*qt[itp][isp];		
			}
			if(oxygen[isp]) psl[isl1][isl2][isp] /= alphat[isp];
 		}
	}
	xmin = 0.;
	ymin = 0.;
	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS-Adobe-2.0\n");
	fprintf(ofp, "%%%%Pages: %i\n",nsp);
	fprintf(ofp, "%%%%EndComments\n");
	for(isp=1; isp<=nsp; isp++){
		fprintf(ofp, "%%%%Page: %i %i\n",isp,isp);
		for(isl1=1; isl1<=nsl1; isl1++)	for(isl2=1; isl2<=nsl2; isl2++) zv[isl1][isl2] = psl[isl1][isl2][isp];
		for(i=1; i<=nl[isp]; i++) cl[i] = clmin[isp] + (i-1)*clint[isp];
		contr_shade(ofp,nsl1,nsl2,scalefac,nl[isp],xmin,xmax,ymin,ymax,cl,zv,1,1,0,0);
//		contr_lines(ofp,nsl1,nsl2,scalefac,nl[isp],xmin,xmax,ymin,ymax,cl,zv);

		fprintf(ofp, "/sl {setlinewidth} def\n");
		fprintf(ofp, "/sc {setrgbcolor} def\n");
		fprintf(ofp, "/s {stroke} def\n");
		fprintf(ofp, "1 setlinecap\n");

//Plot projection of network in contour plane
//plot vessels according to cvseg in order from bottom to top according to z-coordinate
		xzmin = clmin[isp];
		xzmax = clmin[isp] + nl[isp]*clint[isp];
		for(ilevel=1; ilevel<=nlevel; ilevel++){
			zbottom = zmin + (ilevel-1)*(zmax - zmin)/nlevel;
			ztop = zmin + ilevel*(zmax - zmin)/nlevel;
			for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
				zcoord = 0.;
				for(i=1; i<=3; i++)	zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]])/2.*cos[3][i];
				if(zcoord >= zbottom && zcoord < ztop){
					if(xzmin != xzmax) xz = (cvseg[iseg][isp] - xzmin)/(xzmax - xzmin);
					else xz = 0.75;
					xz = FMIN(FMAX(xz,0.),1.);	//restrict xz to [0,1] - added 2/15
					blue = FMIN(FMAX(1.5-4.*fabs(xz-0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
					green= FMIN(FMAX(1.5-4.*fabs(xz-0.5), 0.), 1.);
					red  = FMIN(FMAX(1.5-4.*fabs(xz-0.75), 0.), 1.);

					xs1 = 0.;
					ys1 = 0.;
					xs2 = 0.;
					ys2 = 0.;
					for(i=1; i<=3; i++){
						xs1 += (cnode[i][ista[iseg]]-xsl0[i])*cos[1][i];
						ys1 += (cnode[i][ista[iseg]]-xsl0[i])*cos[2][i];
						xs2 += (cnode[i][iend[iseg]]-xsl0[i])*cos[1][i];
						ys2 += (cnode[i][iend[iseg]]-xsl0[i])*cos[2][i];
					}
					fprintf(ofp,"0 0 0 sc\n");	//Plot vessels slightly larger in black to outline
					fprintf(ofp,"%g sl\n",scalefac*diam[iseg]*diamfac+2.);
					fprintf(ofp, "%g mx %g my m %g mx %g my l s \n", xs1,ys1,xs2,ys2);
					fprintf(ofp,"%f %f %f sc\n",red,green,blue);
					fprintf(ofp,"%g sl\n",scalefac*diam[iseg]*diamfac);//line widths scaled up by diamfac
					fprintf(ofp, "%g mx %g my m %g mx %g my l s \n", xs1,ys1,xs2,ys2);
				}
			}
		}
		fprintf(ofp,"0 0 0 setrgbcolor\n");		//black
		fprintf(ofp, "/Times-Roman findfont 12 scalefont setfont\n");
		fprintf(ofp, "50 30 moveto\n");
		fprintf(ofp, "(Solute %i) show\n",isp);  //show solute number
		fprintf(ofp, "50 10 moveto\n");
		fprintf(ofp, "(Time = %g) show\n",time);  //show solute number
//create a scale bar
		float barlength = 50;
		if(xmax > 250.) barlength = 100.;
		if(xmax > 500.) barlength = 200.;
		if(xmax > 1500.) barlength = 500.;
		fprintf(ofp,"%g sl\n",2.);
		fprintf(ofp,"%g %g m %g %g l stroke\n",120.,25.,120.+scalefac*barlength,25.);
		fprintf(ofp,"%g %g m (%g mm) show\n",120.,30.,barlength/1000.);
		fprintf(ofp, "showpage\n");
	}
	fclose(ofp);
	printf("done\n");
}