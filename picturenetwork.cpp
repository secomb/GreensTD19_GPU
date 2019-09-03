/**********************************************************
picturenetwork.cpp - project network on z = 0 plane.  TWS Dec. 07
Uses parameters from CountourParams.dat
Labels nodes with nodevar and segments with segvar (must be float).
Generates a postscript file.
Version 2.0, May 1, 2010.  Added abbrevations of setlinewidth and stroke, May 09.
Added visualization of tissue points, April 2010.
Version 3.0, May 17, 2011.  Plots in z-order, for better results with 3D networks.
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void picturenetwork(float *nodvar, float *segvar, const char fname[])
{
	extern int max,mxx,myy,mzz,nseg,nnod;
	extern int *segtyp,*ista,*iend;
	extern int ***nbou;//added April 2010
	extern float *axt,*ayt;//added April 2010
	extern float *diam,**cnode,*xsl0,*xsl1,*xsl2;;
	int i,j,k,iseg,inod,showpoint,ilevel,nlevel = 100;
	float xmin,xmax,ymin,ymax,xs,ys,picfac,red,green,blue,xz,xzmin,xzmax;
	float diamfac = 1.,zcoord,zbottom,ztop,zmin,zmax;
	float **cos;
	
	FILE *ofp;

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

	picfac = FMIN(500./xmax,700./ymax);//updated April 2010
	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS-Adobe-2.0\n");
	fprintf(ofp, "%%%%Pages: 1\n");
	fprintf(ofp, "%%%%EndComments\n");
	fprintf(ofp, "%%%%Page: 1 1\n");
	fprintf(ofp, "/mx {%g mul 50 add} def\n",picfac);
	fprintf(ofp, "/my {%g mul 50 add} def\n",picfac);//updated April 2010
	fprintf(ofp, "/cf {closepath fill} def\n");
	fprintf(ofp, "/cs {closepath stroke} def\n");
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/n {newpath} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/sl {setlinewidth} def\n");
	fprintf(ofp, "/sc {setrgbcolor} def\n");
	fprintf(ofp, "/s {stroke} def\n");
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "8 scalefont\n");
	fprintf(ofp, "setfont\n");

	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n",xmin,ymin);
	fprintf(ofp, "%g mx %g my l\n",xmax,ymin);
	fprintf(ofp, "%g mx %g my l\n",xmax,ymax);
	fprintf(ofp, "%g mx %g my l\n",xmin,ymax);
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");
//show tissue points
	fprintf(ofp,"0 0 0 setrgbcolor\n");//black
	for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++){
		showpoint = 0;
		for(k=1; k<=mzz; k++) if(nbou[i][j][k] > 0) showpoint = 1;
		if(showpoint == 1) fprintf(ofp, "%g mx %g my m (.) show\n", axt[i],ayt[j]);
	}
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "12 scalefont\n");
	fprintf(ofp, "setfont\n");
//plot vessels according to segvar in order from bottom to top according to z-coordinate
	xzmin = 1.e6;
	xzmax = -1.e6;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		xzmin = FMIN(xzmin,segvar[iseg]);
		xzmax = FMAX(xzmax,segvar[iseg]);
	}
	for(ilevel=1; ilevel<=nlevel; ilevel++){
		zbottom = zmin + (ilevel-1)*(zmax - zmin)/nlevel;
		ztop = zmin + ilevel*(zmax - zmin)/nlevel;
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
			zcoord = 0.;
			for(i=1; i<=3; i++)	zcoord += (cnode[i][ista[iseg]] + cnode[i][iend[iseg]])/2.*cos[3][i];
			if(zcoord >= zbottom && zcoord < ztop){
				if(xzmin != xzmax) xz = (segvar[iseg] - xzmin)/(xzmax - xzmin);
				else xz = 0.75;
				blue = FMIN(FMAX(1.5-4.*fabs(xz-0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
				green= FMIN(FMAX(1.5-4.*fabs(xz-0.5), 0.), 1.);
				red  = FMIN(FMAX(1.5-4.*fabs(xz-0.75), 0.), 1.);
				fprintf(ofp,"%f %f %f sc\n",red,green,blue);
				fprintf(ofp,"%g sl\n",picfac*diam[iseg]*diamfac);//line widths scaled up by diamfac
				xs = 0.;
				ys = 0.;
				for(i=1; i<=3; i++){
					xs += (cnode[i][ista[iseg]]-xsl0[i])*cos[1][i];
					ys += (cnode[i][ista[iseg]]-xsl0[i])*cos[2][i];
				}
				fprintf(ofp, "%g mx %g my m ", xs,ys);
				xs = 0.;
				ys = 0.;
				for(i=1; i<=3; i++){
					xs += (cnode[i][iend[iseg]]-xsl0[i])*cos[1][i];
					ys += (cnode[i][iend[iseg]]-xsl0[i])*cos[2][i];
				}
				fprintf(ofp, "%g mx %g my l s \n", xs,ys);
			}
		}
	}

//label nodes in black
	fprintf(ofp,"0 0 0 setrgbcolor\n");//black
	for(inod=1; inod<=nnod; inod++){
		xs = 0.;
		ys = 0.;
		for(i=1; i<=3; i++){
			xs += (cnode[i][inod]-xsl0[i])*cos[1][i];
			ys += (cnode[i][inod]-xsl0[i])*cos[2][i];
		}
//comment out next two lines to remove node numbers
		fprintf(ofp, "%g mx %g my m ", xs + 0.5/picfac,ys);
		fprintf(ofp, "(%g) show\n",nodvar[inod]);
	}
//label segments in blue
	fprintf(ofp,"0 0 1 setrgbcolor\n");//blue
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		xs = 0.;
		ys = 0.;
		for(i=1; i<=3; i++){
			xs += ((cnode[i][ista[iseg]]+cnode[i][iend[iseg]])/2.-xsl0[i])*cos[1][i];
			ys += ((cnode[i][ista[iseg]]+cnode[i][iend[iseg]])/2.-xsl0[i])*cos[2][i];
		}
//comment out next two lines to remove segment numbers
		fprintf(ofp, "%g mx %g my m ", xs + 0.5*picfac,ys);
		fprintf(ofp, "(%g) show\n",segvar[iseg]);
	}
//create a color bar
	float c;
	float cbbox = 15.; //size of boxes
	float cbx = 560; //origin of color bar
	float cby = 100;//origin of color bar
	fprintf(ofp, "0.5 setlinewidth\n");
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "8 scalefont\n");
	fprintf(ofp, "setfont\n");
	for(k=0; k<=10; k++){
		xz = k*0.1;
		c = xzmin + (xzmax - xzmin)*xz;
		blue = FMIN(FMAX(1.5-4.*fabs(xz-0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
		green= FMIN(FMAX(1.5-4.*fabs(xz-0.5), 0.), 1.);
		red  = FMIN(FMAX(1.5-4.*fabs(xz-0.75), 0.), 1.);
		fprintf(ofp, "%f %f %f setrgbcolor\n",red,green,blue);
		fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
			cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
		if(k>0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),c);
	}
	fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
		cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*11,cbx,cby+cbbox*11);
	fprintf(ofp, "showpage\n");
	fclose(ofp);
	free_matrix(cos,1,3,1,3);
}