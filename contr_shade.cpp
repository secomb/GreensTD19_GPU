/*******************************************************************
contr_shade - generates contour lines with colored shading between contours
Variables - 
m n nl:  dimensions of array, no. of contour levels
scalefac: determines size of plot
xmin xmax ymin ymax:  boundaries of box
cl(nl):  array of contour levels
zv(m,n):  array of heights
Output to a postscript file.
Version including color shading of region, TWS February 2009.  Updated June 2009.
Rewritten November 2010 to combine polygons with the same color before writing to
the PostScript file.  This yields much smaller files.  They key is that each region
is always enlcosed within a single closed curve, even if it contains 'holes', so
that it is colored correctly.  When polygons are combined, only consecutive sides
merged, so that the curve is not split.  Annular regions contain a cut.
s-polygon: polygon to be shaded within a given rectangle ('square')
r-polygon: polygon to be shaded consisting of multiple s-polygons
flagr - region found
flags - square found
flaga - adjacent square found
Outline of method:
____________________________________________
label entire region with lowest color.
subdivide region into squares
for contour k
	label all squares type 0
	set number of r-polygons nrp = 0
	do
		set flagr = 0
		search for a square of type 0 with >= 1 corner above contour level
		if found
			increment nrp
			do
				set flags = 0
				search for a square of type -nrp 
				if found
					do
						set flaga = 0
						search current then neighboring squares for a square of type = 0
							with >= 1 corner in common with current square above contour level
						if found
							set flagr = 1
							set flags = 1
							set flaga = 1
							define s-polygon for found square
							combine s-polygon with existing r-polygon
							label found square type -nrp
						else label current square type nrp
					while flaga = 1
				endif
			while flags = 1
		endif
	while flagr = 1
end do
*********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl,float xmin,
		   float xmax, float ymin, float ymax, float *cl, float **zv,
		   int showscale, int lowcolor, int hatch, int plotcontour)
{
	fprintf(ofp, "%%!PS\n");
	int i,j,k,iwsp,flagr,flags,flaga,nrp,nspoly,nrpoly,irpoly,ispoly;
	int iin,iint,jjn,jjnt,iseed,jseed,match1,match2,dnrpoly,ispoly1,irpoly1,iipoly;
	int irpolystart,irpolyend,ispolystart,ispolyend;
	int in,in1,in2,in3,in4,inh;
	int *ii,*jj,**corners,**sqtype;
	const int ncmax=100000,nrpolym=10000;
	double xz,xv12,yv12,xv23,yv23,xv34,yv34,xv41,yv41;
	double cx,cy,cbx,cby,cbbox,polytolx,polytoly;
	double *xv,*yv,*red,*blue,*green,*dd,*xvv,*yvv,*xspoly,*yspoly,*xrpoly,*yrpoly;
	double **wsp;

	cx = 50; //origin of contour plot
	cy = 50;//origin of contour plot
	cbbox = 15.; //size of boxes
	cbx = 550; //origin of color bar
	cby = 100;//origin of color bar
	if(showscale == 2){	//horizontal color bar
		cy = 100;
		cbx = 50; //origin of color bar
		cby = 50; //origin of color bar
	}
	polytolx = 1.e-6*(xmax - xmin)/(m-1);
	polytoly = 1.e-6*(ymax - ymin)/(n-1);



	xv = dvector(1,m);
	yv = dvector(1,n);
	xvv = dvector(1,4);
	yvv = dvector(1,4);
	red = dvector(0,nl);
	green = dvector(0,nl);
	blue = dvector(0,nl);
	wsp = dmatrix(1,ncmax,1,4);
	corners = imatrix(1,4,1,2);
	sqtype = imatrix(1,m-1,1,n-1);
	xspoly = dvector(1,6);
	yspoly = dvector(1,6);
	xrpoly = dvector(1,nrpolym);
	yrpoly = dvector(1,nrpolym);
	dd = dvector(1,4);
	ii = ivector(1,4);
	jj = ivector(1,4);
	corners[1][1] = 0;
	corners[1][2] = 0;
	corners[2][1] = 1;
	corners[2][2] = 0;
	corners[3][1] = 1;
	corners[3][2] = 1;
	corners[4][1] = 0;
	corners[4][2] = 1;

	for(i=1; i<=m; i++) xv[i] = xmin + (i - 1)*(xmax - xmin)/(m - 1);
	for(j=1; j<=n; j++) yv[j] = ymin + (j - 1)*(ymax - ymin)/(n - 1);
	fprintf(ofp, "/mx {%g sub %g mul %g add} def\n",xmin,scalefac,cx);
	fprintf(ofp, "/my {%g sub %g mul %g add} def\n",ymin,scalefac,cy);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/n {newpath} def\n");
	fprintf(ofp, "/s {stroke} def\n");
	fprintf(ofp, "/cf {closepath fill} def\n");
	fprintf(ofp, "/cs {closepath stroke} def\n");
	fprintf(ofp, "0.5 setlinewidth\n");
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "8 scalefont\n");
	fprintf(ofp, "setfont\n");
	if(hatch > 0){
		fprintf(ofp, "/diagonals1\n");
		fprintf(ofp, "{ newpath\n");
		fprintf(ofp, "-700 10 550\n");
		fprintf(ofp, "{ 0 moveto\n");
		fprintf(ofp, "700 700 rlineto } for\n");
		fprintf(ofp, "stroke } def\n");

		fprintf(ofp, "/diagonals2\n");
		fprintf(ofp, "{ newpath\n");
		fprintf(ofp, "-700 10 550\n");
		fprintf(ofp, "{ 700 moveto\n");
		fprintf(ofp, "700 -700 rlineto } for\n");
		fprintf(ofp, "stroke } def\n");

		fprintf(ofp, "/verticals\n");
		fprintf(ofp, "{ newpath\n");
		fprintf(ofp, "50 10 550\n");
		fprintf(ofp, "{ 0 moveto\n");
		fprintf(ofp, "0 700 rlineto } for\n");
		fprintf(ofp, "stroke } def\n");

		fprintf(ofp, "/horizontals\n");
		fprintf(ofp, "{ newpath\n");
		fprintf(ofp, "100 10 800\n");
		fprintf(ofp, "{ 0 exch moveto\n");
		fprintf(ofp, "500 0 rlineto } for\n");
		fprintf(ofp, "stroke } def\n");
	}
//Set up colors using Matlab 'jet' scheme
	for(k=0; k<=nl; k++){
		xz = float(k)/float(nl);
		blue[k] = FMIN(FMAX(1.5-4*fabs(xz-0.25), 0.), 1.);
		green[k]= FMIN(FMAX(1.5-4*fabs(xz-0.5), 0.), 1.);
		red[k]  = FMIN(FMAX(1.5-4*fabs(xz-0.75), 0.), 1.);
	}
//Color whole region with lowest color
	if(lowcolor > 0){
		fprintf(ofp, "%f %f %f setrgbcolor\n",red[0],green[0],blue[0]);
		fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cf\n",
			xmin,ymin,xmax,ymin,xmax,ymax,xmin,ymax);
	}
//Analyze each rectangle separately. Overwrite lower colors
	iwsp = 0;
	for(k=1; k<=nl; k++){
		fprintf(ofp, "%f %f %f setrgbcolor\n",red[k],green[k],blue[k]);
		for(i=1; i<m; i++) for(j=1; j<n; j++) sqtype[i][j] = 0;
		nrp = 0;
		do{
//search for a square of type 0 with >= 1 corner above contour level
			flagr = 1;
			for(i=1; i<m; i++) for(j=1; j<n; j++) if(sqtype[i][j] == 0) for(in=1; in<=4; in++)
				if(zv[i + corners[in][1]][j + corners[in][2]] > cl[k]) goto foundit;
			flagr = 0;
			foundit:;			
			if(flagr == 1){
				nrp++;
				iseed = i;
				jseed = j;		//'seed' region nrp
				nrpoly = 0;		//initialize r-polygon
				do{
//search for a square of type -nrp
					flags = 1;
					if(i == iseed && j == jseed && sqtype[i][j] == 0) goto foundit1;
					for(i=1; i<m; i++) for(j=1; j<n; j++) if(sqtype[i][j] == -nrp) goto foundit1;
					flags = 0;
					foundit1:;
					if(flags == 1){
						do{
//search current then neighboring squares for square type 0, with a corner in common with current square above contour level
							flaga = 1;
							for(inh=0; inh<=4; inh++){
								iin = i;
								jjn = j;
								if(inh == 1) iin++;
								if(inh == 2) jjn++;
								if(inh == 3) iin--;
								if(inh == 4) jjn--;
								if(iin > 0 && iin < m && jjn > 0 && jjn < n) if(sqtype[iin][jjn] == 0){
									for(in=1; in<=4; in++){
										iint = iin + corners[in][1];
										jjnt = jjn + corners[in][2];
										if((iint == i || iint == i+1) && (jjnt == j || jjnt == j+1) 
											&& (zv[iint][jjnt] > cl[k])) goto foundit2;
									}
								}
							}
							flaga = 0;
							foundit2:;
							if(flaga == 1){
								sqtype[iin][jjn] = -nrp;
//define s-polygon for found square (iin,jjn) - may have up to 6 sides
								in1 = in;				//this is always in the region
								in2 = in1%4 + 1;
								in3 = in2%4 + 1;
								in4 = in3%4 + 1;
								for(in=1; in<=4; in++){
									ii[in] = iin + corners[in][1];
									jj[in] = jjn + corners[in][2];
									dd[in] = zv[ii[in]][jj[in]] - cl[k];
									xvv[in] = xv[ii[in]];
									yvv[in] = yv[jj[in]];
								}
								if(dd[in1] != dd[in2]){
									xv12 = (dd[in1]*xv[ii[in2]]-dd[in2]*xv[ii[in1]])/(dd[in1]-dd[in2]);
									yv12 = (dd[in1]*yv[jj[in2]]-dd[in2]*yv[jj[in1]])/(dd[in1]-dd[in2]);
								}
								if(dd[in2] != dd[in3]){
									xv23 = (dd[in2]*xv[ii[in3]]-dd[in3]*xv[ii[in2]])/(dd[in2]-dd[in3]);
									yv23 = (dd[in2]*yv[jj[in3]]-dd[in3]*yv[jj[in2]])/(dd[in2]-dd[in3]);
								}
								if(dd[in3] != dd[in4]){
									xv34 = (dd[in3]*xv[ii[in4]]-dd[in4]*xv[ii[in3]])/(dd[in3]-dd[in4]);
									yv34 = (dd[in3]*yv[jj[in4]]-dd[in4]*yv[jj[in3]])/(dd[in3]-dd[in4]);
								}
								if(dd[in4] != dd[in1]){
									xv41 = (dd[in4]*xv[ii[in1]]-dd[in1]*xv[ii[in4]])/(dd[in4]-dd[in1]);
									yv41 = (dd[in4]*yv[jj[in1]]-dd[in1]*yv[jj[in4]])/(dd[in4]-dd[in1]);
								}
								xspoly[1] = xvv[in1];
								yspoly[1] = yvv[in1];
								if(dd[in2] > 0){						//corners 1,2 are this color
									xspoly[2] = xvv[in2];
									yspoly[2] = yvv[in2];
									if(dd[in3] > 0){					//corners 1,2,3 are this color
										xspoly[3] = xvv[in3];
										yspoly[3] = yvv[in3];
										if(dd[in4] > 0){				//corners 1,2,3,4 are this color
											xspoly[4] = xvv[in4];
											yspoly[4] = yvv[in4];
											nspoly = 4;
										}
										else{							//corners 1,2,3,not 4 are this color
											xspoly[4] = xv34;
											yspoly[4] = yv34;										
											xspoly[5] = xv41;
											yspoly[5] = yv41;
											nspoly = 5;
											iwsp++;
											wsp[iwsp][1] = xv34;
											wsp[iwsp][2] = yv34;
											wsp[iwsp][3] = xv41;
											wsp[iwsp][4] = yv41;
										}
									}
									else{								//corners 1,2,not 3 are this color
										xspoly[3] = xv23;
										yspoly[3] = yv23;
										iwsp++;
										wsp[iwsp][1] = xv23;
										wsp[iwsp][2] = yv23;
										if(dd[in4] > 0){				//corners 1,2,not 3,4 are this color
											xspoly[4] = xv34;
											yspoly[4] = yv34;										
											xspoly[5] = xvv[in4];
											yspoly[5] = yvv[in4];
											nspoly = 5;
											wsp[iwsp][3] = xv34;
											wsp[iwsp][4] = yv34;
										}
										else{							//corners 1,2,not 3,not 4 are this color
											xspoly[4] = xv41;
											yspoly[4] = yv41;
											nspoly = 4;
											wsp[iwsp][3] = xv41;
											wsp[iwsp][4] = yv41;
										}
									}
								}
								else{									//corners 1,not 2 are this color
									xspoly[2] = xv12;
									yspoly[2] = yv12;
									iwsp++;
									wsp[iwsp][1] = xv12;
									wsp[iwsp][2] = yv12;
									if(dd[in3] > 0){					//corners 1,not 2,3 are this color
										xspoly[3] = xv23;
										yspoly[3] = yv23;
										xspoly[4] = xvv[in3];
										yspoly[4] = yvv[in3];
										wsp[iwsp][3] = xv23;
										wsp[iwsp][4] = yv23;
										if(dd[in4] > 0){				//corners 1,not 2,3,4 are this color
											xspoly[5] = xvv[in4];
											yspoly[5] = yvv[in4];
											nspoly = 5;
										}
										else{							//corners 1,not 2,3,not 4 are this color
											xspoly[5] = xv34;
											yspoly[5] = yv34;											
											xspoly[6] = xv41;
											yspoly[6] = yv41;
											nspoly = 6;
											iwsp++;
											wsp[iwsp][1] = xv34;
											wsp[iwsp][2] = yv34;
											wsp[iwsp][3] = xv41;
											wsp[iwsp][4] = yv41;
										}
									}
									else{								//corners 1,not 2,not 3 are this color
										if(dd[in4] > 0){				//corners 1,not 2,not 3,4 are this color
											xspoly[3] = xv34;
											yspoly[3] = yv34;											
											xspoly[4] = xvv[in4];
											yspoly[4] = yvv[in4];
											nspoly = 4;
											wsp[iwsp][3] = xv34;
											wsp[iwsp][4] = yv34;
										}
										else{							//corners 1,not 2,not 3,not 4 are this color
											xspoly[3] = xv41;
											yspoly[3] = yv41;
											nspoly = 3;
											wsp[iwsp][3] = xv41;
											wsp[iwsp][4] = yv41;
										}
									}
								}
								if(iwsp > ncmax) printf("*** Error: too many contour points.  Increase ncmax ***\n");
//combine s-polygon with existing r-polygon, eliminating redundant segments
								if(nrpoly == 0){						//initiate r-polygon
									for(ispoly=1; ispoly<=nspoly; ispoly++){
										xrpoly[ispoly] = xspoly[ispoly];
										yrpoly[ispoly] = yspoly[ispoly];
									}
									nrpoly = nspoly;
								}
								else{					//search r-polygon and s-polygon for one side that matches
									for(irpoly=nrpoly; irpoly>=1; irpoly--) for(ispoly=1; ispoly<=nspoly; ispoly++){
										ispoly1 = ispoly%nspoly + 1;
										irpoly1 = irpoly%nrpoly + 1;
										if((fabs(xrpoly[irpoly] - xspoly[ispoly1]) < polytolx)
											&& (fabs(yrpoly[irpoly] - yspoly[ispoly1]) < polytoly)
											&& (fabs(xrpoly[irpoly1] - xspoly[ispoly]) < polytolx)
											&& (fabs(yrpoly[irpoly1] - yspoly[ispoly]) < polytoly)) goto foundit3;
									}
									printf("*** Error: matching segment not found ***\n");
									foundit3:;
									match1 = 0;
									for(iipoly=2; iipoly<nspoly; iipoly++){	//search for further matches nearby 
										irpoly1 = (irpoly - iipoly + nrpoly)%nrpoly + 1;
										ispoly1 = (ispoly + iipoly - 1)%nspoly + 1;
										if((fabs(xrpoly[irpoly1] - xspoly[ispoly1]) < polytolx)
											&& (fabs(yrpoly[irpoly1] - yspoly[ispoly1]) < polytoly)) match1++;
										else goto nomatch1;
									}
									nomatch1:;
									match2 = 0;
									for(iipoly=2; iipoly<nspoly; iipoly++){	//search other way for further matches nearby 
										irpoly1 = (irpoly + iipoly - 1)%nrpoly + 1;
										ispoly1 = (ispoly - iipoly + nspoly)%nspoly + 1;
										if((fabs(xrpoly[irpoly1] - xspoly[ispoly1]) < polytolx)
											&& (fabs(yrpoly[irpoly1] - yspoly[ispoly1]) < polytoly)) match2++;
										else goto nomatch2;
									}
									nomatch2:;
									ispolystart = (ispoly + match1)%nspoly + 1;				//first node of s-polygon to include
									ispolyend = (ispoly - match2 - 1 + nspoly)%nspoly + 1;	//last node of s-polygon to include
									irpolystart = (irpoly + match2)%nrpoly + 1;				//first node of s-polygon to include
									irpolyend = (irpoly - match1 - 1 + nrpoly)%nrpoly + 1;	//last node of s-polygon to include
									if((fabs(xrpoly[irpolystart] - xspoly[ispolyend]) > polytolx)
										|| (fabs(yrpoly[irpolystart] - yspoly[ispolyend]) > polytoly))
										printf("*** Error: r and s nodes do not match ****\n");
									if((fabs(xrpoly[irpolyend] - xspoly[ispolystart]) > polytolx)
										|| (fabs(yrpoly[irpolyend] - yspoly[ispolystart]) > polytoly))
										printf("*** Error: r and s nodes do not match ****\n");
									dnrpoly = nspoly - 2 - 2*match1 - 2*match2;
									if(irpolystart > irpolyend){
										if(dnrpoly > 0){		//expand the arrays xrpoly and yrpoly
											for(iipoly=nrpoly; iipoly>=irpolystart; iipoly--){
												xrpoly[iipoly+dnrpoly] = xrpoly[iipoly];
												yrpoly[iipoly+dnrpoly] = yrpoly[iipoly];
											}
										}
										if(dnrpoly < 0){		//contract the arrays xrpoly and yrpoly
											for(iipoly=irpolystart; iipoly<=nrpoly; iipoly++){
												xrpoly[iipoly+dnrpoly] = xrpoly[iipoly];
												yrpoly[iipoly+dnrpoly] = yrpoly[iipoly];
											}
										}
									}
									nrpoly += dnrpoly;
									if(nrpoly > nrpolym) printf("*** Error: too many polygon points.  Increase nrpolym ***\n");
									if(nrpoly < irpolyend) for(iipoly=irpolyend; iipoly>nrpoly; iipoly--){	//otherwise these values get lost!
										xrpoly[iipoly-nrpoly] = xrpoly[iipoly];
										yrpoly[iipoly-nrpoly] = yrpoly[iipoly];
									}
									for(iipoly=1; iipoly<=nspoly-2-match1-match2; iipoly++){
										irpoly1 = (irpolyend + iipoly - 1)%nrpoly + 1;
										ispoly1 = (ispolystart + iipoly - 1)%nspoly + 1;
										xrpoly[irpoly1] = xspoly[ispoly1];
										yrpoly[irpoly1] = yspoly[ispoly1];
									}
								}
							}
							else sqtype[i][j] = nrp;
						}
						while (flaga == 1);
					}
				}
				while (flags == 1);
//display polygons after combining
				if(hatch > 0) fprintf(ofp, "/clippingpath {\n");
				fprintf(ofp, "n %g mx %g my m ",xrpoly[1],yrpoly[1]);
				for(irpoly=2; irpoly<=nrpoly; irpoly++){
					fprintf(ofp, "%g mx %g my l ",xrpoly[irpoly],yrpoly[irpoly]);
					if(irpoly%5 == 1) fprintf(ofp, "\n");
				}
				if(hatch > 0){
					fprintf(ofp, "} def\n");
					fprintf(ofp, "0 0 0 setrgbcolor\n");//black
					fprintf(ofp, "gsave clippingpath clip\n");
					if(hatch == 1) fprintf(ofp, "diagonals1 grestore\n");
					if(hatch == 2) fprintf(ofp, "diagonals2 grestore\n");
				}
				else fprintf(ofp, "cf\n");
			}
		}
		while (flagr == 1);
	}
//outline contours
	if(plotcontour > 0){
		fprintf(ofp, "0 0 0 setrgbcolor\n");//black
		for(in=1; in<=iwsp; in++) fprintf(ofp, "n %g mx %g my m %g mx %g my l s\n",
				wsp[in][1],wsp[in][2],wsp[in][3],wsp[in][4]);
	}
//draw a box
	fprintf(ofp, "0 0 0 setrgbcolor\n");//black
	fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cs\n",
		xmin,ymin,xmax,ymin,xmax,ymax,xmin,ymax);
//create a vertical color bar
	if(showscale == 1){
		for(k=0; k<=nl; k++){
			fprintf(ofp, "%f %f %f setrgbcolor\n",red[k],green[k],blue[k]);
			fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
				cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
			if(k>0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),cl[k]);
		}
		fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
			cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*(nl+1),cbx,cby+cbbox*(nl+1));
	}
//create a horizontal color bar
	if(showscale == 2){
		for(k=0; k<=nl; k++){
			fprintf(ofp, "%f %f %f setrgbcolor\n",red[k],green[k],blue[k]);
			fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
				cbx+k*cbbox,cby,cbx+k*cbbox,cby+cbbox,cbx+(k+1)*cbbox,cby+cbbox,cbx+(k+1)*cbbox,cby);
			if(k>0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor 90 rotate (%g) show -90 rotate\n",cbx+cbbox*(k+0.2),cby + cbbox*1.2,cl[k]);
		}
		fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
			cbx,cby,cbx+cbbox*(nl+1),cby,cbx+cbbox*(nl+1),cby+cbbox,cbx,cby+cbbox);
	}
	free_dvector(xv,1,m);
	free_dvector(yv,1,n);
	free_dvector(xvv,1,4);
	free_dvector(yvv,1,4);
	free_dvector(red,0,nl);
	free_dvector(green,0,nl);
	free_dvector(blue,0,nl);
	free_dmatrix(wsp,1,ncmax,1,4);
	free_imatrix(corners,1,4,1,2);
	free_imatrix(sqtype,1,m-1,1,n-1);
	free_dvector(xspoly,1,6);
	free_dvector(yspoly,1,6);
	free_dvector(xrpoly,1,nrpolym);
	free_dvector(yrpoly,1,nrpolym);
	free_dvector(dd,1,4);
	free_ivector(ii,1,4);
	free_ivector(jj,1,4);
}