/****************************************************
outboun.cpp
----------------------------------------------------
method = 1: finds the smallest convex region inside the cuboid
which excludes tissue node points that have a distance to the nearest
vessel greater than a value specified by the user (lb).  Any point
outside a region between two planes containing the required points is
excluded.  This is repeated for multiple plane orientations, determined by am.
----------------------------------------------------
method = 2: finds all tissue points within a distance lb of the vessels, but
does not make a convex region.  Fills in 'holes' in interior
----------------------------------------------------
Output:	nnt, total tissue node points inside the region.
nbou > 1 if inside region, value gives tissue point index
TWS 2010
Version 3.0, May 17, 2011.
******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

int outboun(int method)
{
	extern int nseg,mxx,myy,mzz,*segtyp,***nbou;
	extern float totalq,vol,lb;
	extern float *axt,*ayt,*azt,aly,**start,**end,**scos,*lseg,*rseg;
	
	int i,j,k,iseg,jseg,a,b,c,am,nnt;
	int ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax;
	float aa,bb,cc,abc,ddmin,ddmax,dd,t;
	float ds2,de2,disp2,d,x11,x22,x33;

	if(method == 1){
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++) nbou[i][j][k] = 1;
		am = 6;
		for(a=0; a<=am; a++) for(b=-am; b<=am; b++) for(c=-am; c<=am; c++){
			if(a != 0 || b != 0 || c != 0){
				aa = 0.;
				bb = 0.;
				cc = 0.;
				if(a != 0) aa = 1.0/a;
				if(b != 0) bb = 1.0/b;
				if(c != 0) cc = 1.0/c;
				abc = sqrt(SQR(aa) + SQR(bb) + SQR(cc));
				ddmin =  1.e8;
				ddmax = -1.e8;
				for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
					dd = (aa*start[1][iseg] + bb*start[2][iseg] + cc*start[3][iseg])/abc;
					ddmin = FMIN(dd - rseg[iseg] - lb,ddmin);
					ddmax = FMAX(dd + rseg[iseg] + lb,ddmax);
					dd = (aa*end[1][iseg] + bb*end[2][iseg] + cc*end[3][iseg])/abc;
					ddmin = FMIN(dd - rseg[iseg] - lb,ddmin);
					ddmax = FMAX(dd + rseg[iseg] + lb,ddmax);
				}
				for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
					t = (aa*axt[i] + bb*ayt[j] + cc*azt[k])/abc;
					if(t > ddmax || t < ddmin) nbou[i][j][k] = 0;
				}
			}
		}
	}
	if(method == 2){
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
			nbou[i][j][k] = 0;
			for(jseg=1; jseg<=nseg; jseg++) if(segtyp[jseg] == 4 || segtyp[jseg] == 5){
				x11 = (axt[i]-start[1][jseg])*scos[2][jseg]-(ayt[j]-start[2][jseg])*scos[1][jseg];
				x22 = (ayt[j]-start[2][jseg])*scos[3][jseg]-(azt[k]-start[3][jseg])*scos[2][jseg];
				x33 = (azt[k]-start[3][jseg])*scos[1][jseg]-(axt[i]-start[1][jseg])*scos[3][jseg];
				disp2 = SQR(x11) + SQR(x22) + SQR(x33);
				ds2 = SQR(axt[i]-start[1][jseg]) + SQR(ayt[j]-start[2][jseg]) + SQR(azt[k]-start[3][jseg]);
				de2 = SQR(axt[i]-end[1][jseg]) + SQR(ayt[j]-end[2][jseg]) + SQR(azt[k]-end[3][jseg]);
				if(FMAX(ds2,de2)-disp2 > SQR(lseg[jseg])) d = sqrt(FMIN(ds2,de2)) - rseg[jseg];
				else d = sqrt(disp2) - rseg[jseg];
				if(d < lb) nbou[i][j][k] = 1;
			}
		}
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++) if(nbou[i][j][k] == 0){
			imin = mxx;
			imax = 1;
			for(ii=1; ii<=mxx; ii++) if(nbou[ii][j][k] == 1){
				imin = IMIN(ii,imin);
				imax = IMAX(ii,imax);
			}
			jmin = myy;
			jmax = 1;
			for(jj=1; jj<=myy; jj++) if(nbou[i][jj][k] == 1){
				jmin = IMIN(jj,jmin);
				jmax = IMAX(jj,jmax);
			}
			kmin = mzz;
			kmax = 1;
			for(kk=1; kk<=mzz; kk++) if(nbou[i][j][kk] == 1){
				kmin = IMIN(kk,kmin);
				kmax = IMAX(kk,kmax);
			}
			if(i >= imin && i <= imax) if(j >= jmin && j <= jmax) if(k >= kmin && k <= kmax) nbou[i][j][k] = 1;
		}
	}
	nnt = 0;
	for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++) if(nbou[i][j][k] == 1){
		nnt++;
		nbou[i][j][k] = nnt;
	}
	return nnt;
}