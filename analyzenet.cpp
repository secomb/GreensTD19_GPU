/************************************************************************
analyzenet - for Greens07_3D.  TWS October 07
Set up nodtyp, nodseg, nodnod arrays based on flowing segments
Nodes with no flowing segments or with fixed pressure are assigned nodtyp = 0
These values are subsequently overridden by putrank.
Version 3.0, May 17, 2011.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

int outboun(int method);

void analyzenet()
{
	extern int nseg,nnod,nnodbc,nnv,nnt,mxx,myy,mzz,nodsegm;
	extern int *bcnodname,*bcnod,*bctyp,*nodtyp,*segtyp,*nspoint,*nodname,*ista,*iend,*istart;
	extern int **segnodname,**nodnod,**nodseg;
	extern int ***nbou;
	extern float pi1,maxl,tlength,tarea,lb,alx,aly,alz;
	extern float *diam,*lseg,*rseg,*ds,*ss,*q,*qq,*axt,*ayt,*azt,*hd;
	extern float **cnode,**start,**end,**scos;
	extern float ***rsta,***rend;

	int i,j,k,iseg,inod,inod1,inod2,inodbc,m;
	float delx,dely,delz;

//Find node numbers corresponding to segment nodes - for all segments
	for(inod=1; inod<=nnod; inod++) nodtyp[inod] = 0;
	for(iseg=1; iseg<=nseg; iseg++){
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[1][iseg]){
			ista[iseg] = inod;
			goto foundit1;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[1][iseg]);
		foundit1:;
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[2][iseg]){
			iend[iseg] = inod;
			goto foundit2;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[2][iseg]);
		foundit2:;
	}
//Setup nodtyp, nodseg and nodnod
	for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		inod1 = ista[iseg];
		inod2 = iend[iseg];
		nodtyp[inod1]++;
		nodtyp[inod2]++;
		if(nodtyp[inod1] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod1);
		if(nodtyp[inod2] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod2);
		nodseg[nodtyp[inod1]][inod1] = iseg;
		nodseg[nodtyp[inod2]][inod2] = iseg;
		nodnod[nodtyp[inod1]][inod1] = inod2;
		nodnod[nodtyp[inod2]][inod2] = inod1;
	}
//Find node numbers corresponding to boundary nodes
	for (inodbc=1; inodbc<=nnodbc; inodbc++){
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == bcnodname[inodbc]){
			bcnod[inodbc] = inod;
			if(nodtyp[inod] != 1) printf("*** Error: Boundary node %i is not a 1-segment node\n", inod);
			goto foundit;
		}
		printf("*** Error: No matching node found for nodname %i\n", bcnodname[inodbc]);
		foundit:;
	}
//start[k][iseg] = coordinates of starting point of segment i
	tlength = 0.;
	tarea = 0.;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		rseg[iseg] = diam[iseg]/2.0;
		qq[iseg] = fabs(q[iseg]);
		lseg[iseg] = 0.;
		for(k=1; k<=3; k++){
			start[k][iseg] = cnode[k][ista[iseg]];
			end[k][iseg] = cnode[k][iend[iseg]];
			ss[k] = end[k][iseg] - start[k][iseg];
			lseg[iseg] += SQR(ss[k]);
		}
		lseg[iseg] = sqrt(lseg[iseg]);
		if(lseg[iseg] == 0.) printf("*** Error: segment %i has zero length\n",iseg);//added May 2010
		for(j=1; j<=3; j++)	scos[j][iseg] = ss[j]/lseg[iseg];
		tlength += lseg[iseg];
		tarea += pi1*lseg[iseg]*diam[iseg];
	}
//subdivide segments into subsegments as needed, compute coordinates of source points
//nnv = total number of vessel elements
	nnv = 0;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		m = lseg[iseg]/maxl + 1;
		nspoint[iseg] = m;
		istart[iseg] = nnv + 1;
		ds[iseg] = lseg[iseg]/m;
		nnv += m;
	}
	printf("Total vessel points = %i\n", nnv);
//compute coordinates of tissue points
	delx = alx/mxx;
	for(i=1; i<=mxx; i++) axt[i] = (i-0.5f)*delx;
	dely = aly/myy;
	for(i=1; i<=myy; i++) ayt[i] = (i-0.5f)*dely;
	delz = alz/mzz;
	for(i=1; i<=mzz; i++) azt[i] = (i-0.5f)*delz;
//create array of tissue points inside tissue domain boundaries using method 1 or 2
	nnt = outboun(1);
}