/************************************************************************
setuparrays1 - for GreensTD.  TWS June 2012
Set up arrays with dimensions nnod and nseg
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays1(int nseg, int nnod)
{
	extern int nsp,nodsegm;
	extern int *nspoint,*istart,*nk,*nodrank,*nodout,*nodtyp,*ista,*iend;
	extern int **nodnod,**nodseg;
	extern float *segvar,*qq,*rseg,*lseg,*ds,*nodvar,*segc;
	extern float **gamma1,**qvseg,**cvseg,**start,**end,**scos;

	nspoint = ivector(1,nseg);
	istart = ivector(1,nseg);
	nk = ivector(1,nnod); //not nseg - error fixed 20 April 2010
	ista = ivector(1,nseg);
	iend = ivector(1,nseg);
	nodrank = ivector(1,nnod);
	nodout = ivector(1,nnod);
	nodtyp = ivector(1,nnod);

	nodnod = imatrix(1,nodsegm,1,nnod);
	nodseg = imatrix(1,nodsegm,1,nnod);

	segvar = vector(1,nseg);
	qq = vector(1,nseg);
	rseg = vector(1,nseg);
	lseg = vector(1,nseg);
	ds = vector(1,nseg);
	nodvar = vector(1,nnod);
	segc = vector(1,nseg);

	gamma1 = matrix(1,nseg,1,nsp);
	qvseg = matrix(1,nseg,1,nsp);
	cvseg = matrix(1,nseg,1,nsp);
	start = matrix(1,3,1,nseg);
	scos = matrix(1,3,1,nseg);
	end = matrix(1,3,1,nseg);
}