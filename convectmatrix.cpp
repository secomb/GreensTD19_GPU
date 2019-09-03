/************************************************************************
convectmatrix - for greensTD
TWS July 2013
Revised April 16, 2015 with new definitions of ftr and ftf (sign reversed)
Set up convective matrices for convective transport over a finite time deltat
This version tracks the transport of a fluid 'slug' backward through the network 
alphaconv[i][j] = amount of transit time spent in segment j before reaching segment i
betaconv[i][j]  = amount originating in segment j reaching segment i
gammaconv[i][j] = amount originating in boundary node j reaching segment i
zetaconv[i][j]  = amount of transit time spent in segment j before reaching boundary node i
omegaconv[i][j] = amount originating in segment j reaching boundary node i
xiconv[i][j]    = amount originating in boundary node j reaching boundary node i
Boundary nodes are labeled nnv+1 to nnv+nnodbc.
"Slug i" is defined as fluid in segment i at time t=deltat, or fluid entering outflow node
i-nnv during interval [0,deltat].
Quantities used to decribe location of slug i:
ftfseg[j] = arrival time of front of slug to upstream end of current subsegment
ftrseg[j] = arrival time of rear of slug to upstream end of current subsegment
fvseg[j] = volume in slug (may be less than flow rate x slug transit time)
These quantities are also defined at the main nodes
ftfnod[inod] = arrival time of front of slug to node
ftrnod[inod] = arrival time of rear of slug to node
fvnod[ii][inod] = volume in slug in inflow segment ii relative to that node
In particular, ftfseg[i] = deltat - transit[i] and ftrseg[i] = deltat.
For an upstream sefment j connected to segment i with transit time transit[j],
ftfseg[j] = deltat - transit[i] - transit[j] and ftrseg[j] = deltat  - transit[j].
This process is repeated iteratively to show location of segments relative to all segments.
At time t=0, the residence times of the front and rear of a slug are <= -ftfseg[j] and ftrseg[j].
These times are compared with the transit fime transit[j] in subroutines lengthcalc and areacalc,
to compute the contributions to the ij elements of the various matrices. 
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float areacalc(float a, float b, float x0, float x1);
float lengthcalc(float a, float x0, float x1);
 
void convectmatrix()
{
	extern int nnodbc,nnod,nseg,nnodfl,nodsegm,nsp,nnv;
	extern int *mainseg,*nspoint,*istart,*nodtyp,*nodrank,*nodout,*bcnod,*ista,*iend;
	extern int **nodseg;
	extern float *transit,*q,*qq,*ds,*rseg,*volseg;
	extern float **alphaconv,**betaconv,**gammaconv,**omegaconv,**zetaconv,**xiconv;
	extern float deltat,pi1,flowfac,q0fac;

	int i,j,j1,iseg,iseg1,ii,jjj,flag,niter,nodt,nout,nodt1,nout1,iii,in,inod,inod1,inodbc;
	int collisionflag;
	float vel,sum,volout,lratio,aratio;
	float ftrseg1,ftfseg1,fvseg1;	//temporary storage variables 	
	float *ftrseg,*ftfseg,*fvseg,*sumin;	//variables associated with subsegments
	float *ftrnod,*ftfnod,**fvnod;	//variables associated with main nodes

	ftrseg = vector(1,nnv);
	ftfseg = vector(1,nnv);
	fvseg = vector(1,nnv);
	ftrnod = vector(1,nnod);
	ftfnod = vector(1,nnod);
	sumin = vector(1,nnod);
	fvnod = matrix(1,nodsegm,1,nnod);

	printf("Computing convection coefficient matrices...");

	for(i=1; i<=nnv; i++){ 
		for(j=1; j<=nnv; j++) alphaconv[i][j] = 0.;
		for(j=1; j<=nnv; j++) betaconv[i][j] = 0.;
		for(j=1; j<=nnodbc; j++) gammaconv[i][j] = 0.;
		for(j=1; j<=nnodbc; j++) omegaconv[j][i] = 0.;
		for(j=1; j<=nnodbc; j++) zetaconv[j][i] = 0.;
	}
	for(i=1; i<=nnodbc; i++) for(j=1; j<=nnodbc; j++) xiconv[j][i] = 0.;
	if(fabs(q0fac) > 0.){				//skip if all flows are set to zero
		for(i=1; i<=nnv; i++){			//compute transit times, etc.
			iseg = mainseg[i];
			volseg[i] = pi1*SQR(rseg[iseg])*ds[iseg];
			vel = qq[iseg]/pi1/SQR(rseg[iseg])*flowfac;	//velocity in micron/s
			transit[i] = ds[iseg]/vel;	//transit time in seconds
		}
		for(inod=1; inod<=nnod; inod++){	//calculate sum of input flows to each node
			sumin[inod] = 0.;
			nodt = nodtyp[inod];
			nout = nodout[inod];
			for(ii=nout+1; ii<=nodt; ii++){		//inflows to the node
				iseg = nodseg[ii][inod];
				sumin[inod] += qq[iseg];
			}
		}
//-----------------------------------------------------------------------------------------------------------
		for(i=1; i<=nnv+nnodbc; i++){	// Begin main loop: scan over destination subsegments and boundary nodes
			for(j=1; j<=nnv; j++){
				ftrseg[j] = 0.;
				ftfseg[j] = 0.;
				fvseg[j] = 0.;
			}			
			for(inod=1; inod<=nnod; inod++){
				ftrnod[inod] = 0.;
				ftfnod[inod] = 0.;
				for(ii=1; ii<=nodsegm; ii++) fvnod[ii][inod] = 0.;
			}			
			if(i <= nnv){				//initialize segment variables with contents of segment i, propagated upstream by deltat
				ftfseg[i] = deltat - transit[i];
				ftrseg[i] = deltat;
				fvseg[i] = volseg[i];
			}
			else{						//initialize segment variables with node outflow during deltat
				inod = bcnod[i-nnv];
				if(nodout[inod] == 0){	//is an outflow node
					iseg = nodseg[1][inod];
					ftfnod[inod] = 0.;
					ftrnod[inod] = deltat;
					fvnod[1][inod] = qq[iseg]*flowfac*deltat;	//assign slug volume to upstream segment associated with node
																//for boundary nodes, there is only one segment
				}
			}
			flag = 1;
			niter = 0;
			do{
/************************************************************************************************/
				for(in=nnodfl; in>=1; in--){	//scan nodes in upstream order
					inod = nodrank[in];
					nodt = nodtyp[inod];
					nout = nodout[inod];
					if(nodt > nout) for(ii=nout+1; ii<=nodt; ii++){		//process inflow segments to the node, if any
						iseg = nodseg[ii][inod];
						if(q[iseg] > 0.) inod1 = ista[iseg];	//find main node at upstream end of this segment
						else if(q[iseg] < 0.) inod1 = iend[iseg];
						else printf("*** Error: Nonflowing segment included %i\n", iseg);
						if(inod == inod1) printf("*** Error: flow direction mismatch\n");
						nodt1 = nodtyp[inod1];
						nout1 = nodout[inod1];
						collisionflag = 0;
						if(nodt1 > nout1) for(iii=nout1+1; iii<=nodt1; iii++) if(fvnod[iii][inod1] > 0.0001) collisionflag = 1;
						if(collisionflag == 0){	//don't proceed if there is a collision due to diverging and converging pathways.  Try again later!
							for(jjj=1; jjj<=nspoint[iseg]; jjj++){				//propagate between subsegments within segment
								if(q[iseg] > 0.) j = istart[iseg] + nspoint[iseg] - jjj;	//process subsegments in upstream order 
								else j = istart[iseg] + jjj - 1;		//j is the current subsegment that may potentially contribute to slug
								if(jjj == 1 && fvnod[ii][inod] > 0.){	//we have found a node associated with a segment with some slug volume
																		//propagate node values to upstream subsegment
									if(fvseg[j] != 0.) printf("*** Error: this should be zero\n");	//otherwise there is a collision
									ftrseg[j] = ftrnod[inod] - transit[j];	//time label referred to new segment is increased by transit time
									ftfseg[j] = ftfnod[inod] - transit[j];
									fvseg[j] = fvnod[ii][inod];	//transfer slug volume from node to segment
									fvnod[ii][inod] = 0.;	//set node slug volume to zero
								}
								if(fvseg[j] > 0.){							//some slug volume is associated with subsegment
									//residence time at t = 0 is negative of arrival time.
									lratio = lengthcalc(transit[j], -ftrseg[j], -ftfseg[j])/(ftrseg[j] - ftfseg[j]);

									aratio = areacalc(transit[j], deltat, -ftrseg[j], -ftfseg[j])/deltat/(ftrseg[j] - ftfseg[j]);

									if(i <= nnv){
										alphaconv[i][j] += aratio*fvseg[j]/volseg[i];	
										betaconv[i][j] += lratio*fvseg[j]/volseg[i];	//assign lratio*(slug volume) to beta
									}
									else{
										volout = qq[nodseg[1][bcnod[i-nnv]]]*flowfac*deltat;
										zetaconv[i-nnv][j] += aratio*fvseg[j]/volout;
										omegaconv[i-nnv][j] += lratio*fvseg[j]/volout;	//assign lratio*(slug volume) to omega
									}
									fvseg1 = fvseg[j]*(1. - lratio);		//remaining slug volume to be attributed to upstreams subsegment(s)
									ftrseg1 = FMAX(ftrseg[j],0.);
									ftfseg1 = FMAX(ftfseg[j],0.);
									if(jjj<nspoint[iseg]){	//if there is an upstream subsegment in this segment, transfer to it
										if(q[iseg] < 0.) j1 = j + 1;
										else j1 = j - 1;
										fvseg[j1] = fvseg1;		//transfer slug volume to upstream subsegment
										ftfseg[j1] = ftfseg1 - transit[j1];
										ftrseg[j1] = ftrseg1 - transit[j1];
									}
									else{					//otherwise transfer slug volume to next main node upstream
										if(nodt1 > nout1){
											for(iii=nout1+1; iii<=nodt1; iii++){		//segment inflows to the node
												iseg1 = nodseg[iii][inod1];
												fvnod[iii][inod1] = fvseg1*qq[iseg1]/sumin[inod1];	//transfer slug volume from segment to node
																									//ready for transfer to associated segments
											}
											ftfnod[inod1] = ftfseg1;
											ftrnod[inod1] = ftrseg1;
										}
										else{											//must be an inflow node: assign to gamma or xi
											for(inodbc=1; inodbc<=nnodbc; inodbc++){
												if(bcnod[inodbc] == inod1){
													if(i <= nnv) gammaconv[i][inodbc] += fvseg1/volseg[i];
													else{
														volout = qq[nodseg[1][bcnod[i-nnv]]]*flowfac*deltat;
														xiconv[i-nnv][inodbc] += fvseg1/volout;
													}
													goto foundit;
												}
											}
											printf("*** Error: Unable to find boundary node matching node %i\n", inod);
											foundit:;
										}
									}
									fvseg[j] = 0.;		//now slug volume has been disposed of so set to zero
								}
							}
						}
					}					
				}
				flag = 0;
				for(j=1; j<=nnv; j++) if(fvseg[j] > 0.) flag = 1;
				for(ii=1; ii<=nodsegm; ii++) for(j=1; j<=nnod; j++) if(fvnod[ii][j] > 0.) flag = 1;
				niter++;
/************************************************************************************************/
			}
			while(flag == 1 && niter < 1000);
			if(flag == 1) printf("*** Error: failure in convectmatrix\n");
		}
//-----------------------------------------------------------------------------------------------------------
	}
	else{	//non-flowing network for testing purposes
		for(i=1; i<=nnv; i++){
			alphaconv[i][i] = 1.;
			betaconv[i][i] = 1.;
		}
	}
/*
//print out matrices
	FILE *ofp;
	ofp = fopen("ConvectMatrix.out", "w");
	fprintf(ofp,"alphaconv\n");
	for(i=1; i<=nnv; i++){		
		for(j=1; j<=nnv; j++) fprintf(ofp,"%6.3f ",alphaconv[i][j]);
		fprintf(ofp,"\n");
	}
	fprintf(ofp,"betaconv\n");
	for(i=1; i<=nnv; i++){		
		for(j=1; j<=nnv; j++) fprintf(ofp,"%6.3f ",betaconv[i][j]);
		fprintf(ofp,"\n");
	}
	fprintf(ofp,"gammaconv\n");
	for(i=1; i<=nnv; i++){		
		for(j=1; j<=nnodbc; j++) fprintf(ofp,"%6.3f ",gammaconv[i][j]);
		fprintf(ofp,"\n");
	}
	fprintf(ofp,"omegaconv\n");
	for(i=1; i<=nnodbc; i++){		
		for(j=1; j<=nnv; j++) fprintf(ofp,"%6.3f ",omegaconv[i][j]);
		fprintf(ofp,"\n");
	}
	fprintf(ofp,"zetaconv\n");
	for(i=1; i<=nnodbc; i++){		
		for(j=1; j<=nnv; j++) fprintf(ofp,"%6.3f ",zetaconv[i][j]);
		fprintf(ofp,"\n");
	}
	fprintf(ofp,"xiconv\n");
	for(i=1; i<=nnodbc; i++){		
		for(j=1; j<=nnodbc; j++) fprintf(ofp,"%6.3f ",xiconv[i][j]);
		fprintf(ofp,"\n");
	}
	fclose(ofp);
*/
//check input to each subsegment - eq. (33) of supplementary material
	for(i=1; i<=nnv; i++){
		sum = 0.;
		for(j=1; j<=nnv; j++) sum += betaconv[i][j];
		for(inodbc=1; inodbc<=nnodbc; inodbc++) sum += gammaconv[i][inodbc];
		if(fabs(sum - 1.) > 0.001)
			printf("*** Error: input mass conservation violation at subsegment %i\n", i);
	}
//check input to each outflow boundary node - eq. (34) of supplementary material
	for(i=1; i<=nnodbc; i++){
		inod = bcnod[i];
		if(nodout[inod] == 0){
			sum = 0.;
			for(j=1; j<=nnv; j++) sum += omegaconv[i][j];
			for(inodbc=1; inodbc<=nnodbc; inodbc++) sum += xiconv[i][inodbc];
			if(fabs(sum - 1.) > 0.001)
				printf("*** Error: input mass conservation violation at boundary node %i\n", i);
		}
	}
//check transit times to each subsegment do not add up to more than deltat- eq. (35)
	for(i=1; i<=nnv; i++){
		sum = 0.;
		for(j=1; j<=nnv; j++) sum += alphaconv[i][j];
		if(sum > 1.001)
			printf("*** Error: Mass conservation violation (alpha) at subsegment %i\n", i);
	}
//check transit times to node do not add up to more than deltat- eq. (36)
	for(inodbc=1; inodbc<=nnodbc; inodbc++){
		sum = 0.;
		for(j=1; j<=nnv; j++) sum += zetaconv[inodbc][j];
		if(sum > 1.001)
			printf("*** Error: Mass conservation violation (zeta) at subsegment %i\n", i);
	}
//check output from each subsegment - eq. (37) of supplementary material
		for(j=1; j<=nnv; j++){
			sum = 0.;
			for(i=1; i<=nnv; i++) sum += betaconv[i][j]*volseg[i];
			for(inodbc=1; inodbc<=nnodbc; inodbc++) sum += omegaconv[inodbc][j]*qq[nodseg[1][bcnod[inodbc]]]*flowfac*deltat;
			if(fabs(sum/volseg[j] - 1.) > 0.001)
				printf("*** Error: output mass conservation violation at subsegment %i\n", j);
		}
//check output from each inflow boundary node - eq. (38) of supplementary material
	for(j=1; j<=nnodbc; j++){
		inod = bcnod[j];
		if(nodout[inod] == 1){
			sum = 0.;
			for(i=1; i<=nnv; i++) sum += gammaconv[i][j]*volseg[i];
			for(inodbc=1; inodbc<=nnodbc; inodbc++) sum += xiconv[inodbc][j]*qq[nodseg[1][bcnod[inodbc]]]*flowfac*deltat;
			if(fabs(sum/(qq[nodseg[1][inod]]*flowfac*deltat) - 1.) > 0.001)
				printf("*** Error: output mass conservation violation at boundary node %i\n", j);
		}
	}
//check extractions in each segment add up to total extraction  - eq. (39) of supplementary material
	for(j=1; j<=nnv; j++){
		sum = 0.;
		for(i=1; i<=nnv; i++) sum += alphaconv[i][j]*volseg[i];
		for(inodbc=1; inodbc<=nnodbc; inodbc++) sum += zetaconv[inodbc][j]*qq[nodseg[1][bcnod[inodbc]]]*flowfac*deltat;
		if(fabs(sum/volseg[j] - 1.) > 0.001)
			printf("*** Error: total extraction violation at subsegment %i\n", j);
	}
	free_vector(ftrseg,1,nnv);
	free_vector(ftfseg,1,nnv);
	free_vector(fvseg,1,nnv);
	free_vector(ftrnod,1,nnod);
	free_vector(ftfnod,1,nnod);
	free_matrix(fvnod,1,nodsegm,1,nnod);
	free_vector(sumin,1,nnod);
	printf("done\n");
}

float lengthcalc(float a, float x0, float x1)
//compute the length of a line from x0 to x1 contained within the segment from 0 to a
//TWS, July 2013
{
	float length,lplus,lminus;
	if(x1 < x0) printf("*** Error: incorrect arguments for lengthcalc\n");
	length = a;
	if(x0 <= 0.) lminus = 0.;
	else if(x0 <= a) lminus = x0;
	else lminus = a;
	length -= lminus;
	if(x1 <= 0.) lplus = a;
	else if(x1 <= a) lplus = a - x1;
	else lplus = 0.;
	length -= lplus;
	return length;
}

float areacalc(float a, float b, float x0, float x1)
//compute the area of a rectangle contained between two parallel diagonal lines with slope 1
//x intercepts of lines are x0 and x1, where x1 > x0
//TWS, July 2013
{
	float area,aplus,aminus;
	area = a*b;
	if(x1 < x0) printf("*** Error: incorrect arguments for areacalc\n");
	if(a >= b){
		if(x0 <= -b) aplus = 0.;
		else if(x0 <= 0.) aplus = 0.5*SQR(x0 + b);
		else if(x0 <= a-b) aplus = b*(x0 + 0.5*b);
		else if(x0 <= a) aplus = a*b - 0.5*SQR(a - x0);
		else aplus = a*b;
	}
	else{
		if(x0 <= -b) aplus = 0.;
		else if(x0 <= a-b) aplus = 0.5*SQR(x0 + b);
		else if(x0 <= 0.) aplus = a*b - a*(-x0 + 0.5*a);
		else if(x0 <= a) aplus = a*b - 0.5*SQR(a - x0);
		else aplus = a*b;
	}
	area -= aplus;
	if(a >= b){
		if(x1 <= -b) aminus = a*b;
		else if(x1 <= 0.) aminus = a*b - 0.5*SQR(x1 + b);
		else if(x1 <= a-b) aminus = a*b - b*(x1 + 0.5*b);
		else if(x1 <= a) aminus = 0.5*SQR(a - x1);
		else aminus = 0.;
	}
	else{
		if(x1 <= -b) aminus = a*b;
		else if(x1 <= a-b) aminus = a*b - 0.5*SQR(x1 + b);
		else if(x1 <= 0.) aminus = a*(-x1 + 0.5*a);
		else if(x1 <= a) aminus = 0.5*SQR(a - x1);
		else aminus = 0.;
	}
	area -= aminus;
	return area;
}