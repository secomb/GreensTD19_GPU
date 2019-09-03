/******************************************************
input - reads input files.  
Note that input files refer to segname and nodname, but
all arrays use iseg and inod, which are assigned when reading the file, as indices.
Modified to enable multiple cases to be run. January 2015.
Note: Must have irun = 0 on first call
Must not have any changes in number of solutes or number of time points between different runs
*******************************************************/
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void input(int irun, int n_step, int n_diff, int n_perm)
{
	extern int max,nmax,nsl1,nsl2,slsegdiv,ntpts,useGPU;
	extern int mxx,myy,mzz,nnt,nnv,nseg,nnod,nsp,nnodbc,nodsegm,ninterval;
	extern int *permsolute,*bcnodname,*bcnod,*bctyp;
	extern int *segname,*segtyp,*nodname,*nspoint,*nl;
	extern int *oxygen,*diffsolute,*nresis,*intervalnst;
	extern int **tisspoints,**segnodname,**nodseg,**intervalin;

	extern float fn,p50,alphaO2,cs,q0fac,totalq,errfac;
	extern float plow,phigh,clowfac,chighfac,pphighfac;
	extern float lb,maxl,v,vol,req,pi1,alx,aly,alz;
	extern float xmax,ymax,scalefac;
	extern float *intervaldt,*tpts;
	extern float *axt,*ayt,*azt,*ds,*g0,*diff,*alphab,*alphat,*cmin,*cmax,*cmean;
	extern float *crefb,*creft,*cinitb,*cinitt,*lambda,*diam,*q,*hd,*bcprfl,*bchd;
	extern float *x,*xsl0,*xsl1,*xsl2,*clmin,*clint,*p,*cl;
	extern float **start,**scos,**ax,**cnode,**bcp,**bcpbase,**resisdiam,**resis,**tissparam;
	extern float **tconcs,**zv,***psl;

	int i_step,i_diff,i_perm;		//maximum 10

	int i,iseg,isp,nlmax,iinterval,itime,itime1,flag;
	float time;
	FILE *ifp;
	char bb[100],fname1[80],fname2[80],fname3[80];

	//use sequentially numbered files
	i_perm = irun/(n_step*n_diff);
	i_diff = (irun/n_step)%n_diff;
	i_step = irun%n_step;
	
	sprintf(fname1,"SoluteParams%i.dat",i_diff);
	sprintf(fname2,"IntravascRes%i.dat",i_perm);
	sprintf(fname3,"TimeDep%i.dat",i_step);  // use this version for multiple time steps for each case
	//sprintf(fname3,"TimeDep%i.dat",i_diff); // use this version with n_step = 1 for time step matched with diffusivity

	ifp = fopen(fname1, "r");	//"SoluteParams.dat"
	fgets(bb,max,ifp);
	printf("%s\n",bb);
	fscanf(ifp,"%i %i %*[^\n]", &nmax,&useGPU);
	fscanf(ifp,"%f%*[^\n]", &errfac);
	fscanf(ifp,"%f%*[^\n]", &p50);
	fscanf(ifp,"%f%*[^\n]", &fn);
	fscanf(ifp,"%f%*[^\n]", &cs);
	fscanf(ifp,"%f%*[^\n]", &q0fac);
	fscanf(ifp,"%i%*[^\n]", &nsp);
	if(irun == 0){	//only on first pass
		permsolute = ivector(1,nsp);
		diffsolute = ivector(1,nsp);
		oxygen = ivector(1,nsp);
		crefb = vector(1,nsp);
		creft = vector(1,nsp);
		cinitb = vector(1,nsp);
		cinitt = vector(1,nsp);
		diff = vector(1,nsp);
		alphab = vector(1,nsp);
		alphat = vector(1,nsp);
		g0 = vector(1,nsp);
		tissparam = matrix(1,3,1,nsp);
		lambda = vector(1,nsp);
		resisdiam = matrix(1,20,1,nsp);
		resis = matrix(1,20,1,nsp);
		nresis = ivector(1,nsp);
	}
	for(isp=1; isp<=nsp; isp++){
		fgets(bb,max,ifp);
		fgets(bb,max,ifp);
		printf("%s\n",bb);
		fscanf(ifp,"%i %i %i%*[^\n]", &permsolute[isp],&diffsolute[isp],&oxygen[isp]);
		if(diffsolute[isp] != 0 && diffsolute[isp] != 1) printf("*** Error: soluteparams.dat, diffsolute[isp] must be 0 or 1\n");
		if(oxygen[isp] != 0 && oxygen[isp] != 1) printf("*** Error: soluteparams.dat, oxygen[isp] must be 0 or 1\n");
		if(oxygen[isp]) permsolute[isp] = 1; //oxygen is permeable
		if(permsolute[isp]) diffsolute[isp] = 1;  //permeable solutes must be diffusible
		fscanf(ifp,"%f%*[^\n]", &diff[isp]);
		fscanf(ifp,"%f%*[^\n]", &crefb[isp]);
		fscanf(ifp,"%f%*[^\n]", &creft[isp]);
		fscanf(ifp,"%f%*[^\n]", &cinitb[isp]);
		fscanf(ifp,"%f%*[^\n]", &cinitt[isp]);
		fscanf(ifp,"%f%*[^\n]", &alphab[isp]);
		if(oxygen[isp]) alphaO2 = alphab[isp];	//oxygen solubility in blood
		fscanf(ifp,"%f%*[^\n]", &alphat[isp]);
		diff[isp] = diff[isp]*1.e8;
		for(i=1; i<=3; i++)	fscanf(ifp,"%f%*[^\n]", &tissparam[i][isp]);
		fscanf(ifp,"%f%*[^\n]", &lambda[isp]);
	}
	fclose(ifp);
//parameters for blood.  TWS January 2012
	plow = 0.1*p50;
	phigh = 5.*p50;
	clowfac = cs*(1.0 - 1.0/(1.0 + pow((plow/p50),fn)));
	chighfac = cs*(1.0 - 1.0/(1.0 + pow((phigh/p50),fn)));
	pphighfac = cs*fn/p50*pow(phigh/p50,(fn-1))/SQR(1. + pow(phigh/p50,fn));

//intravascular or transvascular resistance data.  Assume zero unless specified in data file.
	ifp = fopen(fname2, "r");	//"IntravascRes.dat"
	for(isp=1; isp<=nsp; isp++){
		fscanf(ifp, "%i", &nresis[isp]);
		if(nresis[isp] == 0 && permsolute[isp]) printf("*** Error: no permeability data in IntravascRes.dat for permeable solute %i\n",isp);
		if(nresis[isp] > 20) printf("*** Error: too many points in IntravascRes.dat, nresis = %i > 20\n",nresis[isp]);
		fgets(bb,max,ifp);
		if(nresis[isp] > 0){
			fgets(bb,max,ifp);
			for(i=1; i<=nresis[isp]; i++) fscanf(ifp,"%f %f", &resisdiam[i][isp],&resis[i][isp]);
		}
	}
	fclose(ifp);
	if(irun == 0){	//only on first pass
		//network data file. Note that only segments of type 4 and 5 are included in computation
		//for compatibility with angiogenesis simulations
		ifp = fopen("Network.dat", "r");
		fgets(bb,max,ifp);
		printf("%s\n",bb);
	//dimensions of box in microns; vertex must be at origin
		fscanf(ifp,"%f %f %f%*[^\n]", &alx,&aly,&alz);
		fscanf(ifp,"%i %i %i%*[^\n]", &mxx,&myy,&mzz);
		fscanf(ifp,"%f%*[^\n]", &lb);
		fscanf(ifp,"%f%*[^\n]", &maxl);
		fscanf(ifp,"%i%*[^\n]", &nodsegm);
	//number of segments in vessel network
		fscanf(ifp,"%i%*[^\n]", &nseg);
		fgets(bb,max,ifp);
		fgets(bb,max,ifp);
	//segment properties: name type nodefrom nodeto diameter flow hematocrit
		segname = ivector(1,nseg);
		segtyp = ivector(1,nseg);
		segnodname = imatrix(1,2,1,nseg);
		diam = vector(1,nseg);
		q = vector(1,nseg);
		hd = vector(1,nseg);
		for(iseg=1; iseg<=nseg; iseg++) fscanf(ifp, "%i %i %i %i %f %f %f%*[^\n]",
			&segname[iseg],&segtyp[iseg],&segnodname[1][iseg],&segnodname[2][iseg],&diam[iseg],&q[iseg],&hd[iseg]);
	//number of nodes in vessel network
		fscanf(ifp,"%i%*[^\n]", &nnod);
		fgets(bb,max,ifp);
		fgets(bb,max,ifp);
	//coordinates of nodes
		nodname = ivector(1,nnod);
		cnode = matrix(1,3,1,nnod);
		for(i=1; i<=nnod; i++) fscanf(ifp, "%i %f %f %f%*[^\n]", &nodname[i],&cnode[1][i],&cnode[2][i],&cnode[3][i]);
	//boundary nodes
		fscanf(ifp,"%i%*[^\n]", &nnodbc);
		fgets(bb,max,ifp);
		fgets(bb,max,ifp);
		bcnodname = ivector(1,nnodbc);
		bcnod = ivector(1,nnodbc);
		bctyp = ivector(1,nnodbc);
		bcprfl = vector(1,nnodbc);
		bchd = vector(1,nnodbc);
		bcp = matrix(1,nnodbc,1,nsp);
		bcpbase = matrix(1,nnodbc,1,nsp);
		totalq = 0.;
		for(i=1; i<=nnodbc; i++){
			fscanf(ifp,"%i %i %f %f%", &bcnodname[i],&bctyp[i],&bcprfl[i],&bchd[i]);
			//read baseline values of solute concentrations at boundary nodes
			for(isp=1; isp<=nsp; isp++) if(permsolute[isp]) fscanf(ifp,"%f",&bcpbase[i][isp]);
			fscanf(ifp,"%*[^\n]");	//ignore any 'extra' solutes in data file
			if(bctyp[i] == 2 && bcprfl[i] > 0.) totalq = totalq + bcprfl[i];
		}
		fclose(ifp);
	
		//scale flows according to factor q0fac.  Modified 2012
		for(iseg=1; iseg<=nseg; iseg++) q[iseg] = q[iseg]*q0fac;
		for(i=1; i<=nnodbc; i++) if(bctyp[i] == 2) bcprfl[i] = bcprfl[i]*q0fac;
		totalq = totalq*q0fac;
		//v = total box volume, vol = volume represented by each tissue point; 
		v = alx*aly*alz;
		vol = v/(mxx*myy*mzz);
		req = pow(vol*0.75/pi1,0.333333);
		//Read parameters for slice on which P is computed for contour plot
		nl = ivector(1,nsp);
		xsl0 = vector(1,3);
		xsl1 = vector(1,3);
		xsl2 = vector(1,3);
		clmin = vector(1,nsp);
		clint = vector(1,nsp);
	
		ifp = fopen("ContourParams.dat", "r");
		fscanf(ifp, "%f %f %f %i%*[^\n]", &xsl0[1],&xsl0[2],&xsl0[3],&slsegdiv);
		fscanf(ifp, "%f %f %f %i%*[^\n]", &xsl1[1],&xsl1[2],&xsl1[3],&nsl1);
		fscanf(ifp, "%f %f %f %i%*[^\n]", &xsl2[1],&xsl2[2],&xsl2[3],&nsl2);
		nlmax = 1;
		for(isp=1; isp<=nsp; isp++){
			fscanf(ifp, "%f %f %i%*[^\n]", &clmin[isp],&clint[isp],&nl[isp]);
			if(nl[isp] > nlmax) nlmax = nl[isp];
		}
		fclose(ifp);
		xmax = sqrt(SQR(xsl1[1]-xsl0[1]) + SQR(xsl1[2]-xsl0[2]) + SQR(xsl1[3]-xsl0[3]));
		ymax = sqrt(SQR(xsl2[1]-xsl0[1]) + SQR(xsl2[2]-xsl0[2]) + SQR(xsl2[3]-xsl0[3]));
		scalefac = FMIN(480./xmax,700./ymax);//updated April 2010
		cl = vector(1,nlmax);
		zv = matrix(1,nsl1,1,nsl2);
		psl = f3tensor(1,nsl1,1,nsl2,1,nsp);
	}

	ifp = fopen(fname3, "r");	//TimeDep.dat
	fgets(bb,max,ifp);
	fscanf(ifp,"%i%*[^\n]", &ninterval);

	intervalnst = ivector(1,ninterval);
	intervalin = imatrix(1,ninterval,1,nsp);
	intervaldt = vector(1,ninterval);
	ntpts = 0;
	for(iinterval=1; iinterval<=ninterval; iinterval++){
			fscanf(ifp,"%f %i ", &intervaldt[iinterval],&intervalnst[iinterval]);
			for(isp=1; isp<=nsp; isp++) if(permsolute[isp]) fscanf(ifp,"%i ", &intervalin[iinterval][isp]);
			fscanf(ifp,"%*[^\n]");
			ntpts += intervalnst[iinterval];
	}
	if(ntpts > 0){
		if(irun > 0){
			free_vector(tpts,1,ntpts);	//computation times may need to be resized
			free_matrix(tconcs,1,ntpts,1,nsp);
		}
		tpts = vector(1,ntpts);	//computation times
		tconcs = matrix(1,ntpts,1,nsp);
		itime = 0;
		time = 0.;
		for(iinterval=1; iinterval<=ninterval; iinterval++)	for(i=1; i<=intervalnst[iinterval]; i++){
			itime++;
			time += intervaldt[iinterval];
			tpts[itime] = time;
			flag = 0;
		for(isp=1; isp<=nsp; isp++) if(permsolute[isp]) {
				if(intervalin[iinterval][isp] == 0) tconcs[itime][isp] = 0.;
				else if(intervalin[iinterval][isp] == 1) tconcs[itime][isp] = 1.;
				else if(intervalin[iinterval][isp] == 2) flag = 1;
				else printf("*** Error in TimeDep.dat\n");
			}
			if(flag){
				fscanf(ifp, "%i", &itime1);	//only read in time-dependent inflow (type 2) solutes
			for(isp=1; isp<=nsp; isp++) if(permsolute[isp]) if(intervalin[iinterval][isp] == 2) fscanf(ifp,"%f",&tconcs[itime][isp]);
				fscanf(ifp,"%*[^\n]");	//ignore any 'extra' solutes in data file
			}
		}
	}
	fclose(ifp);
}