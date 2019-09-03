/************************************************************************
bicgstab - for Greens07_3D.  TWS April 2010
Based on algorithm 12, Parameter-free iterative linear solver by R. Weiss, 1996
system of linear equations:  aa(i,j)*x(j) = b(i)
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

double bicgstab(double **a, double *b, double *x, int n, double eps, int itmax)
{
	double lu,lunew,beta,delta,er,gamma1,t1,t2,err;
	double *r,*rs,*v,*s,*t,*p;
	int i,j,kk;
	r = dvector(1,n);
	rs = dvector(1,n);
	v = dvector(1,n);
	s = dvector(1,n);
	t = dvector(1,n);
	p = dvector(1,n);
	lu = 0.;
	for(i=1; i<=n; i++){
        r[i] = 0.;
		for(j=1; j<=n; j++){
			r[i] += a[i][j]*x[j];
		}
        r[i] -= b[i];
        p[i] = r[i];
 //       rs[i] = 1.;
        rs[i] = r[i];	//this is better if we test first that r is not an exact solution
        lu += r[i]*rs[i];
	}
	kk = 1;
	do
	{
		t1 = 0.;
		for(i=1; i<=n; i++){
			v[i] = 0.;
			for(j=1; j<=n; j++) v[i] += a[i][j]*p[j];
			t1 += v[i]*rs[i];
		}
		delta = -lu/t1;
		for(i=1; i<=n; i++) s[i] = r[i] + delta*v[i];
		for(i=1; i<=n; i++){
			t[i] = 0.;
			for(j=1; j<=n; j++) t[i] += a[i][j]*s[j];
		}
		t1 = 0.;
		t2 = 0.;
		for(i=1; i<=n; i++){
			t1 += s[i]*t[i];
			t2 += t[i]*t[i];
		}
		gamma1 = -t1/t2;
		err = 0.;
		lunew = 0.;
		for(i=1; i<=n; i++){
			r[i] = s[i] + gamma1*t[i];
			er = delta*p[i] + gamma1*s[i];
			x[i] += er;
			if(fabs(er) > err) err = fabs(er);
			lunew += r[i]*rs[i];
		}
		beta = lunew*delta/(lu*gamma1);
		lu = lunew;
		for(i=1; i<=n; i++) p[i] = r[i] + beta*(p[i]+gamma1*v[i]);
		kk += 1;
	}
	while(kk < itmax && err > eps);
	free_dvector(r,1,n);
	free_dvector(rs,1,n);
	free_dvector(v,1,n);
	free_dvector(s,1,n);
	free_dvector(t,1,n);
	free_dvector(p,1,n);
	if(err > eps) printf("*** Warning: linear solution not converged, %i iterations, err = %e\n",kk,err);
	return err;
}