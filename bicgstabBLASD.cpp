/**************************************************************************
BiCGSTABBLASD - double precision
Algorithm 12, Parameter-free iterative linear solver by R. Weiss, 1996
system of linear equations:  aa(i,j)*x(j) = b(i)
Version using BLAS library on GPU
TWS, March 2011
**************************************************************************/
#include <cuda_runtime.h>
#include <cublas.h>
#include <stdio.h>

double bicgstabBLASD(double **a, double *b, double *x, int n, double eps, int itmax)
{
	extern int useGPU;
	extern double *h_x,*h_b,*h_a,*h_rs;
	extern double *d_a, *d_res, *d_x, *d_b;
	extern double *d_r, *d_rs, *d_v, *d_s, *d_t, *d_p, *d_er;
 	const int nmem0 = sizeof(double);	//needed for memcpy
 	const int nmem1 = n*sizeof(double);
	const int nmem2 = n*n*sizeof(double);
	int i,j,kk,ierr;
	double lu,lunew,beta,delta,gamma1,t1,t2,err;

	for(j=0; j<n; j++){
		for(i=0; i<n; i++) h_a[i*n + j] = a[i+1][j+1];
		h_x[j] = x[j+1];
		h_b[j] = b[j+1];
		h_rs[j] = 1.;
	}

	cudaSetDevice( useGPU-1 );

	cudaMemcpy(d_a, h_a, nmem2, cudaMemcpyHostToDevice);//copy variables to GPU
    cudaMemcpy(d_x, h_x, nmem1, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, nmem1, cudaMemcpyHostToDevice);
	cudaMemcpy(d_rs, h_rs, nmem1, cudaMemcpyHostToDevice);

//	r[i] += a[i][j]*x[j];
	cublasDgemv('T', n, n, 1.f, d_a, n, d_x, 1, 0.f, d_r, 1);

//	r[i] -= b[i];
	cublasDaxpy (n, -1.f, d_b, 1, d_r, 1);

//	 p[i] = r[i];
	cublasDcopy (n, d_r, 1, d_p, 1);
	
//	lu += r[i]*rs[i];
	lu = cublasDdot (n, d_r, 1, d_rs, 1);

	kk = 1;
	do{
//		v[i] += a[i][j]*p[j];
		cublasDgemv('T', n, n, 1.f, d_a, n, d_p, 1, 0.f, d_v, 1);

//		t1 += v[i]*rs[i];
		t1 = cublasDdot (n, d_v, 1, d_rs, 1);
		if(t1 == 0){
			printf("t1 = 0, reset to 1e-12\n");
			t1 = 1.e-12;	//added January 2012
		}
		delta = -lu/t1;

//		s[i] = r[i] + delta*v[i];
		cublasDcopy (n, d_r, 1, d_s, 1);
		cublasDaxpy (n, delta, d_v, 1, d_s, 1);
		
//		t[i] += a[i][j]*s[j];
		cublasDgemv('T', n, n, 1.f, d_a, n, d_s, 1, 0.f, d_t, 1);

//		t1 += s[i]*t[i];
		t1 = cublasDdot (n, d_t, 1, d_s, 1);

//		t2 += t[i]*t[i];
		t2 = cublasDdot (n, d_t, 1, d_t, 1);
		if(t2 == 0){
			printf("t2 = 0, reset to 1e-12\n");
			t2 = 1.e-12;	//added January 2012
		}
		gamma1 = -t1/t2;

//		r[i] = s[i] + gamma1*t[i];
		cublasDcopy (n, d_s, 1, d_r, 1);
		cublasDaxpy (n, gamma1, d_t, 1, d_r, 1);

//		er[i] = delta*p[i] + gamma1*s[i];
		cublasDcopy (n, d_s, 1, d_er, 1);
		cublasDscal (n, gamma1, d_er, 1);
		cublasDaxpy (n, delta, d_p, 1, d_er, 1);

//		x[i] += er[i];
		cublasDaxpy (n, 1., d_er, 1, d_x, 1);

//		lunew += r[i]*rs[i];
		lunew = cublasDdot (n, d_r, 1, d_rs, 1);
		if(lunew == 0){
			printf("lunew = 0, reset to 1e-12\n");
			lunew = 1.e-12;	//added January 2012
		}
		if(gamma1 == 0){
			printf("gamma1 = 0, reset to 1e-12\n");
			gamma1 = 1.e-12;	//added January 2012
		}
		beta = lunew*delta/(lu*gamma1);
		lu = lunew;

//		p[i] = r[i] + beta*(p[i]+gamma1*v[i]);
		cublasDaxpy (n, gamma1, d_v, 1, d_p, 1);
		cublasDscal (n, beta, d_p, 1);
		cublasDaxpy (n, 1., d_r, 1, d_p, 1);

		ierr = cublasIdamax (n, d_er, 1);   //Find the maximum value in the array
		cudaMemcpy(&err, d_er+ierr-1, nmem0, cudaMemcpyDeviceToHost);
		kk++;
	}
	while(kk < itmax && fabs(err) > eps);

	cudaMemcpy(h_x, d_x, nmem1, cudaMemcpyDeviceToHost);//bring back x for final result
	for(i=0; i<n; i++) x[i+1] = h_x[i];

	if(fabs(err) > eps) printf("*** Warning: linear solution using BICGSTB not converged, err = %e\n",err);
	return err;	
}
