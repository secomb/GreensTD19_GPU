/**************************************************************************
bicgstabBLASDend - double precision
end bicgstabBLASD
TWS, March 2011
**************************************************************************/
#include <cuda_runtime.h>
#include "nrutil.h"

void bicgstabBLASDend(int nnv)
{
	extern int useGPU;
	extern double *h_x,*h_b,*h_a,*h_rs;
	extern double *d_a, *d_res, *d_x, *d_b;
	extern double *d_r, *d_rs, *d_v, *d_s, *d_t, *d_p, *d_er;
 	const int nmem0 = sizeof(double);	//needed for malloc
 	const int nmem1 = nnv*sizeof(double);	//needed for malloc
	const int nmem2 = nnv*nnv*sizeof(double);

	cudaSetDevice( useGPU-1 );

	cudaFree(d_res);
	cudaFree(d_b);
	cudaFree(d_x);
	cudaFree(d_er);
	cudaFree(d_p);
	cudaFree(d_t);
	cudaFree(d_s);
	cudaFree(d_v);
	cudaFree(d_rs);
	cudaFree(d_r);
	cudaFree(d_a);

	free_dvector(h_rs,0,nnv-1);
	free_dvector(h_b,0,nnv-1);
	free_dvector(h_x,0,nnv-1);
	free_dvector(h_a,0,nnv*nnv-1);
}
