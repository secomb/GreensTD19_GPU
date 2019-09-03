/**************************************************************************
tissueGPUend
end tissueGPU
TWS, December 2011
**************************************************************************/
#include <cuda_runtime.h>
#include "nrutil.h"

void tissueGPUend(int nntGPU, int nnvGPU)
{
	extern int useGPU,nsp,*h_tisspoints,*d_tisspoints;
	extern float *h_ct,*h_ctprev,*h_qt,*h_qtp,*h_pv,*h_qv,*h_gtt,*h_gbartt,*h_tissxyz,*h_vessxyz;
	extern float *d_qt,*d_qtp,*d_ct,*d_ctprev,*d_qv,*d_pv,*d_gtt,*d_gbartt;
	extern float *d_tissxyz,*d_vessxyz;

	free_ivector(h_tisspoints,0,3*nntGPU-1);
	free_vector(h_ct,0,nntGPU-1);
	free_vector(h_ctprev,0,nntGPU-1);
	free_vector(h_qt,0,nntGPU-1);
	free_vector(h_qtp,0,nntGPU-1);
	free_vector(h_pv,0,nnvGPU-1);
	free_vector(h_qv,0,nnvGPU-1);
	free_vector(h_gtt,0,nsp*nntGPU-1);
	free_vector(h_gbartt,0,nsp*nntGPU-1);
	free_vector(h_tissxyz,0,3*nntGPU-1);		//coordinates of tissue points
	free_vector(h_vessxyz,0,3*nnvGPU-1);		//coordinates of vessel points

	cudaSetDevice( useGPU-1 );
	cudaFree(d_tisspoints);
	cudaFree(d_ct);
	cudaFree(d_ctprev);
	cudaFree(d_qt);
	cudaFree(d_qtp);
	cudaFree(d_pv);
	cudaFree(d_qv);
	cudaFree(d_gtt);
	cudaFree(d_gbartt);
	cudaFree(d_tissxyz);
	cudaFree(d_vessxyz);
}
