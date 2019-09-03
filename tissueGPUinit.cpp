/**************************************************************************
tissueGPUinit
initialize tissueGPU
TWS, December 2011
**************************************************************************/
#include <cuda_runtime.h>
#include "nrutil.h"

void tissueGPUinit(int nntGPU, int nnvGPU)
{
	extern int useGPU,nsp;
	extern int *h_tisspoints,*d_tisspoints;
	extern float *h_ct,*h_ctprev,*h_qt,*h_qtp,*h_pv,*h_qv,*h_gtt,*h_gbartt;
	extern float *d_ct,*d_ctprev,*d_qt,*d_qtp,*d_pv,*d_qv,*d_gtt,*d_gbartt;
	extern float *h_tissxyz,*h_vessxyz,*d_tissxyz,*d_vessxyz;

	h_tisspoints = ivector(0,3*nntGPU-1);	//integer coordinates of tissue points
	h_ct = vector(0,nntGPU-1);
	h_ctprev = vector(0,nntGPU-1);
	h_qt = vector(0,nntGPU-1);
	h_qtp = vector(0,nntGPU-1);
	h_pv = vector(0,nnvGPU-1);
	h_qv = vector(0,nnvGPU-1);
	h_gtt = vector(0,nsp*nntGPU-1);			//these must cover all possible distances, depend on solute
	h_gbartt = vector(0,nsp*nntGPU-1);
	h_tissxyz = vector(0,3*nntGPU-1);		//coordinates of tissue points
	h_vessxyz = vector(0,3*nnvGPU-1);		//coordinates of vessel points

	cudaSetDevice( useGPU-1 );
	cudaMalloc((void **)&d_tisspoints, 3*nntGPU*sizeof(int));
	cudaMalloc((void **)&d_ct, nntGPU*sizeof(float));
	cudaMalloc((void **)&d_ctprev, nntGPU*sizeof(float));
	cudaMalloc((void **)&d_qt, nntGPU*sizeof(float));
	cudaMalloc((void **)&d_qtp, nntGPU*sizeof(float));
	cudaMalloc((void **)&d_pv, nnvGPU*sizeof(float));
	cudaMalloc((void **)&d_qv, nnvGPU*sizeof(float));
	cudaMalloc((void **)&d_gtt, nsp*nntGPU*sizeof(float));
	cudaMalloc((void **)&d_gbartt, nsp*nntGPU*sizeof(float));
	cudaMalloc((void **)&d_tissxyz, 3*nntGPU*sizeof(float));
	cudaMalloc((void **)&d_vessxyz, 3*nnvGPU*sizeof(float));
}