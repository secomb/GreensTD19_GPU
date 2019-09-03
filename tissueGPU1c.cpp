/**************************************************************************
tissueGPU1.cpp
program to call tissueGPU1.cu on GPU
TWS, December 2011
**************************************************************************/
#include <cuda_runtime.h>

extern "C" void tissueGPU1(int *d_tisspoints, float *d_gtt, float *d_gbartt, float *d_ct, float *d_ctprev, float *d_qt, int nnt, int nntGPU, int useGPU, int isp);

void tissueGPU1c(int isp)
{
	extern int nnt,nntGPU,useGPU,*d_tisspoints;
	extern float *h_ct,*h_ctprev,*h_qt,*d_qt,*d_ct,*d_ctprev,*d_gtt,*d_gbartt;
	cudaError_t error;

	cudaSetDevice( useGPU-1 );
	error = cudaMemcpy(d_qt, h_qt, nnt*sizeof(float), cudaMemcpyHostToDevice);
	error = cudaMemcpy(d_ctprev, h_ctprev, nnt*sizeof(float), cudaMemcpyHostToDevice);
	tissueGPU1(d_tisspoints,d_gtt,d_gbartt,d_ct,d_ctprev,d_qt,nnt,nntGPU,useGPU,isp);
	error = cudaMemcpy(h_ct, d_ct, nnt*sizeof(float), cudaMemcpyDeviceToHost);
}