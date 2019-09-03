/**************************************************************************
tissueGPUcopy
copy gtt, gbartt and tisspoints matrices to GPU
TWS, August 2014
**************************************************************************/
#include <cuda_runtime.h>

void tissueGPUcopy()
{
	extern int mxx,myy,mzz,nntGPU,nnt,nnv,useGPU,nsp;
	extern int **tisspoints,*h_tisspoints,*d_tisspoints;
	extern float **gtt,*h_gtt,*d_gtt,**gbartt,*h_gbartt,*d_gbartt;
	extern float *axt,*ayt,*azt,**ax;
	extern float *h_tissxyz,*d_tissxyz,*h_vessxyz,*d_vessxyz;

	int jxyz,jtp,ivp,isp;

	for(jtp=0; jtp<nnt; jtp++){
		h_tisspoints[jtp] = tisspoints[1][jtp+1];
		h_tisspoints[jtp+nnt] = tisspoints[2][jtp+1]*mxx;//multiply here to save multiplication in kernel
		h_tisspoints[jtp+2*nnt] = tisspoints[3][jtp+1]*mxx*myy;
		h_tissxyz[jtp] = axt[tisspoints[1][jtp+1]];
		h_tissxyz[jtp+nnt] = ayt[tisspoints[2][jtp+1]];
		h_tissxyz[jtp+2*nnt] = azt[tisspoints[3][jtp+1]];
	}
	for(ivp=0; ivp<nnv; ivp++){
		h_vessxyz[ivp] = ax[1][ivp+1];
		h_vessxyz[ivp+nnv] = ax[2][ivp+1];
		h_vessxyz[ivp+2*nnv] = ax[3][ivp+1];
	}
	for(jxyz=1; jxyz<=nntGPU; jxyz++) for(isp=1; isp<=nsp; isp++){
		h_gtt[jxyz - 1 + (isp - 1)*nntGPU] = gtt[jxyz][isp];
		h_gbartt[jxyz - 1 + (isp - 1)*nntGPU] = gbartt[jxyz][isp];
	}
	cudaSetDevice( useGPU-1 );
	cudaMemcpy(d_tisspoints, h_tisspoints, 3*nnt*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_tissxyz, h_tissxyz, 3*nnt*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_vessxyz, h_vessxyz, 3*nnv*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_gtt, h_gtt, nsp*nntGPU*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_gbartt, h_gbartt, nsp*nntGPU*sizeof(float), cudaMemcpyHostToDevice);
}