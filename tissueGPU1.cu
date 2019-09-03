/***********************************************************
tissueGPU1.cu
GPU kernel to accumulate contributions of tissue source
strengths qt and previoous solute levels ctprev to tissue solute levels ct.
Each tissue point is assigned one or more threads: step is the number of threads
This spreads it over more threads.
TWS September 2014
************************************************************/

__global__ void tissueGPU1Kernel(int *d_tisspoints, float *d_gtt, float *d_gbartt, float *d_ct, float *d_ctprev, float *d_qt, int nnt, int nntGPU, int step, int isp)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
	int itp = i/step;
    int itp1 = i%step;
	int jtp,ixyz,ix,iy,iz,jx,jy,jz,nnt2=2*nnt,istep;
	float p = 0.;
    if(itp < nnt){
		ix = d_tisspoints[itp];
		iy = d_tisspoints[itp+nnt];
		iz = d_tisspoints[itp+nnt2];
		for(jtp=itp1; jtp<nnt; jtp+=step){
			jx = d_tisspoints[jtp];
			jy = d_tisspoints[jtp+nnt];
			jz = d_tisspoints[jtp+nnt2];
			ixyz = abs(jx-ix) + abs(jy-iy) + abs(jz-iz) + (isp - 1)*nntGPU;
			p += d_gtt[ixyz]*d_ctprev[jtp] + d_gbartt[ixyz]*d_qt[jtp];
		}
		if(itp1 == 0) d_ct[itp] = p;
	}
	//The following is apparently needed to assure that d_ct is incremented in sequence from the needed threads
	for(istep=1; istep<step; istep++){
		__syncthreads();
		if(itp1 == istep && itp < nnt) d_ct[itp] += p;
	}
}

extern "C" void tissueGPU1(int *d_tisspoints, float *d_gtt, float *d_gbartt, float *d_ct, float *d_ctprev, float *d_qt, int nnt, int nntGPU, int useGPU, int isp)
{
	int threadsPerBlock = 256;
	int step = 4;//has to be a power of two apparently
	int blocksPerGrid = (step*nnt + threadsPerBlock - 1) / threadsPerBlock;
	tissueGPU1Kernel<<<blocksPerGrid, threadsPerBlock>>>(d_tisspoints,d_gtt,d_gbartt,d_ct,d_ctprev,d_qt,nnt,nntGPU,step,isp);
}
