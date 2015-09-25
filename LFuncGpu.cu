#include "VBW_sc.hh"
#include "LFuncGpu.hh"
#include <cuda.h>
#include <stdio.h>

__global__ void calc(double *x, double *mix_saxs_1d, double alpha_zero, float *fitSaxsMix, int total) {

	float deltamix, smix;

	int i= blockIdx.x * blockDim.x + threadIdx.x;
	int j= blockIdx.y * blockDim.y + threadIdx.y;

	if (i>=total || j>=total) return;

	smix = mix_saxs_1d[i*total+j];
        deltamix = (i!=j) ? -x[i]*x[j] : x[i]*(alpha_zero - x[i]);
	//Most likely atomic add will be needed
	//fitSaxsMix += 1.0;
        fitSaxsMix[i] =smix * deltamix;

}

void checkCuda(cudaError_t result) {
  	if (result != cudaSuccess) {
    		std::cerr<<"CUDA Runtime Error: "<<cudaGetErrorString(result);
		std::exit(0);
  	}
}

void LFuncLoopGpu (block * x, double * mix_saxs_1d, double alpha_zero, double fit_saxs_mix, int total) {

	float *h_f;
	float *d_h; 
	double *d_x, *d_m;
	h_f = ( float * ) malloc( sizeof( float ) * total ) ;
	checkCuda ( cudaMalloc( (void **)&d_h, sizeof( float ) * total));
	checkCuda ( cudaMalloc( (void **)&d_x, sizeof( double ) * total));
	checkCuda ( cudaMalloc( (void **)&d_m, sizeof( double ) * total * total));	
    	checkCuda ( cudaMemcpy( d_h, h_f, sizeof( float ) * total, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy( d_x, x->alphas, sizeof( double ) * total, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy( d_m, mix_saxs_1d, sizeof( double ) * total * total, cudaMemcpyHostToDevice) );
	//int gpu_processor = 1;
	//checkCuda ( cudaSetDevice( gpu_processor ) );
	int gpuBlockSize = 32;
	dim3 GpuBlock( gpuBlockSize,gpuBlockSize );
	dim3 GpuGrid (  (total+gpuBlockSize-1)/gpuBlockSize,  (total+gpuBlockSize-1)/gpuBlockSize  );
	calc<<<GpuGrid,GpuBlock>>>(d_x, d_m, alpha_zero, d_h, total);
	cudaDeviceSynchronize();
	checkCuda( cudaMemcpy(h_f, d_h, sizeof( float ) * total, cudaMemcpyDeviceToHost) );
	//std::cout<<"Copy from Device"<<std::endl;
	for(int i=0; i< total; i++) fit_saxs_mix+=h_f[i];
	free(h_f);
	cudaFree( d_h );
	cudaFree( d_x );
	cudaFree( d_m );
	//cudaDeviceSynchronize();
}
