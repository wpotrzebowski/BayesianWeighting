#include "VBW_sc.hh"
#include "LFuncGpu.hh"
#include <cuda.h>
#include <stdio.h>
#define BLOCKSIZE 1024

__global__ void calc(double *x, double *mix_saxs_1d, double alpha_zero, double *fitSaxsMix, int total) {

	double deltamix, smix;

	int i= blockIdx.x * blockDim.x + threadIdx.x;
	int j= blockIdx.y * blockDim.y + threadIdx.y;

	if (i>=total || j>=total) return;

	smix = mix_saxs_1d[i*total+j];
        deltamix = (i!=j) ? -x[i]*x[j] : x[i]*(alpha_zero - x[i]);
	//Most likely atomic add will be needed
	//fitSaxsMix += 1.0;
        fitSaxsMix[i*total+j] = smix * deltamix;

}

__global__ void strideSum(double *f, int len, int strid){
    	int i = threadIdx.x + blockDim.x * blockIdx.x;
    	if(i+strid<len){
        	f[i]=f[i]+f[i+strid];
    	}
}

__global__ void totalSum(double * input, double * output, int len) {
    //@@ Load a segment of the input vector into shared memory
    __shared__ float partialSum[2 * BLOCKSIZE];
    unsigned int t = threadIdx.x, start = 2 * blockIdx.x * BLOCKSIZE;
    if (start + t < len)
       partialSum[t] = input[start + t];
    else
       partialSum[t] = 0;
    if (start + BLOCKSIZE + t < len)
       partialSum[BLOCKSIZE + t] = input[start + BLOCKSIZE + t];
    else
       partialSum[BLOCKSIZE + t] = 0;
    //@@ Traverse the reduction tree
    for (unsigned int stride = BLOCKSIZE; stride >= 1; stride >>= 1) {
       __syncthreads();
       if (t < stride)
          partialSum[t] += partialSum[t+stride];
    }
    //@@ Write the computed sum of the block to the output vector at the 
    //@@ correct index
    if (t == 0)
       output[blockIdx.x] = partialSum[0];
}


void checkCuda(cudaError_t result) {
  	if (result != cudaSuccess) {
    		std::cerr<<"CUDA Runtime Error: "<<cudaGetErrorString(result);
		std::exit(0);
  	}
}


void LFuncLoopGpu (block * x, double * mix_saxs_1d, double alpha_zero, double *fit_saxs_mix, int total) {

	double *h_f;
	double *d_h; 
	double *d_x, *d_m;
	double * deviceOutput;
	double * hostOutput;
	h_f = ( double * ) malloc( sizeof( double ) * total * total ) ;
	checkCuda ( cudaMalloc( (void **)&d_h, sizeof( double ) * total * total));
	checkCuda ( cudaMalloc( (void **)&d_x, sizeof( double ) * total));
	checkCuda ( cudaMalloc( (void **)&d_m, sizeof( double ) * total * total));
		
    	checkCuda ( cudaMemcpy( d_h, h_f, sizeof( double ) * total * total, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy( d_x, x->alphas, sizeof( double ) * total, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy( d_m, mix_saxs_1d, sizeof( double ) * total * total, cudaMemcpyHostToDevice) );
	//int gpu_processor = 1;
	//checkCuda ( cudaSetDevice( gpu_processor ) );
	int numOutputElements;
	numOutputElements = total*total / (BLOCKSIZE<<1);
    	if (total*total % (BLOCKSIZE<<1)) {
        	numOutputElements++;
    	}
	std::cout<<"Number of output elements: "<<numOutputElements<<std::endl;
	hostOutput = (double*) malloc(numOutputElements * sizeof(double));
	checkCuda ( cudaMalloc( (void **)&deviceOutput, numOutputElements *sizeof( double )));

	int gpuBlockSize = 8;
	dim3 GpuBlock( gpuBlockSize,gpuBlockSize );
	dim3 GpuGrid(( total+gpuBlockSize-1)/gpuBlockSize,  (total+gpuBlockSize-1)/gpuBlockSize  );
	calc<<<GpuGrid,GpuBlock>>>(d_x, d_m, alpha_zero, d_h, total);
	cudaDeviceSynchronize();
	checkCuda( cudaMemcpy( h_f, d_h, sizeof( double ) * total * total, cudaMemcpyDeviceToHost) );
	//for( int i=0; i< total*total; i++ ) *fit_saxs_mix+=h_f[i];
	//std::cout<<"Fit saxs mix"<<*fit_saxs_mix<<std::endl;
	/*int i;	
	for(i=total*total; i>1; i=i/2) {
        	strideSum<<<((i/BLOCKSIZE)+1),BLOCKSIZE>>>(d_h,i,i/2);
        	cudaThreadSynchronize();
   	}
	cudaMemcpy(fit_saxs_mix, d_h, 1*sizeof(double) , cudaMemcpyDeviceToHost);*/
	dim3 dimGrid(numOutputElements, 1, 1);
    	dim3 dimBlock(BLOCKSIZE, 1, 1);

    	totalSum<<<dimGrid, dimBlock>>>(d_h, deviceOutput, total*total);
   	cudaDeviceSynchronize();	
	
	cudaMemcpy(hostOutput, deviceOutput, sizeof(double) * numOutputElements, cudaMemcpyDeviceToHost);
	int ii;	
	for (ii = 1; ii < numOutputElements; ii++) {
        	hostOutput[0] += hostOutput[ii];
    	}
	*fit_saxs_mix = hostOutput[0];
	std::cout<<"Fit saxs mix"<<*fit_saxs_mix<<std::endl;

	free(h_f);
	free(hostOutput);
	cudaFree( d_h );
	cudaFree( d_x );
	cudaFree( d_m );
	cudaFree(deviceOutput);
	//cudaDeviceSynchronize();
}
