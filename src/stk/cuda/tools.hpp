#ifndef __CUDA_CONFIG__
#define __CUDA_CONFIG__

#include <iostream>
#include <cstdlib>
#include <cuda_runtime.h>

void stk_cuInformation()
{
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	
	std::cout << "_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ " << std::endl;
	std::cout << "Device : " << deviceProp.name << std::endl;
	std::cout << "| compute capability = " << deviceProp.major << "." << deviceProp.minor << std::endl;
	std::cout << "| maxGridSize = " << deviceProp.maxGridSize[0] << " x " << deviceProp.maxGridSize[1] << " x " << deviceProp.maxGridSize[2] << std::endl;
	std::cout << "| maxThreadsDim = " << deviceProp.maxThreadsDim[0] << " x " << deviceProp.maxThreadsDim[1] << " x " << deviceProp.maxThreadsDim[2] << std::endl;
	std::cout << "| maxThreadsPerBlock = " << deviceProp.maxThreadsPerBlock << std::endl;
	std::cout << "| multiProcessorCount = " << deviceProp.multiProcessorCount << std::endl;
	std::cout << "| sharedMemPerBlock = " << deviceProp.sharedMemPerBlock << " Octets (" << deviceProp.sharedMemPerBlock/1024 << " Ko)" << std::endl;
	std::cout << "| totalConstMem = " << deviceProp.totalConstMem << " Octets (" << deviceProp.totalConstMem/1024 << " Ko)" << std::endl;
	std::cout << "| totalGlobalMem = " << deviceProp.totalGlobalMem << " Octets (" << deviceProp.totalGlobalMem/(1024*1024) << " Mo)" << std::endl;
	std::cout << "| warpSize = " << deviceProp.warpSize << " threads" << std::endl;
	std::cout << "| _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ " << std::endl << std::endl;
}

void stk_cuGetSizes(dim3& dimGrid, dim3& dimBlock, const int n, const int m)
{
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	
	int used_Gmem = n*16+m*16;
	size_t free_Gmem, tot_Gmem;
	cudaMemGetInfo(&free_Gmem, &tot_Gmem);
	
	if(free_Gmem*0.9 < used_Gmem)
	{
		std::cout << "!! GPU don't have enought global memory" << std::endl;
		exit(EXIT_FAILURE);
	}
	int lim_BxT = deviceProp.maxGridSize[0]*deviceProp.maxThreadsPerBlock;
	if(lim_BxT < n)
	{
		std::cout << "!! GPU reach limit of blocks and threads" << std::endl;
		exit(EXIT_FAILURE);
	}

	//general rules for gpu usage optimization
	int minB = 7*deviceProp.warpSize;
	int maxB = deviceProp.maxThreadsPerBlock;
	int minG = 2*deviceProp.multiProcessorCount;
	int maxG = deviceProp.maxGridSize[0];
	
	//constrain from this particular problem
	int minB_c = (int)ceil(m/float(maxG));
	int maxB_c = (int)ceil(m/float(minG));
	int minG_c = (int)ceil(m/float(maxB));
	int maxG_c = (int)ceil(m/float(minB));

	minB = max(minB, minB_c);
	maxB = min(maxB, maxB_c);
	minG = max(minG, minG_c);
	maxG = min(maxG, maxG_c);
	if(minB>maxB) std::swap(minB, maxB);
	if(minG>maxG) std::swap(minG, maxG);
	
	//~ dimGrid = (minG+maxG)/2;
	//~ dimBlock = (minB+maxB)/2;
	dimGrid = maxG;
	dimBlock = minB;
}

#endif
