#include <stk/cuda/tools.hpp>

#define PI 3.141592653589793238462643

__global__ void stk_cuKrnl_fourierTransform(
	const float2* i_srcPos, const float2* i_srcVal, int i_srcSize,
	const float2* i_dstPos, float2* o_dstVal, int i_dstSize,
	float i_dir, float i_normalization)
{
	int i;
	int u = (blockIdx.x * blockDim.x + threadIdx.x);
	float ph;
	float2 res;
	res.x = 0;
	res.y = 0;
	
	if(u<i_dstSize)
	{
		for(i=0; i<i_srcSize; i++)
		{
			ph = i_srcPos[i].x*i_dstPos[u].x + i_srcPos[i].y*i_dstPos[u].y;
			ph *= i_dir*2.0*PI;
			
			res.x += i_srcVal[i].x * cos(ph) - i_srcVal[i].y * sin(ph);
			res.y += i_srcVal[i].y * cos(ph) + i_srcVal[i].x * sin(ph);
		}
		
		o_dstVal[u].x = res.x / i_normalization;
		o_dstVal[u].y = res.y / i_normalization;
	}
}

void stk_cuFourierTransform(
	const float* i_srcPos, const float* i_srcVal, int i_srcSize,
	const float* i_dstPos, float* o_dstVal, int i_dstSize,
	float i_dir, float i_normalization)
{	
	dim3 dimGrid, dimBlock;
	stk_cuGetSizes(dimGrid, dimBlock, i_srcSize, i_dstSize);

	/* INIT ***********************************************************/
	const int srcByteSz = i_srcSize*sizeof(float2);
	const int destByteSz = i_dstSize*sizeof(float2);
	
	float2* dvcSrcPos;
	cudaMalloc((void**) &dvcSrcPos, srcByteSz);
	float2* dvcSrcVal;
	cudaMalloc((void**) &dvcSrcVal, srcByteSz);
	float2* dvcDestPos;
	cudaMalloc((void**) &dvcDestPos, destByteSz);
	float2* dvcDestVal;
	cudaMalloc((void**) &dvcDestVal, destByteSz);
	
	/* START ***********************************************************/
	cudaMemcpy(dvcSrcPos, i_srcPos, srcByteSz, cudaMemcpyHostToDevice);
	cudaMemcpy(dvcSrcVal, i_srcVal, srcByteSz, cudaMemcpyHostToDevice);
	cudaMemcpy(dvcDestPos, i_dstPos, destByteSz, cudaMemcpyHostToDevice);
	
	stk_cuKrnl_fourierTransform <<< dimGrid, dimBlock >>> (
		dvcSrcPos, dvcSrcVal, i_srcSize,
		dvcDestPos, dvcDestVal, i_dstSize,
		i_dir, i_normalization);
		
	cudaMemcpy(o_dstVal, dvcDestVal, destByteSz, cudaMemcpyDeviceToHost);
	
	cudaFree(dvcSrcPos);
	cudaFree(dvcSrcVal);
	cudaFree(dvcDestPos);
	cudaFree(dvcDestVal);
}
