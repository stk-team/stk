#ifndef __STK_FOURIER_CUDA_TRANSFORM__
#define __STK_FOURIER_CUDA_TRANSFORM__

void stk_cuFourierTransform(
	const float* i_srcPos, const float* i_srcVal, int i_srcSize,
	const float* i_dstPos, float* o_dstVal, int i_dstSize,
	float i_dir, float i_normalization);

#endif
