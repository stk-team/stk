#ifndef __STK_FOURIER_TRANSFORM_CUDA__
#define __STK_FOURIER_TRANSFORM_CUDA__

#ifdef CUDA_ENABLED

#include <stk/fourier/transform-commons.hpp>
#include <stk/fourier/cu-transform.hpp>

namespace stk
{

namespace fourier
{

template<typename VAL>
void _transformCudaVal2Float(const VAL& i_val, float* o_float)
{
	o_float[0] = i_val;
	o_float[1] = 0;
}

template<typename VAL>
void _transformCudaVal2Float(const std::complex<VAL>& i_val, float* o_float)
{
	o_float[0] = real(i_val);
	o_float[1] = imag(i_val);
}

template<typename VAL>
void _transformCudaFloat2Val(const float* o_float, VAL& i_val)
{
	i_val = std::abs(std::complex<float>(o_float[0], o_float[1]));
}

template<typename VAL>
void _transformCudaFloat2Val(const float* o_float, std::complex<VAL>& i_val)
{
	i_val = std::complex<float>(o_float[0], o_float[1]);
}

template<typename POS, typename VALI, typename VALO>
void _transformCuda(
	const PointSet<2, POS, VALI>& i_input,
	PointSet<2, POS, VALO>& _output,
	Direction i_dir)
{
	int srcSize = i_input.size();
	int dstSize = _output.size();

	int blockSize = 1024;
	int nbBlock = dstSize/blockSize;
	int lastBlockSize = 1024;
	
	if(dstSize%blockSize != 0)
	{
		nbBlock++;
		lastBlockSize = dstSize%blockSize;
	}
	
	float* srcPos = new float[srcSize*2];
	float* srcVal = new float[srcSize*2];
	float* dstPos = new float[blockSize*2];
	float* dstVal = new float[blockSize*2];
	
	for(int i=0; i<srcSize; i++)
	{
		srcPos[i*2+0] = i_input[i].pos()[0];
		srcPos[i*2+1] = i_input[i].pos()[1];
		
		_transformCudaVal2Float(i_input[i].val(), srcVal+i*2);
	}

	int currBlockSize = blockSize;
	int arrayShift = 0;
	for(int b=0; b<nbBlock; b++)
	{
		if(b == nbBlock-1) currBlockSize = lastBlockSize;
		
		for(int i=0; i<currBlockSize; i++)
		{
			dstPos[i*2+0] = _output[arrayShift+i].pos()[0];
			dstPos[i*2+1] = _output[arrayShift+i].pos()[1];
		}
		
		stk_cuFourierTransform(
			srcPos, srcVal, srcSize,
			dstPos, dstVal, currBlockSize,
			i_dir, 1.0);
		
		for(int i=0; i<currBlockSize; i++)
		{
			_transformCudaFloat2Val(dstVal+i*2, _output[arrayShift+i].val());
		}

		arrayShift += currBlockSize;
	}
	
	delete[] srcPos;
	delete[] srcVal;
	delete[] dstPos;
	delete[] dstVal;
}

template<typename POS, typename VALI, typename VALO>
void _transformCuda(
	const Histogram<2, POS, VALI>& i_input,
	PointSet<2, POS, VALO>& _output,
	Direction i_dir)
{
	int srcSize = i_input.getArraySize();
	int dstSize = _output.size();

	int blockSize = 1024;
	int nbBlock = dstSize/blockSize;
	int lastBlockSize = 1024;
	
	if(dstSize%blockSize != 0)
	{
		nbBlock++;
		lastBlockSize = dstSize%blockSize;
	}
	
	float* srcPos = new float[srcSize*2];
	float* srcVal = new float[srcSize*2];
	float* dstPos = new float[blockSize*2];
	float* dstVal = new float[blockSize*2];
	
	for(int i=0; i<srcSize; i++)
	{
		srcPos[i*2+0] = i_input.getPosFromIndice(i)[0];
		srcPos[i*2+1] = i_input.getPosFromIndice(i)[1];
		
		_transformCudaVal2Float(i_input.getFromIndice(i), srcVal+i*2);
	}

	int currBlockSize = blockSize;
	int arrayShift = 0;
	for(int b=0; b<nbBlock; b++)
	{		
		if(b == nbBlock-1) currBlockSize = lastBlockSize;
		
		for(int i=0; i<currBlockSize; i++)
		{
			dstPos[i*2+0] = _output[arrayShift+i].pos()[0];
			dstPos[i*2+1] = _output[arrayShift+i].pos()[1];
		}
		
		stk_cuFourierTransform(
			srcPos, srcVal, srcSize,
			dstPos, dstVal, currBlockSize,
			i_dir, 1.0);
		
		for(int i=0; i<currBlockSize; i++)
		{
			_transformCudaFloat2Val(dstVal+i*2, _output[arrayShift+i].val());
		}

		arrayShift += currBlockSize;
	}
	
	delete[] srcPos;
	delete[] srcVal;
	delete[] dstPos;
	delete[] dstVal;
}

template<typename POS, typename VALI, typename VALO>
void _transformCuda(
	const PointSet<2, POS, VALI>& i_input,
	Histogram<2, POS, VALO>& _output,
	Direction i_dir)
{
	int srcSize = i_input.size();
	int dstSize = _output.getArraySize();

	int blockSize = 1024;
	int nbBlock = dstSize/blockSize;
	int lastBlockSize = 1024;
	
	if(dstSize%blockSize != 0)
	{
		nbBlock++;
		lastBlockSize = dstSize%blockSize;
	}
	
	float* srcPos = new float[srcSize*2];
	float* srcVal = new float[srcSize*2];
	float* dstPos = new float[blockSize*2];
	float* dstVal = new float[blockSize*2];
	
	for(int i=0; i<srcSize; i++)
	{
		srcPos[i*2+0] = i_input[i].pos()[0];
		srcPos[i*2+1] = i_input[i].pos()[1];
		
		_transformCudaVal2Float(i_input[i].val(), srcVal+i*2);
	}

	int currBlockSize = blockSize;
	int arrayShift = 0;
	for(int b=0; b<nbBlock; b++)
	{
		if(b == nbBlock-1) currBlockSize = lastBlockSize;
		
		for(int i=0; i<currBlockSize; i++)
		{
			dstPos[i*2+0] = _output.getPosFromIndice(arrayShift+i)[0];
			dstPos[i*2+1] = _output.getPosFromIndice(arrayShift+i)[1];
		}
		
		stk_cuFourierTransform(
			srcPos, srcVal, srcSize,
			dstPos, dstVal, currBlockSize,
			i_dir, 1.0);
		
		for(int i=0; i<currBlockSize; i++)
		{
			_transformCudaFloat2Val(dstVal+i*2, _output.getFromIndice(arrayShift+i));
		}

		arrayShift += currBlockSize;
	}
	
	delete[] srcPos;
	delete[] srcVal;
	delete[] dstPos;
	delete[] dstVal;
}

template<typename POS, typename VALI, typename VALO>
void _transformCuda(
	const Histogram<2, POS, VALI>& i_input,
	Histogram<2, POS, VALO>& _output,
	Direction i_dir)
{
	int srcSize = i_input.getArraySize();
	int dstSize = _output.getArraySize();

	int blockSize = 1024;
	int nbBlock = dstSize/blockSize;
	int lastBlockSize = 1024;
	
	if(dstSize%blockSize != 0)
	{
		nbBlock++;
		lastBlockSize = dstSize%blockSize;
	}
	
	float* srcPos = new float[srcSize*2];
	float* srcVal = new float[srcSize*2];
	float* dstPos = new float[blockSize*2];
	float* dstVal = new float[blockSize*2];
	
	for(int i=0; i<srcSize; i++)
	{
		srcPos[i*2+0] = i_input.getPosFromIndice(i)[0];
		srcPos[i*2+1] = i_input.getPosFromIndice(i)[1];
		
		_transformCudaVal2Float(i_input.getFromIndice(i), srcVal+i*2);
	}

	int currBlockSize = blockSize;
	int arrayShift = 0;
	for(int b=0; b<nbBlock; b++)
	{
		if(b == nbBlock-1) currBlockSize = lastBlockSize;
		
		for(int i=0; i<currBlockSize; i++)
		{
			dstPos[i*2+0] = _output.getPosFromIndice(arrayShift+i)[0];
			dstPos[i*2+1] = _output.getPosFromIndice(arrayShift+i)[1];
		}
		
		stk_cuFourierTransform(
			srcPos, srcVal, srcSize,
			dstPos, dstVal, currBlockSize,
			i_dir, 1.0);
		
		for(int i=0; i<currBlockSize; i++)
		{
			_transformCudaFloat2Val(dstVal+i*2, _output.getFromIndice(arrayShift+i));
		}

		arrayShift += currBlockSize;
	}
	
	delete[] srcPos;
	delete[] srcVal;
	delete[] dstPos;
	delete[] dstVal;
}

}

}

#endif
#endif
