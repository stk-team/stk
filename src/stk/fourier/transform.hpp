#ifndef __STK_FOURIER_TRANSFORM__
#define __STK_FOURIER_TRANSFORM__

#include <stk/fourier/transform-commons.hpp>
#include <stk/fourier/transform-cft.hpp>
#include <stk/fourier/transform-fft.hpp>
#include <stk/fourier/transform-cuda.hpp>

namespace stk
{

namespace fourier
{

const Direction forward = -1.0;
const Direction backward = 1.0;

//Main function
template<typename INPUT, typename OUTPUT>
void transform(
	const INPUT& i_input,
	OUTPUT& o_output,
	Direction i_dir);

//Specialization to CFT : * -> PointSet*
template<typename INPUT, int DIM, typename POS, typename VAL>
void transform(
	const INPUT& i_input,
	PointSet<DIM, POS, VAL>& o_output,
	Direction i_dir)
{
	_transformCft(i_input, o_output, i_dir);
}

//Specialization to CFT : * -> Histogram*
template<typename INPUT, int DIM, typename POS, typename VAL>
void transform(
	const INPUT& i_input,
	Histogram<DIM, POS, VAL>& o_output,
	Direction i_dir)
{
	_transformCft(i_input, o_output, i_dir);
}

#ifdef CUDA_ENABLED

//Specialization to Cuda : PointSet2*c -> PointSet2*c
template<typename POS>
void transform(
	const PointSet<2, POS, Complexd>& i_input,
	PointSet<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}

//Specialization to Cuda : PointSet2*d -> PointSet2*c
template<typename POS>
void transform(
	const PointSet<2, POS, double>& i_input,
	PointSet<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}

//Specialization to Cuda : Histogram2*c -> PointSet2*c
template<typename POS>
void transform(
	const Histogram<2, POS, Complexd>& i_input,
	PointSet<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}

//Specialization to Cuda : Histogram2*d -> PointSet2*c
template<typename POS>
void transform(
	const Histogram<2, POS, double>& i_input,
	PointSet<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}

//Specialization to Cuda : PointSet2*c -> Histogram2*c
template<typename POS>
void transform(
	const PointSet<2, POS, Complexd>& i_input,
	Histogram<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}

//Specialization to Cuda : PointSet2*d -> Histogram2*c
template<typename POS>
void transform(
	const PointSet<2, POS, double>& i_input,
	Histogram<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}

//Specialization to Cuda : Histogram2*c -> Histogram2*c
template<typename POS>
void transform(
	const Histogram<2, POS, Complexd>& i_input,
	Histogram<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}

//Specialization to Cuda : Histogram2*d -> Histogram2*c
template<typename POS>
void transform(
	const Histogram<2, POS, double>& i_input,
	Histogram<2, POS, Complexd>& _output,
	Direction i_dir)
{
	_transformCuda(i_input, _output, i_dir);
}
#endif

//Specialization to FFT : Histogram2dc -> Histogram2dc
static inline void transform(
	const Histogram2dc& i_input,
	Histogram2dc& o_output,
	Direction dir)
{
	if(i_input.getSize() != o_output.getSize())
	{
		#ifdef CUDA_ENABLED
		_transformCuda(i_input, o_output, dir);
		#else
		_transformCft(i_input, o_output, dir);
		#endif
	}
	else _transformFft(i_input, o_output, dir);
}
	
//Specialization to FFT : Histogram2dd -> Histogram2dc
static inline void transform(
	const Histogram2dd& i_input,
	Histogram2dc& o_output,
	Direction dir)
{
	if(i_input.getSize() != o_output.getSize())
	{
		#ifdef CUDA_ENABLED
		_transformCuda(i_input, o_output, dir);
		#else
		_transformCft(i_input, o_output, dir);
		#endif
	}
	else _transformFft(i_input, o_output, dir);
}

/*
 ***********************************************************************
 * FORWARD / BACKWARD
 ***********************************************************************
 */

template<typename INPUT, typename OUTPUT>
void transformForward(
	const INPUT& i_input,
	OUTPUT& o_output)
{
	transform(i_input, o_output, forward);
}

template<typename INPUT, typename OUTPUT>
void transformBackward(
	const INPUT& i_input,
	OUTPUT& o_output)
{
	transform(i_input, o_output, backward);
}

}

}

#endif
