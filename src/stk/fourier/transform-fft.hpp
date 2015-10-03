#ifndef __STK_FOURIER_TRANSFORM_FFT__
#define __STK_FOURIER_TRANSFORM_FFT__

#include <stk/fourier/transform-commons.hpp>

namespace stk
{

namespace fourier
{
 
void _transformFft(
	const Histogram2dc& i_input,
	Histogram2dc& o_output,
	Direction dir);
	
void _transformFft(
	const Histogram2dd& i_input,
	Histogram2dc& o_output,
	Direction dir);

}

}

#endif
