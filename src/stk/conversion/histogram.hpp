#ifndef __STK_CONVERSION_HISTOGRAM__
#define __STK_CONVERSION_HISTOGRAM__

#include <stk/pointset.hpp>
#include <stk/histogram.hpp>

namespace stk
{

namespace conversion
{

template<int DIM, typename POS, typename VAL>
void abs(
	const Histogram< DIM, POS, std::complex<VAL> >& i_input,
	Histogram<DIM, POS, VAL>& o_output)
{
	for(int i=0; i<i_input.getArraySize(); i++)
	{
		o_output.getFromIndice(i) = std::abs(i_input.getFromIndice(i));
	}
	
	o_output.setBoundaries(
		i_input.getMinPosition(),
		i_input.getMaxPosition());
}

template<int DIM, typename POS, typename VAL>
void real(
	const Histogram< DIM, POS, std::complex<VAL> >& i_input,
	Histogram<DIM, POS, VAL>& o_output)
{
	for(int i=0; i<i_input.getArraySize(); i++)
	{
		o_output.getFromIndice(i) = std::real(i_input.getFromIndice(i));
	}
	
	o_output.setBoundaries(
		i_input.getMinPosition(),
		i_input.getMaxPosition());
}

template<int DIM, typename POS, typename VAL>
void imag(
	const Histogram< DIM, POS, std::complex<VAL> >& i_input,
	Histogram<DIM, POS, VAL>& o_output)
{
	for(int i=0; i<i_input.getArraySize(); i++)
	{
		o_output.getFromIndice(i) = std::imag(i_input.getFromIndice(i));
	}
	
	o_output.setBoundaries(
		i_input.getMinPosition(),
		i_input.getMaxPosition());
}

}

}


#endif
