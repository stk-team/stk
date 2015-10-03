#ifndef __STK_FOURIER_NORMALIZE__
#define __STK_FOURIER_NORMALIZE__

#include <stk/histogram.hpp>
#include <stk/pointset.hpp>

namespace stk
{

namespace fourier
{

template<int DIM, typename POS, typename VAL>
void normalize(
	double i_fact,
	PointSet<DIM, POS, VAL>& _output)
{
	for(unsigned int i=0; i<_output.size(); i++)
	{
		_output[i].val() /= i_fact;
	}
}

template<int DIM, typename POS, typename VAL>
void normalize(
	double i_fact,
	Histogram<DIM, POS, VAL>& _output)
{
	for(unsigned int i=0; i<_output.getArraySize(); i++)
	{
		_output.getFromIndice(i) /= i_fact;
	}
}

template<int DIM, typename POS, typename VALI, typename OUTPUT>
void normalize(
	const PointSet<DIM, POS, VALI>& i_src,
	OUTPUT& _output)
{
	normalize(static_cast<double>(i_src.size()), _output);
}

template<int DIM, typename POS, typename VALI, typename OUTPUT>
void normalize(
	const Histogram<DIM, POS, VALI>& i_src,
	OUTPUT& _output)
{
	normalize(static_cast<double>(i_src.getArraySize()), _output);
}

}

}

#endif
