#ifndef __STK_CONVERSION_POINTSET__
#define __STK_CONVERSION_POINTSET__

#include <stk/pointset.hpp>

namespace stk
{

namespace conversion
{

template<int DIM, typename POS, typename VAL>
void abs(
	const PointSet< DIM, POS, std::complex<VAL> >& i_input,
	PointSet<DIM, POS, VAL>& o_output)
{
	o_output.clear();
	for(int i=0; i<i_input.size(); i++)
	{
		o_output.push_back(
			Point<DIM, POS, VAL>(
				i_input[i].pos(),
				std::abs(i_input[i].val())
			)
		);
	}
	
	o_output.setBoundaries(
		i_input.getMinPosition(),
		i_input.getMaxPosition());
}

}

}


#endif
