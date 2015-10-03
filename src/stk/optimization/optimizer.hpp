#ifndef __STK_OPTIMIZATION_OPTIMIZER__
#define __STK_OPTIMIZATION_OPTIMIZER__

#include <stk/pointset.hpp>

namespace stk
{

namespace optimization
{

template<int DIM, typename POS, typename VAL>
class Optimizer
{
	public:
		virtual void optimize(PointSet<DIM, POS, VAL>& _pts) = 0;
		void operator() (PointSet<DIM, POS, VAL>& _pts);
};

template<int DIM, typename POS, typename VAL>
void Optimizer<DIM, POS, VAL>::operator() (PointSet<DIM, POS, VAL>& _pts)
{
	optimize(_pts);
}

}

}

#endif
