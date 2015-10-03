#ifndef __STK_OPTIMIZATION_GRADIENTDESCENT__
#define __STK_OPTIMIZATION_GRADIENTDESCENT__

#include <stk/pointset.hpp>
#include <stk/optimization/optimizer.hpp>

namespace stk
{

namespace optimization
{

template<int DIM, typename POS, typename VAL>
class GradientDescent : public Optimizer<DIM, POS, VAL>
{
	protected:
		POS m_shift;
		POS m_mean;
		POS m_max;
	
	protected:
		virtual void evalFunction(PointSet<DIM, POS, POS>& _data) const = 0;
	
	public:
		GradientDescent(const POS& sh);
		
		void optimize(PointSet<DIM, POS, VAL>& _pts);
		POS getMeanMotion() const;
		POS getMaxMotion() const;
};

template<int DIM, typename POS, typename VAL>
GradientDescent<DIM, POS, VAL>::GradientDescent(const POS& sh) :
	m_shift(sh)
{
	
}

template<int DIM, typename POS, typename VAL>
void GradientDescent<DIM, POS, VAL>::optimize(PointSet<DIM, POS, VAL>& _pts)
{
	PointSet<DIM, POS, POS> data = _pts;
	
	int iSz = _pts.size();
	
	for(int d=0; d<DIM; d++)
	{
		Vector<DIM, POS> shift;
		shift[d] += m_shift;
			
		for(int i=0; i<iSz; i++)
		{
			data.push_back(Point<DIM, POS, POS>(_pts[i].pos() + shift, 0));
		}
	}
	
	evalFunction(data);

	m_max = 0.0;
	m_mean = 0.0;
	if(iSz != 0)
	{
		for(int i=0; i<iSz; i++)
		{
			Vector<DIM, POS> motion;
			
			for(int d=0; d<DIM; d++)
			{
				POS p = (data[i].val() - data[i+iSz*(d+1)].val()) / m_shift;
				_pts[i].pos()[d] -= p;
				motion[d] = p;
			}
			
			POS norm = motion.norm();
			m_mean += norm;
			if(m_max < norm) m_max = norm;
		}
		m_mean /= iSz;
	}
}

template<int DIM, typename POS, typename VAL>
POS GradientDescent<DIM, POS, VAL>::getMeanMotion() const
{
	return m_mean;
}

template<int DIM, typename POS, typename VAL>
POS GradientDescent<DIM, POS, VAL>::getMaxMotion() const
{
	return m_max;
}

}

}

#endif
