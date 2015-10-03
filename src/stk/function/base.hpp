#ifndef __STK_FUNCTION_BASE__
#define __STK_FUNCTION_BASE__

#include <stk/vector.hpp>
#include <stk/function/generic.hpp>

namespace stk
{

namespace function
{

template<int DIM, typename POS, typename VAL>
class ConstantFunction : public GenericFunction<DIM, POS, VAL>
{
	private:
		VAL m_lvl;
		
	public:
		ConstantFunction(const VAL& i_lvl) :
			GenericFunction<DIM, POS, VAL>(),
			m_lvl(i_lvl)
		{ }

		VAL function(const Vector<DIM, POS>& i_pos) const
		{
			return m_lvl;
		}
};

template<int DIM, typename POS, typename VAL>
class GaussianFunction : public GenericFunction<DIM, POS, VAL>
{
	private:
		Vector<DIM, POS> m_pos;
		POS m_dev;
		VAL m_mag;
		
	public:
		GaussianFunction(const Vector<DIM, POS>& i_pos, const POS& i_dev, const VAL& i_mag) :
			GenericFunction<DIM, POS, VAL>(),
			m_pos(i_pos),
			m_dev(i_dev),
			m_mag(i_mag)
		{ }

		VAL function(const Vector<DIM, POS>& i_pos) const
		{
			double x = (i_pos-m_pos).norm();
			return m_mag * exp(-pow(x, 2)/(2*m_dev*m_dev));
		};
};

typedef GaussianFunction<1, double, double> GaussianFunction1dd;
typedef GaussianFunction<2, double, double> GaussianFunction2dd;
typedef GaussianFunction<3, double, double> GaussianFunction3dd;


} /* function */

} /* stk */

#endif
