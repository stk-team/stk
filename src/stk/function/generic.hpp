#ifndef __STK_FUNCTION_GENERIC__
#define __STK_FUNCTION_GENERIC__


#include <stk/vector.hpp>

#include "stk/constant.hpp"

namespace stk {
namespace function {

/**
 * @class GenericFunction
 * @brief Abstract class to derivate specific type of function
 */
template<int DIM, typename POS, typename VAL>
class GenericFunction
{
public:
	/**
	 * @brief Default destructor
	 */
	virtual ~GenericFunction () {};
	
	/**
	 * @brief evaluate function at point/vector
	 * @param i_pos point/vector
	 * @return value at point/vector
	 */
	virtual VAL function(const Vector<DIM, POS>& i_pos) const =0;
	
	/**
	 * @brief evaluate function at point/vector
	 * @param i_pos point/vector
	 * @return value at point/vector
	 */
	virtual VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		this->function(i_pos);
	}
};


/*
template<int DIM, typename POS, typename VAL>
class IntegrableFunction : virtual public GenericFunction<DIM, POS, VAL>
{
public:
	IntegrableFunction()
	: 	GenericFunction<DIM, POS, VAL>()
	{ Integrate(); };
	virtual ~IntegrableFunction ();

	virtual void Integrate() = 0;

protected:
	VAL m_integral;
};
*/

} /* function */
} /* stk */

#endif
