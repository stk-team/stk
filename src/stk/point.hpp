#ifndef __STK_POINT__
#define __STK_POINT__

#include <stk/vector.hpp>
#include <stk/type.hpp>

namespace stk
{
	
/**
 * @class Point
 * @brief A valued sampling point
 */
template<int DIM, typename POS, typename VAL>
class Point
{
	public:
		Vector<DIM, POS> m_pos;
		VAL m_val;

		Point();
		Point(const Vector<DIM, POS>& _pos, const VAL& _val);

		Vector<DIM, POS>& pos();
		const Vector<DIM, POS>& pos() const;

		POS& pos(int i_dim);
		const POS& pos(int i_dim) const;

		VAL& val();
		const VAL& val() const;

		bool operator<(const Point<DIM, POS, VAL>& b) const
		{
			if(m_pos >= b.m_pos) return false;
			if(m_val >= b.m_val) return false;
			return true;
		}

		bool operator>(const Point<DIM, POS, VAL>& b) const
		{
			if(m_pos <= b.m_pos) return false;
			if(m_val <= b.m_val) return false;
			return true;
		}

		bool operator<=(const Point<DIM, POS, VAL>& b) const
		{
			if(m_pos > b.m_pos) return false;
			if(m_val > b.m_val) return false;
			return true;
		}

		bool operator>=(const Point<DIM, POS, VAL>& b) const
		{
			if(m_pos < b.m_pos) return false;
			if(m_val < b.m_val) return false;
			return true;
		}

		bool operator==(const Point<DIM, POS, VAL>& b) const
		{
			return (m_pos == b.m_pos && m_val == b.m_val);
		}

		bool operator!=(const Point<DIM, POS, VAL>& b) const
		{
			return (m_pos != b.m_pos || m_val != b.m_val);
		}
};

typedef Point<2, double, Complexd > Point2dc;
typedef Point<2, double, double > Point2dd;
typedef Point<2, double, int > Point2di;
typedef Point<2, int, int > Point2ii;

template<int DIM, typename POS, typename VAL>
Point<DIM, POS, VAL>::Point()
{
	m_pos = 0;
	m_val = 0;
}

template<int DIM, typename POS, typename VAL>
Point<DIM, POS, VAL>::Point(
	const Vector<DIM, POS>& i_pos,
	const VAL& i_val)
{
	m_pos = i_pos;
	m_val = i_val;
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS>& Point<DIM, POS, VAL>::pos()
{
	return m_pos;
}

template<int DIM, typename POS, typename VAL>
const Vector<DIM, POS>& Point<DIM, POS, VAL>::pos() const
{
	return m_pos;
}

template<int DIM, typename POS, typename VAL>
POS& Point<DIM, POS, VAL>::pos(int i_dim)
{
	return m_pos[i_dim];
}

template<int DIM, typename POS, typename VAL>
const POS& Point<DIM, POS, VAL>::pos(int i_dim) const
{
	return m_pos[i_dim];
}

template<int DIM, typename POS, typename VAL>
VAL& Point<DIM, POS, VAL>::val()
{
	return m_val;
}

template<int DIM, typename POS, typename VAL>
const VAL& Point<DIM, POS, VAL>::val() const
{
	return m_val;
}

}

#endif
