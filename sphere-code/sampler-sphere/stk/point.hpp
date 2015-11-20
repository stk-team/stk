#ifndef __STK_POINT__
#define __STK_POINT__

#include "vector.hpp"

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
};

typedef Point<2, double, double > Point2dd;
typedef Point<2, double, int > Point2di;

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
