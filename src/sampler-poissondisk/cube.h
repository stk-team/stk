// $Id: cube.h 785 2010-09-26 17:50:17Z mag $
#ifndef CUBE_H_
#define CUBE_H_

#include <cassert>

#include "tuple.h"
#include "sample.h"
#include "interval.h"

template <unsigned char N>
class Cube: public Tuple<N,Interval> {
public:

    Cube();
    Cube(const Tuple<N,Interval>& t);

    Sample<N> randomise() const;

    bool overlaps(const Cube<N>& c) const;
    bool contains(const Sample<N>& s) const;

    Cube<N> child(unsigned char l, unsigned char i) const;

    void subdivide(unsigned char l, Cube<N> c[2]) const;

    friend bool operator == (const Cube<N>& c1,
                             const Cube<N>& c2)
    { const Tuple<N,Interval>& t1 = static_cast< const Tuple<N,Interval>& >(c1);
      const Tuple<N,Interval>& t2 = static_cast< const Tuple<N,Interval>& >(c2);
      return t1 == t2; }

};

template <unsigned char N>
inline
Cube<N>::Cube()
{
}

template <unsigned char N>
inline
Cube<N>::Cube(const Tuple<N,Interval>& t) :
    Tuple<N,Interval>(t) {}

template <unsigned char N>
Sample<N>
Cube<N>::randomise() const
{

    Sample<N> s;

    for (unsigned char i = 0; i < N; i++)
        s.get(i) = Tuple<N,Interval>::get(i).randomise();

    return s;

}

template <unsigned char N>
bool
Cube<N>::overlaps(const Cube<N>& c) const
{

    for (unsigned char i = 0; i < N; i++)
       if (!Tuple<N,Interval>::get(i).overlaps(c.get(i)))
          return false;

    return true;

}

template <unsigned char N>
bool
Cube<N>::contains(const Sample<N>& s) const
{

    for (unsigned char i = 0; i < N; i++)
        if (!Tuple<N,Interval>::get(i).contains(s.get(i)))
            return false;

    return true;

}

template <unsigned char N>
Cube<N>
Cube<N>::child(unsigned char l, unsigned char i) const
{

    l %= N;

    Cube<N> c;

    assert(i < 2);

    for (unsigned char j = 0; j < N; j++)
    {

        if (j == l)
        {

           const Interval& k = Tuple<N,Interval>::get(j);
           double m = k.midpoint();

           if (i)
              c.get(j) = Interval(m, k.max());
           else
              c.get(j) = Interval(k.min(), m);

        }
        else
           c.get(j) = Tuple<N,Interval>::get(j);

    }

    return c;

}

template <unsigned char N>
void
Cube<N>::subdivide(unsigned char l, Cube<N> c[2]) const
{

    l %= N;

    for (unsigned char i = 0; i < N; i++) {

        if (i == l)
        {
           Interval j[2];
           Tuple<N,Interval>::get(i).subdivide(j);
           c[0].get(i) = j[0];
           c[1].get(i) = j[1];
        }
        else
        {
           c[0].get(i) = c[1].get(i) = Tuple<N,Interval>::get(i);
        }

    }

}

#endif
