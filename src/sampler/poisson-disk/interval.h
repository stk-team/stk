// $Id: interval.h 785 2010-09-26 17:50:17Z mag $
#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <iostream>
#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace boost;

class Interval {

    double _min;
    double _max;

    static uniform_01<mt19937> rng;

public:

    Interval();
    Interval(double min,
             double max);

    double min() const;
    double max() const;

    double length() const;
    double midpoint() const;

    float randomise() const;

    bool contains(float x) const;
    bool overlaps(const Interval& i) const;

    void subdivide(Interval i[2]) const;

    friend bool operator < (const Interval& i1,
                            const Interval& i2);
    friend bool operator > (const Interval& i1,
                            const Interval& i2);
    friend bool operator != (const Interval& i1,
                             const Interval& i2);

    static void seed(uint32_t s)
    { rng.base().seed(s); }

};

uniform_01<mt19937>
Interval::rng = uniform_01<mt19937>(mt19937());

inline
Interval::Interval() :
    _min(0.0), _max(1.0) {}

inline
Interval::Interval(double min,
                   double max) :
    _min(min), _max(max) {}

inline double
Interval::min() const
{ return _min; }

inline double
Interval::max() const
{ return _max; }

inline double
Interval::length() const
{ return _max - _min; }

inline double
Interval::midpoint() const
{ return 0.5*(_min + _max); }

inline bool
Interval::contains(float x) const
{ return x >= _min && x <= _max; }

inline bool
Interval::overlaps(const Interval& i) const
{
   return i._max >= _min && i._min <= _max;
}

inline void
Interval::subdivide(Interval i[2]) const
{ i[0]._min = _min;
  i[1]._max = _max;
  i[0]._max = i[1]._min = 0.5*(_min + _max); }

inline float
Interval::randomise() const
{ return static_cast<float>(_min + rng()*(_max - _min)); }

inline bool
operator < (const Interval& i1,
            const Interval& i2)
{ return i1._max <= i2._min; }

inline bool
operator > (const Interval& i1,
            const Interval& i2)
{ return i1._min >= i2._max; }

inline bool
operator != (const Interval& i1,
             const Interval& i2)
{ return i1._min != i2._min ||
         i1._max != i2._max; }

inline std::ostream&
operator << (std::ostream& s,
             const Interval& i)
{ return s << '[' << i.min() << ',' << i.max() << ']'; }

#endif /*INTERVAL_H_*/
