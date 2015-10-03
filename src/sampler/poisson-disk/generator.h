// $Id: generator.h 785 2010-09-26 17:50:17Z mag $
#ifndef GENERATOR_H_
#define GENERATOR_H_

#include <cmath>
#include <fstream>
#include <iostream>

#include <stk/stk.hpp>

#include "tree.h"
#include "sample.h"
#include "counter.h"

template <unsigned char N>
class Generator {

    double _v;

    Tree<N>* _t;
    Counter<N>* _c;

    unsigned int _i;

public:

    class Iterator: public Tree<N>::Iterator {

        Iterator(const typename Tree<N>::Iterator& i);

        friend class Generator;

    public:

        Iterator();

        bool operator != (const Iterator& i) const;

    };

    Generator(bool p);
    ~Generator();

    double volume() const;

    Iterator end() const;
    Iterator begin() const;

    bool iterate();

    void spherical_boundary();
    void box_boundary(const Tuple<N,float>& bound);

	int nbSamples() const;
    void output(stk::PointSet2dd& pts) const;
    void output(const std::string& filename) const;

};

template <unsigned char N>
Generator<N>::Generator(bool p) :
    _c(0), _v(0.0), _i(0UL)
{

   _t = new Tree<N>();

   if (p)
      _c = new Counter<N>();

}

template <unsigned char N>
Generator<N>::~Generator()
{
   delete _t;
   delete _c;
}

template <unsigned char N>
void
Generator<N>::spherical_boundary()
{
   if (_t) _t->spherical_boundary();
   delete _c;
   _c = 0;
}

template <unsigned char N>
void
Generator<N>::box_boundary(const Tuple<N,float>& bound)
{
   if (_t) _t->box_boundary(bound);
   delete _c;
   _c = 0;
}

template <unsigned char N>
inline double
Generator<N>::volume() const
{ return _v; }

template <unsigned char N>
inline typename Generator<N>::Iterator
Generator<N>::begin() const
{ return Iterator(_t->begin()); }

template <unsigned char N>
inline typename Generator<N>::Iterator
Generator<N>::end() const
{ return Iterator(_t->end()); }

template <unsigned char N>
inline
Generator<N>::Iterator::Iterator() {}

template <unsigned char N>
inline
Generator<N>::Iterator::Iterator(const typename Tree<N>::Iterator& i) :
    Tree<N>::Iterator(i) {}

template <unsigned char N>
inline bool
Generator<N>::Iterator::operator != (const Iterator& i) const
{ const typename Tree<N>::Iterator& it =
    reinterpret_cast<const typename Tree<N>::Iterator&>(i);
  return it != *this; }

template <unsigned char N>
bool
Generator<N>::iterate()
{
    static uint32_t size0 = 0UL;

    if (!_t->empty()) {

        Sample<N> s(_t->generate());

        assert(_t->valid());

        if (s.valid())
        {

            _t->update(s, Sample<N>::push(s));

            assert(_t->valid());

            if (_c) {

                _c->reset(true);

                while (!_c->finished()) {

                    Sample<N> snew(s);

                    snew += _c->sample_offset();

                    if (!_t->update(snew, Sample<N>::push(snew)))
                        Sample<N>::pop();

                    assert(_t->valid());

                    _c->increment();

                }

            }

        }

        _v = _t->volume();

        uint32_t size1 = static_cast<uint32_t>(100.0*_v);

        if (size1 > size0) {
            size0 = size1;
        }

        ++_i;

        return true;

    }

    _v = 1.0;

    return false;

}

template <unsigned char N>
int
Generator<N>::nbSamples() const
{
    Cube<N> cube;
	int counter=0;
    for (uint32_t i = 0; i < Sample<N>::getNumSamples(); ++i)
    {
       const Sample<N>& sample = Sample<N>::getSample(i);
       if (cube.contains(sample)) counter++;
    }
    
    return counter;
}

template <unsigned char N>
void
Generator<N>::output(stk::PointSet2dd& pts) const
{
    Cube<N> cube;

    for (uint32_t i = 0; i < Sample<N>::getNumSamples(); ++i)
    {
       const Sample<N>& sample = Sample<N>::getSample(i);
       if (cube.contains(sample))
       {
		   pts.push_back(stk::Point2dd(stk::Vector2d(sample.get(0), sample.get(1)), 1));
	   }
    }
}

#endif /*GENERATOR_H_*/
