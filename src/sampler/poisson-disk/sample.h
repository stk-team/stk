// $Id: sample.h 783 2010-09-25 00:26:42Z mag $
#ifndef SAMPLE_H_
#define SAMPLE_H_

#include <cmath>
#include <vector>
#include <limits>
#include <iostream>

#include "tuple.h"

template <unsigned char N>
class Sample: public Tuple<N,float> {
public:

    static float rad;
    static double radsq;
    static uint32_t numsamples;
    static vector< Sample<N> > sample_pool;

    enum Status { In, Out, Over };

    class Iterator: private vector< Sample<N> >::const_iterator
    {

        bool m_periodic;

        Iterator(const typename vector< Sample<N> >::const_iterator& it,
                 bool periodic = false);

        friend class Sample;

    public:

        Iterator();

        const Sample<N>* operator -> () const;

        Iterator& operator ++ (int);

        bool operator != (const Iterator& i) const;

    };

    Sample();
    Sample(const Tuple<N,float>& t);

    Sample& operator += (const Sample& s);

    bool valid() const;
    bool containsRaw(const Tuple<N,float>& t) const;
    bool contains(const Tuple<N,float>& t) const;
    bool contained(const Tuple<N,Interval>& t) const;
    Status intersectsRaw(const Tuple<N,Interval>& t) const;
    Status intersects(const Tuple<N,Interval>& t) const;

    static float getRadius();
    static void setRadius(float r);

    static uint32_t getNumSamples();
    static const Sample<N>& getSample(uint32_t i);
    static uint32_t push(const Sample<N>& s);
    static void pop();

    static Iterator begin(bool periodic = false);
    static Iterator end();

};

template <unsigned char N>
float Sample<N>::rad = 0.0f;

template <unsigned char N>
double Sample<N>::radsq = 0.0;

template <unsigned char N>
uint32_t Sample<N>::numsamples = 0UL;

template <unsigned char N>
vector< Sample<N> > Sample<N>::sample_pool;

template <unsigned char N>
inline
Sample<N>::Sample() :
    Tuple<N,float>(numeric_limits<float>::quiet_NaN()) {}

template <unsigned char N>
inline
Sample<N>::Sample(const Tuple<N,float>& t) :
    Tuple<N,float>(t) {}

template <unsigned char N>
inline Sample<N>&
Sample<N>::operator += (const Sample& s)
{ Tuple<N,float>::operator += (s); return *this; }

template <unsigned char N>
inline ostream&
operator << (ostream& s, const Sample<N>& smp)
{ return s << static_cast< const Tuple<N,float>& >(smp); }

template <unsigned char N>
bool
Sample<N>::valid() const
{ return Tuple<N,float>::get(0) ==
         Tuple<N,float>::get(0); }

template <unsigned char N>
bool
Sample<N>::containsRaw(const Tuple<N,float>& t) const
{

    double d2 = 0.0;

    for (unsigned char i = 0; i < N; i++)
    {
        register double d = Tuple<N,float>::get(i) - t.get(i);
        d2 += d*d;
    }

    return d2 < radsq;

}

template <unsigned char N>
bool
Sample<N>::contains(const Tuple<N,float>& t) const
{
	if(containsRaw(t)) return true;
	
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += -1.0;
		t2.get(1) += -1.0;
		if(containsRaw(t2)) return true;
	}
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += -1.0;
		t2.get(1) += 0.0;
		if(containsRaw(t2)) return true;
	}
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += -1.0;
		t2.get(1) += 1.0;
		if(containsRaw(t2)) return true;
	}
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += 0.0;
		t2.get(1) += -1.0;
		if(containsRaw(t2)) return true;
	}
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += 0.0;
		t2.get(1) += 1.0;
		if(containsRaw(t2)) return true;
	}
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += 1.0;
		t2.get(1) += -1.0;
		if(containsRaw(t2)) return true;
	}
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += 1.0;
		t2.get(1) += 0.0;
		if(containsRaw(t2)) return true;
	}
	{
		Tuple<N,float> t2 = t;
		t2.get(0) += 1.0;
		t2.get(1) += 1.0;
		if(containsRaw(t2)) return true;
	}
	
	return false;
}

template <unsigned char N>
bool
Sample<N>::contained(const Tuple<N,Interval>& t) const
{

    for (unsigned char i = 0; i < N; i++)
        if (!t.get(i).contains(Tuple<N,float>::get(i)))
            return false;

    return true;

}

template <unsigned char N>
typename Sample<N>::Status
Sample<N>::intersectsRaw(const Tuple<N,Interval>& t) const
{
    double dmin = 0.0;
    double dmax = 0.0;

    for (unsigned char i = 0; i < N; i++) {

        float ci = Tuple<N,float>::get(i);

        double mini = t.get(i).min() - ci;

        if (mini > 0.0) {

            if (mini > rad)
                return Out;

            dmin += mini*mini;

            if (dmin > radsq)
                return Out;

            float maxi = static_cast<float>(t.get(i).max() - ci);

            dmax += maxi*maxi;

            continue;

        }

        mini = ci - t.get(i).max();

        if (mini > 0.0) {

            if (mini > rad)
                return Out;

            dmin += mini*mini;

            if (dmin > radsq)
                return Out;

            float maxi = static_cast<float>(ci - t.get(i).min());

            dmax += maxi*maxi;

            continue;

        }

        float maxi = static_cast<float>(max(ci - t.get(i).min(), t.get(i).max() - ci));

        dmax += maxi*maxi;

    }

    return dmax > radsq ? Over : In;

}

template <unsigned char N>
typename Sample<N>::Status
Sample<N>::intersects(const Tuple<N,Interval>& t) const
{
	Sample<N>::Status status = intersectsRaw(t);
	
	if(status == In) return In;
	
	{
		Sample<N> shifted = *this;
		shifted.get(0) += -1.0;
		shifted.get(1) += -1.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	{
		Sample<N> shifted = *this;
		shifted.get(0) += -1.0;
		shifted.get(1) += 0.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	{
		Sample<N> shifted = *this;
		shifted.get(0) += -1.0;
		shifted.get(1) += 1.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	{
		Sample<N> shifted = *this;
		shifted.get(0) += 0.0;
		shifted.get(1) += -1.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	{
		Sample<N> shifted = *this;
		shifted.get(0) += 0.0;
		shifted.get(1) += 1.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	{
		Sample<N> shifted = *this;
		shifted.get(0) += 1.0;
		shifted.get(1) += -1.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	{
		Sample<N> shifted = *this;
		shifted.get(0) += 1.0;
		shifted.get(1) += 0.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	{
		Sample<N> shifted = *this;
		shifted.get(0) += 1.0;
		shifted.get(1) += 1.0;
		Sample<N>::Status status2 = shifted.intersectsRaw(t);
		if(status2 == In) return In;
		else if(status2 == Over) status = Over;
	}
	
	return status;
}

template <unsigned char N>
bool
operator != (const Sample<N>& s1,
             const Sample<N>& s2)
{

    for (unsigned char i = 0; i < N; i++)
        if (s1.get(i) != s2.get(i))
            return true;

    return false;

}

template <unsigned char N>
inline float
Sample<N>::getRadius()
{ return rad; }

template <unsigned char N>
void
Sample<N>::setRadius(float r)
{
    rad = r;
    radsq = static_cast<double>(rad)*
            static_cast<double>(rad);
}

template <unsigned char N>
uint32_t
Sample<N>::getNumSamples()
{ return numsamples; }

template <unsigned char N>
inline const Sample<N>&
Sample<N>::getSample(uint32_t i)
{ return sample_pool.at(i); }

template <unsigned char N>
uint32_t
Sample<N>::push(const Sample<N>& s)
{ 

    static Tuple<N,Interval> cube;

    sample_pool.push_back(s);

    if (s.contained(cube))
        ++numsamples;

    return static_cast<uint32_t>(sample_pool.size() - 1);

}

template <unsigned char N>
inline void
Sample<N>::pop()
{ sample_pool.pop_back(); }

template <unsigned char N>
inline typename Sample<N>::Iterator
Sample<N>::begin(bool periodic)
{ return Iterator(sample_pool.begin(), periodic); }

template <unsigned char N>
inline typename Sample<N>::Iterator
Sample<N>::end()
{ return Iterator(sample_pool.end()); }

template <unsigned char N>
inline
Sample<N>::Iterator::Iterator() :
    m_periodic(false) {}

template <unsigned char N>
inline
Sample<N>::Iterator::Iterator(const typename vector< Sample<N> >::const_iterator& it, bool periodic) :
    vector< Sample<N> >::const_iterator(it), m_periodic(periodic) {}

template <unsigned char N>
inline const Sample<N>*
Sample<N>::Iterator::operator -> () const
{ return vector< Sample<N> >::const_iterator::operator -> (); }

template <unsigned char N>
typename Sample<N>::Iterator&
Sample<N>::Iterator::operator ++(int)
{

    static Tuple<N,Interval> cube;

    vector< Sample<N> >::const_iterator::operator ++ ();

    if (!m_periodic)
        while (*this != sample_pool.end() &&
               !(vector< Sample<N> >::const_iterator::operator * ()).contained(cube))
            vector< Sample<N> >::const_iterator::operator ++ ();

    return *this;

}

template <unsigned char N>
inline bool
Sample<N>::Iterator::operator != (const Iterator& i) const
{ const typename vector< Sample<N> >::const_iterator& it = 
    reinterpret_cast<const typename vector< Sample<N> >::const_iterator&>(i);
  return it != *this; }

#endif /*SAMPLE_H_*/
