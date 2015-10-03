// $Id: tuple.h 790 2011-01-06 19:10:37Z mag $
#ifndef TUPLE_H_
#define TUPLE_H_

#include <vector>
#include <cassert>
#include <iostream>

#include "interval.h"

using namespace std;

template <unsigned char N, typename T>
class Tuple {

    T _v[N];

public:

    Tuple();
    Tuple(const T& x);

    T& get(unsigned char i);
    const T& get(unsigned char i) const;

    Tuple operator - () const;

    template <typename U>
    Tuple& operator += (const Tuple<N,U>& x);

    bool operator == (const Tuple<N,T>& t) const;

};

template <unsigned char N, typename T>
inline
Tuple<N,T>::Tuple() {}

template <unsigned char N, typename T>
inline
Tuple<N,T>::Tuple(const T& x)
{ for (unsigned char i = 0; i < N; i++) _v[i] = x; }

template <unsigned char N, typename T>
inline T&
Tuple<N,T>::get(unsigned char i)
{ assert(i < N); return _v[i]; }

template <unsigned char N, typename T>
inline const T&
Tuple<N,T>::get(unsigned char i) const
{ assert(i < N); return _v[i]; }

template <unsigned char N, typename T>
Tuple<N,T>
Tuple<N,T>::operator - () const
{ Tuple<N,T> t;
  for (unsigned char i = 0; i < N; i++) t._v[i] = -_v[i];
  return t; }

template <unsigned char N, typename T>
template <typename U>
inline Tuple<N,T>&
Tuple<N,T>::operator += (const Tuple<N,U>& x)
{ for (unsigned char i = 0; i < N; i++) _v[i] += x.get(i);
  return *this; }

template <unsigned char N, typename T>
bool
Tuple<N,T>::operator == (const Tuple<N,T>& t) const
{

  for (unsigned char i = 0; i < N; i++)
    if (_v[i] != t._v[i])
      return false;

  return true;

}

template <unsigned char N, typename T>
ostream&
operator << (ostream& s, const Tuple<N,T>& x)
{

    s << '(';

    for (unsigned char i = 0; i < N; i++) {
        s << x.get(i);
        if (i == N - 1)
            break;
        s << ',';
    }

    return s << ')';

}

#endif /*TUPLE_H_*/
