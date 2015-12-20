// $Id: counter.h 790 2011-01-06 19:10:37Z mag $
#ifndef COUNTER_H_
#define COUNTER_H_

#include "tuple.h"
#include "sample.h"

template <unsigned char N>
class Counter {

    unsigned char _i;

    static unsigned char n;
    static Sample<N>* smp;
    static Tuple<N,char>* idx;

public:

    Counter();

    bool finished() const;

    Counter& increment();
    Counter& reset(bool adv);

    const Sample<N>& sample_offset() const;
    const Tuple<N,char>& index_offset() const;

    static void init();

};

template <unsigned char N>
unsigned char
Counter<N>::n = 1U;

template <unsigned char N>
Tuple<N,char>*
Counter<N>::idx = 0;

template <unsigned char N>
Sample<N>*
Counter<N>::smp = 0;

template <unsigned char N>
inline
Counter<N>::Counter() :
    _i(0U) {}

template <unsigned char N>
inline Counter<N>&
Counter<N>::reset(bool adv)
{ _i = (adv) ? 1U : 0U; return *this; }

template <unsigned char N>
inline bool
Counter<N>::finished() const
{ return _i == n; }

template <unsigned char N>
inline const Tuple<N,char>&
Counter<N>::index_offset() const
{ return idx[_i]; }

template <unsigned char N>
inline const Sample<N>&
Counter<N>::sample_offset() const
{ return smp[_i]; }

template <unsigned char N>
inline Counter<N>&
Counter<N>::increment()
{ _i++; return *this; }

template <unsigned char N>
void
Counter<N>::init()
{

    unsigned char i;

    for (i = 0; i < N; i++)
        n *= 3ULL;

    smp = new Sample<N>[n];
    idx = new Tuple<N,char>[n];

    int8_t c[N];

    for (i = 0; i < N; i++)
        c[i] = 0;

    for (unsigned char j = 0; j < n; j++) {

        unsigned char k;

        for (k = 0; k < N; k++) {
            smp[j].get(k) = static_cast<float>((c[k] + 1)%3 - 1);
            idx[j].get(k) = (c[k] + 1)%3 - 1;
        }

        for (k = 0; k < N; k++) {

            c[k] += 1;

            if (c[k] > 2)
                c[k] = 0;
            else
                break;

        }

    }

}

#endif /*COUNTER_H_*/
