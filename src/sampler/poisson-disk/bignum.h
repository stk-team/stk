// $Id: bignum.h 780 2010-09-23 00:43:33Z mag $
#ifndef BIGNUM_H_
#define BIGNUM_H_

#include <boost/random/mersenne_twister.hpp>

template <unsigned char N>
class BigNum{

    unsigned char _n[3*N];

    static boost::mt19937 rng;
    static const unsigned char M = 3*N;

public:

    BigNum();
    BigNum(unsigned char l);
    BigNum(const BigNum<N>& bn);

    operator bool () const;
    operator double () const;

    BigNum<N> random() const;

    BigNum& operator += (const BigNum<N>& bn);
    BigNum& operator -= (const BigNum<N>& bn);

    bool operator < (const BigNum<N>& bn) const;
    bool operator > (const BigNum<N>& bn) const;
    bool operator <= (const BigNum<N>& bn) const;
    bool operator >= (const BigNum<N>& bn) const;
    bool operator == (const BigNum<N>& bn) const;
    bool operator != (const BigNum<N>& bn) const;

    static void seed(uint32_t s) { rng.seed(s); }

};

template <unsigned char N>
boost::mt19937
BigNum<N>::rng;

template <unsigned char N>
inline
BigNum<N>::BigNum()
{ for (unsigned char i = 0; i < M; i++) _n[i] = 0U; }

template <unsigned char N>
BigNum<N>::BigNum(unsigned char l)
{

    l *= N;

    unsigned char b = l/8;

    for (unsigned char i = 0; i < M; i++)
    {
        if (i == b)
            _n[i] = 1U << l%8;
        else
            _n[i] = 0U;
    }

}

template <unsigned char N>
inline
BigNum<N>::BigNum(const BigNum<N>& bn)
{ for (unsigned char i = 0; i < M; i++) _n[i] = bn._n[i]; }

template <unsigned char N>
BigNum<N>::operator bool () const
{ 
    
    for (unsigned char i = 0; i < M; i++)
    {
        if (_n[i])
            return true;
    }

  return false;

}

template <unsigned char N>
BigNum<N>::operator double () const
{

    static const double bi = 256.0;
    static const double b = pow(2.0,-24*N);

    double v = 0.0;
    double p = 1.0;

    for (unsigned char i = 0; i < M; i++) {
        v += _n[i]*p;
        p *= bi;
    }

    return (v != 0.0) ? v*b : 1.0;

}

template <unsigned char N>
BigNum<N>
BigNum<N>::random() const
{

   static const unsigned char M2 = 2*M;

   BigNum<N> bn;

   if (*this)
   {

      unsigned char b[M2];
      unsigned char i, j, k;

      for (i = 0; i < M2; i++)
         b[i] = 0U;

      unsigned char n = 4U;
      uint32_t random = rng();

      for (i = 0; i < M; i++)
      {

         uint16_t carry = 0U;

         for (j = 0, k = i; j < M; j++, k++)
         {
            unsigned char bk = b[k];
            uint16_t prod = static_cast<uint16_t>(random & 0xff)*_n[j] + carry;
            b[k] += static_cast<unsigned char>(prod & 0xff);
            carry = (prod >> 8) + (b[k] < bk);
         }

         for ( ; k < M2; k++)
         {
            unsigned char bk = b[k];
            b[k] += static_cast<unsigned char>(carry);
            carry = b[k] < bk;
            if (carry == 0)
               break;
         }

         if (--n == 0)
         {
            n = 4U;
            random = rng();
         }
         else
            random >>= 8;

      }

      for (i = 0, j = M; i < M; i++, j++)
         bn._n[i] = b[j];

      assert(bn < *this);

   }
   else
   {

      unsigned char n = 4U;
      uint32_t random = rng();

      for (unsigned char i = 0; i < M; i++)
      {

         bn._n[i] = static_cast<unsigned char>(random & 0xff);

         if (--n == 0)
         {
            n = 4U;
            random = rng();
         }
         else
            random >>= 8;

      }

   }

   return bn;

}

template <unsigned char N>
BigNum<N>&
BigNum<N>::operator += (const BigNum<N>& bn)
{

    bool carry = false;

    for (unsigned char i = 0; i < M; i++) {

        unsigned char ni = _n[i];

        _n[i] += bn._n[i] + carry;

        carry = (_n[i] < ni) || (carry && bn._n[i] == 0xff);

    }

    return *this;

}

template <unsigned char N>
BigNum<N>&
BigNum<N>::operator -= (const BigNum<N>& bn)
{

    if (bn) {

        unsigned char i;
        unsigned char n[M];

        for (i = 0; i < M; i++)
            n[i] = ~bn._n[i];

        i = 0;
        bool carry = n[0] == 0xff;

        n[0]++;

        while (carry && ++i < M) {
            carry = n[i] == 0xff;
            n[i]++;
        }

        carry = false;

        for (i = 0; i < M; i++) {

            unsigned char ni = _n[i];

            _n[i] += n[i] + carry;

            carry = (_n[i] < ni) || (carry && n[i] == 0xff);

        }

    }

    return *this;

}

template <unsigned char N>
bool
BigNum<N>::operator < (const BigNum<N>& bn) const
{ 
    
    for (unsigned char i = M-1; i > 0; i--)
        if (_n[i] != bn._n[i])
            return _n[i] < bn._n[i];

    return _n[0] < bn._n[0];

}

template <unsigned char N>
bool
BigNum<N>::operator > (const BigNum<N>& bn) const
{ 

   for (unsigned char i = M-1; i > 0; i--)
      if (_n[i] != bn._n[i])
         return _n[i] > bn._n[i];

   return _n[0] > bn._n[0];

}

template <unsigned char N>
bool
BigNum<N>::operator <= (const BigNum<N>& bn) const
{ 
    
    for (unsigned char i = M-1; i > 0; i--)
        if (_n[i] != bn._n[i])
            return _n[i] < bn._n[i];
  
    return _n[0] <= bn._n[0];

}

template <unsigned char N>
bool
BigNum<N>::operator >= (const BigNum<N>& bn) const
{ 

   for (unsigned char i = M-1; i > 0; i--)
      if (_n[i] != bn._n[i])
         return _n[i] > bn._n[i];

   return _n[0] >= bn._n[0];

}

template <unsigned char N>
bool
BigNum<N>::operator == (const BigNum<N>& bn) const
{ 
    
    for (unsigned char i = 0; i < M; i++)
        if (_n[i] != bn._n[i])
            return false;

    return true;

}

template <unsigned char N>
bool
BigNum<N>::operator != (const BigNum<N>& bn) const
{ 
    
    for (unsigned char i = 0; i < M; i++)
        if (_n[i] != bn._n[i])
            return true;

    return false;

}

#endif /*BIGNUM_H_*/
