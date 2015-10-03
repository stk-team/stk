#ifndef __STK_FOURIER_POWERSPECTRUM__
#define __STK_FOURIER_POWERSPECTRUM__

#include <stk/histogram.hpp>

namespace stk
{

namespace fourier
{

template<typename POS, typename VALI, typename VALO>
void powerSpectrum(
	double i_area,
	const Histogram<2, POS, VALI >& i_freq,
	Histogram<2, POS, VALO >& o_power)
{
	if(i_freq.getSize() != o_power.getSize())
	{
		throw exception::Message("input and output size are different", STK_DBG_INFO);
	}
	
	for(int i=0; i<o_power.getArraySize(); i++)
	{
		o_power[i] = pow(std::abs(i_freq[i]), 2)/i_area;
	}
	
	o_power.setBoundaries(i_freq.getMinPosition(), i_freq.getMaxPosition());
}

template<typename POS, typename VALS, typename VALI, typename VALO>
void powerSpectrum(
	const PointSet<2, POS, VALS >& i_space,
	const Histogram<2, POS, VALI >& i_freq,
	Histogram<2, POS, VALO >& o_power)
{
	powerSpectrum(i_space.size(), i_freq, o_power);
}

template<typename POS, typename VALS, typename VALI, typename VALO>
void powerSpectrum(
	const Histogram<2, POS, VALS >& i_space,
	const Histogram<2, POS, VALI >& i_freq,
	Histogram<2, POS, VALO >& o_power)
{
	powerSpectrum(i_space.getSum(), i_freq, o_power);
}

template<typename POS, typename VAL, typename VALS>
void powerSpectrum(
	const PointSet<2, POS, VALS >& i_space,
	const PointSet<2, POS, std::complex<VAL> >& i_freq,
	PointSet<2, POS, VAL>& o_power)
{
	o_power.clear();
	o_power.setBoundaries(i_freq.getMinPosition(), i_freq.getMaxPosition());
	
	double sum = i_space.size();
	
	for(int i=0; i<i_freq.size(); i++)
	{
		VAL val;
		val = std::pow(std::abs(i_freq[i].val()), 2)/sum;
		o_power.push_back(Point<2, POS, VAL>(i_freq[i].pos(), val));
	}
}

template<typename POS, typename VAL, typename VALS>
void powerSpectrum(
	const PointSet<2, POS, VALS>& i_space,
	const PointSet<2, POS, VAL>& i_freq,
	PointSet<2, POS, VAL>& o_power)
{
	o_power.clear();
	o_power.setBoundaries(i_freq.getMinPosition(), i_freq.getMaxPosition());
	
	double sum = i_space.size();
	
	for(int i=0; i<i_freq.size(); i++)
	{
		VAL val;
		val = pow(i_freq[i].val(), 2)/sum;
		o_power.push_back(Point<2, POS, VAL>(i_freq[i].pos(), val));
	}
	
}

}

}

#endif
