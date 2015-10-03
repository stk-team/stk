#ifndef __STK_CONVERSION_ANISOTROPY__
#define __STK_CONVERSION_ANISOTROPY__

#include <stk/pointset.hpp>
#include <stk/histogram.hpp>

namespace stk
{

namespace conversion
{

template<int DIM, typename POS, typename VAL>
void anisotropy(
	const PointSet<DIM, POS, VAL>& i_spec,
	Histogram<1, POS, VAL>& i_rspec,
	Histogram<1, POS, VAL>& o_aspec)
{
	o_aspec.fill(0);
	
	POS pos;
	int arrayPos;
	
	o_aspec = Histogram<1, POS, VAL>(i_rspec.getSize());
	Histogram<1, POS, int> count(o_aspec.getSize());
	count.fill(0);
	count.setBoundaries(o_aspec.getMinPosition(), o_aspec.getMaxPosition());
	
	for(int i=0; i<i_rspec.size(); i++)
	{
		try
		{
			pos = i_spec[i].pos().norm();
			if( pos < o_aspec.getMaxPosition()[0] )
			{
				count.get(pos)++;
				o_aspec.get(pos) += pow(i_spec[i].val()-i_rspec.get(pos), 2);
			}
		}
		catch(exception::OutOfRange& e) { }
	}
	
	for(int i=0; i<o_aspec.getSize()[0]; i++)
	{
		o_aspec.getData(i) /= (count.getData(i)-1)*pow(i_rspec.getData(i), 2);
	}
}

template<int DIM, typename POS, typename VAL>
void anisotropy(
	const Histogram<DIM, POS, VAL>& i_spec,
	Histogram<1, POS, VAL>& i_rspec,
	Histogram<1, POS, VAL>& o_aspec)
{
	o_aspec.fill(0);
	
	POS pos;
	int arrayPos;
	
	//o_aspec = Histogram<1, POS, VAL>(i_rspec.getSize());
	Histogram<1, POS, int> count(o_aspec.getSize());
	count.fill(0);
	count.setBoundaries(o_aspec.getMinPosition(), o_aspec.getMaxPosition());
	
	for(int i=0; i<i_spec.getArraySize(); i++)
	{
		try
		{
			pos = i_spec.getPosFromIndice(i).norm();
			if( pos < o_aspec.getMaxPosition()[0] )
			{
				count.get(pos)++;
				o_aspec.get(pos) += pow(i_spec.getFromIndice(i)-i_rspec.get(pos), 2);
			}
		}
		catch(exception::OutOfRange& e) { }
	}
	
	for(int i=0; i<o_aspec.getArraySize(); i++)
		o_aspec[i] = 10*log10(o_aspec[i]/((count[i]-1)*pow(i_rspec[i], 2)));
}

}

}

#endif
