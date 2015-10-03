#ifndef __STK_CONVERSION_RADIALMEAN__
#define __STK_CONVERSION_RADIALMEAN__

#include <stk/pointset.hpp>
#include <stk/histogram.hpp>

namespace stk
{

namespace conversion
{

template<int DIM, typename POS, typename VAL>
void radialMean(
	const PointSet<DIM, POS, VAL>& i_spec,
	Histogram<1, POS, VAL>& o_rspec)
{
	o_rspec.fill(0);
	
	POS pos;
	int arrayPos;
	
	Histogram<1, POS, int> count(o_rspec.getSize());
	count.fill(0);
	count.setBoundaries(o_rspec.getMinPosition(), o_rspec.getMaxPosition());
	
	for(int i=0; i<i_spec.size(); i++)
	{
		try
		{
			pos = i_spec[i].pos().norm();
			if( pos < o_rspec.getMaxPosition()[0] )
			{
				count.get(pos)++;
				o_rspec.get(pos) += i_spec[i].val();
			}
		}
		catch(exception::OutOfRange& e) { }
	}
	
	for(int i=0; i<o_rspec.getSize()[0]; i++)
	{
		o_rspec.getData(i) /= count.getData(i);
	}
}

template<int DIM, typename POS, typename VAL>
void radialMean(
	const Histogram<DIM, POS, VAL>& i_spec,
	Histogram<1, POS, VAL>& o_rspec)
{
	o_rspec.fill(0);
	
	POS pos;
	int arrayPos;
	
	Histogram<1, POS, int> count(o_rspec.getSize());
	count.fill(0);
	count.setBoundaries(o_rspec.getMinPosition(), o_rspec.getMaxPosition());
	
	
	for(int i=0; i<i_spec.getArraySize(); i++)
	{
		try
		{
			pos = i_spec.getPosFromIndice(i).norm();
			if( pos < o_rspec.getMaxPosition()[0] )
			{
				count.get(pos)++;
				o_rspec.get(pos) += i_spec.getFromIndice(i);
			}
		}
		catch(exception::OutOfRange& e) { }
	}
	
	for(int i=0; i<o_rspec.getArraySize(); i++)
	{
		o_rspec[i] /= count[i];
	}
}

}

}

#endif
