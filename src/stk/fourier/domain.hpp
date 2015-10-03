#ifndef __STK_FOURIER_DOMAIN__
#define __STK_FOURIER_DOMAIN__

#include <stk/histogram.hpp>
#include <stk/pointset.hpp>
#include <stk/sampler/grid.hpp>

namespace stk
{

namespace fourier
{

template<int DIM, typename POS, typename VAL>
void domain(
	Histogram<DIM, POS, VAL>& o_freq,
	POS i_freqMax)
{
	Vector2d minPos(-i_freqMax);
	Vector2d maxPos(i_freqMax);
	
	o_freq.setBoundaries(minPos, maxPos);
}

template<int DIM, typename POS, typename VAL>
void domain(
	PointSet<DIM, POS, VAL>& o_freq,
	POS i_freqMax)
{
	o_freq.setBoundaries(
		Vector<DIM, POS>(-i_freqMax),
		Vector<DIM, POS>(i_freqMax)
		);
	
	Vector<DIM, int> npts(i_freqMax*2+1);
	
	o_freq.clear();
	
	Vector<DIM, POS> ind;
	Vector<DIM, POS> step;
	
	const Vector<DIM, POS>& m = o_freq.getMinPosition();
	const Vector<DIM, POS>& M = o_freq.getMaxPosition();
	
	for(int i=0; i<DIM; i++)
	{
		step[i] = (M[i]-m[i])/(npts[i]-1);
		ind[i] = m[i];
	}
	
	int d=0;
	while(d < DIM)
	{
		o_freq.push_back(Point<DIM, POS, VAL>(ind, o_freq.getDefaultVal()));
		
		ind[0] += step[d];
		
		d=0;
		while(ind[d] > M[d] && d < DIM)
		{
			ind[d] = m[d];
			d++;
			if(d < DIM) ind[d] += step[d];
		}
	}
}

template<typename OUTPUT, int DIM, typename POS, typename VAL>
void domain(
	const PointSet<DIM, POS, VAL>& i_space,
	OUTPUT& o_freq,
	POS i_freqMaxRel)
{
	POS freqMax = i_freqMaxRel * pow(i_space.size(), 1.0/DIM);
	
	domain(o_freq, freqMax);
}

}

}

#endif
