#ifndef __STK_SAMPLER_GRID__
#define __STK_SAMPLER_GRID__

#include <cmath>

#include <stk/pointset.hpp>

namespace stk
{

namespace sampler
{

template<int DIM, typename POS, typename VAL>
void grid(PointSet<DIM, POS, VAL>& o_pts, const Vector<DIM, int>& i_npts)
{
	o_pts.clear();
	o_pts.reserve(i_npts.volume());

	Vector<DIM, POS> ind;
	Vector<DIM, POS> step;

	const Vector<DIM, POS>& m = o_pts.domain()->boundingBoxMin();
	const Vector<DIM, POS>& M = o_pts.domain()->boundingBoxMax();

	for(int i=0; i<DIM; i++)
	{
		step[i] = (M[i]-m[i])/i_npts[i];
		ind[i] = m[i]+step[i]/2;
	}

	unsigned int d=0;
	while(d < DIM)
	{
		o_pts.push_back(Point<DIM, POS, VAL>(ind, o_pts.getDefaultVal()));

		ind[0] += step[0];
		d=0;
		while(ind[d] > M[d] && d < DIM)
		{
			ind[d] = m[d]+step[d]/2;
			d++;
			if(d < DIM) ind[d] += step[d];
		}
	}
}

template<int DIM, typename POS, typename VAL>
void grid(PointSet<DIM, POS, VAL>& o_pts, const int i_npts)
{
	grid(o_pts, Vector<DIM, int>(pow(i_npts, 1.0/DIM)));
}

}

}

#endif
