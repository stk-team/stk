#ifndef __STK_SAMPLER_WHITENOISE__
#define __STK_SAMPLER_WHITENOISE__

#include <cmath>

#include <stk/pointset.hpp>
#include <stk/vector.hpp>

namespace stk
{

namespace sampler
{

template<int DIM, typename POS, typename VAL>
void whiteNoise(PointSet<DIM, POS, VAL>& o_pts, const unsigned int i_npts)
{
	const Vector<DIM, POS>& s = o_pts.boundingBoxSize();
	const Vector<DIM, POS>& m = o_pts.boundingBoxMin();

	o_pts.clear();
	o_pts.resize(i_npts);

	for(unsigned int i=0; i<i_npts; i++)
	{
		Vector<DIM, POS>& p = o_pts[i].pos();
		
		do
		{
			for(unsigned int j=0; j<DIM; j++)
			{
				p[j] = drand48()*s[j]+m[j];
			}
		}
		while(!o_pts.hitDomainTest(p));
		
		o_pts[i].val() = o_pts.getDefaultVal();
	}
}

}

}

#endif
