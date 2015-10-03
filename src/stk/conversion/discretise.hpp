#ifndef __STK_CONVERSION_DISCRETISE__
#define __STK_CONVERSION_DISCRETISE__

#include <stk/pointset.hpp>
#include <stk/histogram.hpp>

namespace stk
{

namespace conversion
{

/**
 * Compute the DFT resolution according to Schlomer and Deussen (Towards a
 * Standardized Spectral Analysis of Point Sets with Applications in Graphics)
 * which is 2 * floor(4 * sqrt(2 / (sqrt(3) * nbPoints)))
 * We approximate sqrt(2 / sqrt(3)) by 1 (instead of 1.07457).
 * @return The DFT resolution very close to what is recommended by Schlomer and
 * Deussen [2010]
 * @remark For PointSets with 2^(2n) points (256, 1024, 4096...), the result is
 * a power of two (which FFTW likes most)
 */
template<int DIM, typename POS, typename VAL>
int optimalDFTResolution(const PointSet<DIM, POS, VAL>& i_pts)
{
    return 2 * (int)round(4.0 * std::sqrt(i_pts.size()));
}

template<typename T>
int optimalDFTResolution(T i_nPts)
{
    return 2 * (int)round(4.0 * std::sqrt(i_nPts));
}

template<int DIM, typename POS, typename VAL>
void discretise(const PointSet<DIM, POS, VAL>& i_pts, Histogram<DIM, POS, VAL>& o_histo)
{
	o_histo.fill(0);
	o_histo.setBoundaries(i_pts.boundingBoxMin(), i_pts.boundingBoxMax());
	
	for(int i=0; i<i_pts.size(); i++)
	{
		try
		{
			o_histo.get(i_pts[i].pos()) += i_pts[i].val();
		}
		catch(exception::OutOfRange& e)
		{
			
		}
	}
}

}

}


#endif
