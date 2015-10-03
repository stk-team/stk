#include <stk/stk.hpp>

int g_nPts;

extern "C" void module_sampler_init(int nPts)
{
	g_nPts = nPts;
}

template<int BASE>
float generateCoordinate(int n)
{
	float res = 0.0f;
	float dev = 1.0f/BASE;
	
	while(n != 0)
	{
		int d = n % BASE;
		res += d * dev;
		n = (n-d)/BASE;
		dev = dev / BASE;
	}
	
	return res;
}

extern "C" void module_sampler_exec(stk::PointSet2dd& pts)
{
	for(int i=0; i<g_nPts; i++)
	{
		pts.push_back(stk::Point2dd(stk::Vector2d(generateCoordinate<2>(i), generateCoordinate<3>(i)), 1.0));
	}
}

extern "C" void module_sampler_free(int nPts)
{
	
}
