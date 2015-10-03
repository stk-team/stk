#include <stk/stk.hpp>

int g_nPts;

extern "C" void module_sampler_init(int nPts)
{
	g_nPts = nPts;
}

extern "C" void module_sampler_exec(stk::PointSet2dd& pts)
{
	stk::sampler::grid(pts, g_nPts);
}

extern "C" void module_sampler_free(int nPts)
{
	
}
