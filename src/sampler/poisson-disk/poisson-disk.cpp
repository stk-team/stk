// $Id: poisson.cpp 790 2011-01-06 19:10:37Z mag $
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <iomanip>

#include <stk/stk.hpp>

#include "generator.h"

#define NDIM 2

int g_numCriteria;

extern "C" void module_sampler_init(int nPts)
{
	g_numCriteria = nPts;
	
	float radiusFactor;
	if(g_numCriteria <= 1024) radiusFactor = 0.5215;
	else if(g_numCriteria >= 4096) radiusFactor = 0.54701;
	else radiusFactor = 0.5215 + (0.54701-0.5215)*((float)g_numCriteria-1024.0f)/(4096.0f-1024.0f);
    
	uint32_t num = g_numCriteria*1.05;

    float rad = sqrtf(radiusFactor/static_cast<float>(M_PI*num));

    Sample<NDIM>::setRadius(2.0f*rad);
    
    std::cout << "Radius: " << 2.0f*rad << std::endl;
	
    Counter<NDIM>::init();
    Interval::seed(static_cast<uint32_t>(time(NULL)));
    BigNum<NDIM>::seed(static_cast<uint32_t>(time(NULL)));
}

extern "C" void module_sampler_exec(stk::PointSet2dd& pts)
{
	pts.clear();
	Generator<NDIM>* gen = 0;
	
	gen = new Generator<NDIM>(false);

	while (gen->iterate());

	gen->output(pts);

	delete gen;
	Sample<NDIM>::numsamples = 0;
	Sample<NDIM>::sample_pool.clear();
}

extern "C" void module_sampler_free()
{
	
}
