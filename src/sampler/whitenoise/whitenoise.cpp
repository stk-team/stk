#include <stk/stk.hpp>
#include <sys/time.h>

int g_nPts;

extern "C" void module_sampler_init(int nPts)
{
	g_nPts = nPts;
	
	struct timeval time; 
    gettimeofday(&time,NULL);

     // microsecond has 1 000 000
     // Assuming you did not need quite that accuracy
     // Also do not assume the system clock has that accuracy.
    srand48((time.tv_sec * 1000) + (time.tv_usec / 1000));
}

extern "C" void module_sampler_exec(stk::PointSet2dd& pts)
{
	stk::sampler::whiteNoise(pts, g_nPts);
}

extern "C" void module_sampler_free(int nPts)
{
	
}
