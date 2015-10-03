#include <stk/stk.hpp>
#include <sys/time.h>
#include <boost/timer/timer.hpp>

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
	int nPts = g_nPts;
	
	std::vector<int> arrayX;
	std::vector<int> arrayY;
	
	for(unsigned int i=0; i<nPts; ++i)
	{
		arrayX.push_back(i);
		arrayY.push_back(i);
	}
	
	std::random_shuffle(arrayX.begin(), arrayX.end());
	std::random_shuffle(arrayY.begin(), arrayY.end());
	
	pts.clear();
	
	for(unsigned int k=0; k<nPts; ++k)
	{
		pts.push_back(
			stk::Point2dd(
				stk::Vector2d(
					static_cast<double>(arrayX[k])+drand48(),
					static_cast<double>(arrayY[k])+drand48()
				)/static_cast<double>(nPts),
				1.0
			)
		);
	}
}

extern "C" void module_sampler_free(int nPts)
{
	
}
