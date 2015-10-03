#ifndef __STK_HUNGARIAN_ALGORITHM__
#define __STK_HUNGARIAN_ALGORITHM__

#include "stk/hungarian-algorithm-raw.hpp"
#include "stk/vector.hpp"
#include "stk/pointset.hpp"

namespace stk
{

template<int DIM, typename POS, typename VAL>
void hungarianAlgorithm(
	PointSet<DIM, POS, VAL>& i_worker,
	PointSet<DIM, POS, VAL>& i_task,
	std::map<int, int>& i_link
	)
{
	int width = i_worker.size();
	int height = i_task.size();

	POS* matrix = new POS[width*height];

	for(int j=0; j<height; j++)
		for(int i=0; i<width; i++)
		{
			matrix[j*height+i] = std::pow(i_worker.dist(i_worker[i], i_task[j]), 2.0);
		}

	_HungarianAlgorithm<POS> ha(matrix, width, height);
	ha.solve(i_link);

	delete[] matrix;
}

}

#endif
