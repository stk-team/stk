#ifndef __STK_IO_DELAUNAY__
#define __STK_IO_DELAUNAY__

#include <fstream>
#include <sstream>
#include <vector>

#include <stk/delaunay.hpp>
#include <stk/pointset.hpp>
#include <stk/vector.hpp>

namespace stk
{

namespace io
{

void draw(
	std::string i_filename,
	const PointSet2dd& i_pts,
	const std::vector<DelaunayNeighborhood>& i_delaunay,
	const Vector2i& i_size,
	bool colorDistance = false);

}

}

#endif
