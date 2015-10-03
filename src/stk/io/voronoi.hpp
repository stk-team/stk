#ifndef __STK_IO_VORONOI__
#define __STK_IO_VORONOI__

#include <fstream>
#include <sstream>
#include <vector>

#include <stk/voronoi.hpp>
#include <stk/vector.hpp>

namespace stk
{

namespace io
{

void draw(
	std::string i_filename,
	const std::vector<VoronoiCell2d>& i_voronoi,
	const Vector2i& i_size,
	bool connexity = true,
	const stk::Vector2d& i_spaceMin = stk::Vector2d(0.0, 0.0),
	const stk::Vector2d& i_spaceMax = stk::Vector2d(1.0, 1.0));

}

}

#endif
