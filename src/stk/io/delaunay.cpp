#include <stk/io/delaunay.hpp>
#include <stk/exception.hpp>
#include <stk/color.hpp>

#ifdef CAIRO_ENABLED
#include <cairo/cairo.h>
#endif

namespace stk
{

namespace io
{

void draw(
	std::string i_filename,
	const PointSet2dd& i_pts,
	const std::vector<DelaunayNeighborhood>& i_delaunay,
	const Vector2i& i_size,
	bool colorDistance)
{

#ifdef CAIRO_ENABLED
 #ifdef CAIRO_HAS_PNG_FUNCTIONS

	Vector2d minPos = i_pts.boundingBoxMin();
	Vector2d maxPos = i_pts.boundingBoxMax();
	Vector2d spaceSize = i_pts.boundingBoxSize();

	cairo_surface_t* surface;
	cairo_t* context;
	
	surface = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32,
		i_size[0],
		i_size[1]);
	context = cairo_create(surface);

	//Background
	cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
	cairo_paint(context);

	
	float colorR = 166.0/255.0;
	float colorG = 226.0/255.0;
	float colorB = 0.0/255.0;

	//Cells


	//Centroid
	//~ cairo_set_source_rgba(context, 108.0/255.0, 139.0/255.0, 54.0/255.0, 1.0);
	//~ for(int i=0; i<i_voronoi.size(); i++)
	//~ {
		//~ Vector2d center = i_voronoi[i].getCentroid();
		//~ center *= Vector2d(i_size[0], i_size[1]);
		//~ 
		//~ cairo_arc(context, center[0], center[1], 1, 0, 2*M_PI);
		//~ cairo_fill(context);
	//~ }

	//Points
	//~ if(i_pts.isToroidal())
	//~ {
		//~ double meanDist = 1.0/std::sqrt(i_pts.size());
		//~ double meanDistExt = meanDist*5.0;
		//~ spaceSize += Vector2d(meanDistExt, meanDistExt) * 2.0;
		//~ minPos -= Vector2d(meanDistExt, meanDistExt);
		//~ maxPos += Vector2d(meanDistExt, meanDistExt);
//~ 
		//~ cairo_set_line_width (context, 1.0);
		//~ cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);
		//~ cairo_move_to(context,
			//~ 0.5+i_size[0]*(meanDistExt+i_pts.boundingBoxMin()[0])/spaceSize[0],
			//~ 0.5+i_size[1]*(meanDistExt+i_pts.boundingBoxMin()[1])/spaceSize[1]);
		//~ cairo_line_to(context,
			//~ 0.5+i_size[0]*(meanDistExt+i_pts.boundingBoxMax()[0])/spaceSize[0],
			//~ 0.5+i_size[1]*(meanDistExt+i_pts.boundingBoxMin()[1])/spaceSize[1]);
		//~ cairo_line_to(context,
			//~ 0.5+i_size[0]*(meanDistExt+i_pts.boundingBoxMax()[0])/spaceSize[0],
			//~ 0.5+i_size[1]*(meanDistExt+i_pts.boundingBoxMax()[1])/spaceSize[1]);
		//~ cairo_line_to(context,
			//~ 0.5+i_size[0]*(meanDistExt+i_pts.boundingBoxMin()[0])/spaceSize[0],
			//~ 0.5+i_size[1]*(meanDistExt+i_pts.boundingBoxMax()[1])/spaceSize[1]);
		//~ cairo_line_to(context,
			//~ 0.5+i_size[0]*(meanDistExt+i_pts.boundingBoxMin()[0])/spaceSize[0],
			//~ 0.5+i_size[1]*(meanDistExt+i_pts.boundingBoxMin()[1])/spaceSize[1]);
		//~ cairo_stroke(context);
//~ 
		//~ cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
		//~ cairo_set_line_width (context, 0.5);
				//~ 
		//~ for(int i=0; i<i_delaunay.size(); i++)
		//~ {
			//~ Vector2d center = i_pts[i_delaunay[i].id()].pos();
			//~ bool inSpace;
			//~ 
			//~ if(center[0] >= i_pts.boundingBoxMin()[0] &&
				//~ center[0] <= i_pts.boundingBoxMax()[0] &&
				//~ center[1] >= i_pts.boundingBoxMin()[1] &&
				//~ center[1] <= i_pts.boundingBoxMax()[1])
			//~ {
				//~ inSpace = true;
			//~ }
			//~ else inSpace = false;
//~ 
			//~ if(!inSpace) cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);
			//~ else if(colorDistance) cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			//~ cairo_set_line_width (context, 1.0);
//~ 
			//~ center = (center - minPos) / spaceSize;
			//~ center *= Vector2d(i_size[0], i_size[1]);
			//~ 
			//~ const std::vector<Vector2d>& neighborhood = i_delaunay[i].getNeighborhood();
//~ 
			//~ for(int j=0; j<neighborhood.size(); j++)
			//~ {
				//~ Vector2d neighbor = (neighborhood[j] - minPos) / spaceSize;
				//~ 
				//~ if(inSpace && colorDistance)
				//~ {
					//~ double d = (i_delaunay[i].getPoint() - neighborhood[j]).norm();
					//~ double v = std::min(std::max((d-meanDist)/meanDist, -1.0), 1.0)/2.0 + 0.5;
					//~ v = std::min(std::max(v-0.25, -0.5), 0.5);
					//~ Vector3d c = color::hsv2rgb(Vector3d(v, 0.5, 1.0));
					//~ cairo_set_source_rgba(context, c[0], c[1], c[2], 1.0);
				//~ }
				//~ 
				//~ Vector2d neighbor = (neighborhood[j] - minPos) / spaceSize;
				//~ neighbor *= Vector2d(i_size[0], i_size[1]);
				//~ 
				//~ cairo_move_to(context, center[0], center[1]);
				//~ cairo_line_to(context, neighbor[0], neighbor[1]);
				//~ cairo_stroke(context);
			//~ }
		//~ }
	//~ 
		//~ for(int i=0; i<i_pts.size(); i++)
		//~ {
			//~ for(int j=-1; j<=1; j++)
			//~ {
				//~ for(int k=-1; k<=1; k++)
				//~ {
					//~ if(k == 0 && j == 0)
					//~ {
						//~ cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
					//~ }
					//~ else cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);
//~ 
					//~ Vector2d shift(i_pts.boundingBoxSize()[0]*j, i_pts.boundingBoxSize()[1]*k);
					//~ 
					//~ Vector2d center = (i_pts[i].pos() - minPos + shift) / spaceSize;
					//~ center *= Vector2d(i_size[0], i_size[1]);
					//~ 
					//~ cairo_arc(context, center[0], center[1], 2, 0, 2*M_PI);
					//~ cairo_fill(context);
				//~ }
			//~ }
		//~ }
	//~ }
	//~ else
	//~ {
		//~ cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
		//~ cairo_set_line_width (context, 0.5);
				//~ 
		//~ for(int i=0; i<i_delaunay.size(); i++)
		//~ {
			//~ Vector2d center = i_delaunay[i].getPoint();
			//~ center *= Vector2d(i_size[0], i_size[1]);
			//~ 
			//~ const std::vector<Vector2d>& neighborhood = i_delaunay[i].getNeighborhood();
//~ 
			//~ for(int j=0; j<neighborhood.size(); j++)
			//~ {
				//~ Vector2d neighbor = neighborhood[j];
				//~ neighbor *= Vector2d(i_size[0], i_size[1]);
				//~ 
				//~ cairo_move_to(context, center[0], center[1]);
				//~ cairo_line_to(context, neighbor[0], neighbor[1]);
				//~ cairo_stroke(context);
			//~ }
		//~ }
	//~ 
		//~ cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
	//~ 
		//~ for(int i=0; i<i_pts.size(); i++)
		//~ {
			//~ Vector2d center = (i_pts[i].pos() - minPos) / spaceSize;
			//~ center *= Vector2d(i_size[0], i_size[1]);
			//~ 
			//~ cairo_arc(context, center[0], center[1], 2, 0, 2*M_PI);
			//~ cairo_fill(context);
		//~ }
	//~ }
	
	cairo_surface_write_to_png(surface, i_filename.c_str());
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
	
 #else
	throw exception::Message("Could not use stk::io::draw : PNG support for Cairo is not enabled", STK_DBG_INFO);
 #endif
#else
	throw exception::Message("Could not use stk::io::draw : Cairo is not enabled", STK_DBG_INFO);
#endif

}

}

}
