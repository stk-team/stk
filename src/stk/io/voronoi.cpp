#include <stk/io/voronoi.hpp>
#include <stk/exception.hpp>

#ifdef CAIRO_ENABLED
#include <cairo/cairo.h>
#endif

namespace stk
{

namespace io
{

void draw(
	std::string i_filename,
	const std::vector<VoronoiCell2d>& i_voronoi,
	const Vector2i& i_size,
	bool connexity,
	const stk::Vector2d& i_spaceMin,
	const stk::Vector2d& i_spaceMax)
{

#ifdef CAIRO_ENABLED
 #ifdef CAIRO_HAS_PNG_FUNCTIONS

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
	
	stk::Vector2d spaceSize = i_spaceMax - i_spaceMin;

	//Cells
	for(int i=0; i<i_voronoi.size(); i++)
	{
		const std::vector<Vector2d>& polygon = i_voronoi[i].getPolygon();
		
		if(polygon.size() > 1)
		{
			Vector2d center = (polygon[polygon.size()-1] - i_spaceMin)/spaceSize;
			center *= Vector2d(i_size[0], i_size[1]);
				
			cairo_move_to(context, center[0], center[1]);
			
			for(int j=0; j<polygon.size(); j++)
			{
				center = (polygon[j] - i_spaceMin)/spaceSize;
				center *= Vector2d(i_size[0], i_size[1]);
				
				cairo_line_to(context, center[0], center[1]);
			}

			if(connexity)
			{
				float aFact = 1.0-std::min(std::max((double)(polygon.size()-3)/6, 0.0), 1.0);
					
				cairo_set_source_rgba(
					context,
					colorR + (1.0-colorR)*aFact,
					colorG + (1.0-colorG)*aFact,
					colorB + (1.0-colorB)*aFact,
					1.0);
				cairo_fill_preserve(context);
			}
			
			cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			cairo_set_line_width (context, 0.5);
			cairo_stroke(context);
		}
	}

	//Centroid
	cairo_set_source_rgba(context, 108.0/255.0, 139.0/255.0, 54.0/255.0, 1.0);
	for(int i=0; i<i_voronoi.size(); i++)
	{
		Vector2d center = (i_voronoi[i].getCentroid() - i_spaceMin)/spaceSize;
		center *= Vector2d(i_size[0], i_size[1]);
		
		cairo_arc(context, center[0], center[1], 1, 0, 2*M_PI);
		cairo_fill(context);
	}

	//Points
	cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
	for(int i=0; i<i_voronoi.size(); i++)
	{
		Vector2d center = (i_voronoi[i].getPoint() - i_spaceMin)/spaceSize;
		center *= Vector2d(i_size[0], i_size[1]);
		
		cairo_arc(context, center[0], center[1], 1.5, 0, 2*M_PI);
		cairo_fill(context);
	}
	
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
