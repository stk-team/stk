#include <map>
#include <stk/io/array.hpp>
#include <stk/color.hpp>
#include <stk/vector.hpp>

#ifdef CAIRO_ENABLED
#include <cairo/cairo.h>
#endif

namespace stk
{
	
namespace io
{

#ifdef CAIRO_ENABLED
void drawRegion(
	std::string i_filename,
	const Array2i& i_array,
	const Vector2i& i_size)
{
	//Init
	cairo_surface_t* surface;
	cairo_t* context;

	int idMax = i_array.getMax();
	std::map<int, Vector3d> colorScheme;
	for(int i=0; i<i_array.getArraySize(); i++)
	{
		Vector3d c = color::hsv2rgb(Vector3d((double)i_array[i] / (double) idMax, 0.5, 1.0));
		colorScheme[i_array[i]] = c;
	}
	
	surface = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32,
		i_size[0],
		i_size[1]);
	context = cairo_create(surface);
	
	
	//Draw
	int aWidth = i_array.getSize()[0];
	int aHeight = i_array.getSize()[1];

	double stepX = (double)(i_size[0])/(double)(aWidth);
	double stepY = (double)(i_size[1])/(double)(aHeight);
	
	cairo_set_line_width (context, 1.0);
	for(int j=0; j<aHeight; j++)
	{
		for(int i=0; i<aWidth; i++)
		{
			int id = i_array.getData(Vector2i(i, j));
			
			if(id < 0)
			{
				cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			}
			else
			{
				Vector3d mColor = colorScheme[id];
				cairo_set_source_rgba(context, mColor[0], mColor[1], mColor[2], 1.0);
			}
			
			cairo_move_to(context, stepX*i, stepY*j);
			cairo_line_to(context, stepX*(i+1), stepY*j);
			cairo_line_to(context, stepX*(i+1), stepY*(j+1));
			cairo_line_to(context, stepX*i, stepY*(j+1));
			cairo_fill(context);

			cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			
			if(i > 0 && i_array.getData(Vector2i(i-1, j)) != id)
			{
				cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			}
			else cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 0.1);
			
			cairo_move_to(context, stepX*i+0.5, stepY*j);
			cairo_line_to(context, stepX*i+0.5, stepY*(j+1)+1);
			cairo_stroke(context);
			
			if(j > 0 && i_array.getData(Vector2i(i, j-1)) != id)
			{
				cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			}
			else cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 0.1);
			
			cairo_move_to(context, stepX*i, stepY*j+0.5);
			cairo_line_to(context, stepX*(i+1)+1, stepY*j+0.5);
			cairo_stroke(context);
		}
	}
	
	//Export
	cairo_surface_write_to_png(surface, i_filename.c_str());
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
}
#endif

void draw(
	std::string i_filename,
	const Array<2, double>& i_array,
	double i_min,
	double i_max)
{
	int width = i_array.getSize()[0];
	int height = i_array.getSize()[1];
	
	png_bytep* rgb_pixmap = new unsigned char*[height];
	for(int i=0; i<height; i++) rgb_pixmap[i] = new unsigned char[width*3];
	
	double color;
	for(int j=0; j<height; j++)
	{
		for(int i = 0; i<width; i++)
		{
			color = (i_array.getData(Vector2i(i, j))-i_min)/(i_max-i_min);
			if(color < 0.) color = 0.;
			else if(color > 1.) color = 1.;
			
			rgb_pixmap[j][i*3 + 0] = color*255;
			rgb_pixmap[j][i*3 + 1] = color*255;
			rgb_pixmap[j][i*3 + 2] = color*255;
		}
	}
	
	png_structp png_ptr;
	png_infop info_ptr;
	
	FILE *fp = fopen(i_filename.c_str(), "wb");
	if(!fp) throw exception::FileNotWritable(i_filename, STK_DBG_INFO);
	
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr) throw exception::Message("png_create_write_struct failed", STK_DBG_INFO);
	
	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) throw exception::Message("png_create_info_struct failed", STK_DBG_INFO);

	if (setjmp(png_jmpbuf(png_ptr))) throw exception::Message("Error during init_io", STK_DBG_INFO);

	png_init_io(png_ptr, fp);

	/* write header */
	if (setjmp(png_jmpbuf(png_ptr))) throw exception::Message("Error during writing header", STK_DBG_INFO);

	png_set_IHDR(png_ptr, info_ptr, width, height,
		 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);


	/* write bytes */
	if (setjmp(png_jmpbuf(png_ptr))) throw exception::Message("Error during writing bytes", STK_DBG_INFO);

	png_write_image(png_ptr, rgb_pixmap);

	/* end write */
	if (setjmp(png_jmpbuf(png_ptr))) throw exception::Message("Error during end of write", STK_DBG_INFO);

	png_write_end(png_ptr, NULL);
	
	fclose(fp);
	
	for(int i=0; i<height; i++) delete[] rgb_pixmap[i];
	delete[] rgb_pixmap;
}

void writePng(
	std::string i_filename,
	const Array<2, double>& i_array,
	double i_min,
	double i_max)
{
	draw(i_filename, i_array, i_min, i_max);
}

template<int DIM, typename VAL>
void write(
	std::string i_filename,
	const Array<DIM, VAL>& i_array,
	FileType i_type)
{
	switch(i_type)
	{
		case Dat:
			writeDat(i_filename, i_array);
			break;
		default:
			throw exception::Message("Unknown filetype to export array", STK_DBG_INFO);
	}
}

}

}
