#include <fftw3.h>
#include <stk/fourier/transform-fft.hpp>

namespace stk
{

namespace fourier
{

void _transformFft(
	const Histogram2dc& i_input,
	Histogram2dc& o_output,
	Direction direction)
{
	if(o_output.getSize() != i_input.getSize())
	{
		throw exception::Message("FFTW error : image space size and fourier space size are different", STK_DBG_INFO);
	}
	
	int width = i_input.getSize()[0];
	int height = i_input.getSize()[1];
	
	fftw_complex* in = (fftw_complex*) fftw_malloc(width*height*sizeof(fftw_complex));
	fftw_complex* out = (fftw_complex*) fftw_malloc(width*height*sizeof(fftw_complex));
	
	fftw_plan fftp = fftw_plan_dft_2d(
		width, height,
		in, out, direction, FFTW_MEASURE);
	
	for(int j=0; j<height; j++)
	{
		for(int i=0; i<width; i++)
		{
			in[j*width+i][0] = i_input.getData(Vector2i(i, j)).real();
			in[j*width+i][1] = i_input.getData(Vector2i(i, j)).imag();
		}
	}
	
	fftw_execute(fftp);
	
	Vector2i coord;
	
	for(int j=0; j<height; j++)
	{
		for(int i=0; i<width; i++)
		{
			coord[0] = (i + width/2)%width;
			coord[1] = (j + height/2)%height;
            
			o_output.getData(coord) = Complexd(out[j*width+i][0], out[j*width+i][1]);
		}
	}
	
	Vector2d minPos;
	minPos[0] = -i_input.getSize()[0]/2;
	minPos[1] = -i_input.getSize()[1]/2;
	
	Vector2d maxPos;
	maxPos[0] =  i_input.getSize()[0]/2-1;
	maxPos[1] =  i_input.getSize()[1]/2-1;
	
	o_output.setBoundaries(minPos, maxPos);
	
	fftw_destroy_plan(fftp);
	fftw_free(in);
	fftw_free(out);
}

void _transformFft(
	const Histogram2dd& i_input,
	Histogram2dc& o_output,
	Direction direction)
{
	if(o_output.getSize() != i_input.getSize())
	{
		throw exception::Message("FFTW error : image space size and fourier space size are different", STK_DBG_INFO);
	}
	
	int width = i_input.getSize()[0];
	int height = i_input.getSize()[1];
	
	fftw_complex* in = (fftw_complex*) fftw_malloc(width*height*sizeof(fftw_complex));
	fftw_complex* out = (fftw_complex*) fftw_malloc(width*height*sizeof(fftw_complex));
	
	fftw_plan fftp = fftw_plan_dft_2d(
		width, height,
		in, out, direction, FFTW_MEASURE);
	
	for(int j=0; j<height; j++)
	{
		for(int i=0; i<width; i++)
		{
			in[j*width+i][0] = i_input.getData(Vector2i(i, j));
			in[j*width+i][1] = 0;
		}
	}
	
	fftw_execute(fftp);
	
	Vector2i coord;
	
	for(int j=0; j<height; j++)
	{
		for(int i=0; i<width; i++)
		{
			coord[0] = (i + width/2)%width;
			coord[1] = (j + height/2)%height;
			o_output.getData(coord) = Complexd(out[j*width+i][0], out[j*width+i][1]);
		}
	}
	
	Vector2d minPos;
	minPos[0] = -i_input.getSize()[0]/2;
	minPos[1] = -i_input.getSize()[1]/2;
	
	Vector2d maxPos;
	maxPos[0] =  i_input.getSize()[0]/2-1;
	maxPos[1] =  i_input.getSize()[1]/2-1;
	
	o_output.setBoundaries(minPos, maxPos);
	
	fftw_destroy_plan(fftp);
	fftw_free(in);
	fftw_free(out);
}

}

}
