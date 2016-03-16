/*
 * Adrien PILLEBOUE © 2012
 *
 * TODO:
 * - Reintroduce the influence of samples' weights 
 * - rework zoom and resolution parameters to better fit what we need
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>

#include <boost/program_options.hpp>

#include <stk/array.hpp>
#include <stk/pointset.hpp>
#include <stk/io/pointset.hpp>
#include <stk/io/array.hpp>

namespace boostPO = boost::program_options;

class ZonePlate
{
	private:
		stk::Array2d m_data;
		stk::Array2d m_dataCount;
		
		int m_res;
		double m_zoom;
		double m_density;
	
	private:	
		double zoneplate_func(const stk::Vector2d& i_pos)
		{
			//Adrien Pilleboue's formula
			double zoom = 60./512.;
			return 0.5 + 0.5 * sin((i_pos.x()*zoom * i_pos.x()*zoom + i_pos.y()*zoom * i_pos.y()*zoom));
			
			//the function used by abdallah in konstanz
			/*double zoom = 1./512.;
			double factor = 512*3.14159;
			return 0.5 + 0.5 * sin(factor*(i_pos.x()*zoom * i_pos.x()*zoom + i_pos.y()*zoom * i_pos.y()*zoom));*/
			
			
		}

		double mitchell_filter_1D(const double i_x)
		{
			const double x = fabs(i_x);

			if (x < 0.5)
				return (21.0 * x - 18.0) * x * x + 2.0;
			else if (x < 1.0)
				return ((-7.0 * x + 18.0) * x - 15.0) * x + 4.0;
			else
				return 0.0;
		}

		double mitchell_filter(const stk::Vector2d& i_pt)
		{
			return mitchell_filter_1D(i_pt.x()) * mitchell_filter_1D(i_pt.y());
		}
		
		double gaussian_filter(const stk::Vector2d& i_pt)
		{
			double sigma = 0.2;
			double a = sigma*sqrt(2*M_PI);
			a = 1.0/a;
			
			double x_factor = (i_pt[0]*i_pt[0]) / (2*sigma*sigma);
			double y_factor = (i_pt[1]*i_pt[1]) / (2*sigma*sigma);
			
			double f_x = a*exp(-(x_factor+y_factor));
			
			return f_x;
		}

		void apply_mitchell(const stk::Vector2d& i_pt)
		{
			const int mRadius = 2;
			int ix, iy, ind;
			double mFact, val;
			stk::Vector2d mDist;

			for(int j=-mRadius; j<=mRadius; j++)
			{
				for(int i=-mRadius; i<=mRadius; i++)
				{
					ix = (int)round(i_pt[0])+i;
					iy = (int)round(i_pt[1])+j;
					mDist = (stk::Vector2d(0.5 + ix, 0.5 + iy) - i_pt) / mRadius;
					if(ix >= 0 && ix < m_res && iy >= 0 && iy < m_res)
					{
						mFact = mitchell_filter(mDist);
						val = zoneplate_func(i_pt*m_zoom)*mFact;
						
						m_data.getData(stk::Vector2i(ix, iy)) += val;
						m_dataCount.getData(stk::Vector2i(ix, iy)) += mFact;
					}
				}
			}
		}
		
		void apply_box(const stk::Vector2d& i_pt)
		{
			const int mRadius = 0;
			int ix, iy, ind;
			double mFact, val;
			stk::Vector2d mDist;

			for(int j=-mRadius; j<=mRadius; j++)
			{
				for(int i=-mRadius; i<=mRadius; i++)
				{
					ix = (int)round(i_pt[0])+i;
					iy = (int)round(i_pt[1])+j;
					mDist = (stk::Vector2d(0.5 + ix, 0.5 + iy) - i_pt) / mRadius;
					if(ix >= 0 && ix < m_res && iy >= 0 && iy < m_res)
					{
						//mFact = mitchell_filter(mDist);
						val = zoneplate_func(i_pt*m_zoom);
						
						m_data.getData(stk::Vector2i(ix, iy)) += val;
						m_dataCount.getData(stk::Vector2i(ix, iy)) += 1;
					}
				}
			}
		}
		
		void apply_gaussian(const stk::Vector2d& i_pt)
		{
			const int mRadius = 3;
			int ix, iy, ind;
			double mFact, val;
			stk::Vector2d mDist;

			for(int j=-mRadius; j<=mRadius; j++)
			{
				for(int i=-mRadius; i<=mRadius; i++)
				{
					ix = (int)round(i_pt[0])+i;
					iy = (int)round(i_pt[1])+j;
					mDist = (stk::Vector2d(0.5 + ix, 0.5 + iy) - i_pt) / mRadius;
					if(ix >= 0 && ix < m_res && iy >= 0 && iy < m_res)
					{
						mFact = gaussian_filter(mDist);
						val = zoneplate_func(i_pt*m_zoom)*mFact;
						
						m_data.getData(stk::Vector2i(ix, iy)) += val;
						m_dataCount.getData(stk::Vector2i(ix, iy)) += mFact;
					}
				}
			}
		}

	public:
		ZonePlate(int i_res, double i_zoomFact, double i_density) :
			m_res(i_res),
			m_density(i_density),
			m_data(stk::Vector2i(i_res, i_res)),
			m_dataCount(stk::Vector2i(i_res, i_res))
		{
			m_zoom = i_zoomFact;// * 60./512.;
		}
		
		~ZonePlate()
		{
			
		}
		
		void make(const stk::PointSet2dd& i_pts)
		{
			m_data.fill(0.0);
			m_dataCount.fill(0.0);
			
			double scale = sqrt(i_pts.size())/m_density;
			int rep = ceil(m_res/scale);
			
			stk::Vector2d pt, ptShift;
			int ix, iy, ind;
			for(int i=0; i<i_pts.size(); i++)
			{
				pt = i_pts[i].pos();
				
				for(int ry=0; ry<rep; ry++)
				{
					for(int rx=0; rx<rep; rx++)
					{
						ptShift = pt + stk::Vector2d(rx*scale, ry*scale);
						//std::cout << "ptShift " << ptShift.x() << ", " << ptShift.y() << std::endl;
						apply_mitchell(ptShift);
						//apply_box(ptShift);
						//apply_gaussian(ptShift);
					}
				}
			}
			
			double tmp0;
			double tmp1;
			for(int j=0; j<m_res/2; j++)
			{
				for(int i=0; i<m_res; i++)
				{
					int count0 = m_dataCount.getData(stk::Vector2i(i, j));
					int count1 = m_dataCount.getData(stk::Vector2i(i, m_res-j-1));
					
					if(count0 && count1)
					{
						tmp0  = m_data.getData(stk::Vector2i(i, j));
						tmp0 /= count0;
						
						tmp1  = m_data.getData(stk::Vector2i(i, m_res-j-1));
						tmp1 /= count1;
						
						m_data.getData(stk::Vector2i(i, m_res-j-1)) = tmp0;
						m_data.getData(stk::Vector2i(i, j)) = tmp1;
					}
				}
			}
		}

		void gamma_correction()
		{
			for (int i = 0; i < m_data.getArraySize(); i++)
			{
				if (m_data[i] <= 0.0031308)
				{
					m_data[i] *= 12.92;
				}
				else
				{
					m_data[i] = 1.055 * pow(m_data[i], 1./2.4) - 0.055;
				}
			}	
		}
		
		const stk::Array2d& getData() const
		{
			return m_data;
		}
};


int main(int argc, char** argv)
{
	srand48(time(NULL));
	
	/* ARG PARSER *****************************************************/
	std::string fn_input;
	std::string fn_output;
	std::string fn_weight;
	double density;
	double zoom;
	int res;
	
	boostPO::variables_map vm;
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h",
			"produce help message")
		("input-file,i",
			boostPO::value<std::string>(&fn_input)->required(),
			".pts.dat filename for initial distribution")
		("output-file,o",
			boostPO::value<std::string>(&fn_output)->required(),
			".png filename for final image")
		("res,R",
			boostPO::value<int>(&res)->default_value(512),
			"file resolution")
		("zoom,z",
			boostPO::value<double>(&zoom)->default_value(1.0),
			"zoom factor")
		("ppp,p",
			boostPO::value<double>(&density)->default_value(4.0),
			"pts per pixel")
		("srgb,s",
			"write corrected values for sRGB colorspace");
	
	boostPO::positional_options_description p;
	
	try
	{	
		boostPO::store(
			boostPO::command_line_parser(argc, argv).
			  options(desc).positional(p).run(), vm);
		boostPO::notify(vm);
	}
	catch(boost::program_options::error& e)
	{
		std::cout << desc << "\n";
		exit(EXIT_FAILURE);
	}
	
	if(vm.count("help"))
	{
		std::cout << desc << "\n";
		exit(EXIT_SUCCESS);
	}
	
	density = sqrt(density);
	
	/* INIT ***********************************************************/
	stk::PointSet2dd pts;
	stk::io::read(fn_input, pts);
	
	double scale = sqrt(pts.size())/density;
	std::cout << "Scaling " << scale << std::endl;
	int rep = ceil(res/scale);
	
	for(int i=0; i<pts.size(); i++)
	{
		pts[i].pos() *= scale;
	}
	
	/* TEST ***********************************************************/
	ZonePlate zp(res, 1/zoom, density);
	
	zp.make(pts);
	
	if(vm.count("srgb"))
	{
		zp.gamma_correction();
	}
	
	stk::io::writePng(fn_output, zp.getData(), 0, 1);
	
	exit(EXIT_SUCCESS);
}
