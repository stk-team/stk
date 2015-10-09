#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

#include <stk/stk.hpp>

namespace boostPO = boost::program_options;

double* hdrImg;
int hdrImgWidth;
int hdrImgHeight;
std::string g_hdrFilename;

void loadHdrImg()
{
	//Open file
	std::ifstream file(g_hdrFilename.c_str(), std::ios::binary);
	if(!file) throw stk::exception::FileNotFound(g_hdrFilename, STK_DBG_INFO);
	
	//Check signature
	std::string signature;
	std::getline(file, signature);
	if(signature != "#?RADIANCE") throw stk::exception::Message("This file is not a valid .hdr image", STK_DBG_INFO);
	
	//Skip header
	std::string headerLine;
	bool inHeader = true;
	do
	{
		std::getline(file, headerLine);			
		if(headerLine == "") inHeader = false;
	}
	while(inHeader);
	
	//Get image size
	std::getline(file, headerLine);
	std::stringstream headerStream;
	headerStream << headerLine;
	std::string xTag, yTag;
	int width, height;
	headerStream >> yTag >> height >> xTag >> width;
	
	hdrImg = new double[width*height];
	unsigned char rgbe[4];
		
	for(int j=0; j<height; ++j)
	{
		file.read((char*) rgbe, sizeof(char)*4);
		if(!file)
		{
			std::cerr << "error: " << __LINE__ << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		if ((((int)rgbe[2])<<8 | rgbe[3]) != width)
		{
			std::cerr << "error: " << __LINE__ << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		unsigned char* scanline_buffer = new unsigned char[width*4];
		unsigned char* ptr = scanline_buffer;
		unsigned char buf[2];
		  
		for(int c=0;c<4;c++)
		{
			unsigned char* ptr_end = &scanline_buffer[(c+1)*width];
			while(ptr < ptr_end)
			{
				file.read((char*) buf, sizeof(buf[0])*2);
				if(!file)
				{
					std::cerr << "error: " << __LINE__ << std::endl;
					std::exit(EXIT_FAILURE);
				}
				
				if (buf[0] > 128)
				{
					int count = buf[0]-128;
					if ((count == 0)||(count > ptr_end - ptr))
					{
						std::cerr << "error: " << __LINE__ << std::endl;
						std::exit(EXIT_FAILURE);
					}
					while(count-- > 0) *ptr++ = buf[1];
				}
				else
				{
					int count = buf[0];
					if ((count == 0)||(count > ptr_end - ptr))
					{
						std::cerr << "error: " << __LINE__ << std::endl;
						std::exit(EXIT_FAILURE);
					}
					*ptr++ = buf[1];
					if (--count > 0)
					{
						file.read((char*) ptr, sizeof(*ptr)*count);
						if(!file)
						{
							std::cerr << "error: " << __LINE__ << std::endl;
							std::exit(EXIT_FAILURE);
						}
						ptr += count;
					}
				}
			}
		}
		
		for(int i=0;i<width;i++)
		{
			rgbe[0] = scanline_buffer[i];
			rgbe[1] = scanline_buffer[i+width];
			rgbe[2] = scanline_buffer[i+2*width];
			rgbe[3] = scanline_buffer[i+3*width];
			
			double f = std::ldexp(1.0,rgbe[3]-(int)(128+8));
			double r = rgbe[0] * f;
			double g = rgbe[1] * f;
			double b = rgbe[2] * f;
			hdrImg[j*width+i] = (r+g+b)/3.0;
		}
		
		delete[] scanline_buffer;
	}
	
	hdrImgWidth = width;
	hdrImgHeight = height;
	
	std::cout << hdrImgWidth << "x" << hdrImgHeight << std::endl;
	
	//Close file
	file.close();
}

double hdrFunction(double x, double y)
{
	int px = x*hdrImgHeight;
	int py = y*hdrImgHeight;
	
	return hdrImg[py*hdrImgWidth+px];
}

double diskFunction(double x, double y)
{
	if(x*x+y*y < 1.0/16.0) return 4.0/std::sqrt(M_PI);
	else return 0.0;
}

double gaussianFunction(double x, double y)
{
	const double a = 8.0;
	return a*std::exp(-M_PI*(x*x+y*y)*a*a/2.0);
}

double cosFunction(double x, double y)
{
	return std::sqrt(2.0)*std::cos(16.0 * 2.0 * M_PI * x);
}

typedef double (*IntegrandFunction)(double, double);

int main(int argc, char** argv)
{
	/* ARG PARSER *****************************************************/
	std::vector<std::string> fn_input;
	std::string fn_output;
	std::string integrandName;
	int nPts;
	
	boostPO::variables_map vm;
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h",
			"produce help message")
		("input,i",
			boostPO::value(&fn_input)->composing(),
			"list of pointsets: each file must contains multiple pointsets with the same number of points")
		("integrand,I",
			boostPO::value(&integrandName),
			"integrand function: gaussian, disk or a filename (hdr image)")
		("output,o",
			boostPO::value(&fn_output)->required(),
			"output filename");
	
	boostPO::positional_options_description p;
	p.add("input", -1);
	
	try
	{	
		boostPO::store(
			boostPO::command_line_parser(argc, argv).
				options(desc).positional(p).run(), vm);
		boostPO::notify(vm);
	}
	catch(boost::program_options::error& e)
	{
		std::cerr << e.what() << std::endl;
		std::cout << desc << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if(vm.count("help"))
	{
		std::cout << desc << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	/* INIT ***********************************************************/
	
	srand48(time(NULL));
	
	double integrandShiftX = 0.0f;
	double integrandShiftY = 0.0f;
	
	IntegrandFunction integrandFunction;
	if(integrandName == "gaussian")
	{
		integrandFunction = gaussianFunction;
		integrandShiftX = 0.5;
		integrandShiftY = 0.5;
	}
	else if(integrandName == "disk")
	{
		integrandFunction = diskFunction;
		integrandShiftX = 0.5;
		integrandShiftY = 0.5;
	}
	else if(integrandName == "cos")
	{
		integrandFunction = cosFunction;
		integrandShiftX = 0.0;
		integrandShiftY = 0.0;
	}
	else
	{
		g_hdrFilename = integrandName;
		loadHdrImg();
		integrandFunction = hdrFunction;
		integrandShiftX = 0.0;
		integrandShiftY = 0.0;
	}

	std::ofstream file(fn_output.c_str());
	
	for(int n=0; n<fn_input.size(); ++n)
	{
		double nPts = 0;
		int iter = 0;
		double mean = 0;
		double m2 = 0;
		
		stk::io::PointSetInputStream<2, double, double> stream(fn_input[n]);
		do
		{
			//Read pointset
			stk::PointSet2dd pts;
			stream.read(pts);
			
			bool skip = false;
			
			if(!skip)
			{
				nPts += pts.size();
					
				double integration = 0.0f;
				
				stk::Vector2d shift(drand48(), drand48());
				for(int i=0; i<pts.size(); ++i)
				{
					pts[i].pos() += shift;
				}
				pts.normalize();
			
				for(int i=0; i<pts.size(); i++)
				{
					const int r=0;
					for(int y=-r; y<=r; ++y)
					{
						for(int x=-r; x<=r; ++x)
						{
							integration += integrandFunction(
								pts[i].pos()[0]-integrandShiftX+static_cast<double>(x),
								pts[i].pos()[1]-integrandShiftY+static_cast<double>(y)
							);
						}
					}
				}
				integration /= pts.size();
				
				iter++;
				
				double delta = integration - mean;
				mean += delta/iter;
				m2 += delta*(integration - mean);
			}
		}
		while(stream.next());
		
		stream.close();
		
		if(iter > 500)
		{
			nPts /= iter;
			file << nPts << "\t" << mean << "\t" << m2/(iter-1) << std::endl;
			
			std::cout << iter << " x " << nPts << "pts" << std::endl;
		}
	}
	
	file.close();
	
	exit(EXIT_SUCCESS);
}

