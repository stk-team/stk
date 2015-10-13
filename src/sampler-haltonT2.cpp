#include <cstdlib>
#include <iostream>
#include <string>
#include <ctime>
#include <boost/program_options.hpp>
//~ #include <boost/progress.hpp>
#include <boost/timer/timer.hpp>

#include <stk/stk.hpp>

namespace boostPO = boost::program_options;

template<int BASE>
double generateCoordinate(int n)
{
	double res = 0.0f;
	double dev = 1.0f/BASE;
	
	while(n != 0)
	{
		int d = n % BASE;
		res += d * dev;
		n = (n-d)/BASE;
		dev = dev / BASE;
	}
	
	return res;
}

void haltonSampler(stk::PointSet2dd& pts, int nPts)
{
	for(int i=0; i<nPts; i++)
	{
		pts.push_back(stk::Point2dd(stk::Vector2d(generateCoordinate<2>(i), generateCoordinate<3>(i)), 1.0));
	}
}


int main(int argc, char** argv)
{
	/* ARG PARSER *****************************************************/
	std::string fn_output;
	std::string method;
	int nPts;
	int nPatches;
	
	boostPO::variables_map vm;
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h",
			"produce help message")
		("output,o",
			boostPO::value<std::string>(&fn_output),
			"output filename")
		("nPts,n",
			boostPO::value<int>(&nPts)->default_value(4096),
			"number of point per point set")
		("nPatches,m",
			boostPO::value<int>(&nPatches)->default_value(1),
			"number of point sets")
		("binary,b",
			"write in binary mode")
		;
		
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
	
	if(vm.count("output") == 0)
	{
		std::cerr << "the option '--output' is required but missing" << std::endl;
		std::cout << desc << std::endl;
		exit(EXIT_FAILURE);
	}
	
	/* MAIN ***********************************************************/
	
	try
	{		
		//Open Output
		stk::io::PointSetOutputStream<2, double, double> stream;
		stream.setValueType(stk::io::PointSetStream::VAL_NONE);
		stream.setPositionType(stk::io::PointSetStream::POS_DOUBLE);
		if(vm.count("binary"))
		{
			stream.setBinary(true);
		}
		stream.open(fn_output);
		
		double meanTime = 0;
		double meanPts = 0;
		
		std::cout << "Point set " << 1 << "/" << nPatches << "       \r";
		std::cout.flush();
		
		for(int psCur=0; psCur<nPatches; psCur++)
		{
			stk::PointSet2dd pts;
			
			boost::timer::cpu_timer timer;
	
			//Sampling function
			haltonSampler(pts, nPts);
	
			boost::timer::nanosecond_type timeSystem = timer.elapsed().system;
			boost::timer::nanosecond_type timeUser = timer.elapsed().user;
			boost::timer::nanosecond_type timeWall = timer.elapsed().wall;
			meanTime += (static_cast<double>(timeWall)/1000000000.0 - meanTime)/(psCur+1);
	
			meanPts += (static_cast<double>(pts.size()) - meanPts)/(psCur+1);
			
			std::cout << meanTime << " sec, " << meanPts << " pts" << std::endl;
			
			stream.write(pts);
		}
		
		std::cout << nPatches << (nPatches > 1 ? " point sets" : " point set") << " generated                             " << std::endl;
		
		stream.close();
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
	
	exit(EXIT_SUCCESS);
}
