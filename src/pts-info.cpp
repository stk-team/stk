#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include <stk/stk.hpp>

namespace boostPO = boost::program_options;

int main(int argc, char** argv)
{
	/* ARG PARSER *****************************************************/
	std::vector<std::string> fn_input;
	
	int hires, res, nPatches, graph_norm;
	double gaussianSigma;
	double domain;
	
	boostPO::variables_map vm;
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h",
			"produce help message")
		("input,i",
			boostPO::value< std::vector<std::string> >(&fn_input)->composing(),
			".pts filenames of distribution")
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
		std::cerr << e.what() << "\n";
		std::cout << desc << "\n";
		exit(EXIT_FAILURE);
	}
	
	if(vm.count("help") || fn_input.size() == 0)
	{
		std::cout << desc << "\n";
		exit(EXIT_SUCCESS);
	}
	
	/* INIT ***********************************************************/
	
	double numberOfPoint = 0.0;
	int numberOfPointSet = 0;
	
	stk::PointSet2dd pts;
	for(int i=0; i<fn_input.size(); i++)
	{
		stk::io::PointSetInputStream<2, double> stream(fn_input[i]);
		
		switch(stream.positionType())
		{
			case stk::io::PointSetStream::POS_DOUBLE:
				std::cout << "Position type : double" << std::endl;
				break;
			case stk::io::PointSetStream::POS_FLOAT:
				std::cout << "Position type : float" << std::endl;
				break;
			case stk::io::PointSetStream::POS_INT:
				std::cout << "Position type : integer (32 bits)" << std::endl;
				break;
			case stk::io::PointSetStream::POS_INT16:
				std::cout << "Position type : integer (16 bits)" << std::endl;
				break;
			case stk::io::PointSetStream::POS_INT8:
				std::cout << "Position type : integer (8 bits)" << std::endl;
				break;
			default:
				std::cout << "Position type : ???" << std::endl;
		}
		
		switch(stream.valueType())
		{
			case stk::io::PointSetStream::VAL_NONE:
				std::cout << "Value type : none" << std::endl;
				break;
			case stk::io::PointSetStream::VAL_DOUBLE:
				std::cout << "Value type : double" << std::endl;
				break;
			case stk::io::PointSetStream::VAL_FLOAT:
				std::cout << "Value type : float" << std::endl;
				break;
			case stk::io::PointSetStream::VAL_INT:
				std::cout << "Value type : integer (32 bits)" << std::endl;
				break;
			case stk::io::PointSetStream::VAL_UINT:
				std::cout << "Value type : unsigned integer (32 bits)" << std::endl;
				break;
			case stk::io::PointSetStream::VAL_UINT24:
				std::cout << "Value type : unsigned integer (24 bits)" << std::endl;
				break;
			case stk::io::PointSetStream::VAL_COMPLEXD:
				std::cout << "Value type : complex (2 x double)" << std::endl;
				break;
			default:
				std::cout << "Value type : ???" << std::endl;
		}
		
		
		do
		{
			//~ if(numberOfPointSet != 0 && numberOfPointSet%100==0)
			//~ {
				//~ std::cout << "\r" << numberOfPointSet << " pointset (~" << numberOfPoint << " pts)         ";
				//~ std::cout.flush();
			//~ }
			
			stream.read(pts);
			numberOfPointSet++;
			numberOfPoint += (pts.size() - numberOfPoint)/numberOfPointSet;
		}
		while(stream.next());
	}
	
	std::cout << "\rNumber of pointset : " << numberOfPointSet << "         " << std::endl;
	std::cout << "Number of point per pointset : " << numberOfPoint << std::endl;
	
	const stk::UnitDomain<2, double>* unitDomain = dynamic_cast<const stk::UnitDomain<2, double>*>(pts.domain());
	const stk::RectangularDomain<2, double>* rectDomain = dynamic_cast<const stk::RectangularDomain<2, double>*>(pts.domain());
	const stk::BaseDomain<double>* baseDomain = dynamic_cast<const stk::BaseDomain<double>*>(pts.domain());
	const stk::HexagonalDomain<double>* hexaDomain = dynamic_cast<const stk::HexagonalDomain<double>*>(pts.domain());
	
	if(unitDomain) std::cout << "Domain : Unit";
	if(rectDomain) std::cout << "Domain : Rectangular";
	if(baseDomain) std::cout << "Domain : Basis";
	if(hexaDomain) std::cout << "Domain : Hexagonal";
	
	if(pts.isToroidal()) std::cout << " (toroidal)";
	std::cout << std::endl;
	
	exit(EXIT_SUCCESS);
}
