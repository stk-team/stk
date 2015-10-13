#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <map>

#include <boost/program_options.hpp>
#include <stk/stk.hpp>

namespace boostPO = boost::program_options;

namespace stk
{

void validate(
	boost::any& v, 
	const std::vector<std::string>& values,
	Vector2d*, int)
{
	boostPO::validators::check_first_occurrence(v);
	const std::string& s = boostPO::validators::get_single_string(values);

	std::stringstream inString(s);
	double a, b;
	if(inString >> a >> b)
	{
		v = boost::any(stk::Vector2d(a, b));
	}
	else
	{
		throw boostPO::validation_error(boostPO::validation_error::invalid_option_value);
	}
}

}

namespace boostPO = boost::program_options;

class drawPtsOptions
{
	public:
		boostPO::variables_map vm;
		
		std::string fn_input;
		std::vector<std::string> fn_previous;
		std::vector<unsigned int> flag_hot;
		std::vector<unsigned int> flag_cold;
		std::vector<unsigned int> flag_wrong;
		std::vector<unsigned int> flag_fixed;
		std::string fn_output;
		int res;
		int subPatch;
		int zoom_id;
		double zoom_scale;
		stk::Vector2d zoom_coord;
		double radius;
};

template<typename VAL>
void genericPtsDraw(
	stk::io::PointSetInputStream<2, double>& stream,
	const drawPtsOptions& options
)
{
	int subPatchCounter=0;
	stk::PointSet<2, double, VAL> pts;
	do
	{
		if(options.subPatch == subPatchCounter)
		{
			stream.read(pts);
			break;
		}
		else
		{
			stream.ignore();
			subPatchCounter++;
		}
	}
	while(stream.next());
	
	if(options.vm.count("normalize")) pts.normalize();
	
	stk::io::PointSetGraphics<VAL> ptsGraph;
	ptsGraph.pointset(pts);
	
	for(int i=0; i<options.flag_hot.size(); i++)
	{
		ptsGraph.addFlag(options.flag_hot[i], stk::io::PointSetGraphics<VAL>::Hot);
	}
	for(int i=0; i<options.flag_cold.size(); i++)
	{
		ptsGraph.addFlag(options.flag_cold[i], stk::io::PointSetGraphics<VAL>::Cold);
	}
	for(int i=0; i<options.flag_wrong.size(); i++)
	{
		ptsGraph.addFlag(options.flag_wrong[i], stk::io::PointSetGraphics<VAL>::Wrong);
	}
	for(int i=0; i<options.flag_fixed.size(); i++)
	{
		ptsGraph.addFlag(options.flag_fixed[i], stk::io::PointSetGraphics<VAL>::Fixed);
	}
	
	if(options.vm.count("no-periodic")) ptsGraph.borderSize(0);
	int reducedRes = options.res - 2.0*ptsGraph.borderSize();
	int resX = reducedRes * pts.boundingBoxSize()[0]/pts.boundingBoxSize()[1] + 2.0*ptsGraph.borderSize();
	int resY = reducedRes + 2.0*ptsGraph.borderSize();
	
	std::vector<stk::PointSet<2, double, VAL> > previousPts;
	if(options.vm.count("previous-file"))
	{
		previousPts.resize(options.fn_previous.size());
		for(int i=0; i<options.fn_previous.size(); i++)
		{
			stk::io::read(options.fn_previous[i], previousPts[i]);
			ptsGraph.timeTrace(previousPts[i], i);
		}
		ptsGraph.timeTrace(pts, options.fn_previous.size());
	}
	
	if(options.vm.count("delaunay")) ptsGraph.enableDelaunay();
	
	if(options.vm.count("voronoi-connexity")) ptsGraph.enableVoronoi(stk::io::PointSetGraphics<VAL>::VOR_CONNEXITY);
	else if(options.vm.count("voronoi-area")) ptsGraph.enableVoronoi(stk::io::PointSetGraphics<VAL>::VOR_AREA);
	else if(options.vm.count("voronoi-class-bound")) ptsGraph.enableVoronoi(stk::io::PointSetGraphics<VAL>::VOR_CLASS_BOUND);
	else if(options.vm.count("voronoi-class-color")) ptsGraph.enableVoronoi(stk::io::PointSetGraphics<VAL>::VOR_CLASS_COLOR);
	else if(options.vm.count("voronoi")) ptsGraph.enableVoronoi();
	
	if(options.vm.count("radius")) ptsGraph.radiusSize(options.radius);
	
	if(options.vm.count("previous-dot")) ptsGraph.enableTimeTracePoint();
	
	if(options.vm.count("label-id")) ptsGraph.enableLabel(stk::io::PointSetGraphics<VAL>::Id);
	else if(options.vm.count("label-value") || options.vm.count("label-valuei") || options.vm.count("label-valuef"))
	{
		ptsGraph.enableLabel(stk::io::PointSetGraphics<VAL>::Value);
	}
	
	if(options.vm.count("zoom-coord")) ptsGraph.enableZoom(options.zoom_coord, options.zoom_scale);
	else if(options.vm.count("zoom-point") && options.zoom_id < pts.size() && options.zoom_id >= 0)
		ptsGraph.enableZoom(pts[options.zoom_id].pos(), options.zoom_scale);
	
	ptsGraph.draw(options.fn_output, stk::Vector2i(resX, resY));
}

int main(int argc, char** argv)
{
	srand48(time(NULL));
	drawPtsOptions options;
	
	/* ARG PARSER *****************************************************/
	
	//Option list :
	// ...d...hi....nop.....v....
	// .............N.P.R...V....
	
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h",
			"produce help message")
		("input-file,i",
			boostPO::value<std::string>(&options.fn_input)->required(),
			".pts.dat filename for initial distribution")
		("previous-file,p",
			boostPO::value< std::vector<std::string> >(&options.fn_previous)->multitoken(),
			".pts.dat filename for initial distribution")
		("flag-hot",
			boostPO::value< std::vector<unsigned int> >(&options.flag_hot)->multitoken(),
			"List of hot point")
		("flag-cold",
			boostPO::value< std::vector<unsigned int> >(&options.flag_cold)->multitoken(),
			"List of cold point")
		("flag-fixed",
			boostPO::value< std::vector<unsigned int> >(&options.flag_fixed)->multitoken(),
			"List of fixed point")
		("flag-wrong",
			boostPO::value< std::vector<unsigned int> >(&options.flag_wrong)->multitoken(),
			"List of wrong point")
		("output-file,o",
			boostPO::value<std::string>(&options.fn_output)->required(),
			".png filename for final image")
		("res,R",
			boostPO::value<int>(&options.res)->default_value(512),
			"file resolution")
		("subpatch,N",
			boostPO::value<int>(&options.subPatch)->default_value(0),
			"")
		("voronoi,v",
			"draw voronoi cells")
		("delaunay,d",
			"draw delaunay triangulation")
		("voronoi-connexity,V",
			"draw voronoi cells with connexity color code")
		("voronoi-area",
			"draw voronoi cells with area")
		("voronoi-class-bound",
			"draw voronoi cells with class color")
		("voronoi-class-color",
			"draw voronoi cells with boundaries between class")
		("normalize,n",
			"normalize pointset according to the domain")
		("no-periodic,P",
			"dont show pointset periodic replication")
		("previous-dot",
			"draw previous point")
		("label-id",
			"Show id of each point")
		("label-value",
			"Show value of each point")
		("zoom-point",
			boostPO::value<int>(&options.zoom_id),
			"Zoom on a specific point")
		("zoom-scale",
			boostPO::value<double>(&options.zoom_scale)->default_value(1.0),
			"Zoom scale factor")
		("zoom-coord",
			boostPO::value<stk::Vector2d>(&options.zoom_coord),
			"Zoom on a specific coordinate. Format : \"x y\"")
		("radius",
			boostPO::value<double>(&options.radius)->default_value(1.0),
			"Radius of each point");
	
	boostPO::positional_options_description p;
	
	try
	{	
		boostPO::store(
			boostPO::command_line_parser(argc, argv).
			  options(desc).positional(p).run(), options.vm);
		boostPO::notify(options.vm);
	}
	catch(boost::program_options::error& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << desc << "\n";
		exit(EXIT_FAILURE);
	}
	
	if(options.vm.count("help"))
	{
		std::cout << desc << "\n";
		exit(EXIT_SUCCESS);
	}
	
	/* INIT ***********************************************************/
	try
	{
		stk::io::PointSetInputStream<2, double> stream(options.fn_input);
		
		switch(stream.valueType())
		{
			case stk::io::PointSetStream::VAL_INT:
				genericPtsDraw<int>(stream, options);
				break;
			case stk::io::PointSetStream::VAL_UINT:
			case stk::io::PointSetStream::VAL_UINT24:
				genericPtsDraw<unsigned int>(stream, options);
				break;
			case stk::io::PointSetStream::VAL_DOUBLE:
			case stk::io::PointSetStream::VAL_FLOAT:
			case stk::io::PointSetStream::VAL_NONE:
			default:
				genericPtsDraw<double>(stream, options);
		}
	}
	catch(std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
	
	exit(EXIT_SUCCESS);
}
