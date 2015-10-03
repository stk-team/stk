#include <cstdlib>
#include <iostream>
#include <string>
#include <ctime>
#include <boost/program_options.hpp>
//~ #include <boost/progress.hpp>
#include <boost/timer/timer.hpp>

#include <stk/stk.hpp>

#include <dlfcn.h>

namespace boostPO = boost::program_options;

typedef void (*ModuleSamplerInit)(int);
typedef void (*ModuleSamplerExec)(stk::PointSet2dd&);
typedef void (*ModuleSamplerFree)();
		
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
		("sampler,s",
			boostPO::value<std::string>(&method),
			"sampling method")
		("nPts,n",
			boostPO::value<int>(&nPts)->default_value(4096),
			"number of point per point set")
		("nPatches,m",
			boostPO::value<int>(&nPatches)->default_value(1),
			"number of point sets")
		("sampler-list,l",
			"list of all available sampler")
		("binary,b",
			"write in binary")
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
	
	if(vm.count("sampler-list"))
	{
		std::cout << "Sampler list :" << std::endl;
		std::cout << " * whitenoise" << std::endl;
		std::cout << " * grid" << std::endl;
		std::cout << " * stratified" << std::endl;
		std::cout << " * poisson-disk" << std::endl;
		std::cout << " * ccvt" << std::endl;
		std::cout << " * fpo" << std::endl;
		std::cout << " * sobol" << std::endl;
		std::cout << " * halton" << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	if(vm.count("output") == 0)
	{
		std::cerr << "the option '--output' is required but missing" << std::endl;
		std::cout << desc << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if(vm.count("sampler") == 0)
	{
		std::cerr << "the option '--sampler' is required but missing" << std::endl;
		std::cout << desc << std::endl;
		exit(EXIT_FAILURE);
	}
	
	/* MAIN ***********************************************************/
	
	try
	{
		//Load Sampler
		std::stringstream fn_module;
		fn_module << DIR_MODULE << "samplers/" << method << ".so";
		void* handle = dlopen(fn_module.str().c_str(), RTLD_LAZY);
		
		if (!handle)
		{
			std::stringstream errorMsg;
			errorMsg << "Cannot open library: " << dlerror();
			throw stk::exception::Message(errorMsg.str(), STK_DBG_INFO);
		}
		
		dlerror();
		ModuleSamplerInit module_sampler_init = (ModuleSamplerInit) dlsym(handle, "module_sampler_init");
		const char *dlsym_error = dlerror();
		{
			const char *dlsym_error = dlerror();
			if (dlsym_error)
			{
				std::stringstream errorMsg;
				errorMsg << "Cannot open library " << fn_module.str() << " : " << dlsym_error;
				dlclose(handle);
				throw stk::exception::Message(errorMsg.str(), STK_DBG_INFO);
			}
		}
		
		dlerror();
		ModuleSamplerExec module_sampler_exec = (ModuleSamplerExec) dlsym(handle, "module_sampler_exec");
		{
			const char *dlsym_error = dlerror();
			if (dlsym_error)
			{
				std::stringstream errorMsg;
				errorMsg << "Cannot open library " << fn_module.str() << " : " << dlsym_error;
				dlclose(handle);
				throw stk::exception::Message(errorMsg.str(), STK_DBG_INFO);
			}
		}
		
		dlerror();
		ModuleSamplerFree module_sampler_free = (ModuleSamplerFree) dlsym(handle, "module_sampler_free");
		{
			const char *dlsym_error = dlerror();
			if (dlsym_error)
			{
				std::stringstream errorMsg;
				errorMsg << "Cannot open library " << fn_module.str() << " : " << dlsym_error;
				dlclose(handle);
				throw stk::exception::Message(errorMsg.str(), STK_DBG_INFO);
			}
		}
		
		//Open Output
		stk::io::PointSetOutputStream<2, double, double> stream;
		stream.setValueType(stk::io::PointSetStream::VAL_NONE);
		stream.setPositionType(stk::io::PointSetStream::POS_DOUBLE);
		if(vm.count("binary")) stream.setBinary(true);
		stream.open(fn_output);
		
		module_sampler_init(nPts);
		
		double meanTime = 0;
		double meanPts = 0;
		
		std::cout << "Point set " << 1 << "/" << nPatches << "       \r";
		std::cout.flush();
			
		//~ boost::timer samplerTimer;
		
		for(int psCur=0; psCur<nPatches; psCur++)
		{
			stk::PointSet2dd pts;
			
			boost::timer::cpu_timer timer;
	
			module_sampler_exec(pts);
	
			boost::timer::nanosecond_type timeSystem = timer.elapsed().system;
			boost::timer::nanosecond_type timeUser = timer.elapsed().user;
			boost::timer::nanosecond_type timeWall = timer.elapsed().wall;
			meanTime += (static_cast<double>(timeWall)/1000000000.0 - meanTime)/(psCur+1);
	
			meanPts += (static_cast<double>(pts.size()) - meanPts)/(psCur+1);
			
			std::cout << meanTime << " sec, " << meanPts << " pts" << std::endl;
			
	
			stream.write(pts);

			//~ stk::PointSet2dd pts;
			//~ 
			//~ double t0 = samplerTimer.elapsed();
			//~ module_sampler_exec(pts);
			//~ double t1 = samplerTimer.elapsed();
			//~ 
			//~ meanTime += ((t1 - t0) - meanTime)/(psCur+1);
			//~ float elt = (double)(nPatches-psCur-1)*meanTime;
			//~ int eltS = elt;
			//~ int eltM = eltS/60; eltS %= 60;
			//~ int eltH = eltM/60; eltM %= 60;
			//~ int eltD = eltH/24; eltH %= 24;
			//~ 
			//~ std::cout
				//~ << "Point set " << (psCur+1) << "/" << nPatches << ", "
				//~ << "estimated left time : ";
			//~ 
			//~ if(eltD > 1) std::cout << eltD << " days and ";
			//~ else if(eltD == 1) std::cout << eltD << " day and ";
			//~ 
			//~ std::cout << eltH << ":" << eltM << ":" << eltS;
			//~ std::cout << "          \r";
			//~ std::cout.flush();
			//~ 
			//~ stream.write(pts);
		}
		
		module_sampler_free();
		
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
