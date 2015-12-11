#include <cstdlib>
#include <iostream>
#include <string>
#include <ctime>
#include <boost/program_options.hpp>
//~ #include <boost/progress.hpp>
#include <boost/timer/timer.hpp>

#include <stk/stk.hpp>

namespace boostPO = boost::program_options;

void sobolSampler(stk::PointSet2dd& pts, int nPts, bool shift)
{
	const int D = 2;
	
	// L = max number of bits needed 
	unsigned L = (unsigned)std::ceil(std::log((double)nPts)/log(2.0)); 

	// C[i] = index from the right of the first zero bit of i
	unsigned *C = new unsigned [nPts];
	C[0] = 1;
	for (unsigned i=1;i<=nPts-1;i++)
	{
		C[i] = 1;
		unsigned value = i;
		while (value & 1)
		{
			value >>= 1;
			C[i]++;
		}
	}

	// POINTS[i][j] = the jth component of the ith point
	//                with i indexed from 0 to N-1 and j indexed from 0 to D-1
	double *POINTS = new double[nPts*D];
	POINTS[0] = 0.0;
	POINTS[1] = 0.0;

	// ----- Compute the first dimension -----

	// Compute direction numbers V[1] to V[L], scaled by pow(2,32)
	unsigned *V = new unsigned [L+1]; 
	for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

	// Evalulate X[0] to X[N-1], scaled by pow(2,32)
	unsigned *X = new unsigned [nPts];
	X[0] = 0;
	for (unsigned i=1;i<=nPts-1;i++)
	{
		X[i] = X[i-1] ^ V[C[i-1]];
		POINTS[i*D+0] = (double)X[i]/pow(2.0,32); // *** the actual points
		//        ^ 0 for first dimension
	}

	// Clean up
	delete [] V;
	delete [] X;

	//Computing remaning dimensions (just the second one in fact)
	{
		const int j=1;
		const unsigned d=2;
		const unsigned s=1;
		const unsigned a=0;

		// Read in parameters from file
		unsigned *m = new unsigned [s+1];
		m[1] = 1;

		// Compute direction numbers V[1] to V[L], scaled by pow(2,32)
		unsigned *V = new unsigned [L+1];
		if (L <= s)
		{
			for (unsigned i=1;i<=L;i++)
			{
				V[i] = m[i] << (32-i); 
			}
		}
		else
		{
			for (unsigned i=1;i<=s;i++)
			{
				V[i] = m[i] << (32-i);
			}
			for (unsigned i=s+1;i<=L;i++)
			{
				V[i] = V[i-s] ^ (V[i-s] >> s); 
				
				for (unsigned k=1;k<=s-1;k++) 
				{
					V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
				}
			}
		}

		// Evalulate X[0] to X[N-1], scaled by pow(2,32)
		unsigned *X = new unsigned [nPts];
		X[0] = 0;
		for (unsigned i=1;i<=nPts-1;i++)
		{
			X[i] = X[i-1] ^ V[C[i-1]];
			POINTS[i*D+j] = (double)X[i]/pow(2.0,32); // *** the actual points
			//        ^ j for dimension (j+1)
		}

		// Clean up
		delete [] m;
		delete [] V;
		delete [] X;
	}
	
	delete [] C;
	
	double x = 0;
	double y = 0;
	if (shift)
	{
		 x = drand48();
		 y = drand48();
	}
	
	for(int i=0; i<nPts; i++)
	{
		double px = POINTS[i*D+0]+x;
		if (px > 1)
			px -= 1;
		double py = POINTS[i*D+1]+y;
		if (py > 1)
			py -= 1;
		
		pts.push_back(stk::Point2dd(stk::Vector2d(px, py), 1.0));
	}
	
	delete [] POINTS;
}

int main(int argc, char** argv)
{
	/* ARG PARSER *****************************************************/
	std::string fn_output;
	std::string method;
	int nPts;
	int nPatches;
	
	srand48(time(NULL));
	
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
		("shift,s",
			"random shifts the pointsets")
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
			sobolSampler(pts, nPts, vm.count("shift"));
	
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
