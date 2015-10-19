#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/timer.hpp>

#include <stk/stk.hpp>

namespace boostPO = boost::program_options;

stk::PointSet2dd pts;
std::vector<int> histoX;
std::vector<double> histoY;

enum
{
	STAR_DISCREPANCY,
	RANDOM_STAR_DISCREPANCY,
	DISK_DISCREPANCY
};

bool sortXfunc(int i, int j)
{
	return (pts[i].pos()[0] < pts[j].pos()[0]);
}

bool sortYfunc(int i, int j)
{
	return (pts[i].pos()[1] < pts[j].pos()[1]);
}

double computeDiskDiscrepancy(int nPts)
{
	double discrepancy = 0.0;
	if(nPts > pts.size())
	{
		std::cerr << "The pointset is too small : " << pts.size() << std::endl;
		return discrepancy;
	}
	
	//Complexity N log(N)
	std::vector<double> sortedPts;
	sortedPts.resize(nPts);
	for(int i=0; i<nPts; ++i) sortedPts[i] = pts[i].pos().norm();
	std::sort(sortedPts.begin(), sortedPts.end());
	
	for(int i=0; i<nPts; ++i)
	{
		double radius = sortedPts[i];
		
		if(radius >= 1.0) break;
		
		double d0 = std::abs(static_cast<double>(i)/static_cast<double>(nPts) - M_PI*radius*radius/4.0);
		double d1 = std::abs(static_cast<double>(i+1)/static_cast<double>(nPts) - M_PI*radius*radius/4.0);
		
		if(d0 > discrepancy) discrepancy = d0;
		if(d1 > discrepancy) discrepancy = d1;
	}
	
	return discrepancy;
}

double computeStarDiscrepancy(int nPts)
{
	double discrepancy = 0.0;
	if(nPts > pts.size())
	{
		std::cerr << "The pointset is too small : " << pts.size() << std::endl;
		return discrepancy;
	}
	
	//Complexity N log(N)
	std::vector<int> sortedX;
	sortedX.resize(nPts);
	for(int i=0; i<nPts; ++i) sortedX[i] = i;
	std::sort(sortedX.begin(), sortedX.end(), sortXfunc);
	
	//Complexity N log(N)
	std::vector<int> sortedY;
	sortedY.resize(nPts);
	for(int i=0; i<nPts; ++i) sortedY[i] = i;
	std::sort(sortedY.begin(), sortedY.end(), sortYfunc);
	
	int* gridCurrent = new int[nPts+1];
	int* gridBlock = new int[nPts+1];
	int* gridX = new int[nPts+1];
	int* gridY = new int[nPts+1];
	
	//Complexity N
	for(int i=0; i<=nPts; i++)
	{
		gridCurrent[i] = 0;
		gridBlock[i] = 0;
		gridX[i] = 0;
	}
	
	//Complexity N*N
	{
		double boxX0;
		double boxY0;
		double boxX1;
		double boxY1;
		
		for(int j=0; j<=nPts; j++)
		{
			gridY[0] = 0;
			
			int ptsIndice = -1;
			
			//~ #pragma omp parallel for reduction(max : ptsIndice)
			//~ #pragma omp parallel for
			for(int i=0; i<=nPts; i++)
			{
				int v;
				if(i>0 && j>0 && (sortedX[i-1] == sortedY[j-1]))
				{
					v = i;
				}
				else v = -1;
				
				//~ #pragma omp flush(ptsIndice)
				if (v > ptsIndice)
				{
					//~ #pragma omp critical
					{
						if (v > ptsIndice) ptsIndice = v;
					}
				}
			}
			if(ptsIndice < 0) ptsIndice = std::numeric_limits<int>::max();
			
			//~ #pragma omp parallel for 
			for(int i=0; i<=nPts; i++)
			{
				if(i<ptsIndice) gridY[i] = 0;
				else gridY[i] = 1;
			}
			
			//~ #pragma omp parallel for 
			for(int i=0; i<=nPts; i++)
			{
				gridCurrent[i] = gridX[i];
				if(i>0 && j>0)
				{
					gridCurrent[i] += gridY[i-1] + gridBlock[i-1];
				
					if(sortedX[i-1] == sortedY[j-1])
					{
						if(ptsIndice != i) std::cerr << "Ahhhhg! " << i << " != " << ptsIndice << std::endl;
						++gridX[i];
						++gridCurrent[i];
					}
				}
			}
			
			//~ #pragma omp parallel for
			for(int i=0; i<=nPts; i++)
			{
				gridBlock[i] = gridCurrent[i];
			}
			
			if(j>0) boxY0 = pts[sortedY[j-1]].pos()[1];
			else boxY0 = 0.0;
			
			if(j<nPts) boxY1 = pts[sortedY[j]].pos()[1];
			else boxY1 = 1.0;
			
			//~ #pragma omp parallel for reduction(max : discrepancy) 
			//~ #pragma omp parallel for
			for(int i=0; i<=nPts; i++)
			{
				if(i>0) boxX0 = pts[sortedX[i-1]].pos()[0];
				else boxX0 = 0.0;
				
				if(i<nPts) boxX1 = pts[sortedX[i]].pos()[0];
				else boxX1 = 1.0;
			
				double d = discrepancy;
				double error0 = std::abs(static_cast<double>(gridCurrent[i])/nPts - boxX0*boxY0);
				if(error0 > d)
				{
					d = error0;
				}
				
				double error1 = std::abs(static_cast<double>(gridCurrent[i])/nPts - boxX1*boxY1);
				if(error1 > d)
				{
					d = error1;
				}
				
				double error2 = std::abs(static_cast<double>(gridCurrent[i])/nPts - boxX0*boxY1);
				if(error2 > d)
				{
					d = error2;
				}
				
				double error3 = std::abs(static_cast<double>(gridCurrent[i])/nPts - boxX1*boxY0);
				if(error3 > d)
				{
					d = error3;
				}
				
				//~ #pragma omp flush(discrepancy)
				if (d > discrepancy)
				{
					//~ #pragma omp critical
					{
						if (d > discrepancy) discrepancy = d;
					}
				}
			}
		}
	}
	
	delete[] gridCurrent;
	delete[] gridBlock;
	delete[] gridX;
	
	return discrepancy;
}

template<int BASE>
float generateCoordinate(int n)
{
	float res = 0.0f;
	float dev = 1.0f/BASE;
	
	while(n != 0)
	{
		int d = n % BASE;
		res += d * dev;
		n = (n-d)/BASE;
		dev = dev / BASE;
	}
	
	return res;
}

double computeRandomStarDiscrepancy(int nPts, int iterMax)
{
	double discrepancy = 0.0;
	if(nPts > pts.size())
	{
		std::cerr << "The pointset is too small : " << pts.size() << std::endl;
		return discrepancy;
	}
	
	for(int i=0; i<iterMax; ++i)
	{
		double x = drand48();
		double y = drand48();
		
		int counter = 0;
		
		for(int j=0; j<nPts; ++j)
		{
			if(pts[j].pos()[0] < x && pts[j].pos()[1] < y) ++counter;
		}
		
		double error = std::abs(static_cast<double>(counter)/static_cast<double>(nPts) - x*y);
		
		if(error > discrepancy) discrepancy = error;
	}
	
	return discrepancy;
}

int main(int argc, char** argv)
{
	/* ARG PARSER *****************************************************/
	std::string fn_input;
	std::string fn_output;
	std::string fn_histo;
	std::string method;
	int iterMax;
    srand48( time(NULL) + 100.*(double)clock() );
	
	boostPO::variables_map vm;
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h",
			"produce help message")
		("input,i",
			boostPO::value(&fn_input)->required(),
			"Filename of the input point set")
		("output,o",
			boostPO::value(&fn_output)->required(),
			"output filename")
		("iteration,I",
			boostPO::value(&iterMax)->required(),
			"number of iteration for random star discrepancy")
		("method,m",
			boostPO::value(&method)->required(),
			"discrepancy method : star, disk, random-star")
		("histo-file,H",
			boostPO::value(&fn_histo),
			"histogram to load (at least a file with one bin per line)");
	
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
	
	stk::io::read(fn_input, pts);
	
	//Histogram
	if(vm.count("histo-file"))
	{
		std::cout << "Loading bins descriptor: ";
		std::cout.flush();
	
		std::ifstream fileHisto(fn_histo.c_str());
		int trial;
		int x;
		double y;
		double vx;
		double vy;
		
		std::string line;
		while(!fileHisto.eof())
		{
			std::getline(fileHisto, line);
			std::stringstream lineStream;
			lineStream << line;
			
			int numberOfPoints;
			lineStream >> numberOfPoints;
			histoX.push_back(numberOfPoints);
			histoY.push_back(0.0);
		}
		
		fileHisto.close();
		
		std::cout << "done! (" << histoX.size() << " bins created)" << std::endl;
	}
	else
	{			
		histoX.push_back(pts.size());
		histoY.push_back(0.0);
		
		std::cout << "Loading bins descriptor: done! (one bin created)" << std::endl;
	}
	
	int methodId;
	if(method == "star") methodId = STAR_DISCREPANCY;
	else if(method == "random-star") methodId = RANDOM_STAR_DISCREPANCY;
	else if(method == "disk") methodId = DISK_DISCREPANCY;
	else std::exit(EXIT_FAILURE);
	
	for(int i=0; i<histoX.size(); ++i)
	{
		double discrepancy;
		if(methodId == STAR_DISCREPANCY) discrepancy = computeStarDiscrepancy(histoX[i]);
		else if(methodId == RANDOM_STAR_DISCREPANCY) discrepancy = computeRandomStarDiscrepancy(histoX[i], iterMax);
		else if(methodId == DISK_DISCREPANCY) discrepancy = computeDiskDiscrepancy(histoX[i]);
		
		histoY[i] = discrepancy;
		
		std::ofstream file(fn_output.c_str());
		if(!file) std::cerr << "Can't open ouptut" << std::endl;
		
		file << std::setprecision(9);
		for(int j=0; j<i+1; ++j)
		{
			file << histoX[j] << "\t" << histoY[j] << std::endl;
		}
		
		file.close();
		
		std::cout << "Discrepancy computed for N = " << histoX[i] << std::endl;
	}
	
	exit(EXIT_SUCCESS);
}
