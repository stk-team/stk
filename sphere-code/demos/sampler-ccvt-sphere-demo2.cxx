#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <queue>
#include <exception>
#include <boost/program_options.hpp>
#include <boost/random.hpp>

#include <iomanip>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <sys/stat.h>

#include "./../sampler-sphere/ccvt-sphere.h"
#include "./../sampler-sphere/stk/pointset.hpp"
#include "./../domaintransform/domaintransform.h"
#include "./../core/utils-sphere.h"
#include "./../core/utils.h"


int main(int argc, char** argv)
{
    if(argc != 5){
        std::cerr << "Parameters: nsites nptsPerSite samplingpattern numTrials[1 or more]" << std::endl;
        return 0;
    }
    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
#ifdef __APPLE__
    create_folders(resultFolder, datafiles, images, graphs);
#elif __linux__
    create_folders(resultFolder, datafiles, images, graphs);
#endif
    //#####################Declaring Variables############################
    srand48(time(NULL));
    int nsites = atoi(argv[1]);
    int npts = atoi(argv[2]);
    std::string samplingpattern = argv[3];
    int numTrials = atoi(argv[4]);
    std::stringstream oss;
    /* INIT ***********************************************************/

    for(int trial=1; trial <= numTrials; trial++){
        //fprintf(stderr,"\rTrial: %d", trial);
        stk::PointSet2dd sites;
        stk::PointSet2di points;
        std::vector< std::vector<int> > ptsPerSites;

        //std::cerr << trial << std::endl;
        clock_t t1,t2;
        t1 = clock();
        coherent_initialization(sites, points, ptsPerSites, nsites, npts, samplingpattern);
        t1 = clock() - t1;

        std::cerr <<  "time for intialization: " << ((float)t1)/CLOCKS_PER_SEC << std::endl;

        std::cerr << sites.size() << " "<< points.size() << std::endl;

        for(int i=0; i < sites.size(); i++)
            std::cout << i << " " <<  ptsPerSites[i].size() << std::endl;

        exit(-2);
        t2 = clock();
        ccvt_sphere(sites, points, ptsPerSites);
        t2 = clock() - t2;
        std::cerr <<  "time for ccvt-sphere: " << ((float)t2)/CLOCKS_PER_SEC << std::endl;

        oss.str(std::string());
        oss << trial;
        std::string s1 = oss.str();
        paddedzerosN(s1, trial);

        std::ofstream file;
        oss.str(std::string());
        oss << datafiles << "ccvt-sphere-" << samplingpattern << "-p" << npts << "-n" << nsites << "-" << s1 << ".txt";
        //std::cerr << oss.str() << std::endl;

        file.open(oss.str().c_str());

        for(int i=0; i<sites.size(); i++){
            double sxyz[] = {0,0,0};
            thetaphi2xyz(sxyz, sites[i].pos()[0], sites[i].pos()[1]);
            file << sxyz[0] <<" " << sxyz[1] <<" " << sxyz[2] << std::endl;
        }
        file.close();
    }
    std::cerr << std::endl;
    exit(EXIT_SUCCESS);
}

