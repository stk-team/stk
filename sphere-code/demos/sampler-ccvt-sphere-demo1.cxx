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
#include <sys/stat.h>

#include "./../sampler-sphere/ccvt-sphere.h"
#include "./../sampler-sphere/stk/pointset.hpp"
#include "./../domaintransform/domaintransform.h"
#include "./../core/utils-sphere.h"
#include "./../core/utils.h"


int main(int argc, char** argv)
{
    if(argc != 6){
        std::cerr << "Parameters: nsites nptsPerSite InitialSamplingpattern trialBegin trialEnd" << std::endl;
        return 0;
    }

    //#####################Declaring Variables############################
    srand48(time(NULL));
    int nsites = atoi(argv[1]);
    int npts = atoi(argv[2]);
    std::string samplingpattern = argv[3];
    int trialBegin = atoi(argv[4]);
    int trialEnd = atoi(argv[5]);

    std::stringstream oss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    oss.str(std::string());
    oss << resultFolder << "ccvt-n" << nsites << "-" << samplingpattern <<  "-trial-" << trialBegin << "-" << trialEnd << "/";
    mkdir(oss.str().c_str(),0755);

#ifdef __APPLE__
    create_folders(oss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(oss.str(), datafiles, images, graphs);
#endif

    /* INIT ***********************************************************/

    for(int trial=trialBegin; trial <= trialEnd; trial++){
        fprintf(stderr,"\rTrial: %d", trial);
        stk::PointSet2dd sites;
        stk::PointSet2di points;
        std::vector< std::vector<int> > ptsPerSites;

        //std::cerr << trial << std::endl;
        //clock_t t1,t2;
        //t1 = clock();
        generate_polar_sample_space(sites, points, ptsPerSites, nsites, npts, samplingpattern);
        //exit(-2);
        //t1 = clock() - t1;
        //std::cerr <<  "time for intialization: " << ((float)t1)/CLOCKS_PER_SEC << std::endl;

        ccvt_sphere(sites, points, ptsPerSites);
        //t2 = clock() - t2;
        //std::cerr <<  "time for ccvt-sphere: " << ((float)t2)/CLOCKS_PER_SEC << std::endl;
        //exit(-2);
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
