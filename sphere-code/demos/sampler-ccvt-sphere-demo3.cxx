#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <queue>
#include <exception>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <cstdlib>      // std::rand, std::srand

#include <iomanip>
#include <sys/stat.h>

#include "./../sampler-sphere/ccvt-sphere.h"
#include "./../sampler-sphere/stk/pointset.hpp"
#include "./../domaintransform/domaintransform.h"
#include "./../core/utils-sphere.h"
#include "./../core/utils.h"
#include "./../io/read-pointset.h"

int main(int argc, char** argv)
{
    if(argc != 5){
        std::cerr << "Parameters: nsites nptsPerSite samplingpattern filename" << std::endl;
        return 0;
    }

    //#####################Declaring Variables############################

    srand48(time(NULL));
    int nsites = atoi(argv[1]);
    int npts = atoi(argv[2]);
    std::string samplingpattern = argv[3];
    std::string filename = argv[4];

    //##########################################################
    std::stringstream oss;
    std::vector<double> sphereSamples;
    read_pointsetnD(filename, sphereSamples, 3);
    int N = sphereSamples.size()/3.0;
    std::srand ( unsigned ( std::time(NULL) ) );
    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    oss.str(std::string());
    oss << resultFolder << "pointset-sphere-ccvt-with-" << samplingpattern << "/";
    mkdir(oss.str().c_str(),0755);

#ifdef __APPLE__
    create_folders(oss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(oss.str(), datafiles, images, graphs);
#endif

    /* INIT ***********************************************************/

    stk::PointSet2dd sites;
    stk::PointSet2di points;
    std::vector< std::vector<int> > ptsPerSites;

    //generate_polar_sample_space(sites, points, ptsPerSites, nsites, npts, samplingpattern);
    generate_sample_from_file(sites, points, ptsPerSites, nsites, npts, sphereSamples);


    for(int i=0; i<sites.size(); i++){
        for(int j=0; j < ptsPerSites[i].size(); j++){
            double pxyz[] = {0,0,0};
            thetaphi2xyz(pxyz, ptsPerSites[i][j], ptsPerSites[i][j]);
            std::cout << pxyz[0] << " " << pxyz[1] << " " << pxyz[2] << " ";
        }
        std::cout << std::endl;
    }

    //exit(-2);

    ccvt_sphere(sites, points, ptsPerSites);
    int trial = 1;
    oss.str(std::string());
    oss << trial;
    std::string s1 = oss.str();
    paddedzerosN(s1, trial);

    std::ofstream file, fvictor;
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

    if(1){
        oss.str(std::string());
        oss << datafiles << "victor-ccvt-sphere-" << samplingpattern << "-p" << npts << "-n" << nsites << "-" << s1 << ".txt";
        //std::cerr << oss.str() << std::endl;

        fvictor.open(oss.str().c_str());
        for(int i=0; i<sites.size(); i++){
            for(int j=0; j < ptsPerSites[i].size(); j++){
                double pxyz[] = {0,0,0};
                thetaphi2xyz(pxyz, ptsPerSites[i][j], ptsPerSites[i][j]);
                fvictor << pxyz[0] << " " << pxyz[1] << " " << pxyz[2] << " ";
            }
            if(i < sites.size()-1)
                fvictor << std::endl;
        }
        fvictor.close();
    }

    exit(EXIT_SUCCESS);
}

