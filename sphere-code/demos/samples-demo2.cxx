/*** Created by Gurprit Singh ***/

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../core/utils.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc != 3){
        std::cerr << "Parameters: samplingPattern N" << std::endl;
        exit(-2);
    }

    //#####################Declaring Variables############################

    srand48(time(NULL));
    std::string samplingpattern = argv[1];
    int N = atoi(argv[2]);
    std::stringstream ss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    ss.str(std::string());
    ss << resultFolder << "pointsets-sphere-" << samplingpattern << "/";
    mkdir(ss.str().c_str(),0755);

#ifdef __APPLE__
    create_folders(ss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(ss.str(), datafiles, images, graphs);
#endif


    //##########################################################
    std::vector<double> sphereSamples;//(N, inival);

    if(samplingpattern == "whitenoise")
        sphereSamples = spherical_whitenoise_samples(N);
    else if(samplingpattern == "regular")
        sphereSamples = spherical_regular_samples(N);//,-1,true);
    else if(samplingpattern == "stratified")
        sphereSamples = spherical_stratified_samples(N);
    else if(samplingpattern == "dartthrowing")
        sphereSamples = spherical_dart_throwing_samples(N);
    else if(samplingpattern == "poissondisk"){
        int initN = 1048576; ///1M samples for initialzation of the Sample Space S_i

        std::vector<double> initSamples;
        initSamples = spherical_stratified_samples(initN);
        double* sampleSpace = new double[initSamples.size()] ();

        int sampleSpaceSize = initSamples.size()/3.0;

        for(int i=0; i<sampleSpaceSize;i++){
            sampleSpace[3*i+0] = initSamples[3*i+0];
            sampleSpace[3*i+1] = initSamples[3*i+1];
            sampleSpace[3*i+2] = initSamples[3*i+2];
        }
        sphereSamples = spherical_poisson_disk_samples(N, sampleSpace, sampleSpaceSize);
        delete [] sampleSpace;
    }
    else{
        std::cerr << "sampling pattern do not match !!!" << std::endl;
        std::cerr << "Patterns: whitenoise, stratified, poissondisk" << std::endl;
        return 0;
    }

    std::cerr << "Numpts: " << sphereSamples.size()/3.0 << std::endl;

    //Dump into a file
    //##########################################################
    //Normalize powSpecAccumulator
    std::ofstream outfile;
    ss.str(std::string());
    ss << datafiles << "samples-"  << samplingpattern << "-n" << N << ".txt";
    outfile.open(ss.str().c_str());

    for(int i=0;i<sphereSamples.size(); i+=3)
        outfile << sphereSamples[i+0] << " " << sphereSamples[i+1] << " " << sphereSamples[i+2] << std::endl;

    outfile.close();

    //########################################################33

    return 0;
}

