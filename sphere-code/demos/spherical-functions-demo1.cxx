/*** Created by Gurprit Singh ***/

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../sphere-analysis/spherical-functions.h"
#include "./../core/utils.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

int main(int argc, char* argv[]){

    if(argc != 4){
        std::cerr << "Parameters: samplingPattern sphericalFunction N" << std::endl;
        exit(-2);
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

    std::string samplingpattern = argv[1];
    std::string sphericalFunction = argv[2];
    int N = atoi(argv[3]);

    //##########################################################

    std::stringstream ss;

    std::vector<double> sphereSamples;//(N, inival);

    if(samplingpattern == "whitenoise")
        sphereSamples = spherical_whitenoise_samples(N);
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

    //for(int i = 0; i < sphereSamples.size(); i+=3)
    //    std::cout << sphereSamples[i+0] << " "<< sphereSamples[i+1] <<" " << sphereSamples[i+2] << std::endl;
    //exit(-2);

    //##########################################################
    double ref[] = {0.0,1.0,0.0};
    std::vector<double> outSamples, colorValues;
    for(int i=0; i < sphereSamples.size(); i+=3){
        double xyz[] = {sphereSamples[i+0],sphereSamples[i+1],sphereSamples[i+2]};
        if(sphericalFunction == "spherical-cap"){
            //std::cout << "here " << std::endl;
            if( spherical_cap_function(ref, xyz, 30.0) == 1){
                outSamples.push_back(sphereSamples[i+0]);
                outSamples.push_back(sphereSamples[i+1]);
                outSamples.push_back(sphereSamples[i+2]);
            }
        }
        else if(sphericalFunction == "spherical-cap-notused"){
            if( spherical_disc_northpole(xyz, 30.0) == 1){
                outSamples.push_back(sphereSamples[i+0]);
                outSamples.push_back(sphereSamples[i+1]);
                outSamples.push_back(sphereSamples[i+2]);
            }
        }
        else if(sphericalFunction == "spherical-gaussian"){
            double ref[] = {0,0,1};
            double value = spherical_gaussian_function(ref, xyz, 2.0);
            colorValues.push_back(value);
            outSamples.push_back(sphereSamples[i+0]);
            outSamples.push_back(sphereSamples[i+1]);
            outSamples.push_back(sphereSamples[i+2]);
        }
        else if(sphericalFunction == "spherical-harmonic"){
            double value = spherical_harmonic_function(xyz, 4,2, "real");
            //double x = value *
        }
    }

    std::cerr << std::endl;
    //##########################################################
    //Normalize powSpecAccumulator
    std::ofstream outfile, outval;
    ss.str(std::string());
    ss << datafiles << "fsample-" << sphericalFunction << "-" << samplingpattern << "-n" << N << ".txt";
    outfile.open(ss.str().c_str());

    ss.str(std::string());
    ss << datafiles << "fvalue-" << sphericalFunction << "-" << samplingpattern << "-n" << N << ".txt";
    outval.open(ss.str().c_str());

    for(int i=0;i<outSamples.size(); i+=3)
        outfile << outSamples[i+0] << " " << outSamples[i+1] << " " << outSamples[i+2] << std::endl;

    for(int i=0;i<colorValues.size(); i++){
        outval << colorValues[i] << " "<< colorValues[i] <<" " << colorValues[i] << std::endl;
    }
    std::cerr << outSamples.size() <<" " << colorValues.size() << std::endl;
    outfile.close();
    outval.close();
    //########################################################33

    return 0;
}



