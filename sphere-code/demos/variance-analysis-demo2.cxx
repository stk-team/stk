/* Variance analysis of Spherical Harmonic Function
 * using poissondisk sampler
 * */
#include "./../core/utils.h"
#include "./../core/utils-sphere.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sphere-analysis/spherical-functions.h"

#include <sstream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc != 4){
        std::cerr << "Parameters: samplingpattern L M" << std::endl;
        exit(-1);
    }


    //#####################Declaring Variables############################

    srand48(time(NULL));
    std::string samplingpattern = argv[1];
    int L = atoi(argv[2]);
    int M = atoi(argv[3]);
    std::stringstream ss;
    int totalPts = 105000;
    int numTrials = 200;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./results/";
    ss.str(std::string());
    ss << resultFolder << "variance-analysis-sphericalharmonic-" << samplingpattern << "-trials-" << numTrials  << "/";
    mkdir(ss.str().c_str(),0755);

#ifdef __APPLE__
    create_folders(ss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(ss.str(), datafiles, images, graphs);
#endif

    //##########################################################
    double *sampleSpace;
    int sampleSpaceSize = 0;
    if(samplingpattern == "poissondisk"){
        int initN = 1048576; ///1M samples for initialzation of the Sample Space S_i

        std::vector<double> initSamples;
        initSamples = spherical_stratified_samples(initN);
        sampleSpace = new double[initSamples.size()] ();
        sampleSpaceSize = initSamples.size()/3.0;
        for(int i=0; i<sampleSpaceSize;i++){
            sampleSpace[3*i+0] = initSamples[3*i+0];
            sampleSpace[3*i+1] = initSamples[3*i+1];
            sampleSpace[3*i+2] = initSamples[3*i+2];
        }
        std::cerr << "PoissonDisk initial discrete sampleSpace generated !" << std::endl;
    }

    //##########################################################

    std::ofstream fvar, fmean, fvarReN;
    ss.str(std::string());
    ss << datafiles << "variance-spherical-harmonic-l" << L << "-m" << M << "-" << samplingpattern << ".txt";
    fvar.open(ss.str().c_str(), std::ios::app);

    ss.str(std::string());
    ss << datafiles << "variance-spherical-harmonic-withrealN-l" << L << "-m" << M << "-" << samplingpattern << ".txt";
    fvarReN.open(ss.str().c_str(), std::ios::app);

    ss.str(std::string());
    ss << datafiles << "mean-spherical-harmonic-l" << L << "-m" << M << "-" << samplingpattern << ".txt";
    fmean.open(ss.str().c_str(), std::ios::app);

    //##############################################################################

    //##############################################################################

    for(int N = 1; N <= totalPts;){
        //fprintf(stderr,"\r Numpts (%d) ",N);

        double mean = 0;
        double variance = 0;
        int realN = N;
        for(int trial = 1; trial <= numTrials; trial++){
            fprintf(stderr,"\r trial/Numpts (%d/%d) ",trial, N);

            //##############################################################################
            std::vector<double> sphereSamples;
            if(samplingpattern == "dartthrowing")
                sphereSamples = spherical_dart_throwing_samples(N);
            else if(samplingpattern == "poissondisk"){
                sphereSamples = spherical_poisson_disk_samples(N, sampleSpace, sampleSpaceSize);
            }
            else{
                std::cerr << "sampling pattern do not match !!!" << std::endl;
                std::cerr << "Patterns: dartthrowing, poissondisk" << std::endl;
                return 0;
            }
            int tempN = sphereSamples.size() / 3.0;
            //##############################################################################

            double integral = 0;
            for(int i=0; i < tempN; i++){
                double xyz[] = {sphereSamples[3*i+0],sphereSamples[3*i+1],sphereSamples[3*i+2]};
                std::complex<double> sphcoeff = spherical_harmonic_function(xyz, L,M, "complex");
                integral += sphcoeff.real();
            }
            double factor = (4*PI) / double(tempN);
            integral *= factor;

            //##############################################################################

            //Incremental Mean and Variance computation
            mean = ((trial-1)/double(trial)) * mean + (1/double(trial)) * integral;
            if(trial < 2){
                variance = 0;
            }
            else{
                variance = ((trial-1)/double(trial)) * variance + (1/double(trial-1)) * (integral - mean) * (integral - mean);
            }
            realN = sphereSamples.size()/3.0;
            //##############################################################################
        }
        //varianceRef /= double(numTrials);
        fmean << std::fixed << N << std::setprecision(15) << " " << mean << std::endl;
        fvar << std::fixed << N << std::setprecision(15) << " " << variance << std::endl;
        fvarReN << std::fixed << realN << std::setprecision(15) << " " << variance << std::endl;

        //##############################################################################
        if(N > 10000)
            N += 10000;
        else if(N >= 1000 && N <= 10000)
            N += 5000;
        else if(N >= 500 && N < 1000){
            N+=250;
        }
        else if(N >= 150 && N < 500){
            N+=50;
        }
        else if(N < 150){
            N += 10;
        }
        //##############################################################################
    }
    fvar.close();
    fmean.close();
    fvarReN.close();

    delete [] sampleSpace;
    return 0;
}



