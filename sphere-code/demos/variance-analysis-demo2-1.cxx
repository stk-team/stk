/* Variance analysis of Spherical Cap: CCVT Sampler
 * Reads all files for a given N, and compute variance for that N.
 * */

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "./../core/utils.h"
#include "./../io/read-pointset.h"
#include "./../core/utils-sphere.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sphere-analysis/spherical-functions.h"


int main(int argc, char* argv[]){

    if(argc != 6){
        std::cerr << "Parameters: filename samplingPattern N L M" << std::endl;
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
    //########################################################33
    srand48(time(NULL));

    std::string filename = argv[1];
    std::string samplingpattern = argv[2];
    int approxN = atoi(argv[3]);
    int L = atoi(argv[4]);
    int M = atoi(argv[5]);

    //########################################################33
    std::vector<double> vec;
    read_pointsetnD(filename, vec, 3);
    int N = vec.size()/3.0;
    std::stringstream ss;

    //#####################Declaring Variables############################
    double ref[] = {0.,0.,1.};
    float theta = 30.0;
    std::ofstream fvar, fmean;//,fvarRef;
    ss.str(std::string());
    ss << datafiles << "variance-spherical-harmonic-l" << L << "-m" << M << "-" << samplingpattern << ".txt";
    fvar.open(ss.str().c_str(), std::ios::app);

    ss.str(std::string());
    ss << datafiles << "mean-spherical-harmonic-l" << L << "-m" << M << "-" << samplingpattern << ".txt";
    fmean.open(ss.str().c_str(), std::ios::app);

    //########################################################33

    //std::cerr << "Number of points: " << N << std::endl;

    int numTrials = 500;
    double mean = 0, variance = 0;
    for(int trial=1; trial<= numTrials; trial++){
        //fprintf(stderr,"\rTrials (%d/%d) ",trial,numTrials);

        //Read all files containing sample points from the directory
        filename = argv[1];
        ss.str(std::string());
        filename.erase(filename.end()-10, filename.end());
        if(trial < 10)
            ss << "00000" << trial << ".txt";
        else if( trial > 9 && trial < 100)
            ss << "0000" << trial << ".txt";
        else if( trial > 99 && trial < 1000)
            ss << "000" << trial << ".txt";
        else if( trial > 999 && trial < 9999)
            ss << "00" << trial << ".txt";

        filename += ss.str();
        //std::cerr << filename << std::endl;
        //continue;

        std::vector<double> sphereSamples;
        read_pointsetnD(filename, sphereSamples, 3);
        N = sphereSamples.size() / 3.0;

        if(N==0){
            std::cerr << "No Data available !!!" << std::endl;
            return 0;
        }
        float thetaInDeg = theta;
        //##############################################################################
        //std::vector<double> rotatedSamples;
        //rotatedSamples = randomshift_on_sphere(sphereSamples, sphereSamples.size());
        //##############################################################################

        double integral = 0;
        for(int i=0; i < sphereSamples.size(); i+=3){
            double xyz[] = {sphereSamples[i+0],sphereSamples[i+1],sphereSamples[i+2]};
            std::complex<double> sphcoeff = spherical_harmonic_function(xyz, L,M, "complex");
            integral += sphcoeff.real();
        }
        double factor = (4*PI) / double(N);
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
        //##############################################################################
    }
    //std::cerr << std::endl;
    fmean << std::fixed << approxN << std::setprecision(15) << " " << mean << std::endl;
    fvar << std::fixed << approxN << std::setprecision(15) << " " << variance << std::endl;
    std::cerr << std::fixed << approxN << " " << std::setprecision(15) << " " << variance << std::endl;
    //##############################################################################
    fvar.close();
    fmean.close();
    return 0;
}



