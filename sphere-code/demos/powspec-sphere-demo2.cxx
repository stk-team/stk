/** Read one file, of points : x y z
 * Generate power spectrum data file as output
 * */

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../io/read-pointset.h"
#include "./../core/utils.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

int main(int argc, char* argv[]){

    if(argc != 4){
        std::cerr << "Parameters: filename samplingPattern N" << std::endl;
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

    std::string filename = argv[1];
    std::string samplingpattern = argv[2];
    int testN = atoi(argv[3]);

    //##########################################################
    std::stringstream ss;
    std::vector<double> sphereSamples;
    read_pointsetnD(filename, sphereSamples, 3);
    int N = sphereSamples.size()/3.0;
    //##########################################################

    //##########################################################

    int nlmax = 256;
    //Compute Radial Distribution Function of spherical harmonic coefficients
    double* lm = new double[nlmax*nlmax] ();
    double* powSpec = new double[nlmax] ();
    double* lmAccum = new double[nlmax*nlmax] ();

    powspec_spherical_harmonics_parallel(lm, powSpec, sphereSamples, N, nlmax);
    std::cerr << std::endl;
    //########################################################33


    //########################################################33
    //Normalize powSpecAccumulator
    std::ofstream fpowspec;
    ss.str(std::string());
    ss << datafiles << "powspec-sphere-" << samplingpattern << "-n" << testN << ".txt";
    fpowspec.open(ss.str().c_str());

    for(int i=0;i<nlmax; i++)
        fpowspec << std::fixed << std::setprecision(8) << i <<" " << powSpec[i] << std::endl;
    fpowspec.close();
    //#########################################################

    delete [] lm;
    delete [] lmAccum;
    delete [] powSpec;

    return 0;
}





