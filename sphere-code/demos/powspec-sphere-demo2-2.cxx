/* Read n files, with each file containing point samples: x y z
 * and generate corresponding power spectrum 1D plot.
 * *

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../io/read-pointset.h"
#include "./../core/utils.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc != 4){
        std::cerr << "Parameters: filename samplingPattern N" << std::endl;
        exit(-2);
    }

    //########################################################33
    srand48(time(NULL));
    std::string filename = argv[1];
    std::string samplingpattern = argv[2];
    int approxN = atoi(argv[3]);

    //########################################################33
    std::vector<double> vec;
    read_pointsetnD(filename, vec, 3);
    int N = vec.size()/3.0;
    std::stringstream ss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    ss.str(std::string());
    ss << resultFolder << "power-spectrum-radial-" << samplingpattern << "-n" << approxN  << "/";
    mkdir(ss.str().c_str(),0755);

#ifdef __APPLE__
    create_folders(ss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(ss.str(), datafiles, images, graphs);
#endif

    //#####################Declaring Variables############################
    int nlmax = 356;
    double* lmAccum = new double[nlmax*nlmax] ();
    double* powSpecAccum = new double[nlmax] ();
    double* t_lm = new double[nlmax*nlmax] ();
    double* t_powspec = new double[nlmax] ();

    //########################################################33

    std::cerr << "Number of points: " << N << std::endl;

    int nIters = 1000;
    for(int iter=1; iter<= nIters; iter++){
        fprintf(stderr,"\rIterations (%d) ",iter);

        //Read all files containing sample points from the directory
        filename = argv[1];
        ss.str(std::string());
        filename.erase(filename.end()-10, filename.end());
        if(iter < 10)
            ss << "00000" << iter << ".txt";
        else if( iter > 9 && iter < 100)
            ss << "0000" << iter << ".txt";
        else if( iter > 99 && iter < 1000)
            ss << "000" << iter << ".txt";
        else if( iter > 999 && iter < 9999)
            ss << "00" << iter << ".txt";

        filename += ss.str();
        //std::cerr << filename << std::endl;
        //continue;
        //exit(-2);

        std::vector<double> sphereSamples;
        read_pointsetnD(filename, sphereSamples, 3);
        N = sphereSamples.size() / 3.0;


        if(N==0){
            std::cerr << "No Data available !!!" << std::endl;
            return 0;
        }

        for(int r=0;r<nlmax;r++){
            t_powspec[r] = 0;
            for(int c=0;c<nlmax;c++){
                t_lm[r*nlmax+c] = 0;
            }
        }
        //std::cerr << "N : " << N << std::endl;
        //Compute SH coeffs for a given l and m and corresponding Power Spectra from these values
        //harmonics_powspec(lm, powSpec, sphereSamples, N, nlmax);
        powspec_spherical_harmonics_parallel(t_lm, t_powspec, sphereSamples, N, nlmax);

        for(int r=0; r< nlmax; r++)
            for(int c = 0; c < nlmax; c++)
                lmAccum[r*nlmax+c] += t_lm[r*nlmax+c];
        for(int l=0; l< nlmax; l++)
            powSpecAccum[l] += t_powspec[l];

        //########################################################33

        //Normalize powSpecAccumulator
        //if(iter%9 == 0 || iter%10 == 0 || iter == 1 || iter == 0)
        {
            std::ofstream fpowspec;
            ss.str(std::string());
            ss << datafiles << "powspec-sphere-" << samplingpattern << "-n" << approxN << "-t" << iter << ".txt";
            fpowspec.open(ss.str().c_str());

            for(int i=0;i<nlmax; i++)
                fpowspec << std::fixed << std::setprecision(8) << i <<" " << powSpecAccum[i] * 1./double(iter) << std::endl;
            fpowspec.close();
        }
    }

    //##########################################################

    delete [] lmAccum;
    delete [] t_lm;
    delete [] t_powspec;
    delete [] powSpecAccum;

    return 0;
}

/* */

