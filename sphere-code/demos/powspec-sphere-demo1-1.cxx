/* Compute power spectrum from built-in sampling methods.
 * Created by Gurprit Singh
 * */

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

    //##########################################################

    std::string samplingpattern = argv[1];
    int N = atoi(argv[2]);

    std::stringstream ss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    ss.str(std::string());
    ss << resultFolder << "powspec-sphere-" << samplingpattern  << "-n" << N << "/";
    mkdir(ss.str().c_str(),0755);

#ifdef __APPLE__
    create_folders(ss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(ss.str(), datafiles, images, graphs);
#endif
    //#####################Declaring Variables############################

    srand48(time(NULL));

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

    int nlmax = 256;
    double* t_lm = new double[nlmax*nlmax] ();
    double* t_powSpec = new double[nlmax] ();
    double* lmAccum = new double[nlmax*nlmax] ();
    double* powSpecAccum = new double[nlmax] ();

    int nTrials=1000;
    //Compute Radial Distribution Function of the above histoAccum
    for(int trial=1; trial<= nTrials; trial++){
        fprintf(stderr," Trials (%7d%) ",trial);


        //##############################################################################
        std::vector<double> sphereSamples;

        if(samplingpattern == "whitenoise")
            sphereSamples = spherical_whitenoise_samples(N);
        else if(samplingpattern == "stratified")
            sphereSamples = spherical_stratified_samples(N);
        else if(samplingpattern == "dartthrowing")
            sphereSamples = spherical_dart_throwing_samples(N);
        else if(samplingpattern == "poissondisk"){
            sphereSamples = spherical_poisson_disk_samples(N, sampleSpace, sampleSpaceSize);
        }
        else{
            std::cerr << "sampling pattern do not match !!!" << std::endl;
            std::cerr << "Patterns: dartthrowing, poissondisk" << std::endl;
            return 0;
        }
        int newN = sphereSamples.size()/3.0;
        //##############################################################################

        //Initialize powspec and lm arrays to zero for the new iteration
        for(int r=0;r<nlmax;r++){
            t_powSpec[r] = 0;
            for(int c=0;c<nlmax;c++){
                t_lm[r*nlmax+c] = 0;
            }
        }

        powspec_spherical_harmonics_parallel(t_lm, t_powSpec, sphereSamples, newN, nlmax);

        for(int r=0; r< nlmax; r++)
            for(int c = 0; c < nlmax; c++)
                lmAccum[r*nlmax+c] += t_lm[r*nlmax+c];
        for(int l=0; l< nlmax; l++)
            powSpecAccum[l] += t_powSpec[l];

        //##########################################################

        //Normalize powSpecAccumulator

        if(trial%10 == 0 || trial%9 == 0 || trial==1)
        {
            std::ofstream fpowspec;
            ss.str(std::string());
            ss << trial;
            std::string s1 = ss.str();
            paddedzerosN(s1, nTrials);
            ss.str(std::string());
            ss << datafiles << "powspec-sphere-" << samplingpattern << "-l" << nlmax << "-n" << N << "-trial-" << s1 << ".txt";
            fpowspec.open(ss.str().c_str());

            for(int i=0;i<nlmax; i++)
                fpowspec << std::fixed << std::setprecision(8) << i <<" " << (powSpecAccum[i] * 1./double(trial)) << std::endl;
            fpowspec.close();
        }
    }
    std::cerr << std::endl;
    //##########################################################

    //##########################################################

    delete [] lmAccum;
    delete [] t_lm;
    delete [] t_powSpec;
    delete [] powSpecAccum;
    delete [] sampleSpace;
}

