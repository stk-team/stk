/* Variance Analysis demo 1 spherical cap!
 * Here our function the spherical cap,
 * sampled with given sampler.
 * Integration done using Monte Carlo estimation.
 * */


//#include "./../../sphere-harmonics/sphere-analysis.h"
//#include "./../../sphere-harmonics/spherical-functions.h"
//#include "./../../sampler-sphere/sampler-sphere.h"
//#include "./../../domaintransform/domaintransform.h"
//#include "./../../core/utils-sphere.h"
//#include "./../../core/utils.h"

#include "./../core/utils.h"
#include "./../core/utils-sphere.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sphere-analysis/spherical-functions.h"
#include "./../domaintransform/domaintransform.h"

#include <fstream>
#include <iomanip>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc != 4){
        std::cerr << "Parameters: samplingpattern theta numTrials" << std::endl;
        exit(-1);
    }

    //#####################Declaring Variables############################

    std::string samplingpattern = argv[1];
    double theta = atof(argv[2]);
    int numTrials = atoi(argv[3]);

    std::stringstream ss;
    int totalPts = 200000;

    srand48(time(NULL));

   //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./results/";
    ss.str(std::string());
    ss << resultFolder << "variance-analysis-sphericalcap-" << samplingpattern << "-trials" << numTrials <<"/";
    mkdir(ss.str().c_str(), 0755);

#ifdef __APPLE__
    create_folders(ss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(ss.str(), datafiles, images, graphs);
#endif

    //##########################################################
    //Global Sample Space for PoissonDisk Sampler
    //std::vector<double> sampleSpace;
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


//    if(samplingpattern == "poissondisk"){
//        int initN = 1048576; ///1M samples for initialzation of the Sample Space S_i
//        sampleSpace = spherical_stratified_samples(initN);
//        sampleSpaceSize = sampleSpace.size() / 3.0;
//    }
    //##########################################################

    std::ofstream fvar, fmean;
    ss.str(std::string());
    ss << datafiles << "variance-spherical-cap-theta-" << theta << "-" << samplingpattern << ".txt";
    fvar.open(ss.str().c_str(), std::ios::app);
    //file << "#Numpts \t variance \t mean " << std::endl;

    ss.str(std::string());
    ss << datafiles << "mean-spherical-cap-theta-" << theta << "-" << samplingpattern << ".txt";
    fmean.open(ss.str().c_str(), std::ios::app);
    //fmean << "#Numpts \t variance \t mean " << std::endl;


    //##########################################################

    double rx=0,ry=0,rz=0;
    uniformSampleSphere(drand48(), drand48(), &rx, &ry, &rz);
    //center of the spherical disc computed randomly on the sphere
    //vec3 ref(rx, ry, rz);
    double ref[] = {0,0,1};

    //##########################################################

    for(int N = 1; N <= totalPts;){
        double mean = 0;
        double variance = 0;
        for(int trial = 1; trial <= numTrials; trial++){
            fprintf(stderr,"\r Numpts/Trial (%7d/%7d)",N, trial);
            //arr<vec3> sphereSamples;
            std::vector<double> sphereSamples;

            //##########################################################

            if(samplingpattern == "whitenoise")
                sphereSamples = spherical_whitenoise_samples(N);
            else if(samplingpattern == "dartthrowing")
                sphereSamples = spherical_dart_throwing_samples(N);
            else if(samplingpattern == "poissondisk"){
                sphereSamples = spherical_poisson_disk_samples(N, sampleSpace, sampleSpaceSize);
            }
            else{
                std::cerr << "sampling pattern do not match !!!" << std::endl;
                std::cerr << "Patterns: whitenoise, stratified, poissondisk" << std::endl;
                return 0;
            }

            //##########################################################
            int tempN = sphereSamples.size() / 3.0;
            std::vector<double> rotatedSamples(3*tempN, 0.0);
            random_rotate_spherical_samples(rotatedSamples, sphereSamples, tempN);

            //##########################################################

            double integral = 0;
            for(int i=0; i < tempN; i++){
                //pointing sample;
                double point [] = {rotatedSamples[3*i+0],rotatedSamples[3*i+1],rotatedSamples[3*i+2]};
                double value = spherical_cap_function(ref, point, theta);
                integral += value;
            }
            double factor = ((4*PI) / double(tempN));
            integral *= factor;

            //##########################################################

            //Incremental Mean and Variance computation
            mean = ((trial-1)/double(trial)) * mean + (1/double(trial)) * integral;
            if(trial < 2)
                variance = 0;
            else
                variance = ((trial-1)/double(trial)) * variance + (1/double(trial-1)) * (integral - mean) * (integral - mean);

            //##########################################################
        }
        fvar << std::fixed << std::setprecision(15) << N << " " << variance << std::endl;
        fmean << std::fixed << std::setprecision(15) << N << " " << mean << std::endl;

        //##########################################################

        if(N >= 20000)
            N += 10000;
        else if(N >= 1000 && N < 20000){
            N+=5000;
        }
        else if(N >= 500 && N < 1000)
            N += 250;
        else if(N >= 150 && N < 500){
            N+=10;
        }
        else if( N < 150){
            N += 5;
        }
        //##########################################################
    }

    fvar.close();
    fmean.close();

    delete [] sampleSpace;

    return 0;
}

