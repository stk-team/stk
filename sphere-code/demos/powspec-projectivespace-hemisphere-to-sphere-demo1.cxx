/* Power spectrum in Projective 2-space, of a hemispherical sampling pattern
 * computed via spherical harmonics after taking the
 * hemisphere to the Projective 2-space.
 * */

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../core/utils.h"

#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <time.h>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc != 4){
        std::cerr << "Parameters: samplingPattern N nlmax" << std::endl;
        exit(-2);
    }

    //#####################Declaring Variables############################

    srand48(time(NULL));

    std::string samplingpattern = argv[1];
    int N = atoi(argv[2]);
    int nlmax = atoi(argv[3]);

    int bandlength = 2*nlmax + 1;
    std::stringstream ss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./ results/";
    ss.str(std::string());
    ss << resultFolder << "powspec-projectivespace-hemisphere-to-sphere-" << samplingpattern  << "-n" << N << "/";
    mkdir(ss.str().c_str(),0755);

    create_folders(ss.str(), datafiles, images, graphs);

    //##########################################################

    double* lm = new double[bandlength*bandlength] ();
    double* lmAccum = new double[bandlength*bandlength] ();
    double* shMagnitude = new double[bandlength*bandlength] ();

    int powSpeclength = 0.5*nlmax;  //since all odd l valued frequencies are null!
    std::vector<double> powSpec(powSpeclength, 0.0);

    //##########################################################

    std::vector<double> sphereSamples;

    //Generating white noise on a sphere
    if(samplingpattern == "whitenoise")
        sphereSamples = hemisphere_to_sphere_mapping_whitenoise_samples(N);
    else if(samplingpattern == "stratified")
        sphereSamples = hemisphere_to_sphere_mapping_stratified_samples(N);
    else if(samplingpattern == "dartthrowing")
        sphereSamples = hemisphere_to_sphere_mapping_dartthrowing_samples(N);
    else
        exit(-2);

    int tempN = sphereSamples.size() / 3.0;
    std::cerr << "Numpts: " << tempN << std::endl;

    //##########################################################

    std::vector<double> rotatedSamples(sphereSamples.size(), 0.0);
    random_rotate_spherical_samples(rotatedSamples, sphereSamples, tempN);

    //##########################################################

    //powspec_spherical_harmonics_parallel(lm, powSpec, rotatedSamples, tempN, nlmax);
    powspec_projectivespace_sharmonics_parallel(lm, powSpec, rotatedSamples, tempN, nlmax);

    std::cerr << std::endl;
    //##########################################################
    //Normalize powSpecAccumulator
    std::ofstream fpowspec;
    ss.str(std::string());
    ss << datafiles << "powspec-projectivespace-hemisphere-to-sphere-" << samplingpattern << "-n" << N << ".txt";
    fpowspec.open(ss.str().c_str());

    for(int i=0;i<powSpec.size(); i++){
        if(i == 0)
            fpowspec << std::fixed << std::setprecision(8) << i <<" " <<  powSpec[i] << std::endl;
        else
            fpowspec << std::fixed << std::setprecision(8) << i <<" " <<  0.5*powSpec[i] << std::endl;
    }
    fpowspec.close();

    //########################################################33

    delete [] lm;
    delete [] lmAccum;
    return 0;
}
