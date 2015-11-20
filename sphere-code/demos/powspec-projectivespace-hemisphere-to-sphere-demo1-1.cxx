/* Power spectrum in Projective 2-space, of a hemispherical sampling pattern
 * computed via spherical harmonics after taking the
 * hemisphere to the Projective 2-space.
 * Averaged over multiple realizations!
 * */

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../core/utils.h"
#include "./../core/utils-samples.h"

#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <time.h>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc != 5){
        std::cerr << "Parameters: samplingPattern N nlmax numTrials" << std::endl;
        exit(-2);
    }

    //#####################Declaring Variables############################

    srand48(time(NULL));

    std::string samplingpattern = argv[1];
    int N = atoi(argv[2]);
    int nlmax = atoi(argv[3]);
    int numTrials = atoi(argv[4]);

    int bandlength = 2*nlmax + 1;
    std::stringstream ss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    ss.str(std::string());
    ss << resultFolder << "powspec-projectivespace-hemisphere-to-sphere-" << samplingpattern  << "-n" << N << "/";
    mkdir(ss.str().c_str(),0755);

    create_folders(ss.str(), datafiles, images, graphs);

    //##########################################################

    double* lm = new double[bandlength*bandlength] ();
    double* lmAccum = new double[bandlength*bandlength] ();

    int powSpeclength = 0.5*nlmax;  //since all odd l valued frequencies are null!
    std::vector<double> powSpec(powSpeclength, 0.0);
    std::vector<double> powSpecAccum(powSpeclength, 0.0);

    //##########################################################

    for(int trial=1; trial <= numTrials; trial++){
        fprintf(stderr,"\r trial : (%5d) ", trial);
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

        //##########################################################

        std::vector<double> rotatedSamples(sphereSamples.size(), 0.0);
        random_rotate_spherical_samples(rotatedSamples, sphereSamples, tempN);

        //##########################################################

        for(int r=0; r< bandlength; r++)
            for(int c = 0; c < bandlength; c++){
                lm[r*bandlength+c] = 0.0;
            }

        for(int i=0; i < powSpeclength; i++)
            powSpec[i] = 0.0;

        //##########################################################

        powspec_projectivespace_sharmonics_parallel(lm, powSpec, rotatedSamples, tempN, nlmax);

        //##########################################################

        for(int i=0; i < powSpeclength; i++)
            powSpecAccum[i] += powSpec[i];

        if(trial%10 == 0 || trial == 1){
            std::ofstream fpowspec;
            ss.str(std::string());
            ss << trial;
            std::string s1 = ss.str();
            paddedzerosN(s1, numTrials);
            ss.str(std::string());
            ss << datafiles << "powspec-projectivespace-hemisphere-to-sphere-" << samplingpattern << "-l" << nlmax << "-n" << N << "-trial-" << s1 << ".txt";
            fpowspec.open(ss.str().c_str());

            for(int i=0;i<powSpec.size(); i++){
                if(i == 0)
                    fpowspec << std::fixed << std::setprecision(8) << i <<" " <<  powSpecAccum[i]/double(trial) << std::endl;
                else
                    fpowspec << std::fixed << std::setprecision(8) << i <<" " <<  0.5*(powSpecAccum[i]/double(trial)) << std::endl;
            }
        }
    }
    std::cerr << std::endl;
    //########################################################33

    delete [] lm;
    delete [] lmAccum;
    return 0;
}
