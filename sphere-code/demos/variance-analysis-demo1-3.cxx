/* Variance Analysis Spherical cap function Stratified and Regular Sampler!
 * Both Samplers are deduced from Victor Ostromukhov's Healpix Sampler adaptation.
 * Spherical cap function sampled with stratified sampling.
 * Integration done using Monte Carlo approach.
 * */

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sphere-analysis/spherical-functions.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../domaintransform/domaintransform.h"
#include "./../core/utils.h"
#include "./../core/utils-sphere.h"
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc !=2){
        std::cerr << "Paramters : samplingpattern" << std::endl;
        std::cerr << "stratified \n" <<
                     "regular"   << std::endl;
        return 0;
    }

    //#####################Declaring Variables############################

    srand48(time(NULL));
    std::string samplingpattern = argv[1];
    std::stringstream ss;
    int numTrials = 1000;

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
    double rx=0,ry=0,rz=0;
    uniformSampleSphere(drand48(), drand48(), &rx, &ry, &rz);
    //center of the spherical disc computed randomly on the sphere
    //vec3 ref(rx, ry, rz);
    double ref[] = {0.0,0.0,1.0};
    double theta = 60.0;
    //##########################################################

    //for(double theta = 30.0; theta <= 150.0; theta+=30.0)
    {

        std::ofstream fvar, fvarRef, fmean;
        ss.str(std::string());
        ss << datafiles << "variance-spherical-cap-theta-" << theta << "-" << samplingpattern << ".txt";
        fvar.open(ss.str().c_str(), std::ios::app);

        ss.str(std::string());
        ss << datafiles << "variance-withrealmean-spherical-cap-theta-" << theta << "-" << samplingpattern << ".txt";
        fvarRef.open(ss.str().c_str(), std::ios::app);

        ss.str(std::string());
        ss << datafiles << "mean-spherical-cap-theta-" << theta << "-" << samplingpattern << ".txt";
        fmean.open(ss.str().c_str(), std::ios::app);

        //##########################################################

        //Compute Radial Distribution Function of the spherical harmonic coeffs
        double surfaceAreaMean = 2*PI*(1-cos(theta*deg2rad));

        for(int nside = 1; nside <= 100; nside++){
            //fprintf(stderr,"\r Numpts (%d) ",N);

            //##########################################################

            double mean = 0;
            double variance = 0, varianceRef = 0;
            int N = 12 * nside * nside ;

            //##########################################################

            for(int trial = 1; trial <= numTrials; trial++){
                fprintf(stderr,"\r trial/Numpts/theta (%6d/%10d/%3.1f) ",trial, N, theta);

                //arr<vec3> sphereSamples;
                std::vector<double> sphereSamples;

                if(samplingpattern == "regular"){
                    sphereSamples = spherical_regular_samples(N, nside);
                }
                else if(samplingpattern == "stratified")
                    sphereSamples = spherical_stratified_samples(N);

                int tempN = sphereSamples.size() / 3.0;
                //##########################################################

                std::vector<double> rotatedSamples(3*tempN, 0.0);
                random_rotate_spherical_samples(rotatedSamples, sphereSamples, tempN);

                //##########################################################

                double integral = 0;
                for(int i=0; i < tempN; i++){
                    double point [] = {rotatedSamples[3*i+0],rotatedSamples[3*i+1],rotatedSamples[3*i+2]};
                    double value = spherical_cap_function(ref, point, theta);
                    integral += value;
                }
                double factor = (4*PI) / double(tempN);
                integral *= factor;

                //##########################################################

                //Incremental Mean and Variance computation
                mean = ((trial-1)/double(trial)) * mean + (1/double(trial)) * integral;
                if(trial < 2){
                    variance = 0;
                    varianceRef = 0;
                }
                else{
                    variance = ((trial-1)/double(trial)) * variance + (1/double(trial-1)) * (integral - mean) * (integral - mean);
                }

                //##########################################################
                varianceRef +=  (integral - surfaceAreaMean) * (integral - surfaceAreaMean);
            }
            varianceRef /= double(numTrials);

            fmean << std::fixed << N << std::setprecision(15) << " " << mean << std::endl;
            fvar << std::fixed << N << std::setprecision(15) << " " << variance << std::endl;
            fvarRef <<  std::fixed << N << std::setprecision(15) << " " << varianceRef << std::endl;

            //##########################################################
        }
        fvar.close();
        fmean.close();
        fvarRef.close();
    }
    std::cerr << std::endl;
    return 0;
}


