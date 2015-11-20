/* Read n files, with each file containing point samples: x y z
 * and generate corresponding power spectrum 1D plot.
 * */

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../io/read-pointset.h"
#include "./../core/utils.h"
#include "./../core/utils-samples.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc != 6){
        std::cerr << "Parameters: filename samplingPattern N nlmax numTrials" << std::endl;
        exit(-2);
    }

    //########################################################33
    srand48(time(NULL));
    std::string filename = argv[1];
    std::string samplingpattern = argv[2];
    int approxN = atoi(argv[3]);
    int nlmax = atoi(argv[4]);
    int numTrials = atoi(argv[5]);

    //########################################################33
    std::vector<double> vec;
    read_pointsetnD(filename, vec, 3);
    int N = vec.size()/3.0;
    std::stringstream ss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    ss.str(std::string());
    ss << resultFolder << "powspec-projectivespace-hemisphere-to-sphere-" << samplingpattern << "-n" << approxN  << "/";
    mkdir(ss.str().c_str(),0755);

    create_folders(ss.str(), datafiles, images, graphs);

    //#####################Declaring Variables############################

    double* lmAccum = new double[nlmax*nlmax] ();
    double* t_lm = new double[nlmax*nlmax] ();

    int powSpeclength = 0.5*nlmax;  //since all odd l valued frequencies are null!
    std::vector<double> powSpec(powSpeclength, 0.0);
    std::vector<double> powSpecAccum(powSpeclength, 0.0);


    //########################################################33

    std::cerr << "Number of points: " << N << std::endl;

    for(int trial=1; trial<= numTrials; trial++){
        fprintf(stderr,"\rIterations (%d) ",trial);

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
//        std::cerr << filename << std::endl;
//        continue;
        //exit(-2);
        //########################################################

        std::vector<double> sphereSamples, projectiveSpaceSamples, hemiSphereSamples;
        read_pointsetnD(filename, sphereSamples, 3);
        //########################################################
        N = sphereSamples.size() / 3.0;

        for(int i=0; i< N; i++){
            if(sphereSamples[3*i+2] >= 0){
                hemiSphereSamples.push_back(sphereSamples[3*i+0]);
                hemiSphereSamples.push_back(sphereSamples[3*i+1]);
                hemiSphereSamples.push_back(sphereSamples[3*i+2]);
            }
        }
        N = hemiSphereSamples.size() / 3.0;

        projectiveSpaceSamples = projective_space_hemisphere_to_sphere_samples(hemiSphereSamples, N);
        //write_samples(projectiveSpaceSamples,3);
        N = projectiveSpaceSamples.size() / 3.0;
        //exit(-2);
        //########################################################

        if(N==0){
            std::cerr << "No Data available !!!" << std::endl;
            return 0;
        }

        for(int i=0; i < powSpeclength; i++)
            powSpec[i] = 0.0;

        //########################################################
        //std::cerr << "N : " << N << std::endl;
        //Compute SH coeffs for a given l and m and corresponding Power Spectra from these values
        //harmonics_powspec(lm, powSpec, sphereSamples, N, nlmax);
        powspec_projectivespace_sharmonics_parallel(t_lm, powSpec, projectiveSpaceSamples, N, nlmax);

        //########################################################33

        for(int i=0; i < powSpeclength; i++)
            powSpecAccum[i] += powSpec[i];

        //##########################################################

        //Normalize powSpecAccumulator
        //if(iter%9 == 0 || iter%10 == 0 || iter == 1 || iter == 0)
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
    //##########################################################

    delete [] lmAccum;
    delete [] t_lm;
    return 0;
}



