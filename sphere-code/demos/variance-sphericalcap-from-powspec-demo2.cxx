/* Compute power spectrum from built-in sampling methods.
 * Created by Gurprit Singh
 * */

#include "./../sphere-analysis/sphere-analysis.h"
#include "./../sampler-sphere/sampler-sphere.h"
#include "./../core/utils.h"
#include "./../core/constants.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <sys/stat.h>

int main(int argc, char* argv[]){

    if(argc < 5){
        std::cerr << "Parameters: profile samplingpattern limit alpha" << std::endl;
        return 0;
    }

    //##########################################################

    std::string profile = argv[1];
    std::string samplingpattern = argv[2];
    double limit = atof(argv[3]);
    double alpha = atof(argv[4]);
    int nlmax = 256;
    double theta = 30.0;
    std::stringstream ss;

    //###############Creating Folders###################################
    std::string datafiles, images, graphs;
    std::string resultFolder = "./../results/";
    ss.str(std::string());
    ss << resultFolder << "variance-sphericalcap-from-powspec" << "/";
    mkdir(ss.str().c_str(),0755);

#ifdef __APPLE__
    create_folders(ss.str(), datafiles, images, graphs);
#elif __linux__
    create_folders(ss.str(), datafiles, images, graphs);
#endif
    //#####################Declaring Variables############################

    srand48(time(NULL));

    std::ofstream fprofile, fpowspec, fvarpowspec;
    ss.str(std::string());
    ss << datafiles << profile << "-" << samplingpattern << "-limit-" << limit << "-alpha-" << alpha << "-l" << nlmax << ".txt";
    fprofile.open(ss.str().c_str());

    ss.str(std::string());
    ss << datafiles << "powspec-scap-"  << profile << "-" << samplingpattern << "-limit-" << limit << "-alpha-" << alpha <<  "-l" << nlmax << ".txt";
    fpowspec.open(ss.str().c_str());

    ss.str(std::string());
    ss << datafiles << "variance-scap-"  << profile << "-" << samplingpattern << "-limit-" << limit << "-alpha-" << alpha <<  "-l" << nlmax << ".txt";
    fvarpowspec.open(ss.str().c_str());


    for(int l=0; l < nlmax; l++){
        double value = 0;
        if(profile == "constant")
            value = constantProfile(l, limit);
        else if(profile == "linear")
            value = linearProfile(l, limit,alpha,4096);
        else if(profile == "quadratic")
            value = quadraticProfile(l, limit,alpha,4096);

        fprofile << l << " "<< value << std::endl;
        fpowspec << l <<" " << powspec_spherical_cap(theta,l) << std::endl;
    }
    fprofile.close();
    fpowspec.close();
    for(int N=1; N < 200000;){
        fprintf(stderr,"\r Numpts (%d) ",N);
        double variance = 0;
        for(int l=0; l < nlmax; l++){
            if(profile == "linear")
                variance += (1.0/(2*l+1)) * linearProfile(l, limit,alpha,N) * powspec_spherical_cap(theta, l);

            else if(profile == "quadratic")
                variance += (1.0/(2*l+1)) * quadraticProfile(l, limit,alpha,N) * powspec_spherical_cap(theta, l);
        }
        variance *= 16*M_PI*M_PI/N;
        fvarpowspec << N << " " << variance << std::endl;

        //##############################################################################
        if(N >= 1000)
            N += 2000;
        else if(N >= 500 && N < 1000){
            N+=250;
        }
        else if(N >= 150 && N < 500){
            N+=10;
        }
        else if(N < 150){
            N += 5;
        }
        //##############################################################################
    }

    fprofile.close();
    fpowspec.close();
    fvarpowspec.close();

    return 0;
}
