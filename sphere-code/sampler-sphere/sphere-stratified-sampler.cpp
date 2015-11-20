/*** Created by Gurprit Singh ***/

#include "sampler-sphere.h"
#include "sphere-subdivision.h"
#include "./../core/utils.h"

#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>

boost::mt19937 mtrng( time( NULL ) );

void getSphericalStratifiedPatchPoints(std::vector<double>& vec, int face, int nsubdivisions, int nsteps){
    for(int uu = 1; uu <= nsubdivisions; uu++){
        for(int vv = 1; vv <= nsubdivisions; vv++){

            boost::random::uniform_real_distribution<> urange((uu-1)/double(nsubdivisions), uu/double(nsubdivisions) );
            boost::random::uniform_real_distribution<> vrange((vv-1)/double(nsubdivisions), vv/double(nsubdivisions) );

            //boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution <> > urandom( mtrng, urange);
            //boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution <> > vrandom( mtrng, vrange);
            double u = urange(mtrng);
            double v = vrange(mtrng);

            double x=0,y=0,z=0;
            double theta=0,phi=0;

            uv2xy(x, y, u, v, face, nsteps);
            xy2phitheta(phi, theta, x, y);
            phitheta2xyz(x, y, z, phi, theta);
           {
                vec.push_back(x);
                vec.push_back(y);
                vec.push_back(z);
            }
        }
    }
}


void getHemiSphericalStratifiedPatchPoints(std::vector<double>& vec, int face, int nsubdivisions, int nsteps){
    for(int uu = 1; uu <= nsubdivisions; uu++){
        for(int vv = 1; vv <= nsubdivisions; vv++){

            boost::random::uniform_real_distribution<> urange((uu-1)/double(nsubdivisions), uu/double(nsubdivisions) );
            boost::random::uniform_real_distribution<> vrange((vv-1)/double(nsubdivisions), vv/double(nsubdivisions) );

            //boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution <> > urandom( mtrng, urange);
            //boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution <> > vrandom( mtrng, vrange);
            double u = urange(mtrng);
            double v = vrange(mtrng);

            double x=0,y=0,z=0;
            double theta=0,phi=0;

            uv2xy(x, y, u, v, face, nsteps);
            xy2phitheta(phi, theta, x, y);
            phitheta2xyz(x, y, z, phi, theta);
           if(z >= 0){
                vec.push_back(x);
                vec.push_back(y);
                vec.push_back(z);
            }
        }
    }
}

//###########################################################################

std::vector<double> spherical_stratified_samples(const int &N, int Nside){
    int nsteps = 1;

    int nsubdivisions = Nside;
    if(Nside < 0)
        nsubdivisions = numpts_to_nsubdivsions(N);

    std::vector<double> samples;
    //Generate numpts
    for(int face = 1; face <= 12; face++){
        getSphericalStratifiedPatchPoints(samples, face, nsubdivisions, nsteps);
    }
    //int numpts = samples.size()/3.0;
    return samples;
}

std::vector<double> hemispherical_stratified_samples(const int &N, std::string domain){

    int nsteps = 1;
    int nsubdivisions = numpts_to_nsubdivsions(N,domain);
    std::vector<double> samples;
    for(int face = 1; face <= 12; face++){
        getHemiSphericalStratifiedPatchPoints(samples, face, nsubdivisions, nsteps);
    }

    return samples;
}

std::vector<double> hemisphere_to_sphere_mapping_stratified_samples(const int &N){
    std::vector<double> samples = hemispherical_stratified_samples(N, "sphere");
    int tempN = samples.size() / 3.0;

    std::vector<double> results;
    for(int i=0; i< tempN; i++){
        //std::cout << samples[3*i+0] << " "<< samples[3*i+1] << " "<< samples[3*i+2] << std::endl;
        double projective_space_sample[3] = {0.0,0.0,0.0};

        //Add hemispherical samples and their copies f(x) = f(-x) in results.
        for(int k=0; k<3; k++){
            projective_space_sample[k] = -samples[3*i+k];
            results.push_back(projective_space_sample[k]);
        }
        for(int k=0; k<3; k++){
            results.push_back(samples[3*i+k]);
        }
    }
    return results;
}




