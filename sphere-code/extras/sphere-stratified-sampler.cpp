#include "./../sphere-harmonics/sphere-analysis.h"
#include "./../core/utils.h"
#include "./../core/constants.h"
#include "sphere-subdivision.h"

#include <boost/random.hpp>

boost::mt19937 mtrng( time( NULL ) );

void getSphericalStratifiedPatchPoints(std::vector<double>& vec, int face, int nsubdivisions, int nsteps);
void getHemiSphericalStratifiedPatchPoints(std::vector<double>& vec, int face, int nsubdivisions, int nsteps);


arr<vec3> spherical_stratified_samples(const int &N, int Nside){
    int nsteps = 1;
    int nsubdivisions = numpts_to_nsubdivsions(N);
    std::vector<double> samples;
    //Generate numpts
    for(int face = 1; face <= 12; face++){
        if(Nside >= 1)
            getSphericalStratifiedPatchPoints(samples, face, Nside, nsteps);
        else
            getSphericalStratifiedPatchPoints(samples, face, nsubdivisions, nsteps);
    }
    int numpts = samples.size()/3.0;

    vec3 initval(0,0,0);
    arr<vec3> sphereSamples(numpts, initval);

    for(int i=0;i<numpts;i++){
        sphereSamples[i].x = samples[3*i];
        sphereSamples[i].y = samples[3*i+1];
        sphereSamples[i].z = samples[3*i+2];
        //std::cout << sphereSamples[i].x <<" " << sphereSamples[i].y << " " << sphereSamples[i].z << " " << 1 << std::endl;
    }
    return sphereSamples;
}

void hemispherical_stratified_samples(arr<vec3> &sphereSamples, const int &N){
    int nsteps = 1;
    int nsubdivisions = numpts_to_nsubdivsions(N,"hemisphere");
    std::vector<double> samples;
    for(int face = 1; face <= 12; face++){
        getHemiSphericalStratifiedPatchPoints(samples, face, nsubdivisions, nsteps);
    }
    int numpts = samples.size()/3.0;

    sphereSamples.resize(numpts);

    for(int i=0;i<numpts;i++){
        sphereSamples[i].x = samples[3*i];
        sphereSamples[i].y = samples[3*i+1];
        sphereSamples[i].z = samples[3*i+2];
    }
}

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

