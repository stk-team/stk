#include "sphere-subdivision.h"
#include "sampler-sphere.h"
#include <vector>

void getSphericalRegularPatchPoints(std::vector<double>& vec, int face, int nsubdivisions,int nsteps);
void getHemiSphericalRegularPatchPoints(std::vector<double>& vec, int face, int nsubdivisions,int nsteps);

std::vector<double> spherical_regular_samples(const int &N, int Nside, bool randomshift){
    int nsteps = 1;
    int nsubdivisions = Nside;

    if(Nside < 0)
        nsubdivisions = numpts_to_nsubdivsions(N);

    std::vector<double> samples;

    //Generate numpts

    for(int face = 1; face <= 12; face++){
        getSphericalRegularPatchPoints(samples, face, nsubdivisions, nsteps);
    }

    int numpts = samples.size()/3.0;
    //Preparing pointset for output
    if(randomshift){
        std::vector<double> sphereSamples = std::vector<double>();
        sphereSamples = randomshift_on_sphere(samples, numpts);
        return sphereSamples;
    }
    else{
        return samples;
    }
}

std::vector<double> hemispherical_regular_samples(const int &N, int Nside, bool randomshift){
    int nsteps = 1;
    int nsubdivisions = Nside;

    if(Nside < 0)
        nsubdivisions = numpts_to_nsubdivsions(N);

    std::vector<double> samples;

    //Generate numpts
    for(int face = 1; face <= 12; face++){
        getHemiSphericalRegularPatchPoints(samples, face, nsubdivisions, nsteps);
    }
    //int numpts = samples.size()/3.0;
    return samples;
}


void getSphericalRegularPatchPoints(std::vector<double>& vec, int face, int nsubdivisions,int nsteps){
    for(int uu = 1; uu <= nsubdivisions; uu++){
        for(int vv = 1; vv <= nsubdivisions; vv++){
            double u=0, v=0;

            u = (uu-0.5)/double(nsubdivisions);
            v = (vv-0.5)/double(nsubdivisions);

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

void getHemiSphericalRegularPatchPoints(std::vector<double>& vec, int face, int nsubdivisions,int nsteps){
    for(int uu = 1; uu <= nsubdivisions; uu++){
        for(int vv = 1; vv <= nsubdivisions; vv++){
            double u=0, v=0;

            u = (uu-0.5)/double(nsubdivisions);
            v = (vv-0.5)/double(nsubdivisions);

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

