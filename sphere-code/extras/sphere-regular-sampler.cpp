#include "sphere-subdivision.h"
#include "sampler-sphere.h"
#include <healpix_base.h>
#include <vector>

void getSphericalRegularPatchPoints(std::vector<double>& vec, int face, int nsubdivisions,int nsteps);
void getHemiSphericalRegularPatchPoints(std::vector<double>& vec, int face, int nsubdivisions,int nsteps);

arr<vec3> spherical_regular_samples(const int &N, int Nside, bool randomshift){
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

    vec3 initval(0,0,0);
    arr<vec3> sphereSamples(numpts, initval);

    //Preparing pointset for output

    if(randomshift){
        sphereSamples = randomshift_on_sphere(samples, numpts);
    }
    else{
        for(int i=0;i<numpts;i++){
            sphereSamples[i].x = samples[3*i+0];
            sphereSamples[i].y = samples[3*i+1];
            sphereSamples[i].z = samples[3*i+2];
        }
    }
    return sphereSamples;
}

arr<vec3> hemispherical_regular_samples(const int &N, int Nside, bool randomshift){
    int nsteps = 1;
    int nsubdivisions = Nside;

    if(Nside < 0)
        nsubdivisions = numpts_to_nsubdivsions(N);

    std::vector<double> samples;

    //Generate numpts
    for(int face = 1; face <= 12; face++){
        getHemiSphericalRegularPatchPoints(samples, face, nsubdivisions, nsteps);
    }
    int numpts = samples.size()/3.0;

    vec3 initval(0,0,0);
    arr<vec3> sphereSamples(numpts, initval);

    //Preparing pointset for output
    for(int i=0;i<numpts;i++){
        sphereSamples[i].x = samples[3*i];
        sphereSamples[i].y = samples[3*i+1];
        sphereSamples[i].z = samples[3*i+2];
        //std::cout << sphereSamples[i].x <<" " << sphereSamples[i].y << " " << sphereSamples[i].z << " " << 1 << std::endl;
    }
    return sphereSamples;
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

