/*** Created by Gurprit Singh ***/

#include "sampler-sphere.h"
#include "./../domaintransform/domaintransform.h"
#include <stdlib.h>

std::vector<double> spherical_whitenoise_samples(const int &N){
    //vec3 initval(0,0,0);
    //arr<vec3>  sphsamples(N, initval);
    std::vector<double> sphsamples;
    for(int i=0; i<N; i++){
        double sx=0,sy=0,sz=0;
        uniformSampleSphere(drand48(), drand48(), &sx,&sy,&sz);
        sphsamples.push_back(sx);
        sphsamples.push_back(sy);
        sphsamples.push_back(sz);
    }
    //std::cerr << "Whitenoise sphere samples: " << N << std::endl;
    return sphsamples;
}

std::vector<double> hemispherical_whitenoise_samples(const int &N){
    std::vector<double> sphsamples(3*N,0.0);
    for(int i=0; i<N; i++){
        double sx=0,sy=0,sz=0;
        uniformSampleHemiSphere(drand48(), drand48(), &sx,&sy,&sz);
        sphsamples[3*i+0] = sx;
        sphsamples[3*i+1] = sy;
        sphsamples[3*i+2] = sz;
    }
    //std::cerr << "Whitenoise sphere samples: " << N << std::endl;
    return sphsamples;
}

std::vector<double> hemisphere_to_sphere_mapping_whitenoise_samples(const int &N){
    int tempN = 0.5*N;
    std::vector<double> samples = hemispherical_whitenoise_samples(tempN);
    std::vector<double> results;
    for(int i=0; i< tempN; i++){
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

std::vector<double> projective_space_hemisphere_to_sphere_samples(const std::vector<double> &hemiSphereSamples, int tempN){
    std::vector<double> results;
    for(int i=0; i< tempN; i++){
        double projective_space_sample[3] = {0.0,0.0,0.0};
        //Add hemispherical samples and their copies f(x) = f(-x) in results.
        for(int k=0; k<3; k++){
            projective_space_sample[k] = -hemiSphereSamples[3*i+k];
            results.push_back(projective_space_sample[k]);
        }
        for(int k=0; k<3; k++){
            results.push_back(hemiSphereSamples[3*i+k]);
        }
    }
    return results;
}

