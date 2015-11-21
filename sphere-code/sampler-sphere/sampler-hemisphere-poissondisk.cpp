#include "sampler-sphere.h"
#include "./../core/utils-sphere.h"
#include "./../domaintransform/domaintransform.h"

#include <cmath>
#include <iostream>
#include <set>
#include <unordered_set>

//Slow Dart -throwing version of Poisson Disk Sampling.
//For Fast version, look the implementation after this one!
std::vector<double> hemispherical_dart_throwing_samples(const int &N){
    std::cerr << N << std::endl;
    double radius = 110.0/sqrt(N);
    std::vector<double> buffer;
    double sx,sy,sz;
    uniformSampleHemiSphere(drand48(), drand48(), &sx,&sy,&sz);
    buffer.push_back(sx);
    buffer.push_back(sy);
    buffer.push_back(sz);
    for(int i = 0; i < N*500; i++){
        double sx=0,sy=0,sz=0;
        uniformSampleHemiSphere(drand48(), drand48(), &sx,&sy,&sz);
        bool conflict = false;
        for(int j = 0; j < buffer.size(); j+=3){
            double ax[] =  {buffer[j+0],buffer[j+1],buffer[j+2]};
            double bx[] = {sx,sy,sz};

            normalize_vector(ax);
            normalize_vector(bx);

            double adotb = 0;
            for(int k=0; k<3;k++)
                adotb += ax[k]*bx[k];

            double geodesic = acos(adotb) * rad2deg;    //geodesic is in degrees

            if(geodesic < radius){
                conflict = true;
                break;
            }
        }
        if(!conflict){
            buffer.push_back(sx);
            buffer.push_back(sy);
            buffer.push_back(sz);
        }
        if(buffer.size() == (3*N))
            break;
    }
    //std::cerr << "# darts: " << buffer.size() / 3.0 << std::endl;

    int tempN = N;
    if(buffer.size() != (3*N)){
        tempN = buffer.size()/3.0;
        //sphereSamples.resize(tempN);
       std::cerr << "# samples generated: " << tempN << std::endl;
    }
    //vec3 initval(0,0,0);
    //arr<vec3> sphereSamples(tempN, initval);
    std::vector<double> dartSamples;
    for(int i=0; i < tempN; i++){
        dartSamples.push_back(buffer[3*i+0]);
        dartSamples.push_back(buffer[3*i+1]);
        dartSamples.push_back(buffer[3*i+2]);
        //std::cout << buffer[3*i] << " " << buffer[3*i+1] << " " << buffer[3*i+2] << " " << 1 << std::endl;
    }
    return dartSamples;
}

std::vector<double> hemisphere_to_sphere_mapping_dartthrowing_samples(const int &N){
    int tempN = 0.5 * N;

    std::vector<double> samples = hemispherical_dart_throwing_samples(tempN);

    tempN = samples.size() / 3.0;

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
