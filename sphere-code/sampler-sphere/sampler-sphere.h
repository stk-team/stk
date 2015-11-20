/*** Created by Gurprit Singh ***/
#ifndef SAMPLE_SPHERE_H
#define SAMPLE_SPHERE_H

#include <vector>
#include <string>

struct cell{
    bool valid, visited;
    int cellid,r,c,s;
    std::vector<double> samples;
    cell(){
        valid = false;
        visited = false;
        cellid = 0;
        r = 0;
        c = 0;
        s = 0;
    }
};

std::vector<double> spherical_whitenoise_samples(const int &N);
std::vector<double> spherical_regular_samples(const int &N, int Nside=-1, bool randomshift=false);
std::vector<double> spherical_stratified_samples(const int &N, int Nside=-1);
std::vector<double> spherical_dart_throwing_samples(const int &N); //Slow version
std::vector<double> spherical_poisson_disk_samples(const int &N, const double *sampleSpace, int sampleSpaceSize); //Light Speed Version

//Hemispherical Sampler
std::vector<double> hemispherical_regular_samples(const int &N, int Nside=-1, bool randomshift=false);
std::vector<double> hemispherical_whitenoise_samples(const int &N);
std::vector<double> hemispherical_stratified_samples(const int &N, std::string domain="hemisphere");
std::vector<double> hemispherical_dart_throwing_samples(const int &N);

//Hemispherical to spherical Projective 2-space sampling
std::vector<double>hemisphere_to_sphere_mapping_regular_samples_(const int &N, int Nside=-1, bool randomshift=false);
std::vector<double> hemisphere_to_sphere_mapping_whitenoise_samples(const int &N);
std::vector<double> hemisphere_to_sphere_mapping_stratified_samples(const int &N);
std::vector<double> hemisphere_to_sphere_mapping_dartthrowing_samples(const int &N);
std::vector<double> projective_space_hemisphere_to_sphere_samples(const std::vector<double> &sphereSamples, int tempN);

//utils
void random_rotate_spherical_samples(std::vector<double> &outsamples, const std::vector<double> &insamples, int N);
std::vector<double> randomshift_on_sphere(const std::vector<double> &insamples, int num_pts);
#endif
