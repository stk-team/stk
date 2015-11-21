#ifndef PROJECT_SPHERE
#define PROJECT_SPHERE

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <complex>

//void powspec_spherical_harmonics(arr<double> &powbuffer, const arr<vec3> &data, int N, int nlmax);
void powspec_spherical_harmonics_parallel(double* t_lm, double* t_powspec, std::vector<double> &data, int N, int nlmax);
double powspec_spherical_cap(double theta0, int l);
void powspec_projectivespace_sharmonics_parallel(double* t_lm, std::vector<double> &t_powspec,
                                                 const std::vector<double> &data, int N, int nlmax);


//profiles
double constantProfile(int l, double a);
double linearProfile(int l, double limit, double alpha,int n);
double quadraticProfile(int l, double limit, double alpha,int n);
#endif
