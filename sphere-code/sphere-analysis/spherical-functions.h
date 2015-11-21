#ifndef SPHERICAL_FUNCTIONS_H
#define SPHERICAL_FUNCTIONS_H

#include <string>
#include <complex>

double spherical_disc_northpole(double* coords, double thetaInDeg);
double spherical_cap_function(double* ref, double* p, double thetaInDeg);
double spherical_gaussian_function(double* ref, double* v, double bwidth);
std::complex<double> spherical_harmonic_function(double *v, int l, int m, std::string type="real");

#endif
