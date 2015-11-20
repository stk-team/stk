#include "spherical-functions.h"
#include "./../core/utils-sphere.h"
#include "./../core/utils.h"

#include <complex>      // std::complex, std::abs
#include <boost/math/special_functions/spherical_harmonic.hpp>


double spherical_disc_northpole(double* coords, double thetaInDeg){
    double theta=0, phi=0;
    xyz2thetaphi(theta, phi, coords[0], coords[1],coords[2]);
    if((theta * rad2deg) <= thetaInDeg)
        return 1.0;
    else
        return 0.0;
}

double spherical_cap_function(double* ref, double* p, double thetaInDeg){
    double geodesic = geodesic_sphere(ref,p,3);
    if(geodesic * rad2deg <= thetaInDeg){
        //std::cout << geodesic * rad2deg << std::endl;
        return 1.0;
    }
    else
        return 0.0;
}

double spherical_gaussian_function(double *ref, double *v, double bwidth){

    double geodesic = geodesic_sphere(ref, v, 3);
    //double sphgauss = exp(-bwidth)*exp(bwidth*geodesic);
    double sphgauss = exp(-bwidth*geodesic);
    return sphgauss;
}

std::complex<double> spherical_harmonic_function(double *v, int l, int m, std::string type){

    double theta=0,phi=0;
    xyz2thetaphi(theta, phi, v[0],v[1],v[2]);

    std::complex<double> sphcoeffs;

    if(type == "real")
        sphcoeffs = boost::math::spherical_harmonic_r(l, m, theta, phi);
    else if(type == "imag")
        sphcoeffs = boost::math::spherical_harmonic_i(l, m, theta, phi);
    else
        sphcoeffs = boost::math::spherical_harmonic(l, m, theta, phi);

    return sphcoeffs;
}

