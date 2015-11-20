/*** Created by Gurprit Singh ***/

#include "sphere-analysis.h"
#include "./../core/constants.h"
#include "./../core/utils-sphere.h"

#include <cmath>
#include <complex>
#include <iomanip>
#include <omp.h>
#include <iostream>

#include <boost/math/special_functions/spherical_harmonic.hpp>


double compute_omp(double* t_lm, int l, const std::vector<double> &data, int N, int nlmax);

void powspec_spherical_harmonics_parallel(double* t_lm, double* t_powspec, std::vector<double>& data, int N, int nlmax){

    for(int l = 0; l < nlmax; l++){
        //fprintf(stderr,"\rAt l: (%d) ", l);
        double myCl = compute_omp(t_lm, l, data, N, nlmax);

        //Normalization by the total number of rotations (m) for each l
        myCl /= double(2*l+1);
        ///Normalization to value 1 for any number of points N
        t_powspec[l] = (N/double(4*PI))*myCl;

        ///
        /// This step is necessary if we want to emphasize the SH definition for spectral analysis.
        /// alm *= sqrt(4*PI);
        /// In that case the power spectrum needs to be normalized with N/(4*pi)^2 .
        ///
        //t_powspec[l] = (N/double(16*PI*PI))*myCl;
    }
}

void powspec_projectivespace_sharmonics_parallel(double* t_lm, std::vector<double> &t_powspec,
                                                 const std::vector<double> &data, int N, int nlmax){
    int index = 0;
    for(int l = 0; l < nlmax; l+=2){
        //fprintf(stderr,"\r l: (%5d ) ", l);
        double myCl = compute_omp(t_lm, l, data, N, nlmax);

        //Normalization by the total number of rotations (m) for each l
        myCl /= double(2*l+1);
        ///Normalization to value 1 for any number of points N
        t_powspec[index] = (N/double(4*PI))*myCl;

        ///
        /// This step is necessary if we want to emphasize the SH definition for spectral analysis.
        /// Check above function: powspec_spherical_harmonics_spectral_normalization_parallel
        /// alm *= sqrt(4*PI);
        /// In that case the power spectrum needs to be normalized with N/(4*pi)^2 .
        ///
        //t_powspec[l] = (N/double(16*PI*PI))*myCl;
        index++;
    }
}


double compute_omp(double* t_lm, int l, const std::vector<double> &data, int N, int nlmax){
    double Cl = 0;

//    int height = 2*nlmax+1;
//    int width = 2*nlmax+1;
//    int half_height = 0.5*height;
//    int half_width = 0.5*width;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:Cl) num_threads(6)
    for(int m=-l; m<=l; m++){
        std::complex<double> alm(0.0,0.0);
        for(int i=0; i < N; i++){
            double phi=0, theta=0;
            xyz2thetaphi(theta, phi, data[3*i+0],data[3*i+1],data[3*i+2]);

            std::complex<double> sphcoeffs = boost::math::spherical_harmonic(l,m,theta, phi);
            alm += sphcoeffs;
        }
        alm *= ((4*PI)/double(N));

        ///
        /// This step is necessary if we want to emphasize the SH definition for spectral analysis.
        /// In that case the power spectrum needs to be normalized with N/(4*pi)^2 .
        ///alm *= sqrt(4*PI);
        ///
        //Compute power spectra = magnitude_square
        double pow_alm = std::norm(alm);

        Cl += pow_alm;
    }
    return Cl;
}

