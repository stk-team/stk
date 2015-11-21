
#include <boost/math/special_functions/legendre.hpp>
#include <cmath>
#include <complex>

#include "./../core/constants.h"
#include "sphere-analysis.h"
#include <iostream>

double powspec_spherical_cap(double theta0, int l){

    double fcap_coeff = (boost::math::legendre_p(l+1,0,std::cos(theta0 * deg2rad)) -
                                         boost::math::legendre_p(l-1,0,std::cos(theta0 *deg2rad)));

    double powspec = (4.0*PI*PI*std::abs(fcap_coeff)*std::abs(fcap_coeff))/(static_cast<double>(2*l+1)*static_cast<double>(2*l+1));
    return powspec;
}

