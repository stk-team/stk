#ifndef UTILSSPHERE_H
#define UTILSSPHERE_H

#include <cmath>        // std::atan2
#include "./../core/constants.h"

inline double fmodulo (double v1, double v2){
    using namespace std;
    if (v1>=0)
        return (v1<v2) ? v1 : fmod(v1,v2);
    double tmp=fmod(v1,v2)+v2;
    return (tmp==v2) ? 0. : tmp;
    return (v1>=0) ? ((v1<v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2);
}

inline void normalize_theta(double& theta, double& phi){
    theta=fmodulo(theta,twopi);
    if (theta>M_PI)
    {
        phi+=M_PI;
        theta=twopi-theta;
    }
}

inline void normalize_theta_phi(double& theta, double &phi){
    normalize_theta(theta, phi);
    phi=fmodulo(phi,twopi);
}

inline double safe_atan2 (double y, double x){
    using namespace std;
    return ((x==0.) && (y==0.)) ? 0.0 : atan2(y,x);
}

inline void xyz2thetaphi (double& theta, double& phi, const double& x, const double& y, const double& z){
    theta = atan2(sqrt(x*x + y*y),z);
    phi = safe_atan2 (y,x);
    if (phi<0.) phi += twopi;
}

inline void thetaphi2xyz(double *x, const double &theta, const double &phi){
    //here theta[0,180] ,
    //where phi [0, 360];
    x[0] = cos(phi) * sin(theta);
    x[1] = sin(phi) * sin(theta);
    x[2] = cos(theta);
}

inline void normalize_vector(double* a){
    double adota=0;
    for(int k=0;k<3;k++)
        adota += a[k]*a[k];
    for(int k=0;k<3;k++){
        a[k] /= sqrt(adota);
    }
}

inline double length(double *a){
    double adota=0;
    for(int k=0;k<3;k++)
        adota += a[k]*a[k];
    return sqrt(adota);
}

inline double geodesic_sphere(double* a, double* b, const int ndims){

    double adotb = 0;
    for(int k=0; k<ndims;k++)
        adotb += a[k]*b[k];

    return acos(adotb);
}

inline void phi_theta_from_vec3 (double& phi, double& theta, const double& x, const double& y, const double& z){
    theta = std::atan2(sqrt(x*x + y*y),z);
    phi = safe_atan2 (y,x);
    if (phi<0.) phi += twopi;
}


#endif // UTILSSPHERE_H
