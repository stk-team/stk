/*** Created by Gurprit Singh ***/

#include <cmath>
#include "domaintransform.h"
#include <algorithm>    // std::max

void uniformSampleHemiSphere(double u1, double u2, double* dx, double* dy, double* dz){
    *dz = u1;
    double r = sqrtf(std::max(0.0, 1-(*dz)*(*dz)));
    double phi = 2 * M_PI * u2;
    *dx = r * cos(phi);
    *dy = r * sin(phi);
}

void uniformSampleSphere(double u1, double u2, double* dx, double* dy, double* dz){
    *dz = (1 - 2 *u1);
    double r = sqrtf(std::max(0.0, 1-(*dz)*(*dz)));
    double phi = 2 * M_PI * u2;
    *dx = r * cos(phi);
    *dy = r * sin(phi);
}

void cosineSampleHemiSphere(double u1, double u2, double* dx, double* dy, double* dz){
    concentricSampleDisk(u1, u2, dx, dy);
    *dz = sqrtf(std::max(0., 1.0 - (*dx) * (*dx) - (*dy)*(*dy) ));
}
