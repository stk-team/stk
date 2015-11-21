/*** Created by Gurprit Singh ***/

#include <cmath>
#include "domaintransform.h"
void concentricSampleDisk(double u1, double u2, double *dx, double *dy){
    double r, theta;

    //Map uniform random numbers to [-1,1]^2
    double sx = 2 * u1 - 1;
    double sy = 2 * u2 - 1;

    //Map square to (r,theta)
        // Handle degeneracy at the origin
    if (sx == 0.0 && sy == 0.0) {
        *dx = 0.0;
        *dy = 0.0;
        return;
    }
    if (sx >= -sy) {
        if (sx > sy) {
            // Handle first region of disk
            r = sx;
            if (sy > 0.0) theta = sy/r;
            else          theta = 8.0f + sy/r;
        }
        else {
            // Handle second region of disk
            r = sy;
            theta = 2.0f - sx/r;
        }
    }
    else {
        if (sx <= sy) {
            // Handle third region of disk
            r = -sx;
            theta = 4.0f - sy/r;
        }
        else {
            // Handle fourth region of disk
            r = -sy;
            theta = 6.0f + sx/r;
        }
    }
    theta *= M_PI / 4.f;
    *dx = r * cosf(theta);
    *dy = r * sinf(theta);
}
