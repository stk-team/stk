#include "sampler-sphere.h"
#include <vector>
#include <stdlib.h>
#include <iostream>
#include "./../core/utils-sphere.h"
#include "./../core/constants.h"
#include "./../quaternions/quaternion.h"
#include "./../domaintransform/domaintransform.h"

void random_rotate_spherical_samples(std::vector<double> &outsamples, const std::vector<double> &insamples, int N){

    //First rotation randomly anywhere on the sphere
    double x=0, y=0, z=0;
    uniformSampleSphere(drand48(), drand48(), &x, &y, &z);
    double randvec[] = {x,y,z};
    normalize_vector(randvec);

    double vecref[] = {0.0,0.0,1.0};
    double quat[] = {0,0,0,1};
    QuaternionFromTwoVectors(vecref, randvec, quat);

    for(int i = 0; i< N; i++){
        double qmat[16], qvout[4];
        double qvec[] = {insamples[3*i+0], insamples[3*i+1], insamples[3*i+2],1};
        ConvertQuaternionToMatrix(quat, qmat);
        MultiplyVecMat(qmat, qvec, qvout);

        outsamples[3*i+0] =  qvout[0];
        outsamples[3*i+1] =  qvout[1];
        outsamples[3*i+2] =  qvout[2];
    }

    //Second rotation along phi[0,2*PI], can be called as Zonal Rotation.
    quat[0] = 0, quat[1]=0, quat[2]=0,quat[3]=0;
    double axis[] = {randvec[0]-0, randvec[1]-0, randvec[2]-0, 1.0};
    double angle = drand48() * 2 * PI;       //between [0,2*PI)
    SetQuaternionFromAxisAngle(axis, angle, quat);

    for(int i = 0; i< N; i++){
        double qmat[16], qvout[4];
        double qvec[] = {outsamples[3*i+0], outsamples[3*i+1], outsamples[3*i+2],1};
        ConvertQuaternionToMatrix(quat, qmat);
        MultiplyVecMat(qmat, qvec, qvout);

        outsamples[3*i+0] =  qvout[0];
        outsamples[3*i+1] =  qvout[1];
        outsamples[3*i+2] =  qvout[2];
    }
}

