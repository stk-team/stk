#include "./../domaintransform/domaintransform.h"
#include "./../quaternions/quaternion.h"
#include "sampler-sphere.h"
#include <stdlib.h>
//#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

std::vector<double> randomshift_on_sphere(const std::vector<double> &insamples, int num_pts){

    std::vector<double> output(3*num_pts, 0.0);

    //Find two random vectors on the sphere
    double rx=0,ry=0,rz=0;
    uniformSampleSphere(drand48(), drand48(), &rx, &ry, &rz);
    double temp1[] = {rx, ry, rz};
    uniformSampleSphere(drand48(), drand48(), &rx, &ry, &rz);
    double temp2[]  = {rx, ry, rz};

    //Compute a quaternion from these vectors
    double quat[] = {0,0,0,1};
    QuaternionFromTwoVectors(temp1, temp2, quat);

    //Use the above quaternion to rotate all samples.
    for(int b = 0; b < num_pts; b++){
        double qmat[16], qvout[4];
        double qvec[] = {insamples[3*b+0], insamples[3*b+1], insamples[3*b+2],1};
        ConvertQuaternionToMatrix(quat, qmat);
        MultiplyVecMat(qmat, qvec, qvout);

        //rotated samples are here
        output[3*b+0] = (qvout[0]);
        output[3*b+1] = (qvout[1]);
        output[3*b+2] = (qvout[2]);
    }

    return output;
}

