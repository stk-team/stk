#include "quaternion.h"
#include <math.h>

void quatNormalize(double *quat){
    double mag = sqrt(quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3]);
    for(int i=0;i<4;i++)
        quat[i] /= mag;
}

double quatLength(double *quat){
    return sqrt(quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3]);

}

// Routine to set a quaternion from a rotation axis and angle
// ( input axis = double[3] angle = double  output: quat = double[4] )
void SetQuaternionFromAxisAngle(const double *axis, double angle, double *quat)
{
    double sina2, norm;
    sina2 = (double)sin(0.5f * angle);
    norm = (double)sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    quat[0] = sina2 * axis[0] / norm;
    quat[1] = sina2 * axis[1] / norm;
    quat[2] = sina2 * axis[2] / norm;
    quat[3] = (double)cos(0.5f * angle);
}

void SetQuaternionFromScalarVector(const double& scalar, const double* vec, double *quat){
    quat[0] = vec[0];
    quat[1] = vec[1];
    quat[2] = vec[2];
    quat[3] = scalar;
}

// Routine to convert a quaternion to a 4x4 matrix
// ( input: quat = double[4]  output: mat = double[4*4] )
void ConvertQuaternionToMatrix(const double *quat, double *mat)
{
    double yy2 = 2.0f * quat[1] * quat[1];
    double xy2 = 2.0f * quat[0] * quat[1];
    double xz2 = 2.0f * quat[0] * quat[2];
    double yz2 = 2.0f * quat[1] * quat[2];
    double zz2 = 2.0f * quat[2] * quat[2];
    double wz2 = 2.0f * quat[3] * quat[2];
    double wy2 = 2.0f * quat[3] * quat[1];
    double wx2 = 2.0f * quat[3] * quat[0];
    double xx2 = 2.0f * quat[0] * quat[0];

    mat[0*4+0] = - yy2 - zz2 + 1.0f;
    mat[0*4+1] = xy2 + wz2;
    mat[0*4+2] = xz2 - wy2;
    mat[0*4+3] = 0;
    mat[1*4+0] = xy2 - wz2;
    mat[1*4+1] = - xx2 - zz2 + 1.0f;
    mat[1*4+2] = yz2 + wx2;
    mat[1*4+3] = 0;
    mat[2*4+0] = xz2 + wy2;
    mat[2*4+1] = yz2 - wx2;
    mat[2*4+2] = - xx2 - yy2 + 1.0f;
    mat[2*4+3] = 0;
    mat[3*4+0] = mat[3*4+1] = mat[3*4+2] = 0;
    mat[3*4+3] = 1;
}



// Routine to multiply 2 quaternions (ie, compose rotations)
// ( input q1 = double[4] q2 = double[4]  output: qout = double[4] )
void MultiplyQuaternions(const double *q1, const double *q2, double *qout)
{
    double qr[4];
    qr[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
    qr[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
    qr[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
    qr[3]  = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
    qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}

void MultiplyVecMat(const double *qmat, const double *qvec, double *qvout){
    //double qr[4];
    qvout[0] = qmat[0]*qvec[0]  + qmat[4]*qvec[1]  + qmat[8] *qvec[2]  + qmat[12]*qvec[3];
    qvout[1] = qmat[1]*qvec[0]  + qmat[5]*qvec[1]  + qmat[9] *qvec[2]  + qmat[13]*qvec[3];
    qvout[2] = qmat[2]*qvec[0]  + qmat[6]*qvec[1]  + qmat[10]*qvec[2]  + qmat[14]*qvec[3];
    qvout[3] = qmat[3]*qvec[0]  + qmat[7]*qvec[1]  + qmat[11]*qvec[2]  + qmat[15]*qvec[3];
    for(int i=0;i<4;i++)
        if(fabs(qvout[i])<1e-5)
            qvout[i] = 0;
    //qvout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}

//return quaternion quat that rotates vector u to v.
void QuaternionFromTwoVectors(double* u, double* v, double *quat){
    double dotproduu = 0.0, dotprodvv = 0.0, dotproduv = 0.0;
    for(int k=0;k<3;k++){
        dotproduu += u[k]*u[k];
        dotproduu += v[k]*v[k];
        dotproduv += u[k]*v[k];
    }
    double norm_u_norm_v = sqrt(dotproduu * dotprodvv);
    double real_part = norm_u_norm_v + dotproduv;
    double w[] = {0.,0.,0.};
    if (real_part < 1.e-6f * norm_u_norm_v)
    {
        /* If u and v are exactly opposite, rotate 180 degrees
         * around an arbitrary orthogonal axis. Axis normalisation
         * can happen later, when we normalise the quaternion. */
        real_part = 0.0f;
        //w = fabs(u[0]) > fabs(u[2]) ? vec3(-u[1], u[0], 0.f) : vec3(0.f, -u[2], u[1]);
        if(fabs(u[0]) > fabs(u[2])){
            w[0] = -u[1];
            w[1] = u[0];
            w[2] = 0.0;
        }
        else{
            w[0] = 0.0;
            w[1] = -u[2];
            w[2] = u[1];
        }
    }
    else{
        /* Otherwise, build quaternion the standard way. */
        //w = crossprod(u, v);
        w[0] = u[1]*v[2] - v[1]*u[2];
        w[1] = v[0]*u[2] - u[0]*v[2];
        w[2] = u[0]*v[1] - v[0]*u[1];
    }

    quat[0] = w[0];
    quat[1] = w[1];
    quat[2] = w[2];
    quat[3] = real_part;

    //Normalize the quaternion
    double mag = sqrt(quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3]);
    for(int i=0;i<4;i++)
        quat[i] /= mag;
    //return normalize(quat(real_part, w.x, w.y, w.z));
}
