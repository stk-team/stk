#ifndef QUATERNION_H
#define QUATERNION_H

void quatNormalize(double *quat);
double quatLength(double *quat);

void SetQuaternionFromAxisAngle(const double *axis, double angle, double *quat);
void SetQuaternionFromScalarVector(const double& scalar, const double* vec, double *quat);
void ConvertQuaternionToMatrix(const double *quat, double *mat);
void MultiplyQuaternions(const double *q1, const double *q2, double *qout);
void MultiplyVecMat(const double *qmat, const double *qvec, double *qvout);
void QuaternionFromTwoVectors(double* u, double* v, double *quat);

#endif // QUATERNION_H
