#ifndef SAMPLING_H
#define SAMPLING_H

void concentricSampleDisk(double u1, double u2, double *dx, double *dy);
void uniformSampleHemiSphere(double u1, double u2, double* dx, double* dy, double* dz);
void uniformSampleSphere(double u1, double u2, double* dx, double* dy, double* dz);
void cosineSampleHemiSphere(double u1, double u2, double* dx, double* dy, double* dz);

#endif // CORE_H
