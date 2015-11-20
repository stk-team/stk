#ifndef SPHERESUBDIVISION_H
#define SPHERESUBDIVISION_H

#include <string>

int numpts_to_nsubdivsions(const int &N, std::string domain="sphere");
void xy2phitheta(double &phi, double &theta, const double &x, const double& y);
void phitheta2xyz(double &x, double &y, double &z, const double &phi, const double &theta);
void uv2xy(double &x, double &y, const double &u, const double &v, int face, int nstep=1);

#endif // SPHERESUBDIVISION_H
