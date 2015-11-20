#ifndef CCVTSPHERE_H
#define CCVTSPHERE_H

#include <vector>
#include <list>
#include <cmath>

#include "stk/pointset.hpp"

template<int DIM, typename POS>
double distance(const stk::Vector<DIM, POS>& a, const stk::Vector<DIM, POS>& b);

void generate_polar_sample_space(stk::PointSet2dd& sites, stk::PointSet2di& points,
                                 std::vector< std::vector<int> >& ptsPerSites,
                                 int nsites, int npts, std::string samplingpattern="whitenoise");

void generate_sample_from_file(stk::PointSet2dd& sites, stk::PointSet2di& points,
                                 std::vector< std::vector<int> >& ptsPerSites,
                                 int nsites, int npts, std::vector<double> &allpoints);


void coherent_initialization(stk::PointSet2dd& sites, stk::PointSet2di& points,
                                 std::vector< std::vector<int> >& ptsPerSites, int nsites,
                                 int npts, std::string samplingpattern);

template<int DIM, typename POS>
void ccvt_sphere(stk::PointSet<DIM, POS, double>& sites, //Final point position
                 const stk::PointSet<DIM, POS, int>& points, //Discretization
                 std::vector< std::vector<int> >& ptsPerSites);

#endif // CCVTSPHERE_H
