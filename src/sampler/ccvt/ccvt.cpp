/*
  Copyright (C) 2009 Michael Balzer (michael.balzer@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdio>
#include <cstring>
#include <fstream>
#include <ctime>
#include "ccvt_metric.h"
#include "ccvt_optimizer.h"
#include "ccvt_point.h"
#include "ccvt_site.h"

#include <stk/stk.hpp>

using namespace ccvt;

// discrete space with constant density;
// the points form a regular grid
void constant_regular_density(Point2::List& points, const int numberOfPoints, const double torusSize) {
  double n = sqrt((double)(numberOfPoints));
  for (int x = 0; x < n; ++x) {
    for (int y = 0; y < n; ++y) {
      double dx = x / n * torusSize;
      double dy = y / n * torusSize;
      points.push_back(Point2(dx, dy));
    }
  }
}

int g_nPts;

extern "C" void module_sampler_init(int nPts)
{
	g_nPts = nPts;
	
	srand(time(NULL));
}

extern "C" void module_sampler_exec(stk::PointSet2dd& pts)
{
	int nbSamples = g_nPts;
	int nbLocations = 1024;
	
    // intializing the underlying discrete space
    Point2::List points;
    constant_regular_density(points, nbSamples * nbLocations, 1.0);

    // initializing the Voronoi sites with equal capacity
    unsigned int overallCapacity = (int)(points.size());
    Site<Point2>::List sites;
    for (int i = 0; i < nbSamples; ++i) {
      double x = (double)(rand() % RAND_MAX) / RAND_MAX;
      double y = (double)(rand() % RAND_MAX) / RAND_MAX;
      int capacity = overallCapacity / (nbSamples - i);
      overallCapacity -= capacity;
      sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
    }

    // initializing the CCVT
    Optimizer<Site<Point2>, Point2, MetricToroidalEuclidean2> optimizer;
    MetricToroidalEuclidean2 metric(Point2(1.0, 1.0));
    optimizer.initialize(sites, points, metric);

    // optimization
    int iteration = 0;
    bool stable;
    do
    {
      stable = optimizer.optimize(true);
    } while (!stable);

    // write output
    const Site<Point2>::Vector& result = optimizer.sites();
    for (unsigned int i = 0; i < result.size(); i++)
    {
		pts.push_back(stk::Point2dd(stk::Vector2d(result[i].location.x, result[i].location.y), 1.0));
	}
}

extern "C" void module_sampler_free(int nPts)
{
	
}
