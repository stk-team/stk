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
#include <vector>

#include <stk/stk.hpp>

#include "DT.h"
#include "util.h"

int g_nPts;
MTRand rng(time(NULL));

enum Strategy {
    GLOBAL, LOCAL, HYBRID
};

static bool interruptFlag = false;

int mkdir(const char *name) {
    char command[64];
#ifdef _WIN32
    sprintf(command, "%s %s", "md"   , name);
#else
    sprintf(command, "%s %s", "mkdir", name);
#endif
    return system(command);
}

void generate_random(vector<Point_2> &points, unsigned npoints, MTRand &rng)
{
    points.reserve(npoints);
    for (unsigned i = 0; i < npoints; ++i) {
        Point_2 p(rng.randExc(), rng.randExc());
        points.push_back(p);
    }
}

void generate_darts(vector<Point_2> &points, unsigned npoints, MTRand &rng)
{
    const double rnorm  = 0.725;
    const double mdnorm = sqrt(SQRT3 * 0.5 * npoints);
    const double md     = rnorm / mdnorm;
    const double sqr_md = md * md;
    
    points.reserve(npoints);
    while (points.size() < npoints) {
        while (true) {
            Point_2 cand(rng.randExc(), rng.randExc());
            bool hit = true;
            for (unsigned i = 0; i < points.size(); ++i) {
                if (sqr_dist_unit_torus(cand, points[i]) < sqr_md) {
                    hit = false;
                    break;
                }
            }
            if (hit) {
                points.push_back(cand);
                break;
            }
        }
    }
}

void update_statistics(const DT &dt, unsigned it, double &global_md,
                       double &avg_md, bool output = true)
{
    dt.get_statistics(global_md, avg_md);
    
    const double mdnorm = sqrt(SQRT3 * 0.5 * dt.number_of_vertices());
    global_md *= mdnorm;
    avg_md    *= mdnorm;
    
    //~ if (output)
        //~ std::cout << std::fixed << std::setprecision(6) << std::setw(5)
        //~ << it << " " << global_md << " " << avg_md << "         \r";
}

static void sighandler(int sig)
{
    std::cout << "Aborting...\n";
    interruptFlag = true;
}

extern "C" void module_sampler_init(int nPts)
{
	g_nPts = nPts;
	
	srand(time(NULL));
}

extern "C" void module_sampler_exec(stk::PointSet2dd& pts)
{
    unsigned max_iter  = 10000;
    double max_md      = 0.925;
    Strategy strategy  = GLOBAL;
    bool silent        = false;
    unsigned npoints   = g_nPts;
    unsigned npointsets = 1;
    unsigned long seed = 0;
    
	double global_md = 0;
	double avg_md    = 0;
	double old_avg   = avg_md;
    
	// Initial points
	vector<Point_2> points;
	generate_random(points, npoints, rng);
	
	// Set up initial triangulation
	DT dt(points, true);
	
	update_statistics(dt, 0, global_md, avg_md, !silent);
	
	unsigned it = 0;
	vector<VH> neighbors;
	neighbors.reserve(10);
	
	// Set up vertex processing order
	// We choose a random order to ensure there is no correlation between the
	// original order of the point set and farthest-point optimization
	unsigned order[npoints];
	for (unsigned i = 0; i < npoints; ++i) order[i] = i;
	shuffle(order, npoints, rng);
	
	while (!interruptFlag)
	{
		// Main loop that moves each point to the farthest point, i.e. the
		// center of the largest circumcircle of a triangle in the DT
		for (unsigned i = 0; i < npoints; ++i)
		{
			// Pick removal candidate
			unsigned r = order[i];
			Point_2 cand = dt.get_vertex(r)->point();
			
			neighbors.clear();                 // Candidate's neighborhood
			dt.incident_vertices(r, neighbors);
			
			double cand_md = DBL_MAX;          // Candidate's local mindist
			for (unsigned i = 0; i < neighbors.size(); ++i) {
				cand_md = std::min(cand_md, 
						 CGAL::squared_distance(cand, neighbors[i]->point()));
			}
			
			Circle_2 c(cand, cand_md);         // Empty circle about candidate

			// Remove candidate
			dt.clear_vertex(r); 
			
			// Search for largest circumcircle in triangulation
			FH face;
			Circle_2 l = Circle_2(Point_2(0, 0), 0);
			switch (strategy) {
				case GLOBAL:
					l = dt.global_largest_circumcircle(face);
					break;
				case LOCAL:
					l = dt.local_largest_circumcircle(neighbors, face);
					break;
				case HYBRID:
					if (old_avg < 0.930)
						l = dt.global_largest_circumcircle(face);
					else
						l = dt.local_largest_circumcircle(neighbors, face);
					break;
			}

			// Set center of largest circumcircle as new vertex
			if (l.squared_radius() > c.squared_radius()) {
				Point_2 l_center = wrap_unit_torus(l.center());
				dt.set_vertex(r, l_center, face);
			} else
				dt.set_vertex(r, c.center(), face);
		}
		++it;
		update_statistics(dt, it, global_md, avg_md, !silent);
		
		if (it >= max_iter || global_md >= max_md)
			break;

		if (avg_md - old_avg == 0.0) break;
		old_avg = avg_md;
	}
	
	dt.save_vertices(pts);
}

extern "C" void module_sampler_free(int nPts)
{
	
}
