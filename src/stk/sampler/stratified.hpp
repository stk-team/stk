#ifndef __STK_SAMPLER_STRATIFIED__
#define __STK_SAMPLER_STRATIFIED__

#include <cmath>

#include <stk/pointset.hpp>
#include <stk/sampler/grid.hpp>
#include <stk/voronoi.hpp>

namespace stk
{

namespace sampler
{

template<int DIM, typename POS, typename VAL>
void stratified(PointSet<DIM, POS, VAL>& o_pts, const Vector<DIM, int>& i_npts)
{
	grid(o_pts, i_npts);

	Vector<DIM, POS> step;
	const Vector<DIM, POS>& s = o_pts.boundingBoxSize();
	int npts = i_npts.volume();

	for(int i=0; i<DIM; i++)
		step[i] = s[i]/i_npts[i];

	for(int i=0; i<npts; i++)
	{
		for(int j=0; j<DIM; j++)
		{
			o_pts[i].pos(j) += step[j]*(drand48()-0.5);
		}
	}
}

template<int DIM, typename POS, typename VAL>
void stratified(PointSet<DIM, POS, VAL>& o_pts, const int i_npts)
{
	stratified(o_pts, Vector<DIM, int>(pow(i_npts, 1.0/DIM)));
}

template<typename POS>
bool voronoiStratifiedTest(
	const std::vector< stk::Vector<2, POS> >& polygon,
	const stk::Vector<2, POS>& pt)
{
	int sz = polygon.size();
	for(int j=0; j<sz; j++)
	{
		stk::Vector<2, POS> a = polygon[j];
		stk::Vector<2, POS> b = polygon[(j+1)%sz];
		stk::Vector<2, POS> seg = b-a;
		if(seg.norm() < 1e-10) continue;
		stk::Vector<2, POS> rPt = pt-a;

		if(seg[0]*rPt[1]-seg[1]*rPt[0] < 0) return false;
	}
	return true;
}

template<typename POS, typename VAL>
void voronoiStratified(PointSet<2, POS, VAL>& o_pts, const PointSet<2, POS, VAL>& i_pts, double strengh = 1.0)
{
	o_pts.resize(i_pts.size());
	
	std::vector<stk::VoronoiCell<POS, VAL> > voronoi;
	stk::voronoi(i_pts, voronoi);

	for(int i=0; i<i_pts.size(); i++)
	{
		const std::vector< stk::Vector<2, POS> >& polygon = voronoi[i].getPolygon();
		
		stk::Vector<2, POS> bbMin = i_pts[i].pos();
		stk::Vector<2, POS> bbMax = i_pts[i].pos();
		for(int j=0; j<polygon.size(); j++)
		{
			if(bbMin[0] > polygon[j][0]) bbMin[0] = polygon[j][0];
			if(bbMax[0] < polygon[j][0]) bbMax[0] = polygon[j][0];
			if(bbMin[1] > polygon[j][1]) bbMin[1] = polygon[j][1];
			if(bbMax[1] < polygon[j][1]) bbMax[1] = polygon[j][1];
		}
		stk::Vector<2, POS> bbSize = bbMax - bbMin;
		
		stk::Vector<2, POS> r;
		int ok = 1;
		do
		{
			if(ok > 100) break;
			
			r = stk::Vector<2, POS>(bbMin[0]+drand48()*bbSize[0], bbMin[1]+drand48()*bbSize[1]);

			if(voronoiStratifiedTest(polygon, r)) ok = 0;
			else ok++;
		}
		while(ok);
		if(ok) r = i_pts[i].pos();
		
		stk::Vector<2, POS> vorCenter = voronoi[i].getCentroid();
		o_pts[i].pos() = ( r-vorCenter ) * strengh + vorCenter;
		o_pts[i].val() = o_pts.getDefaultVal();
	}
	
	o_pts.normalize();
}

}

}

#endif
