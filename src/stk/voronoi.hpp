#ifndef __STK_VORONOI__
#define __STK_VORONOI__

#include <vector>

#include <stk/pointset.hpp>

namespace stk
{

template<typename POS, typename VAL>
class VoronoiCell
{
	private:
		std::vector<Vector<2, POS> > m_polygon;
		Vector<2, POS> m_pt;
		Vector<2, POS> m_centroid;
	
	public:
		VoronoiCell(const Vector<2, POS>& i_pt);
		
		void addVertex(const Vector<2, POS>& i_vertex);
		void close();
		
		const Vector<2, POS>& getCentroid() const;
		POS getArea() const;
		const Vector<2, POS>& getPoint() const;
		const std::vector<Vector<2, POS> >& getPolygon() const;
};

typedef VoronoiCell<double, double> VoronoiCell2d;
typedef VoronoiCell<double, double> VoronoiCell2dd;
typedef VoronoiCell<double, int> VoronoiCell2di;

template<typename POS, typename VAL>
VoronoiCell<POS, VAL>::VoronoiCell(const Vector<2, POS>& i_pt) :
	m_pt(i_pt)
{
	
}

template<typename POS, typename VAL>
void VoronoiCell<POS, VAL>::addVertex(const Vector<2, POS>& i_vertex)
{
	m_polygon.push_back(i_vertex);
}

template<typename POS, typename VAL>
void VoronoiCell<POS, VAL>::close()
{
	double a = 0;
	double m;
	
	if(m_polygon.size() > 1)
	{
		m_centroid = Vector<2, POS>(0, 0);
		
		for(int i=1; i<m_polygon.size(); i++)
		{
			m  = m_polygon[i-1][0] * m_polygon[i][1];
			m -= m_polygon[i-1][1] * m_polygon[i][0];
			
			a += 0.5 * m;
			
			m_centroid += (m_polygon[i-1] + m_polygon[i]) * m;
		}
		
		int last = m_polygon.size()-1;
		
		m  = m_polygon[last][0] * m_polygon[0][1];
		m -= m_polygon[0][0]   * m_polygon[last][1];
		
		a += 0.5 * m;
		
		m_centroid[0] += m * (m_polygon[last][0] + m_polygon[0][0]);
		m_centroid[1] += m * (m_polygon[last][1] + m_polygon[0][1]);
	
		if (a != 0.0)
		{
			m_centroid[0] /= 6*a;
			m_centroid[1] /= 6*a;
		}
	}
	else m_centroid = m_pt;
}

template<typename POS, typename VAL>
const Vector<2, POS>& VoronoiCell<POS, VAL>::getCentroid() const
{
	return m_centroid;
}

template<typename POS, typename VAL>
const Vector<2, POS>& VoronoiCell<POS, VAL>::getPoint() const
{
	return m_pt;
}

template<typename POS, typename VAL>
POS VoronoiCell<POS, VAL>::getArea() const
{
	POS area = 0;
	
	for(int i=1; i<m_polygon.size(); i++)
	{
		area += m_polygon[i-1][0] * m_polygon[i][1];
		area -= m_polygon[i][0] * m_polygon[i-1][1];
	}
	
	area += m_polygon[m_polygon.size()-1][0] * m_polygon[0][1];
	area -= m_polygon[0][0] * m_polygon[m_polygon.size()-1][1];
	
	return area/2.0;
}

template<typename POS, typename VAL>
const std::vector<Vector<2, POS> >& VoronoiCell<POS, VAL>::getPolygon() const
{
	return m_polygon;
}

namespace doNotUse
{

class _cgalWeightedVoronoi
{
	public:
		void* m_tri;
		void* m_vertices;
		int m_counter;
		
	public:
		_cgalWeightedVoronoi(int i_size);
		~_cgalWeightedVoronoi();
		void add(double i_x, double i_y, double i_w = 1.0);
		void get(std::vector<std::vector<Vector2d> >& i_polygons);
};

}

template<typename POS, typename VAL>
void weightedVoronoi(
	const PointSet<2, POS, VAL>& i_pts,
	std::vector<VoronoiCell<POS, VAL> >& o_polygons)
{
	std::vector< std::vector<Vector2d> > polygons;
	
	if(i_pts.isToroidal())
	{
		std::vector< Vector<2, POS> > replication;
		std::vector< int > replicationId;
		typename PointSet<2, POS, VAL>::replicationIterator iter;
		for(iter = i_pts.replicationBegin(true, 7.0); iter != i_pts.replicationEnd(); iter++)
		{
			replication.push_back(iter.pos());
			replicationId.push_back(iter.id());
		}
		
		polygons.resize(replication.size());
		doNotUse::_cgalWeightedVoronoi cgalVor(replication.size());
		for(int i=0; i<replication.size(); i++)
		{
			cgalVor.add(replication[i][0], replication[i][1], i_pts[replicationId[i]].val());
		}
		cgalVor.get(polygons);
	}
	else
	{
		polygons.resize(i_pts.size());
		doNotUse::_cgalWeightedVoronoi cgalVor(i_pts.size());
		for(int i=0; i<i_pts.size(); i++)
		{
			cgalVor.add(i_pts[i].pos()[0], i_pts[i].pos()[1], i_pts[i].val());
		}
		cgalVor.get(polygons);
	}

	for(int i=0; i<i_pts.size(); i++)
	{
		VoronoiCell<POS, VAL> cell(i_pts[i].pos());
		for(int j=0; j<polygons[i].size(); j++)
		{
			cell.addVertex(
				Vector<2, POS>(
					(POS) polygons[i][j][0],
					(POS) polygons[i][j][1]));
		}
		
		cell.close();
		o_polygons.push_back(cell);
	}
}

template<typename POS, typename VAL>
void voronoi(
	const PointSet<2, POS, VAL>& i_pts,
	std::vector<VoronoiCell<POS, VAL> >& o_polygons)
{
	std::vector< std::vector<Vector2d> > polygons;
	
	if(i_pts.isToroidal())
	{
		std::vector< Vector<2, POS> > replication;
		typename PointSet<2, POS, VAL>::replicationIterator iter;
		for(iter = i_pts.replicationBegin(true, 7.0); iter != i_pts.replicationEnd(); iter++)
		{
			replication.push_back(iter.pos());
		}
		
		polygons.resize(replication.size());
		doNotUse::_cgalWeightedVoronoi cgalVor(replication.size());
		for(int i=0; i<replication.size(); i++)
		{
			cgalVor.add(replication[i][0], replication[i][1], 1.0);
		}
		cgalVor.get(polygons);
	}
	else
	{
		polygons.resize(i_pts.size());
		doNotUse::_cgalWeightedVoronoi cgalVor(i_pts.size());
		for(int i=0; i<i_pts.size(); i++)
		{
			cgalVor.add(i_pts[i].pos()[0], i_pts[i].pos()[1], 1.0);
		}
		cgalVor.get(polygons);
	}

	for(int i=0; i<i_pts.size(); i++)
	{
		VoronoiCell<POS, VAL> cell(i_pts[i].pos());
		for(int j=0; j<polygons[i].size(); j++)
		{
			cell.addVertex(
				Vector<2, POS>(
					(POS) polygons[i][j][0],
					(POS) polygons[i][j][1]));
		}
		
		cell.close();
		o_polygons.push_back(cell);
	}
}

}

#endif
