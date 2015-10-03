#ifndef __STK_DELAUNAY__
#define __STK_DELAUNAY__

#include <vector>

#include <stk/pointset.hpp>

namespace stk
{

/* DELAUNAY : COMMONS *************************************************/

class DelaunayNeighborhood;

class DelaunayTriangulator
{
	public:
		void* m_tri;
		void* m_vertices;
		
	public:
		DelaunayTriangulator();
		~DelaunayTriangulator();
		void add(double i_x, double i_y, int i_id);
		bool getNeighborhood(int i_index, DelaunayNeighborhood& i_nei);
		bool getAll(int i_index, std::vector<int>& neighborhood, std::vector<stk::Vector2d>& polygon);
};

/* DELAUNAY : NEIGHBORHOOD ********************************************/

class DelaunayNeighborhood
{
	private:
		int m_id;
		std::vector<int> m_neighborhood;
	
	public:
		DelaunayNeighborhood();
		DelaunayNeighborhood(int i_id);
		
		void id(int i);
		void addNeighbor(int i_nei);
		
		int id() const;
		const std::vector<int>& neighborhood() const;
};

template<typename POS, typename VAL>
void delaunay(
	const PointSet<2, POS, VAL>& i_pts,
	std::vector<DelaunayNeighborhood>& o_result)
{
	DelaunayTriangulator delTri;
		
	if(i_pts.isToroidal())
	{
		typename PointSet<2, POS, VAL>::replicationIterator iter;
		for(iter = i_pts.replicationBegin(true, 7.0); iter != i_pts.replicationEnd(); iter++)
		{
			delTri.add(iter.pos()[0], iter.pos()[1], iter.id());
		}
	}
	else
	{
		for(int i=0; i<i_pts.size(); i++)
		{
			delTri.add(i_pts[i].pos()[0], i_pts[i].pos()[1], i);
		}
	}
	
	
	for(int i=0; i<i_pts.size(); i++)
	{
		DelaunayNeighborhood n;
			
		if(delTri.getNeighborhood(i, n))
		{
			o_result.push_back(n);
		}
	}
}

/* DELAUNAY : NEIGHBORHOOD & VORONOI POLYGON **************************/

template<typename POS>
class DelaunayElement
{
	private:
		int m_id;
		
		std::vector<int> m_neighborhood;
		std::vector< Vector<2, POS> > m_polygon;
		Vector<2, POS> m_centroid;
		Vector<2, POS> m_pt;
		double m_area;
	
	public:
		DelaunayElement()
		{
			m_id = -1;
		}
		
		DelaunayElement(int i_id)
		{
			id(i_id);
		}
		
		void id(int i)
		{
			m_id = i;
		}
		
		int id() const
		{
			return m_id;
		}
		
		const std::vector<int>& neighborhood() const
		{
			return m_neighborhood;
		}
		
		const std::vector< Vector<2, POS> >& polygon() const
		{
			return m_polygon;
		}
		
		std::vector<int>& neighborhood()
		{
			return m_neighborhood;
		}
		
		std::vector< Vector<2, POS> >& polygon()
		{
			return m_polygon;
		}

		void vertex(const Vector<2, POS>& i_vertex)
		{
			m_pt = i_vertex;
		}

		const Vector<2, POS>& vertex() const
		{
			return m_pt;
		}
		
		const double& area() const
		{
			return m_area;
		}

		const Vector<2, POS> centroid() const
		{
			return m_centroid;
		}

		void closePolygon()
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
			
			m_area = a;
		}
};

template<typename POS, typename VAL>
void delaunay(
	const PointSet<2, POS, VAL>& i_pts,
	std::vector< DelaunayElement<POS> >& o_result)
{
	DelaunayTriangulator delTri;
		
	if(i_pts.isToroidal())
	{
		typename PointSet<2, POS, VAL>::replicationIterator iter;
		for(iter = i_pts.replicationBegin(true, 7.0); iter != i_pts.replicationEnd(); iter++)
		{
			delTri.add(iter.pos()[0], iter.pos()[1], iter.id());
		}
	}
	else
	{
		for(int i=0; i<i_pts.size(); i++)
		{
			delTri.add(i_pts[i].pos()[0], i_pts[i].pos()[1], i);
		}
	}
	
	
	for(int i=0; i<i_pts.size(); i++)
	{
		DelaunayElement<POS> delElem(i);
		std::vector<Vector2d> polygon;
			
		if(delTri.getAll(i, delElem.neighborhood(), polygon))
		{
			for(int j=0; j<polygon.size(); j++)
			{
				delElem.polygon().push_back(Vector<2, POS>(polygon[j][0], polygon[j][1]));
			}
			delElem.closePolygon();
			o_result.push_back(delElem);
		}
	}
}

}

#endif
