#define CGAL_NO_PRECONDITIONS
#define CGAL_NO_POSTCONDITIONS
#define CGAL_NO_ASSERTIONS

#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <stk/delaunay.hpp>

namespace CGAL
{
	typedef Exact_predicates_inexact_constructions_kernel K;
	typedef Triangulation_euclidean_traits_2<K> Traits;
	
	typedef Triangulation_face_base_2<Traits> Fb;
	typedef Triangulation_vertex_base_2<Traits> VbOld;
	typedef Triangulation_vertex_base_with_info_2<int,Traits,VbOld> Vb;
	typedef Triangulation_data_structure_2<Vb,Fb> Tds;
	
	typedef Delaunay_triangulation_2<Traits,Tds> DelaunayTri;
	typedef DelaunayTri::Point Point;
	typedef DelaunayTri::Vertex_handle VertexHandle;
	typedef DelaunayTri::Vertex_circulator VertexCirculator;
	typedef DelaunayTri::Face_circulator FaceCirculator;
}

namespace stk
{

DelaunayNeighborhood::DelaunayNeighborhood()
	: m_id(0)
{
	
}

DelaunayNeighborhood::DelaunayNeighborhood(int i_id)
	: m_id(i_id)
{
	
}

void DelaunayNeighborhood::addNeighbor(int i_nei)
{
	m_neighborhood.push_back(i_nei);
}

void DelaunayNeighborhood::id(int i)
{
	m_id = i;
}

int DelaunayNeighborhood::id() const
{
	return m_id;
}

const std::vector<int>& DelaunayNeighborhood::neighborhood() const
{
	return m_neighborhood;
}

DelaunayTriangulator::DelaunayTriangulator()
{
	m_tri = (void*) new CGAL::DelaunayTri;
	m_vertices = (void*) new std::vector<CGAL::VertexHandle>;
}
	
DelaunayTriangulator::~DelaunayTriangulator()
{
	delete (CGAL::DelaunayTri*) m_tri;
	delete (std::vector<CGAL::VertexHandle>*) m_vertices;
}

void DelaunayTriangulator::add(double i_x, double i_y, int i_ptsId)
{
	CGAL::DelaunayTri* tri = (CGAL::DelaunayTri*) m_tri;
	std::vector<CGAL::VertexHandle>* vertices = (std::vector<CGAL::VertexHandle>*) m_vertices;
	
	CGAL::VertexHandle vertHand = tri->insert(CGAL::Point(i_x, i_y));
	vertHand->info() = i_ptsId;
	vertices->push_back(vertHand);
}

bool DelaunayTriangulator::getNeighborhood(int i_index, DelaunayNeighborhood& i_nei)
{
	CGAL::Traits::Segment_2 s;

	CGAL::DelaunayTri* tri = (CGAL::DelaunayTri*) m_tri;
	std::vector<CGAL::VertexHandle>* vertices = (std::vector<CGAL::VertexHandle>*) m_vertices;
	
	if(tri->is_edge(vertices->at(i_index), tri->infinite_vertex())) return false;
	
	i_nei.id(vertices->at(i_index)->info());
	
	//List all incident vertices in delaunay triangulation
	{
		CGAL::VertexCirculator ec, start;
		ec = start = tri->incident_vertices(vertices->at(i_index));
		do
		{
			i_nei.addNeighbor(ec->info());
		}
		while ( ++ec != start );
	}

	return true;
}

bool DelaunayTriangulator::getAll(int i_index, std::vector<int>& neighborhood, std::vector<stk::Vector2d>& polygon)
{
	CGAL::Traits::Segment_2 s;

	CGAL::DelaunayTri* tri = (CGAL::DelaunayTri*) m_tri;
	std::vector<CGAL::VertexHandle>* vertices = (std::vector<CGAL::VertexHandle>*) m_vertices;
	
	if(tri->is_edge(vertices->at(i_index), tri->infinite_vertex())) return false;
	
	//List all incident vertices in delaunay triangulation
	{
		CGAL::VertexCirculator ec, start;
		ec = start = tri->incident_vertices(vertices->at(i_index));
		do
		{
			neighborhood.push_back(ec->info());
		}
		while ( ++ec != start );
	}
	
	//List all incident faces in delaunay triangulation
	{
		CGAL::FaceCirculator fc, start;
		fc = start = tri->incident_faces(vertices->at(i_index));
		if(fc != 0)
		{
			do
			{
				if(!tri->is_infinite(fc))
				{
					//The dual is a point in the voronoi polygon
					CGAL::Point currVert = tri->dual(fc);
					polygon.push_back(Vector2d(currVert.x(), currVert.y()));
				}
			}
			while (++fc != start);
		}
	}

	return true;
}

}
