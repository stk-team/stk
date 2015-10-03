#define CGAL_NO_PRECONDITIONS
#define CGAL_NO_POSTCONDITIONS
#define CGAL_NO_ASSERTIONS

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>

#include <stk/voronoi.hpp>

namespace CGAL
{
	typedef Exact_predicates_inexact_constructions_kernel K;
	typedef Regular_triangulation_filtered_traits_2<K> Traits;
	typedef Regular_triangulation_face_base_2<Traits> Fb;
	typedef Regular_triangulation_vertex_base_2<Traits> VbOld;
	typedef Triangulation_vertex_base_with_info_2<int,Traits,VbOld> Vb;
	typedef Triangulation_data_structure_2<Vb,Fb> Tds;
	
	typedef Regular_triangulation_2<Traits,Tds> Triangulation;
	typedef Triangulation::Vertex_handle VertexHandle;
	typedef Triangulation::Weighted_point Point;
	typedef Triangulation::Face_circulator FaceCirculator;
}

namespace stk
{

namespace doNotUse
{

_cgalWeightedVoronoi::_cgalWeightedVoronoi(int i_size)
{
	CGAL::Triangulation* tri = new CGAL::Triangulation;
	CGAL::VertexHandle* vertices = new CGAL::VertexHandle[i_size];

	m_tri = (void*) tri;
	m_vertices = (void*) vertices;
	m_counter = 0;
}
	
_cgalWeightedVoronoi::~_cgalWeightedVoronoi()
{
	delete[] ((CGAL::VertexHandle*) m_vertices);
	delete ((CGAL::Triangulation*) m_tri);
}

void _cgalWeightedVoronoi::add(double i_x, double i_y, double i_w)
{
	i_w = std::sqrt(i_w/M_PI);
	
	CGAL::Triangulation* tri = (CGAL::Triangulation*) m_tri;
	CGAL::VertexHandle* vertices = (CGAL::VertexHandle*) m_vertices;
	
	vertices[m_counter] = tri->insert(CGAL::Point(CGAL::Point(i_x, i_y), i_w));
	vertices[m_counter]->info() = m_counter;

	m_counter++;
}

void _cgalWeightedVoronoi::get(std::vector<std::vector<Vector2d> >& i_polygons)
{
	CGAL::Triangulation* tri = (CGAL::Triangulation*) m_tri;
	CGAL::VertexHandle* vertices = (CGAL::VertexHandle*) m_vertices;

	CGAL::Triangulation::Finite_vertices_iterator iter;
	CGAL::Triangulation::Finite_vertices_iterator iterEnd;
	iter = tri->finite_vertices_begin();
	iterEnd = tri->finite_vertices_end();

	while(iter != iterEnd)
	{
		CGAL::VertexHandle vh(iter);
		int id = std::abs(vh->info());

		if(!tri->is_infinite(iter))
		{
			CGAL::FaceCirculator fc = tri->incident_faces(iter);
			CGAL::FaceCirculator fend = fc;
			if(fc != 0)
			{
				do
				{
					if(!tri->is_infinite(fc))
					{
						CGAL::Triangulation::Point currVert = tri->dual(fc);
						Vector2d pt = Vector2d(currVert.x(), currVert.y());
						
						i_polygons[id].push_back(pt);
					}
					++fc;
				}
				while(fc != fend);
			}
		}
		
		++iter;
	}
}

}

}
