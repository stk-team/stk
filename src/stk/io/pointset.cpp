#include <map>

#include <stk/io/pointset.hpp>
#include <stk/color.hpp>
#include <stk/type.hpp>
#include <stk/voronoi.hpp>
#include <stk/delaunay.hpp>

#ifdef CAIRO_ENABLED
#include <cairo/cairo.h>
#endif

namespace stk
{

namespace io
{

#ifdef CAIRO_ENABLED

template<typename VAL>
PointSetGraphics<VAL>::DomainToImage::DomainToImage(const Vector2d& size)
{
	imgSize = size;
}

template<typename VAL>
PointSetGraphics<VAL>::DomainToImageGlobal::DomainToImageGlobal(const Vector2d& size, const Domain<2, double>* domain, int bSize)
	: DomainToImage(size)
{
	Vector2d boxMin = domain->boundingBoxMin();
	Vector2d boxMax = domain->boundingBoxMax();
	Vector2d boxSize = domain->boundingBoxSize();
	
	double boxRatio = (double)boxSize[0]/(double)boxSize[1];
	
	double borderWidth;
	if(this->imgSize[0] < this->imgSize[1])
	{
		borderWidth = boxSize[0]/(((double)this->imgSize[0]/bSize)-2);
	}
	else
	{
		borderWidth = boxSize[1]/(((double)this->imgSize[1]/bSize)-2);
	}
	
	Vector2d borderSize = Vector2d(boxRatio, 1.0)*borderWidth;
	viewMin = boxMin - borderSize;
	viewSize = boxSize + borderSize*2.0;
}

template<typename VAL>
Vector2d PointSetGraphics<VAL>::DomainToImageGlobal::convert(const Vector2d& p) const
{
	return ((p-viewMin)/viewSize)*this->imgSize;
}

template<typename VAL>
PointSetGraphics<VAL>::DomainToImageZoom::DomainToImageZoom(const Vector2d& size, const Domain<2, double>* domain, int bSize, const Vector2d& c, double s)
	: DomainToImage(size)
{
	Vector2d boxMin = domain->boundingBoxMin();
	Vector2d boxMax = domain->boundingBoxMax();
	Vector2d boxSize = domain->boundingBoxSize();
	
	double boxRatio = (double)boxSize[0]/(double)boxSize[1];
	
	double borderWidth;
	if(this->imgSize[0] < this->imgSize[1])
	{
		borderWidth = boxSize[0]/(((double)this->imgSize[0]/bSize)-2);
	}
	else
	{
		borderWidth = boxSize[1]/(((double)this->imgSize[1]/bSize)-2);
	}
	
	Vector2d borderSize = Vector2d(boxRatio, 1.0)*borderWidth;
	viewMin = boxMin - borderSize;
	viewSize = boxSize + borderSize*2.0;
	
	center = c;
	scale = s;
}

template<typename VAL>
Vector2d PointSetGraphics<VAL>::DomainToImageZoom::convert(const Vector2d& p) const
{
	return ((p-center)/viewSize)*scale*this->imgSize + this->imgSize/2;
}

template<typename VAL>
void PointSetGraphics<VAL>::draw(const std::string& i_filename, const Vector2i& i_size) const
{
	cairo_surface_t* surface;
	cairo_t* context;
	
	surface = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32,
		i_size[0],
		i_size[1]);
	context = cairo_create(surface);
	
	//Background
	if(m_pointset->isToroidal())
	{
		//~ cairo_set_source_rgba(context, 0.95, 0.95, 0.95, 1.0);
		cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
	}
	else
	{
		cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
	}
	cairo_paint(context);
	
	//Define space transformation (domain to image)
	Vector2d imgSize = Vector2d(i_size[0], i_size[1]);
	DomainToImage* domainToImage;
	if(m_zoom) domainToImage = new DomainToImageZoom(imgSize, m_pointset->domain(), m_borderSize, m_zoomCenter, m_zoomScale);
	else domainToImage = new DomainToImageGlobal(imgSize, m_pointset->domain(), m_borderSize);
	
	//Compute delaunay triangulation
	std::vector< stk::DelaunayElement<double> > delaunayTri;
	if(m_computeDelaunay) stk::delaunay(*m_pointset, delaunayTri);
	
	if(m_voronoi)
	{
		double vorAreaMean = m_pointset->domain()->volume() / m_pointset->size();
				
		//Precompute class color
		std::map<VAL, double> vorClassColor;
		if(m_voronoiColorType == VOR_CLASS_COLOR)
		{
			for(int i=0; i<m_pointset->size(); i++)
			{
				vorClassColor[m_pointset->at(i).val()] = 0;
			}
			
			for(
				typename std::map<VAL, double>::iterator iter = vorClassColor.begin();
				iter != vorClassColor.end();
				iter++)
			{
				iter->second = drand48();
			}
		}
		
		if(m_voronoiColorType != VOR_CLASS_BOUND)
		{
			for(int i=0; i<delaunayTri.size(); i++)
			{			
				const std::vector<Vector2d>& polygon = delaunayTri[i].polygon();
				
				if(polygon.size() <= 2) continue;
				
				//Avoid bad cells (ex: 14084)
				if(!m_pointset->isToroidal())
				{
					bool badCell = false;
					for(int j=0; j<polygon.size(); j++)
					{
						stk::Vector2d p = polygon[j];
						if(!m_pointset->domain()->hitTest(p))
						{
							badCell = true;
							break;
						}
					}
					if(badCell) continue;
				}		
				
				Vector2d center = domainToImage->convert(polygon[polygon.size()-1]);
				cairo_move_to(context, center[0], center[1]);
				
				for(int j=0; j<polygon.size(); j++)
				{
					center = domainToImage->convert(polygon[j]);
					cairo_line_to(context, center[0], center[1]);
				}
				
				if(m_voronoiColorType == VOR_CONNEXITY)
				{
					if(polygon.size() == 5)
					{
						cairo_set_source_rgba(context, 211.0/255.0, 1.0, 1.0, 1.0);
					}
					else if(polygon.size() < 5)
					{
						cairo_set_source_rgba(context, 178.0/255.0, 223.0/255.0, 1.0, 1.0);
					}
					else if(polygon.size() == 7)
					{
						cairo_set_source_rgba(context, 1.0, 249.0/255.0, 203.0/255.0, 1.0);
					}
					else if(polygon.size() > 7)
					{
						cairo_set_source_rgba(context, 1.0, 232.0/255.0, 187.0/255.0, 1.0);
					}
					else cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
					
					cairo_fill_preserve(context);
					cairo_set_source_rgba(context, 0.5, 0.5, 0.5, 1.0);
					cairo_set_line_width (context, 0.6);
					cairo_stroke(context);
				}
				else if(m_voronoiColorType == VOR_AREA)
				{
					double area = std::min(std::max(delaunayTri[i].area()/vorAreaMean, 0.5), 1.5)-0.5;
					Vector3d cHsv = color::hsv2rgb(Vector3d(area/2.0, 0.6, 1.0));
					
					cairo_set_source_rgba(context, cHsv[0], cHsv[1], cHsv[2], 1.0);
					cairo_fill_preserve(context);
					cairo_set_source_rgba(context, 0.5, 0.5, 0.5, 1.0);
					cairo_set_line_width (context, 0.6);
					cairo_stroke(context);
				}
				else if(m_voronoiColorType == VOR_CLASS_COLOR)
				{
					int idClass = m_pointset->at(delaunayTri[i].id()).val();
					Vector3d cHsv = color::hsv2rgb(Vector3d(vorClassColor[idClass], 0.3, 1.0));
					
					cairo_set_source_rgba(context, cHsv[0], cHsv[1], cHsv[2], 1.0);
					cairo_fill_preserve(context);
					cairo_set_line_width (context, 0.6);
					cairo_stroke(context);
				}
				else
				{
					cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
					cairo_fill_preserve(context);
					cairo_set_source_rgba(context, 0.5, 0.5, 0.5, 1.0);
					cairo_set_line_width (context, 0.6);
					cairo_stroke(context);
				}
			}
		}
		
		if(m_voronoiColorType == VOR_CLASS_COLOR || m_voronoiColorType == VOR_CLASS_BOUND)
		{
			for(int i=0; i<delaunayTri.size(); i++)
			{
				const std::vector<Vector2d>& polygon = delaunayTri[i].polygon();
				
				if(polygon.size() <= 2) continue;
				if(!m_pointset->isToroidal())
				{
					bool badCell = false;
					for(int j=0; j<polygon.size(); j++)
					{
						stk::Vector2d p = polygon[j];
						if(!m_pointset->domain()->hitTest(p))
						{
							badCell = true;
							break;
						}
					}
					if(badCell) continue;
				}
				
				int idClassCenter = m_pointset->at(delaunayTri[i].id()).val();
				const std::vector<int>& neighborhood = delaunayTri[i].neighborhood();
				int sz = polygon.size();
				if(m_voronoiColorType == VOR_CLASS_COLOR)
				{
					cairo_set_source_rgba(context, 0.2, 0.2, 0.2, 1.0);
				}
				else
				{
					cairo_set_source_rgba(context, 0.5, 0.5, 0.5, 1.0);
				}
				for(int j=0; j<sz; j++)
				{
					int idClassNei = m_pointset->at(neighborhood[j]).val();
					if(idClassNei != idClassCenter)
					{
						Vector2d center = domainToImage->convert(polygon[(sz+j-1)%sz]);
						cairo_move_to(context, center[0], center[1]);
						
						cairo_set_line_width (context, 0.7);
						center = domainToImage->convert(polygon[(sz+j)%sz]);
						cairo_line_to(context, center[0], center[1]);
						
						cairo_stroke(context);
					}
				}
			}
		}
	}
	
	//Draw timeTrace
	{
		if(m_timeTrace.size() >= 2)
		{
			cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			for(int i=0; i<m_pointset->size(); i++)
			{
				bool first = true;
				typename std::map<float, const PointSet<2, double, VAL>*>::const_iterator iter;
				for(iter = m_timeTrace.begin(); iter != m_timeTrace.end(); iter++)
				{
					Vector2d refPos = m_pointset->at(i).pos();
					Vector2d curPos = iter->second->at(i).pos();
					Vector2d pos =  refPos + m_pointset->diff(refPos, curPos);
					Vector2d center = domainToImage->convert(pos);
					
					if(first)
					{
						cairo_move_to(context, 0.5+center[0], 0.5+center[1]);
						first = false;
					}
					else cairo_line_to(context, 0.5+center[0], 0.5+center[1]);
				}
					
				cairo_stroke(context);
			}
			
			if(m_timeTracePoint)
			{
				for(int i=0; i<m_pointset->size(); i++)
				{
					typename std::map<unsigned int, PointSetGraphicsFlag>::const_iterator iter;
					iter = m_flags.find(i);
					PointSetGraphicsFlag flag = (iter != m_flags.end() ? iter->second : Default);
					
					if(flag != Fixed)
					{
						typename std::map<float, const PointSet<2, double, VAL>*>::const_iterator iter;
						for(iter = m_timeTrace.begin(); iter != m_timeTrace.end(); iter++)
						{
							Vector2d refPos = m_pointset->at(i).pos();
							Vector2d curPos = iter->second->at(i).pos();
							Vector2d pos =  refPos + m_pointset->diff(refPos, curPos);
							Vector2d center = domainToImage->convert(pos);
							
							cairo_set_source_rgba(context, 172.0/255.0, 10.0/255.0, 128.0/255.0, 1.0);
							cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
							cairo_fill(context);
						}
					}
				}
			}
		}
	}
	
	//Draw delaunay
	if(m_delaunay)
	{
		cairo_set_source_rgba(context, 6.0/255.0, 46.0/255.0, 108.0/255.0, 1.0);
		
		for(int i=0; i<delaunayTri.size(); i++)
		{
			const std::vector<int>& nei = delaunayTri[i].neighborhood();
			int ptId = delaunayTri[i].id();
			
			Vector2d ptPos = domainToImage->convert(m_pointset->at(ptId).pos());
			
			for(int j=0; j<nei.size(); j++)
			{
				int neiId = nei[j];
				
				Vector2d neiPos = m_pointset->at(ptId).pos() + m_pointset->diff(m_pointset->at(ptId).pos(), m_pointset->at(neiId).pos());
				neiPos = domainToImage->convert(neiPos);
				
				cairo_move_to(context, 0.5+ptPos[0], 0.5+ptPos[1]);
				cairo_line_to(context, 0.5+neiPos[0], 0.5+neiPos[1]);
				cairo_stroke(context);
			}
		}
	}
	
	//Draw replication
	if(m_pointset->isToroidal())
	{
		cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);
		
		typename PointSet<2, double, VAL>::replicationIterator iter;
		for(iter = m_pointset->replicationBegin(false); iter != m_pointset->replicationEnd(); iter++)
		{
			Vector2d center = domainToImage->convert(iter.pos());
			
			cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
			cairo_fill(context);
		}
	}
	
	//Draw points
	{
		for(int i=0; i<m_pointset->size(); i++)
		{
			Vector2d center = domainToImage->convert(m_pointset->at(i).pos());
			
			typename std::map<unsigned int, PointSetGraphicsFlag>::const_iterator iter;
			iter = m_flags.find(i);
			PointSetGraphicsFlag flag = (iter != m_flags.end() ? iter->second : Default);
			if(flag == Default)
			{
				if(!m_pointset->domain()->hitTest(m_pointset->at(i).pos()))
				{
					flag = Wrong;
				}
			}
			
			switch(flag)
			{
				case Wrong:
					cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 1.0);
					cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
					cairo_fill(context);
					break;
				case Hot:
					cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 1.0);
					cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
					cairo_fill(context);
					
					cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 0.6);
					cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
					cairo_stroke(context);
					
					cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 0.3);
					cairo_arc(context, center[0], center[1], 3*m_pointRadius, 0, 2*M_PI);
					cairo_stroke(context);
					break;
				case Cold:
					cairo_set_source_rgba(context, 0.129411765, 0.462745098, 0.792156863, 1.0);
					cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
					cairo_fill(context);
					
					cairo_set_source_rgba(context, 0.129411765, 0.462745098, 0.792156863, 0.6);
					cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
					cairo_stroke(context);
					
					cairo_set_source_rgba(context, 0.129411765, 0.462745098, 0.792156863, 0.3);
					cairo_arc(context, center[0], center[1], 3*m_pointRadius, 0, 2*M_PI);
					cairo_stroke(context);
					break;
				case Fixed:	
					cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
					cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
					cairo_stroke(context);
					break;
				case Default:				
				default:
					cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
					cairo_arc(context, center[0], center[1], 2*m_pointRadius, 0, 2*M_PI);
					cairo_fill(context);
			};
		}
	}
	
	//Draw boundaries
	{
		const BaseDomain<double>* baseDomain = dynamic_cast<const BaseDomain<double>*>(m_pointset->domain());
		const HexagonalDomain<double>* hexaDomain = dynamic_cast<const HexagonalDomain<double>*>(m_pointset->domain());
		
		if(baseDomain != NULL)
		{
			Vector2d v0 = domainToImage->convert(baseDomain->referencePoint());
			Vector2d v1 = domainToImage->convert(baseDomain->referencePoint() + baseDomain->vector(0));
			Vector2d v2 = domainToImage->convert(baseDomain->referencePoint() + baseDomain->vector(0) + baseDomain->vector(1));
			Vector2d v3 = domainToImage->convert(baseDomain->referencePoint() + baseDomain->vector(1));
			
			cairo_set_line_width (context, 1.0);
			cairo_move_to(context, 0.5+v0[0], 0.5+v0[1]);
			cairo_line_to(context, 0.5+v1[0], 0.5+v1[1]);
			cairo_line_to(context, 0.5+v2[0], 0.5+v2[1]);
			cairo_line_to(context, 0.5+v3[0], 0.5+v3[1]);
			cairo_line_to(context, 0.5+v0[0], 0.5+v0[1]);
			
			//~ if(m_pointset->isToroidal() && !m_voronoi)
			//~ {
				//~ cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
				//~ cairo_fill_preserve(context);
			//~ }
			cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 0.5);
			cairo_stroke(context);
		}
		else if(hexaDomain != NULL)
		{
			Vector2d center = hexaDomain->referencePoint();
			double r = hexaDomain->incircleRadius();
			double R = hexaDomain->circumcircleRadius();
			
			Vector2d v0 = center+Vector2d(r, R/2);
			Vector2d v1 = center+Vector2d(0, R);
			Vector2d v2 = center+Vector2d(-r, R/2);
			Vector2d v3 = center+Vector2d(-r, -R/2);
			Vector2d v4 = center+Vector2d(0, -R);
			Vector2d v5 = center+Vector2d(r, -R/2);
			
			cairo_set_line_width (context, 1.0);
			cairo_move_to(context, 0.5+v0[0], 0.5+v0[1]);
			cairo_line_to(context, 0.5+v1[0], 0.5+v1[1]);
			cairo_line_to(context, 0.5+v2[0], 0.5+v2[1]);
			cairo_line_to(context, 0.5+v3[0], 0.5+v3[1]);
			cairo_line_to(context, 0.5+v4[0], 0.5+v4[1]);
			cairo_line_to(context, 0.5+v5[0], 0.5+v5[1]);
			cairo_line_to(context, 0.5+v0[0], 0.5+v0[1]);
			
			if(m_pointset->isToroidal() && !m_voronoi)
			{
				cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
				cairo_fill_preserve(context);
			}
			cairo_set_source_rgba(context, 0.7, 0.7, 0.7, 1.0);
			cairo_stroke(context);
		}
		else
		{
			Vector2d v0 = domainToImage->convert(Vector2d(m_pointset->boundingBoxMin()[0], m_pointset->boundingBoxMin()[1]));
			Vector2d v1 = domainToImage->convert(Vector2d(m_pointset->boundingBoxMax()[0], m_pointset->boundingBoxMin()[1]));
			Vector2d v2 = domainToImage->convert(Vector2d(m_pointset->boundingBoxMax()[0], m_pointset->boundingBoxMax()[1]));
			Vector2d v3 = domainToImage->convert(Vector2d(m_pointset->boundingBoxMin()[0], m_pointset->boundingBoxMax()[1]));
			
			cairo_set_line_width (context, 1.0);
			cairo_move_to(context, 0.5+v0[0], 0.5+v0[1]);
			cairo_line_to(context, 0.5+v1[0], 0.5+v1[1]);
			cairo_line_to(context, 0.5+v2[0], 0.5+v2[1]);
			cairo_line_to(context, 0.5+v3[0], 0.5+v3[1]);
			cairo_line_to(context, 0.5+v0[0], 0.5+v0[1]);
			
			//~ if(m_pointset->isToroidal() && !m_voronoi)
			//~ {
				//~ cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
				//~ cairo_fill_preserve(context);
			//~ }
			cairo_set_source_rgba(context, 0.7, 0.7, 0.7, 1.0);
			cairo_stroke(context);
		}
	}
	
	if(m_labelType != None)
	{
		cairo_select_font_face(context, "Ubuntu",
			CAIRO_FONT_SLANT_NORMAL,
			CAIRO_FONT_WEIGHT_NORMAL);
		cairo_set_font_size(context, 10);
		
		cairo_text_extents_t extents;
		
		for(int i=0; i<m_pointset->size(); i++)
		{
			Vector2d center = domainToImage->convert(m_pointset->at(i).pos());
			
			std::stringstream text;
			switch(m_labelType)
			{
				case Value:
					text << m_pointset->at(i).val();
					break;
				default:
					text << i;
			}
			cairo_text_extents(context, text.str().c_str(), &extents);
			
			//Flag
			typename std::map<unsigned int, PointSetGraphicsFlag>::const_iterator iter;
			iter = m_flags.find(i);
			PointSetGraphicsFlag flag = (iter != m_flags.end() ? iter->second : Default);
			if(flag == Default)
			{
				if(!m_pointset->domain()->hitTest(m_pointset->at(i).pos()))
				{
					flag = Wrong;
				}
			}
			
			if(flag != Default)
			{
				cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
				cairo_rectangle(context, center[0]-extents.width/2-3, center[1]-extents.height/2-3, extents.width+6, extents.height+6);
				cairo_fill_preserve(context);
				switch(flag)
				{
					case Wrong:
					case Hot:
						cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 0.75);
						break;
					case Cold:
						cairo_set_source_rgba(context, 0.129411765, 0.462745098, 0.792156863, 0.75);
						break;
					case Fixed:				
					default:
						cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 0.75);
				};
				cairo_rectangle(context, center[0]-extents.width/2-3, center[1]-extents.height/2-3, extents.width+6, extents.height+6);
				cairo_stroke(context);
			}
			else
			{
				cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 0.75);
				cairo_rectangle(context, center[0]-extents.width/2-3, center[1]-extents.height/2-3, extents.width+6, extents.height+6);
				cairo_fill(context);
			}
			
			switch(flag)
			{
				case Wrong:
				case Hot:
					cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 1.0);
					break;
				case Cold:
					cairo_set_source_rgba(context, 0.129411765, 0.462745098, 0.792156863, 1.0);
					break;
				case Default:				
				default:
					cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
			};
			
			cairo_move_to(context, center[0], center[1]);
			cairo_rel_move_to(context, -extents.width/2, extents.height/2);			
			cairo_show_text(context, text.str().c_str());
		}
	}
	
	cairo_status_t status;
	status = cairo_surface_write_to_png(surface, i_filename.c_str());
	if(status != CAIRO_STATUS_SUCCESS)
	{
		throw stk::exception::CairoError(status, STK_DBG_INFO)
			.addVar("filename", i_filename)
			.addVar("size", i_size)
			.addVar("pointset", *m_pointset);
	}
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
	
	delete domainToImage;
}

template class PointSetGraphics<double>;
template class PointSetGraphics<int>;
template class PointSetGraphics<unsigned int>;

void draw(
	std::string i_filename,
	const PointSet2dd& i_pts,
	const Vector2i& i_size)
{
	PointSetGraphics<double> graph;
	graph.pointset(i_pts);
	graph.draw(i_filename, i_size);
}

void drawWeight(
	std::string i_filename,
	const PointSet2dd& i_pts,
	const Vector2i& i_size)
{
	cairo_surface_t* surface;
	cairo_t* context;
	
	surface = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32,
		i_size[0],
		i_size[1]);
	context = cairo_create(surface);

	//Background
	cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
	cairo_paint(context);

	Vector2d minPos = i_pts.boundingBoxMin();
	Vector2d maxPos = i_pts.boundingBoxMax();
	Vector2d spaceSize = i_pts.boundingBoxSize();

	if(i_pts.isToroidal())
	{
		double meanDist = 5.0*i_pts.meanDistance();
		spaceSize += Vector2d(meanDist, meanDist) * 2.0;
		minPos -= Vector2d(meanDist, meanDist);
		maxPos += Vector2d(meanDist, meanDist);

		cairo_set_line_width (context, 1.0);
		cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);
		cairo_move_to(context,
			0.5+i_size[0]*(meanDist+i_pts.boundingBoxMin()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_pts.boundingBoxMin()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_pts.boundingBoxMax()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_pts.boundingBoxMin()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_pts.boundingBoxMax()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_pts.boundingBoxMax()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_pts.boundingBoxMin()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_pts.boundingBoxMax()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_pts.boundingBoxMin()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_pts.boundingBoxMin()[1])/spaceSize[1]);
		cairo_stroke(context);
	
		for(int i=0; i<i_pts.size(); i++)
		{
			double ptRadius = 4.0*std::sqrt(i_pts[i].val());
			
			for(int j=-1; j<=1; j++)
			{
				for(int k=-1; k<=1; k++)
				{
					if(k == 0 && j == 0)
					{
						cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
					}
					else cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);

					Vector2d shift(i_pts.boundingBoxSize()[0]*j, i_pts.boundingBoxSize()[1]*k);
					
					Vector2d center = (i_pts[i].pos() - minPos + shift) / spaceSize;
					center *= Vector2d(i_size[0], i_size[1]);
					
					cairo_arc(context, center[0], center[1], ptRadius, 0, 2*M_PI);
					cairo_fill(context);
				}
			}
		}
	}
	else
	{
		cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
	
		for(int i=0; i<i_pts.size(); i++)
		{
			double ptRadius = 4.0*std::sqrt(i_pts[i].val());
			
			Vector2d center = (i_pts[i].pos() - minPos) / spaceSize;
			center *= Vector2d(i_size[0], i_size[1]);
			
			cairo_arc(context, center[0], center[1], ptRadius, 0, 2*M_PI);
			cairo_fill(context);
		}
	}
	
	cairo_surface_write_to_png(surface, i_filename.c_str());
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
}

void draw(
	std::string i_filename,
	const PointSet2di& i_pts,
	const Vector2i& i_size)
{
	cairo_surface_t* surface;
	cairo_t* context;

	int idMax = 0;
	for(int i=0; i<i_pts.size(); i++)
	{
		if(idMax < i_pts[i].val()) idMax = i_pts[i].val();
	}
	
	std::map<int, Vector3d> colorScheme;
	for(int i=0; i<i_pts.size(); i++)
	{
		Vector3d c = color::hsv2rgb(Vector3d((double)i_pts[i].val() / (double) idMax, 0.5, 1.0));
		colorScheme[i_pts[i].val()] = c;
	}
	
	surface = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32,
		i_size[0],
		i_size[1]);
	context = cairo_create(surface);

	//Background
	cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
	cairo_paint(context);
	
	for(int i=0; i<i_pts.size(); i++)
	{
		Vector3d mColor = colorScheme[i_pts[i].val()];
		cairo_set_source_rgba(context, mColor[0], mColor[1], mColor[2], 1.0);
		
		Vector2d center = (i_pts[i].pos() - i_pts.boundingBoxMin()) / i_pts.boundingBoxSize();
		center *= Vector2d(i_size[0], i_size[1]);
		
		cairo_arc(context, center[0], center[1], 2, 0, 2*M_PI);
		cairo_fill(context);
	}
	
	cairo_surface_write_to_png(surface, i_filename.c_str());
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
}

void drawMotion(
	std::string i_filename,
	const PointSet2dd& i_before,
	const PointSet2dd& i_after,
	const Vector2i& i_size)
{
	cairo_surface_t* surface;
	cairo_t* context;
	
	surface = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32,
		i_size[0],
		i_size[1]);
	context = cairo_create(surface);

	Vector2d minPos = i_after.boundingBoxMin();
	Vector2d maxPos = i_after.boundingBoxMax();
	Vector2d spaceSize = i_after.boundingBoxSize();

	//Background
	cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
	cairo_paint(context);
	
	if(i_after.isToroidal())
	{
		double meanDist = 5.0/std::sqrt(i_after.size());
		spaceSize += Vector2d(meanDist, meanDist) * 2.0;
		minPos -= Vector2d(meanDist, meanDist);
		maxPos += Vector2d(meanDist, meanDist);

		cairo_set_line_width (context, 1.0);
		cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);
		cairo_move_to(context,
			0.5+i_size[0]*(meanDist+i_after.boundingBoxMin()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_after.boundingBoxMin()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_after.boundingBoxMax()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_after.boundingBoxMin()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_after.boundingBoxMax()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_after.boundingBoxMax()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_after.boundingBoxMin()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_after.boundingBoxMax()[1])/spaceSize[1]);
		cairo_line_to(context,
			0.5+i_size[0]*(meanDist+i_after.boundingBoxMin()[0])/spaceSize[0],
			0.5+i_size[1]*(meanDist+i_after.boundingBoxMin()[1])/spaceSize[1]);
		cairo_stroke(context);
	}

	//Before
	//~ cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 1.0);
	//~ for(int i=0; i<i_before.size(); i++)
	//~ {
		//~ Vector2d center = (i_before[i].pos() - minPos) / spaceSize;
		//~ center *= Vector2d(i_size[0], i_size[1]);
		//~ 
		//~ cairo_arc(context, center[0], center[1], 2, 0, 2*M_PI);
		//~ cairo_fill(context);
	//~ }

	//Line
	cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 1.0);
	cairo_set_line_width (context, 1.0);
	for(int i=0; i<i_before.size(); i++)
	{
		Vector2d dist = i_after.diff(i_after[i].pos(), i_before[i].pos());
		Vector2d center0 = (i_after[i].pos() + dist - minPos) / spaceSize;
		Vector2d center1 = (i_after[i].pos() - minPos) / spaceSize;

		center0 *= Vector2d(i_size[0], i_size[1]);
		center1 *= Vector2d(i_size[0], i_size[1]);

		cairo_move_to(context, center0[0], center0[1]);
		cairo_line_to(context, center1[0], center1[1]);
		cairo_stroke(context);
	}
	
	//After
	if(i_after.isToroidal())
	{
		for(int i=0; i<i_after.size(); i++)
		{
			for(int j=-1; j<=1; j++)
			{
				for(int k=-1; k<=1; k++)
				{
					if(k == 0 && j == 0)
					{
						cairo_set_source_rgba(context, 0.0, 1.0, 0.0, 1.0);
					}
					else cairo_set_source_rgba(context, 0.8, 0.8, 0.8, 1.0);

					Vector2d shift(i_after.boundingBoxSize()[0]*j, i_after.boundingBoxSize()[1]*k);
					
					Vector2d center = (i_after[i].pos() - minPos + shift) / spaceSize;
					center *= Vector2d(i_size[0], i_size[1]);
					
					cairo_arc(context, center[0], center[1], 2, 0, 2*M_PI);
					cairo_fill(context);
				}
			}
		}
	}
	else
	{
		cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
	
		for(int i=0; i<i_after.size(); i++)
		{
			Vector2d center = (i_after[i].pos() - minPos) / spaceSize;
			center *= Vector2d(i_size[0], i_size[1]);
			
			cairo_arc(context, center[0], center[1], 2, 0, 2*M_PI);
			cairo_fill(context);
		}
	}
	
	cairo_surface_write_to_png(surface, i_filename.c_str());
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
}

#endif

}

}
