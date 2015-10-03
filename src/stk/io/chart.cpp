#include <stk/exception.hpp>
#include <stk/io/chart.hpp>
#include <stk/color.hpp>

#ifdef CAIRO_ENABLED
#include <cairo/cairo.h>
#endif

namespace stk
{
	
namespace io
{

#ifdef CAIRO_ENABLED

static const double dashed3[] = {3.0};

Chart::LineStyle::LineStyle() :
	m_strokeEnabled(true),
	m_defaultStrokeColor(true),
	m_fillingEnabled(false),
	m_dashEnabled(false)
{
	m_strokeWidth = 1.0;
	m_strokeColor = Vector4d(1.0, 0.0, 0.0, 1.0);
	m_fillColor = m_strokeColor;
	m_fillColor[3] = m_strokeColor[3]*0.25;
}

Chart::LineStyle& Chart::LineStyle::stroke(bool enable)
{
	m_strokeEnabled = enable;
	return *this;
}

bool Chart::LineStyle::strokeEnabled() const
{
	return m_strokeEnabled;
}

bool Chart::LineStyle::defaultStrokeColor() const
{
	return m_defaultStrokeColor;
}

Chart::LineStyle& Chart::LineStyle::strokeWidth(float width)
{
	m_strokeWidth = width;
	return *this;
}

float Chart::LineStyle::getStrokeWidth() const
{
	return m_strokeWidth;
}

Chart::LineStyle& Chart::LineStyle::strokeColor(const Vector4d& color)
{
	m_defaultStrokeColor = false;
	m_strokeColor = color;
	m_fillColor = m_strokeColor;
	m_fillColor[3] = m_strokeColor[3]*0.25;
	return *this;
}

const Vector4d& Chart::LineStyle::getStrokeColor() const
{
	return m_strokeColor;
}

Chart::LineStyle& Chart::LineStyle::filling(bool enable)
{
	m_fillingEnabled = enable;
	return *this;
}

bool Chart::LineStyle::fillingEnabled() const
{
	return m_fillingEnabled;
}

const Vector4d& Chart::LineStyle::getFillColor() const
{
	return m_fillColor;
}

Chart::LineStyle& Chart::LineStyle::dash(bool enable)
{
	m_dashEnabled = enable;
	return *this;
}

bool Chart::LineStyle::dashEnabled() const
{
	return m_dashEnabled;
}

Chart::Chart()
{
	m_yminAuto = true;
	m_ymaxAuto = true;
	m_xminAuto = true;
	m_xmaxAuto = true;
	m_xmin = 0.0;
	m_xmax = 1.0;
	m_leftBorder = 42;
	m_bottomBorder = 20;
	m_topBorder = 10;
	m_rightBorder = 10;
}

void Chart::add(const Histogram1dd* i_data, const LineStyle& style)
{
	m_data.push_back(i_data);
	m_dataStyle.push_back(style);
}

void Chart::add(const Histogram1dd& i_data, const LineStyle& style)
{
	m_data.push_back(&i_data);
	m_dataStyle.push_back(style);
}

void Chart::addHLine(double y, const LineStyle& style)
{
	m_hlines.push_back(std::pair<double,LineStyle>(y,style));
}

void Chart::addVLine(double x, const LineStyle& style)
{
	m_vlines.push_back(std::pair<double,LineStyle>(x,style));
}

void Chart::getXAxisLimit(double& xmin, double& xmax) const
{
	double tmp;
	double rmin;
	double rmax;
	
	if(m_xminAuto)
	{
		std::vector<const Histogram1dd*>::const_iterator iter;
		
		iter = m_data.begin();
		if(iter != m_data.end())
		{
			rmin = (*iter)->getMinPosition()[0];
			rmax = (*iter)->getMaxPosition()[0];
			
			iter++;
			while(iter != m_data.end())
			{
				tmp = (*iter)->getMinPosition()[0];
				if(tmp < rmin) rmin = tmp;
				
				tmp = (*iter)->getMaxPosition()[0];
				if(tmp > rmax) rmax = tmp;
				
				iter++;
			}

			xmin = rmin;
			xmax = rmax;
		}
		else throw stk::exception::Message("Chart : No data found to find xmin/xmax", STK_DBG_INFO);
	}
	else
	{
		xmin = m_xmin;
		xmax = m_xmax;
	}
}

void Chart::getYAxisLimit(double& ymin, double& ymax) const
{
	double tmp;
	double rmin;
	double rmax;
	
	if(m_yminAuto)
	{
		std::vector<const Histogram1dd*>::const_iterator iter;
		
		iter = m_data.begin();
		if(iter != m_data.end())
		{
			rmin = (*iter)->getMin();
			rmax = (*iter)->getMax();
			
			iter++;
			while(iter != m_data.end())
			{
				tmp = (*iter)->getMin();
				if(tmp < rmin) rmin = tmp;
				
				tmp = (*iter)->getMax();
				if(tmp > rmax) rmax = tmp;
				
				iter++;
			}

			ymin = rmin;
			ymax = rmax;
		}
		else throw stk::exception::Message("Chart : No data found to find ymin/ymax", STK_DBG_INFO);
	}
	else
	{
		ymin = m_ymin;
		ymax = m_ymax;
	}
}

void Chart::setXMin(double xmin)
{
	m_xmin = xmin;
	m_xminAuto = false;
}

void Chart::setXMax(double xmax)
{
	m_xmax = xmax;
	m_xmaxAuto = false;
}

void Chart::setYMin(double ymin)
{
	m_ymin = ymin;
	m_yminAuto = false;
}

void Chart::setYMax(double ymax)
{
	m_ymax = ymax;
	m_ymaxAuto = false;
}

void Chart::setXLabel(const std::string& label)
{
	m_xLabel = label;
}

void Chart::setYLabel(const std::string& label)
{
	m_yLabel = label;
}

Vector2d Chart::pos2Image(
	const Vector2d& i_pos,
	const Vector2d& i_size,
	const Vector2d& i_min,
	const Vector2d& i_max) const
{
	Vector2d size = i_size;
	size[0] -= m_leftBorder + m_rightBorder;
	size[1] -= m_topBorder + m_bottomBorder;
	
	Vector2d res = (i_pos-i_min)/(i_max-i_min);
	res *= size;

	res[0] += m_leftBorder;
	res[1] += m_bottomBorder;

	//Reverse y
	res[1] = (i_size[1]-1.0)-res[1];
	
	return res;
}

void Chart::getAxisDivision(std::vector<double>& o_axis, double i_min, double i_max) const
{
	double diff = i_max - i_min;
	double unit = pow(10.0, floor(log10(diff)));
	double round_min = floor(i_min/unit)*unit;
	double round_max = ceil(i_max/unit)*unit;
	double round_diff = round_max - round_min;

	int type = round_diff/unit;

	double nbPart;
	switch(type)
	{
		case 1:
			nbPart = 5.0;
			break;
		case 2:
			nbPart = 8.0;
			break;
		case 3:
			nbPart = 6.0;
			break;
		case 4:
			nbPart = 8.0;
			break;
		case 5:
			nbPart = 5.0;
			break;
		case 6:
			nbPart = 6.0;
			break;
		case 7:
			nbPart = 7.0;
			break;
		case 8:
			nbPart = 8.0;
			break;
		case 9:
			nbPart = 9.0;
			break;
		case 10:
			nbPart = 5.0;
			break;
		default:
			nbPart = 1.0;
	}	
	double step = round_diff/nbPart;
	for(double i=round_min; i<=round_max; i += step)
	{
		//if(abs(i) < step/1000.0) i = 0.0;
		if(i >= i_min && i <= i_max)
		{
			o_axis.push_back(i);
		}
	}
}

void Chart::draw(const std::string& i_filename, const Vector2i& i_size) const
{
	double dashedLine[] = {3.0};
	
	//Auto window
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	getXAxisLimit(xmin, xmax);
	getYAxisLimit(ymin, ymax);

	Vector2d imgSize(i_size[0], i_size[1]);
	Vector2d valMin(xmin, ymin);
	Vector2d valMax(xmax, ymax);

	cairo_surface_t* surface;
	cairo_t* context;
	
	surface = cairo_image_surface_create(
		CAIRO_FORMAT_ARGB32,
		i_size[0],
		i_size[1]);
	context = cairo_create(surface);

	/* Background *****************************************************/
	cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 1.0);
	cairo_paint(context);

	/* Lines **********************************************************/
	cairo_set_line_width (context, 1.0);
		
	std::vector< std::pair<double,LineStyle> >::const_iterator liter;
	liter = m_vlines.begin();
	while(liter != m_vlines.end())
	{
		//Style
		const LineStyle& style = liter->second;
		
		cairo_set_line_width (context, style.getStrokeWidth());
	
		if(style.dashEnabled()) cairo_set_dash(context, dashedLine, 1, 0);
		else cairo_set_dash(context, NULL, 0, 0);
		
		Vector4d strokeColor(0.8, 0.8, 0.8, 1.0);
		if(!style.defaultStrokeColor()) strokeColor = style.getStrokeColor();
		cairo_set_source_rgba(context, strokeColor[0], strokeColor[1], strokeColor[2], strokeColor[3]);
		
		//Geometry
		Vector2d ptImg;

		ptImg = pos2Image(Vector2d(liter->first, valMin[1]), imgSize, valMin, valMax);
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		
		ptImg = pos2Image(Vector2d(liter->first, valMax[1]), imgSize, valMin, valMax);
		cairo_line_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		liter++;

		cairo_stroke(context);
	}
		
	liter = m_hlines.begin();
	while(liter != m_hlines.end())
	{
		//Style
		const LineStyle& style = liter->second;
		
		cairo_set_line_width (context, style.getStrokeWidth());
		
		if(style.dashEnabled()) cairo_set_dash(context, dashedLine, 1, 0);
		else cairo_set_dash(context, NULL, 0, 0);
		
		Vector4d strokeColor(0.8, 0.8, 0.8, 1.0);
		if(!style.defaultStrokeColor()) strokeColor = style.getStrokeColor();
		cairo_set_source_rgba(context, strokeColor[0], strokeColor[1], strokeColor[2], strokeColor[3]);
		
		//Geometry
		Vector2d ptImg;

		ptImg = pos2Image(Vector2d(valMin[0], liter->first), imgSize, valMin, valMax);
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		
		ptImg = pos2Image(Vector2d(valMax[0], liter->first), imgSize, valMin, valMax);
		cairo_line_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		liter++;

		cairo_stroke(context);
	}
	cairo_set_dash(context, NULL, 0, 0);
	cairo_set_line_width (context, 1.0);

	/* Data ***********************************************************/
	for(int i=0; i<m_data.size(); i++)
	{
		const Histogram1dd* histo = m_data[i];
		const LineStyle& style = m_dataStyle[i];
		
		//Filling
		if(style.fillingEnabled())
		{
			Vector4d color = style.getFillColor();
			cairo_set_source_rgba(context, color[0], color[1], color[2], color[3]);
			
			if(histo->getArraySize() > 0)
			{
				{
					Vector2d ptHisto(histo->getPosFromIndice(0)[0], 0.0);
					Vector2d ptImg = pos2Image(ptHisto, imgSize, valMin, valMax);
					cairo_move_to(context, ptImg[0], ptImg[1]-0.5);
				}
			
				int sz = histo->getArraySize();
				for(int i=0; i<sz; i++)
				{
					Vector2d ptHisto(histo->getPosFromIndice(i)[0], histo->getFromIndice(i));
					Vector2d ptImg = pos2Image(ptHisto, imgSize, valMin, valMax);

					cairo_line_to(context, ptImg[0], ptImg[1]-0.5);
				}
				
				{
					Vector2d ptHisto(histo->getPosFromIndice(sz-1)[0], 0.0);
					Vector2d ptImg = pos2Image(ptHisto, imgSize, valMin, valMax);
					cairo_line_to(context, ptImg[0], ptImg[1]-0.5);
				}
			}
			
			cairo_fill(context);
		}
		
		if(style.strokeEnabled())
		{
			cairo_set_line_width (context, style.getStrokeWidth());
			
			Vector4d color(1.0, 0.0, 0.0, 1.0);
			if(!style.defaultStrokeColor()) color = style.getStrokeColor();
			cairo_set_source_rgba(context, color[0], color[1], color[2], color[3]);
			
			if(style.dashEnabled()) cairo_set_dash(context, dashedLine, 1, 0);
			else cairo_set_dash(context, NULL, 0, 0);
			
			if(histo->getArraySize() > 0)
			{
				{
					Vector2d ptHisto(histo->getPosFromIndice(0)[0], histo->getFromIndice(0));
					Vector2d ptImg = pos2Image(ptHisto, imgSize, valMin, valMax);
						
					cairo_move_to(context, ptImg[0], ptImg[1]-0.5);
				}
			
				for(int i=1; i<histo->getArraySize(); i++)
				{
					Vector2d ptHisto(histo->getPosFromIndice(i)[0], histo->getFromIndice(i));
					Vector2d ptImg = pos2Image(ptHisto, imgSize, valMin, valMax);

					cairo_line_to(context, ptImg[0], ptImg[1]-0.5);
				}
			}
			
			cairo_stroke(context);
		}
		cairo_set_dash(context, NULL, 0, 0);
		cairo_set_line_width (context, 1.0);
	}

	cairo_set_source_rgb(context, 0.1, 0.1, 0.1);
	//Axis
	{
		Vector2d ptImg;
		ptImg = pos2Image(
			Vector2d(0.0, 1.0),
			imgSize,
			Vector2d(0.0, 0.0),
			Vector2d(1.0, 1.0));
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		
		ptImg = pos2Image(
			Vector2d(0.0, 0.0),
			imgSize,
			Vector2d(0.0, 0.0),
			Vector2d(1.0, 1.0));
		cairo_line_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		
		ptImg = pos2Image(
			Vector2d(1.0, 0.0),
			imgSize,
			Vector2d(0.0, 0.0),
			Vector2d(1.0, 1.0));
		cairo_line_to(context, ptImg[0]-0.5, ptImg[1]-0.5);

		cairo_stroke(context);
	}

	//Numbers and labels
	cairo_select_font_face(context, "Ubuntu",
		CAIRO_FONT_SLANT_NORMAL,
		CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(context, 11);
	
	cairo_text_extents_t extents;

	std::vector<double> yAxis;
	getAxisDivision(yAxis, valMin[1], valMax[1]);
	for(int i=0; i<yAxis.size(); i++)
	{
		double y = yAxis[i];
		Vector2d ptImg = pos2Image(
			Vector2d(valMin[0], y),
			imgSize,
			valMin,
			valMax);

		//Line
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		cairo_rel_line_to(context, -3.0, 0.0);
		cairo_stroke(context);
		
		std::stringstream text;
		text << y;

		//get text size
		cairo_text_extents(context, text.str().c_str(), &extents);
		
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		cairo_rel_move_to(context, -6.0-extents.width, extents.height/2);
		cairo_show_text(context, text.str().c_str());
	}

	std::vector<double> xAxis;
	getAxisDivision(xAxis, valMin[0], valMax[0]);
	for(int i=0; i<xAxis.size(); i++)
	{
		double x = xAxis[i];
		Vector2d ptImg = pos2Image(
			Vector2d(x, valMin[1]),
			imgSize,
			valMin,
			valMax);

		//Line
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		cairo_rel_line_to(context, 0.0, 3.0);
		cairo_stroke(context);
		
		std::stringstream text;
		text << x;

		//get text size
		cairo_text_extents(context, text.str().c_str(), &extents);
		
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		cairo_rel_move_to(context, -extents.width/2, 6.0+extents.height);
		cairo_show_text(context, text.str().c_str());
	}
	
	if(m_xLabel != "")
	{
		Vector2d ptImg = pos2Image(
			Vector2d(1.0, 0.0),
			imgSize,
			Vector2d(0.0, 0.0),
			Vector2d(1.0, 1.0));
			
		cairo_text_extents(context, m_xLabel.c_str(), &extents);
		
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		cairo_rel_move_to(context, -extents.width-6.0, -6.0);
		cairo_show_text(context, m_xLabel.c_str());
	}
	
	if(m_yLabel != "")
	{
		Vector2d ptImg = pos2Image(
			Vector2d(0.0, 1.0),
			imgSize,
			Vector2d(0.0, 0.0),
			Vector2d(1.0, 1.0));
			
		cairo_text_extents(context, m_yLabel.c_str(), &extents);
		
		cairo_move_to(context, ptImg[0]-0.5, ptImg[1]-0.5);
		cairo_save(context);
		cairo_rotate(context, -M_PI/2);
		cairo_rel_move_to(context, -extents.width, extents.height+6.0);
		cairo_show_text(context, m_yLabel.c_str());
		cairo_restore(context);
	}
	
	//Export
	cairo_surface_write_to_png(surface, i_filename.c_str());
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
}

Chart2d::Chart2d()
{
	m_zminAuto = true;
	m_zmaxAuto = true;
	m_zmin = 0.0;
	m_zmax = 1.0;
	m_autoz_mask = NULL;
	m_show_min = false;
	m_show_max = false;
	m_enable_showPeakCircle = false;
	m_enable_hideData = false;
	m_colorFunction = DEFAULT;
	m_colorFunctionCrop = true;
}

void Chart2d::add(const Histogram2dd& i_data)
{
	m_data = &i_data;
}

void Chart2d::add(const Histogram2dd* i_data)
{
	m_data = i_data;
}

void Chart2d::addCircle(const Vector2d& i_center, double i_radius)
{
	m_circle_center.push_back(i_center);
	m_circle_radius.push_back(i_radius);
}

void Chart2d::addLabel(const std::string& i_msg)
{
	m_labels.push_back(i_msg);
}

void Chart2d::setColorFunction(const ColorFunction& func, bool crop)
{
	m_colorFunction = func;
	m_colorFunctionCrop = crop;
}

void Chart2d::setZMin(double zmin)
{
	m_zmin = zmin;
	m_zminAuto = false;
}

void Chart2d::setZMax(double zmax)
{
	m_zmax = zmax;
	m_zmaxAuto = false;
}

void Chart2d::setAutoZMask(const Histogram2di& i_mask)
{
	m_autoz_mask = &i_mask;
}

void Chart2d::setAutoZMask(const Histogram2di* i_mask)
{
	m_autoz_mask = i_mask;
}

void Chart2d::showMaximumWithColor()
{
	m_show_max = true;
}

void Chart2d::showMinimumWithColor()
{
	m_show_min = true;
}

void Chart2d::showPeaksWithCircle(
	double minPeak,
	double maxPeak,
	double minRadius,
	double maxRadius,
	bool onlyPeak)
{
	m_enable_showPeakCircle = true;
	m_enable_hideData = onlyPeak;
	m_peakCircleZMin = minPeak;
	m_peakCircleZMax = maxPeak;
	m_peakCircleRMin = minRadius;
	m_peakCircleRMax = maxRadius;
}

void Chart2d::generateColorFromValueHue(double c, unsigned char* p) const
{
	Vector3d rgb = color::hsv2rgb(stk::Vector3d(c, 1.0, 1.0));
	p[0] = rgb[0]*255;
	p[1] = rgb[1]*255;
	p[2] = rgb[2]*255;
	p[3] = 255;
}

void Chart2d::generateColorFromValueGrey(double c, unsigned char* p) const
{
	p[0] = c*255;
	p[1] = c*255;
	p[2] = c*255;
	p[3] = 255;
}

void Chart2d::generateColorFromValue(double c, unsigned char* p) const
{
	if(c > 1 && m_show_max)
	{
		p[0] = 1;
		p[1] = 1;
		p[2] = 221;
		p[3] = 255;
	}
	else if(c < 0 && m_show_min)
	{
		p[0] = 211;
		p[1] = 132;
		p[2] = 0;
		p[3] = 255;
	}
	else
	{
		if(m_colorFunctionCrop) c = std::min(std::max(c, 0.0), 1.0);
		else c = c - std::floor(c);
		
		switch(m_colorFunction)
		{
			case HUE:
				generateColorFromValueHue(c, p);
				break;
			default:
				generateColorFromValueGrey(c, p);
		}
	}
}

void Chart2d::draw(const std::string& i_filename) const
{
	//Compute rescaling
	double data_min = m_data->getFromIndice(0);
	double data_max = m_data->getFromIndice(0);
	
	if(m_autoz_mask != NULL)
	{
		for(int i=0; i<m_data->getArraySize(); i++)
		{
			if(m_autoz_mask->getFromIndice(i) > 0)
			{
				const double& v = m_data->getFromIndice(i);
				if(data_min > v) data_min = v;
			}
		}
	}
	else
	{
		for(int i=0; i<m_data->getArraySize(); i++)
		{
			const double& v = m_data->getFromIndice(i);
			if(data_min > v) data_min = v;
		}
	}
	
	if(m_autoz_mask != NULL)
	{
		for(int i=0; i<m_data->getArraySize(); i++)
		{
			if(m_autoz_mask->getFromIndice(i) > 0)
			{
				const double& v = m_data->getFromIndice(i);
				if(data_max < v) data_max = v;
			}
		}
	}
	else
	{
		for(int i=0; i<m_data->getArraySize(); i++)
		{
			const double& v = m_data->getFromIndice(i);
			if(data_max < v) data_max = v;
		}
	}

	double zmin = (m_zminAuto) ? data_min : m_zmin;
	double zmax = (m_zmaxAuto) ? data_max : m_zmax;
			
	//Generate data
	int width = m_data->getSize()[0];
	int height = m_data->getSize()[1];
	int stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, width);
	
	unsigned char* rgb_pixmap = new unsigned char[stride*height];
	if(m_enable_hideData)
	{
		for(int j=0; j<height; j++)
		{
			for(int i = 0; i<width; i++)
			{
				unsigned char* p = rgb_pixmap + j*stride + i*4;
				p[0] = 255;
				p[1] = 255;
				p[2] = 255;
				p[3] = 255;
			}
		}
	}
	else
	{
		double c;
		for(int j=0; j<height; j++)
		{
			for(int i = 0; i<width; i++)
			{
				c = (m_data->getData(Vector2i(i, j))-zmin)/(zmax-zmin);
				unsigned char* p = rgb_pixmap + j*stride + i*4;
				
				generateColorFromValue(c, p);
			}
		}
	}

	//Search peak
	std::vector<Vector2i> peakPos;
	std::vector<double> peakRadius;
	if(m_enable_showPeakCircle)
	{
		for(int i=0; i<m_data->getArraySize(); i++)
		{
			if(m_autoz_mask == NULL || m_autoz_mask->getFromIndice(i) > 0)
			{
				double v = m_data->getFromIndice(i);
				if(v > m_peakCircleZMin)
				{
					double radius = m_peakCircleZMin + (v - m_peakCircleZMin)/(m_peakCircleZMax - m_peakCircleZMin);

					radius = m_peakCircleRMin + radius*(m_peakCircleRMax - m_peakCircleRMin);
					
					peakPos.push_back(m_data->getCoordFromIndice(i));
					peakRadius.push_back(radius);
				}
			}
		}
	}

	//Import in cairo
	cairo_surface_t* surface;
	cairo_t* context;

	surface = cairo_image_surface_create_for_data(
		rgb_pixmap,
		CAIRO_FORMAT_ARGB32,
		width,
		height,
		stride);
	context = cairo_create(surface);

	//Draw circles
	cairo_set_line_width (context, 1.0);
	Vector2d minPos = m_data->getMinPosition();
	Vector2d maxPos = m_data->getMaxPosition();
	Vector2d spaceSize = maxPos - minPos;

	cairo_set_source_rgba(context, 0.0, 1.0, 0.0, 1.0);
	
	//~ cairo_set_dash(context, dashed3, 1, 0);
	for(int i=0; i<m_circle_radius.size(); i++)
	{
		Vector2d center = (m_circle_center[i] - minPos) / spaceSize;
		center *= Vector2d(width, height);
		
		cairo_arc(context, center[0], center[1], m_circle_radius[i], 0, 2*M_PI);
		cairo_stroke(context);
	}
	cairo_set_dash(context, dashed3, 0, 0);

	//Draw Peaks
	cairo_set_line_width (context, 1.0);
	for(int i=0; i<peakPos.size(); i++)
	{
		double radius = peakRadius[i];
		if(radius > m_peakCircleRMax)
		{
			cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 0.5);
			cairo_arc(context,
				peakPos[i][0]+0.5, peakPos[i][1]+0.5,
				m_peakCircleRMax, 0, 2*M_PI);
			cairo_fill(context);
			
			cairo_set_source_rgba(context, 1.0, 0.0, 0.0, 1.0);
			cairo_arc(context,
				peakPos[i][0]+0.5, peakPos[i][1]+0.5,
				m_peakCircleRMax, 0, 2*M_PI);
			cairo_stroke(context);
		}
		else
		{
			cairo_set_source_rgba(context, 1.0, 1.0, 0.0, 0.5);
			cairo_arc(context,
				peakPos[i][0]+0.5, peakPos[i][1]+0.5,
				radius, 0, 2*M_PI);
			cairo_fill(context);
			
			cairo_set_source_rgba(context, 0.8, 0.8, 0.0, 1.0);
			cairo_arc(context,
				peakPos[i][0]+0.5, peakPos[i][1]+0.5,
				radius, 0, 2*M_PI);
			cairo_stroke(context);
		}
	}

	//Draw labels
	cairo_select_font_face(context, "Ubuntu",
		CAIRO_FONT_SLANT_NORMAL,
		CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(context, 13);
	
	cairo_text_extents_t extents;

	cairo_move_to(context, 10.5, 10.5);
	double labelWidth;
	int labelHeight = 0;
	for(int i=0; i<m_labels.size(); i++)
	{
		cairo_text_extents(context, m_labels[i].c_str(), &extents);
		labelWidth = std::max(labelWidth, extents.width);
		labelHeight += extents.height+6;
	}
	
	cairo_set_source_rgba(context, 1.0, 1.0, 1.0, 0.75);
	cairo_rectangle(context, 10.5, 10.5, labelWidth+8, (double) labelHeight);
	cairo_fill(context);
	
	cairo_set_source_rgba(context, 0.0, 0.0, 0.0, 1.0);
	int hOffset = 0;
	for(int i=0; i<m_labels.size(); i++)
	{
		cairo_text_extents(context, m_labels[i].c_str(), &extents);
		hOffset += extents.height;
		
		cairo_move_to(context, 12.5, 10.5 + hOffset);
		cairo_show_text(context, m_labels[i].c_str());
		
		hOffset += 6;
	}
	
	
	//Export
	cairo_surface_write_to_png(surface, i_filename.c_str());
	
	cairo_destroy(context);
	cairo_surface_destroy(surface);
}

#endif

}

}
