#ifndef __STK_IO_CHART__
#define __STK_IO_CHART__

#include <vector>
#include <string>

#include <stk/histogram.hpp>

namespace stk
{

namespace io
{

#ifdef CAIRO_ENABLED
class Chart
{
	public:
		enum DashStyle
		{
			LINESTYLE_DEFAULT,
			LINESTYLE_DASHED
		};
	
		class LineStyle
		{
			protected:
			
			public:
				LineStyle();
				
				//Stroke enabled
				bool m_strokeEnabled;
				LineStyle& stroke(bool enabled);
				bool strokeEnabled() const;
				
				//Stroke color
				bool m_defaultStrokeColor;
				bool defaultStrokeColor() const;
				Vector4d m_strokeColor;
				LineStyle& strokeColor(const Vector4d& color);
				const Vector4d& getStrokeColor() const;
				
				//Fill color
				float m_strokeWidth;
				LineStyle& strokeWidth(float width);
				float getStrokeWidth() const;
				
				//Filling enabled
				bool m_fillingEnabled;
				LineStyle& filling(bool enabled);
				bool fillingEnabled() const;
				
				//Fill color
				Vector4d m_fillColor;
				const Vector4d& getFillColor() const;
				
				//Dash enabled
				bool m_dashEnabled;
				LineStyle& dash(bool enabled);
				bool dashEnabled() const;
		};
		
	protected:
		std::vector<const Histogram1dd*> m_data;
		std::vector<LineStyle> m_dataStyle;
		std::vector< std::pair<double,LineStyle> > m_vlines;
		std::vector< std::pair<double,LineStyle> > m_hlines;

		bool m_xminAuto;
		bool m_xmaxAuto;
		bool m_yminAuto;
		bool m_ymaxAuto;

		double m_xmin;
		double m_xmax;
		double m_ymin;
		double m_ymax;
		
		std::string m_xLabel;
		std::string m_yLabel;

		int m_leftBorder;
		int m_bottomBorder;
		int m_topBorder;
		int m_rightBorder;

	private:
		void getXAxisLimit(double& xmin, double& xmax) const;
		void getYAxisLimit(double& ymin, double& ymax) const;
		
		Vector2d pos2Image(
			const Vector2d& i_pos,
			const Vector2d& i_size,
			const Vector2d& i_min,
			const Vector2d& i_max) const;
	
	public:	
		Chart();
		
		void add(const Histogram1dd* i_data, const LineStyle& style = LineStyle());
		void add(const Histogram1dd& i_data, const LineStyle& style = LineStyle());
		void addHLine(double y, const LineStyle& style = LineStyle());
		void addVLine(double x, const LineStyle& style = LineStyle());

		void setXMin(double xmin);
		void setXMax(double xmax);
		void setYMin(double ymin);
		void setYMax(double ymax);
		void setXLabel(const std::string& label);
		void setYLabel(const std::string& label);
	
		void draw(const std::string& i_filename, const Vector2i& i_size) const;
		
		void getAxisDivision(std::vector<double>& o_axis, double i_min, double i_max) const;
};

class Chart2d
{
	public:
		enum ColorFunction
		{
			DEFAULT,
			GREY,
			HUE
		};
	
	protected:
		const Histogram2dd* m_data;
		const Histogram2di* m_autoz_mask;
		std::vector<Vector2d> m_circle_center;
		std::vector<double> m_circle_radius;
		std::vector<std::string> m_labels;

		bool m_zminAuto;
		bool m_zmaxAuto;

		double m_zmin;
		double m_zmax;

		double m_show_min;
		double m_show_max;

		bool m_enable_showPeakCircle;
		bool m_enable_hideData;
		double m_peakCircleZMin;
		double m_peakCircleZMax;
		double m_peakCircleRMin;
		double m_peakCircleRMax;
		
		ColorFunction m_colorFunction;
		bool m_colorFunctionCrop;

	protected:
		void generateColorFromValue(double c, unsigned char* p) const;
		void generateColorFromValueGrey(double c, unsigned char* p) const;
		void generateColorFromValueHue(double c, unsigned char* p) const;

	public:
		Chart2d();

		void add(const Histogram2dd& i_data);
		void add(const Histogram2dd* i_data);
		void addCircle(const Vector2d& i_center, double i_radius);
		void addLabel(const std::string& i_msg);

		void setZMax(double zmax);
		void setZMin(double zmin);
		void setAutoZMask(const Histogram2di& i_mask);
		void setAutoZMask(const Histogram2di* i_mask);
		void showMaximumWithColor();
		void showMinimumWithColor();
		void showPeaksWithCircle(
			double minPeak,
			double maxPeak,
			double minRadius,
			double maxRadius,
			bool onlyPeak = true);
		void setColorFunction(const ColorFunction& func, bool crop = true);

		void draw(const std::string& i_filename) const;
};

#endif

}

}

#endif
