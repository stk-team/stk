#ifndef __STK_FUNCTION_WINDOW__
#define __STK_FUNCTION_WINDOW__


#include <stdlib.h>
#include <math.h>
#include <sstream>

#include "stk/function/generic.hpp"

using namespace std;

namespace stk {
namespace function {

/**
 * @class Window
 * @brief Abstract class defining additionnal parameter specific to windowing function
 */
template<int DIM, typename POS>
class Window
{
public:
	/**
	 * @brief Constructor
	 * @param i_metric metric to use in evaluation (see Vector::norm)
	 * @param i_bound space bound of window
	 * @param i_pos center position of window
	 */
	Window(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	m_metric(i_metric),
		m_bound(i_bound),
		m_pos(i_pos)
	{};

/**
 * @brief Accessor to m_metric attribute
 */
	unsigned int& metric() { return m_metric; };
/**
 * @brief Accessor to m_metric attribute (const version)
 */
	const unsigned int& metric() const { return m_metric; };
/**
 * @brief Mutator for m_metric attribute
 * @param i_metric value for m_metric attribute
 */
	void metric(unsigned int& i_metric) { m_metric = i_metric; };

/**
 * @brief Accessor to m_bound attribute
 */
	POS& bound() { return m_bound; };
/**
 * @brief Accessor to m_bound attribute (const version)
 */
	const POS& bound() const { return m_bound; };
/**
 * @brief Mutator for m_bound attribute
 * @param i_bound value for m_bound attribute
 */
	void bound(POS& i_bound) { m_bound = i_bound; };

/**
 * @brief Accessor to m_pos attribute
 */
	Vector<DIM, POS>& pos() { return m_pos; };
/**
 * @brief Accessor to m_pos attribute (const version)
 */
	const Vector<DIM, POS>& pos() const { return m_pos; };
/**
 * @brief Mutator for m_pos attribute
 * @param i_pos value for m_pos attribute
 */
	void pos(Vector<DIM, POS>& i_pos) { m_pos = i_pos; };

protected:
	unsigned int m_metric; /**< @brief space metric to use (see Vector::norm)*/
	POS m_bound; /**< @brief limit of window (-m_bound < ||x|| < +m_bound)*/
	Vector<DIM, POS> m_pos; /**< @brief position of window center*/
};

/**
 * @class WindowFunction
 * @brief Abstract class defining windowing function
 * @sa GenericFunction
 * @sa Window
 */
template<int DIM, typename POS, typename VAL>
class WindowFunction : public GenericFunction<DIM, POS, VAL>, public Window<DIM, POS>
{
public:
/**
 * @copydoc GenericFunction::GenericFunction
 * @copydoc Window::Window
 */
	WindowFunction(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	GenericFunction<DIM, POS, VAL>(),
		Window<DIM, POS>(i_metric, i_bound, i_pos)
	{};
};


/**
 * @class ConstantWindow
 * @brief window with constant value (box window)
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class ConstantWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 * @param i_lvl value of the constant function
 */
	ConstantWindow(const VAL& i_lvl = 1, const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos),
		m_lvl(i_lvl)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		if(i_pos.norm(this->m_metric) < this->m_bound)
			return m_lvl;
		else
			return 0;
	};

/**
 * @brief Accessor to m_lvl attribute
 */
	VAL& lvl() { return m_lvl; };
/**
 * @brief Accessor to m_lvl attribute (const version)
 */
	const VAL& lvl() const { return m_lvl; };
/**
 * @brief Mutator for m_lvl attribute
 * @param i_lvl value for m_lvl attribute
 */
	void lvl(VAL& i_lvl) { m_lvl = i_lvl; };

protected:
	VAL m_lvl; /**< @brief level(value) of constant function*/
};

/**
 * @class HannWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class HannWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 */
	HannWindow(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos < this->m_bound)
			return (1. + cos(PI*pos/this->m_bound))/2.;
		else
			return 0;
	};
};

/**
 * @class HammingWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class HammingWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 */
	HammingWindow(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos < this->m_bound)
			return .54 + .46*cos(PI*pos/this->m_bound);
		else
			return 0;
	};
};

/**
 * @class TukeyWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class TukeyWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 * @param i_alpha half-size of flat top
 */
	TukeyWindow(const POS& i_alpha = 0.25, const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos),
		m_alpha(i_alpha)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos <= m_alpha)
			return 1;
		else if(pos < this->m_bound)
			return (1. + cos( PI*(pos-m_alpha) / (this->m_bound-m_alpha) ))/2.;
		else
			return 0;
	};

/**
 * @brief Accessor to m_alpha attribute
 */
	POS& alpha() { return m_alpha; };
/**
 * @brief Accessor to m_alpha attribute (const version)
 */
	const POS& alpha() const { return m_alpha; };
/**
 * @brief Mutator for m_alpha attribute
 * @param i_alpha value for m_alpha attribute
 */
	void alpha(POS& i_alpha) { m_alpha = i_alpha; };

protected:
	POS m_alpha; /**< @brief half-size of Tukey window flat top*/
};

/**
 * @class CosineWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class CosineWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 */
	CosineWindow(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos < this->m_bound)
			return cos( PI*pos/(2*this->m_bound) );
		else
			return 0;
	};
};

/**
 * @class LanczosWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class LanczosWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 */
	LanczosWindow(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos > this->m_bound)
			return 0;
		else if(pos == 0)
			return 1;
		else
		{
			POS tmp = PI*pos/this->m_bound;
			return sin(tmp)/tmp;
		}
	};
};

/**
 * @class TriangleWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class TriangleWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 */
	TriangleWindow(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos < this->m_bound)
			return (this->m_bound-pos)/this->m_bound;
		else
			return 0;
	};
};

/**
 * @class GaussWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class GaussWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 */
	GaussWindow(const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos < this->m_bound)
			return exp(-pow(pos/(.4*this->m_bound), 2));
		else
			return 0;
	};
};

/**
 * @class BlackmanWindow
 * @sa WindowFunction
 */
template<int DIM, typename POS, typename VAL>
class BlackmanWindow : public WindowFunction<DIM, POS, VAL>
{
public:
/**
 * @copydoc WindowFunction::WindowFunction
 * @param i_alpha parameter in Blackman window equation
 */
	BlackmanWindow(const POS& i_alpha = 0.16, const unsigned int& i_metric = 2, const POS& i_bound = 0.5, const Vector<DIM, POS>& i_pos = Vector<DIM, POS>(0))
	:	WindowFunction<DIM, POS, VAL>(i_metric, i_bound, i_pos),
		m_alpha(i_alpha)
	{};

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		POS pos = i_pos.norm(this->m_metric);
		if(pos < this->m_bound)
			return ((1.-m_alpha)/2. - cos(PI*(pos-this->m_bound)/this->m_bound)/2. + m_alpha*cos(2*PI*(pos-this->m_bound)/this->m_bound)/2.);
		else
			return 0;
	};

/**
 * @brief Accessor to m_alpha attribute
 */
	POS& alpha() { return m_alpha; };
/**
 * @brief Accessor to m_alpha attribute (const version)
 */
	const POS& alpha() const { return m_alpha; };
/**
 * @brief Mutator for m_alpha attribute
 * @param i_alpha value for m_alpha attribute
 */
	void alpha(POS& i_alpha) { m_alpha = i_alpha; };

protected:
	POS m_alpha; /**< @brief alpha parameter in Blackman parametrization*/
};




////////////////////////////////////////////////////////////////////////////////////
/*
inline double hammingWindow(double x, const double sz)
{
	double h=0.1, w=0.5;
	if(std::abs(x)>sz) return h;
	else return ( h + PI*(1-h*4*sz*sz) * (1 + cos(PI*x/sz)) / (w*w*(PI*PI-4)) );
}
*/

/*
void testWindows(std::string dir)
{
	double sz=1, szP=sz+.2, step=0.01;
	vector<double> data;
	std::stringstream out;

	cout << "Generate points for rectangular window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(rectWindow(i,sz));
	out.str("");
	out << dir << "/winRect.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for hann window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(hannWindow(i,sz));
	out.str("");
	out << dir << "/winHann.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for hamming window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(hammingWindow(i,sz));
	out.str("");
	out << dir << "/winHamming.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Tukey (alpha=0.25) window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(tukeyWindow(i,sz,0.25));
	out.str("");
	out << dir << "/winTukey_25.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Tukey (alpha=0.50) window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(tukeyWindow(i,sz,0.5));
	out.str("");
	out << dir << "/winTukey_50.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Tukey (alpha=0.75) window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(tukeyWindow(i,sz,0.75));
	out.str("");
	out << dir << "/winTukey_75.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Cosine window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(cosineWindow(i,sz));
	out.str("");
	out << dir << "/winCosine.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Lanczos window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(lanczosWindow(i,sz));
	out.str("");
	out << dir << "/winLanczos.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Triangle window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(triangWindow(i,sz));
	out.str("");
	out << dir << "/winTriang.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Gaussian window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(gaussWindow(i,sz));
	out.str("");
	out << dir << "/winGauss.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());

	cout << "Generate points for Blackman window" << endl;
	data.clear();
	for(double i=-szP; i<=szP; i+=step)
		data.push_back(blackmanWindow(i,sz));
	out.str("");
	out << dir << "/winBlackman.dat";
	cout << "   |--> save points in \"" << out.str() << "\"" << endl;
	saveData1D(data, out.str());
}
*/

} /* function */
} /* stk */

#endif
