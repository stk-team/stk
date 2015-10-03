#ifndef __STK_FUNCTION_SPECTRUM__
#define __STK_FUNCTION_SPECTRUM__

#include <iostream>
#include <cassert>
#include <cmath>

#include "stk/vector.hpp"
#include "stk/function/generic.hpp"
#include "stk/function/base.hpp"


namespace stk {
namespace function {

/**
 * @class Spectrum
 * @brief Abstract class defining additionnal parameter specific to spectrum function
 */
template<int DIM, typename POS, typename VAL>
class Spectrum : public GenericFunction<DIM, POS, VAL>
{
public:
	/**
	 * @brief Constructor for absolute spectrum
	 */
	Spectrum()
	:	GenericFunction<DIM, POS, VAL>(),
		m_nbSamplingPt(0),
		m_isRelative(false)
	{};
	/**
	 * @brief Constructor for relative spectrum
	 * @param i_nbSamplingPt number of sampling point
	 */
	Spectrum(const unsigned int& i_nbSamplingPt)
	:	GenericFunction<DIM, POS, VAL>(),
		m_nbSamplingPt(i_nbSamplingPt),
		m_isRelative(true)
	{};

	VAL getRelVal(const Vector<DIM, POS>& i_pos)
	{
		return this->function( toRelativeFreq(i_pos) );
	}

	VAL operator() (const Vector<DIM, POS>& i_pos) const
	{
		if(m_isRelative)
			return toAbsoluteMag( this->function( toRelativeFreq(i_pos) ) );
		else
			return this->function(i_pos);
	};

protected:

	/**
	 * @brief convert absolute frequency to relative frequency coordinate
	 * @param i_pos point/vector to convert
	 * @return point/vector converted
	 */
	Vector<DIM, POS> toRelativeFreq(const Vector<DIM, POS>& i_pos) const
	{ return i_pos / pow(m_nbSamplingPt, 1./DIM); };

	/**
	 * @brief convert relative power magnitude to absolute magnitude
	 * @param i_mag magnitude value to convert
	 * @return magnitude value converted
	 */
	VAL toAbsoluteMag(const VAL& i_mag) const
	{ return sqrt( m_nbSamplingPt * i_mag); };

public:
/**
 * @brief Accessor to m_nbSamplingPt attribute
 */
	unsigned int& nbSamplingPt() { return m_nbSamplingPt; };
/**
 * @brief Accessor to m_nbSamplingPt attribute (const version)
 */
	const unsigned int& nbSamplingPt() const { return m_nbSamplingPt; };

/**
 * @brief Accessor to m_isRelative attribute (const version)
 */
	const bool& isRelative() const { return m_isRelative; };
/**
 * @brief Mutator for m_isRelative and m_nbSamplingPt attributes (set)
 * @param i_nbSamplingPt value for m_nbSamplingPt attribute
 */
	void setRelative(unsigned int& i_nbSamplingPt)
	{
		m_nbSamplingPt = i_nbSamplingPt;
		m_isRelative = true;
	};
/**
 * @brief Mutator for m_isRelative attribute (unset)
 */
	void unsetRelative()
	{
		m_isRelative = false;
		m_nbSamplingPt = 0;
	};

protected:
	unsigned int m_nbSamplingPt; /**< @brief number of sampling point from which other spectrum is compute*/
	bool m_isRelative; /**< @brief does the function have to be relative in input frequency and output magnitude*/
};

/**
 * @class WhiteNoise
 * @brief WhiteNoise spectrum
 */
template<int DIM, typename POS, typename VAL>
class WhiteNoise : public Spectrum<DIM, POS, VAL>
{
public:
/*
 * @copydoc Spectrum::Spectrum(i_nbSamplingPt)
 */
	WhiteNoise(const unsigned int& i_nbSamplingPt)
	:	Spectrum<DIM, POS, VAL>(i_nbSamplingPt)
	{};

protected:
	VAL function(const Vector<DIM, POS>& i_pos) const
	{ return 1; };
};

/**
 * @class StepNoise
 * @brief StepNoise spectrum
 */
template<int DIM, typename POS, typename VAL>
class StepNoise : public Spectrum<DIM, POS, VAL>
{
protected:
	POS m_step0;
	POS m_step1;
	VAL m_leftVal;
	VAL m_rightVal;

public:

	/*
	 * @copydoc Spectrum::Spectrum(i_nbSamplingPt)
	 */
	StepNoise(const unsigned int& i_nbSamplingPt, POS i_step)
		: Spectrum<DIM, POS, VAL>(i_nbSamplingPt),
		m_step0(0),
		m_step1(i_step),
		m_leftVal(0.0),
		m_rightVal(1.0)
	{}

	/*
	 * @copydoc Spectrum::Spectrum(i_nbSamplingPt)
	 */
	StepNoise(const unsigned int& i_nbSamplingPt, POS i_step0, POS i_step1)
		: Spectrum<DIM, POS, VAL>(i_nbSamplingPt),
		m_step0(i_step0),
		m_step1(i_step1),
		m_leftVal(0.0),
		m_rightVal(1.0)
	{}

	/*
	 * @copydoc Spectrum::Spectrum(i_nbSamplingPt)
	 */
	StepNoise(
		const unsigned int& i_nbSamplingPt,
		POS i_step,
		VAL i_leftVal,
		VAL i_rightVal
	)
		: Spectrum<DIM, POS, VAL>(i_nbSamplingPt),
		m_step0(0.0),
		m_step1(i_step),
		m_leftVal(i_leftVal),
		m_rightVal(i_rightVal)
	{}

	/*
	 * @copydoc Spectrum::Spectrum(i_nbSamplingPt)
	 */
	StepNoise(
		const unsigned int& i_nbSamplingPt,
		POS i_step0,
		POS i_step1,
		VAL i_leftVal,
		VAL i_rightVal
	)
		: Spectrum<DIM, POS, VAL>(i_nbSamplingPt),
		m_step0(i_step0),
		m_step1(i_step1),
		m_leftVal(i_leftVal),
		m_rightVal(i_rightVal)
	{}
	

protected:
	VAL function(const Vector<DIM, POS>& i_pos) const
	{
		if(i_pos.norm() <= m_step1 && i_pos.norm() >= m_step0) return m_leftVal;
		else return m_rightVal;
	}
};

/**
 * @class AnisoGaussianNoise
 * @brief AnisoGaussianNoise spectrum
 */
template<int DIM, typename POS, typename VAL>
class AnisoGaussianNoise : public Spectrum<DIM, POS, VAL>
{
protected:
	GaussianFunction<DIM, POS, double> m_gauss;

public:

	/*
	 * @copydoc Spectrum::Spectrum(i_nbSamplingPt)
	 */
	AnisoGaussianNoise(
		const unsigned int& i_nbSamplingPt,
		const Vector<DIM, POS>& i_pos,
		POS i_dev,
		VAL i_mag
	)
		: Spectrum<DIM, POS, VAL>(i_nbSamplingPt),
		m_gauss(i_pos, i_dev, i_mag)
	{}
	

protected:
	VAL function(const Vector<DIM, POS>& i_pos) const
	{
		return m_gauss(i_pos) + m_gauss(i_pos*-1.0);
	}
};




/*************************
Default spectrum class use

//! @class DC
//! @brief DC spectrum
template<int DIM, typename POS, typename VAL>
class DC : public Spectrum<DIM, POS, VAL>
{
public:
//! @copydoc Spectrum::Spectrum(i_nbSamplingPt)
	DC(const unsigned int& i_nbSamplingPt)
	:	Spectrum(i_nbSamplingPt)
	{};

protected:
	VAL function(const Vector<DIM, POS>& i_pos) const
	{};
};
*************************/


/////////////////////////////////////////////////////////////////////////

template<unsigned int DIM, class T>
T func_whiteNoise(const Vector<DIM,T> x)
{ return 1.; }

template<unsigned int DIM, class T>
T func_null(const Vector<DIM,T> x)
{ return 0.; }

template<unsigned int DIM, class T>
T func_ulichney(const Vector<DIM,T> x)
{
	double _aGauss = 0.2,
		_kDelta = 1.3,
		_deltaWidth = 0.00000000001,
		_sqrtPi = sqrt(PI),
		_X = x.norm();
	return ( _aGauss*_kDelta*
			( erf(_sqrtPi*(_deltaWidth - 2.*(_X-1.))/(2.*_aGauss))
			+ erf(_sqrtPi*(_deltaWidth + 2.*(_X-1.))/(2.*_aGauss)) )
			/_deltaWidth + erfc(_sqrtPi*(1.-_X)/_aGauss) ) /2.;
}

template<unsigned int DIM, class T>
T func_ulichney2(const Vector<DIM,T> x)
{
	int n = 1024;
	double height = 1.3;
	double aGauss = 0.2;
	double pPos = 2;
	
	double dWidth = 0.00000000001;
	double sqrtPi = sqrt(PI);
	double X = x.norm();
		
	double kFourier = 0.875*sqrt(n); //Normalisation à la con
	
	double diracFunc = aGauss/dWidth * (
			  erf(sqrtPi*( 2. + dWidth - (2.*X)/pPos)/(2.*aGauss))
			+ erf(sqrtPi*(-2. + dWidth + (2.*X)/pPos)/(2.*aGauss)) );
	double stepFunc = erfc(sqrtPi*(1.-X/pPos)/aGauss);
		
	return kFourier*(height * diracFunc + stepFunc )/2.;
}
	
template<unsigned int DIM, class T>
T func_freqWeight(const Vector<DIM, T> x)
{
	std::cout << "Error : no weight definition for dimensions other than 1d and 2d." << std::endl;
	exit(EXIT_FAILURE);
}
template<class T>
T func_freqWeight(const Vector<1, T> x)
{
	if(x[0]!=0.) return 1./(2e5 * pow(x[0], 2));
	else return 1;
}
template<class T>
T func_freqWeight(const Vector<2, T> x)
{
	if(x.norm()!=0.) return 1./(75 * pow(x.norm(), 2));
	else return 1;
}



// BASE functions
template<class T>
T func_gauss_norm(const T x, const T pos, const T dev)
{
	return exp(-pow(x-pos, 2)/(2*dev*dev))/(dev*sqrt(2*PI));
}
template<class T>
T func_gauss(const T x, const T pos, const T dev, const T a)
{
	return a * exp(-pow(x-pos, 2)/(2*dev*dev));
}
template<class T>
T func_box(const T x, const T pos, const T rad, const T a)
{
	return ((pos-rad)<x && x<(pos+rad)) ? a : 0;
}


//WEI PAPER functions
template<unsigned int DIM, class T>
T wei_pink(const Vector<DIM,T> x)
{
	T _x=x.norm(), pos=0, dev=0.6, a=3;
	return 1 + func_gauss<T>(_x, pos, dev, a);
}
template<unsigned int DIM, class T>
T wei_green(const Vector<DIM,T> x)
{
	T _x=x.norm(), pos=0.546875, dev=0.0625, a=5;
	return 1 + func_gauss<T>(_x, pos, dev, a);
}
template<unsigned int DIM, class T>
T wei_magenta1(const Vector<DIM,T> x)
{
	T _x=x.norm(), pos=0.546875, dev=0.0625, a=1;
	return 1 - func_gauss<T>(_x, pos, dev, a);
}
template<unsigned int DIM, class T>
T wei_magenta2(const Vector<DIM,T> x)
{
	T _x=x.norm(), pos=0.546875, rad=0.1, a=1;
	return 1 - func_box<T>(_x, pos, rad, a);
}
template<unsigned int DIM, class T>
T wei_blue2(const Vector<DIM,T> x)
{
	T _x=x.norm(), pos=0.3, rad=pos, a=1;
	return 1 - func_box<T>(_x, pos, rad, a);
}




#define BWidth 0.1
#define Freq_1_m 1.2
#define Freq_2_m 0.7
#define Freq_3_m 0.9
#define Freq_1_M Freq_1_m + BWidth
#define Freq_2_M Freq_2_m + BWidth
#define Freq_3_M Freq_3_m + BWidth

template<unsigned int DIM, class T>
T func_oneBump(const Vector<DIM,T> x)
{
	double _X = x.norm();
	if(Freq_1_m<_X && _X<Freq_1_M) return 2.;
	else return 1.;
}

template<unsigned int DIM, class T>
T func_oneHole(const Vector<DIM,T> x)
{
	double _X = x.norm();
	if(Freq_1_m<_X && _X<Freq_1_M) return 0.;
	else return 1.;
}
template<unsigned int DIM, class T>
T func_twoHoles(const Vector<DIM,T> x)
{
	double _X = x.norm();
	if((Freq_1_m<_X && _X<Freq_1_M) || (Freq_2_m<_X && _X<Freq_2_M)) return 0.;
	else return 1.;
}
template<unsigned int DIM, class T>
T func_threeHoles(const Vector<DIM,T> x)
{
	double _X = x.norm();
	if((Freq_1_m<_X && _X<Freq_1_M) || (Freq_2_m<_X && _X<Freq_2_M) || (Freq_3_m<_X && _X<Freq_3_M)) return 0.;
	else return 1.;
}
template<unsigned int DIM, class T>
T func_oneHill(const Vector<DIM,T> x)
{
	double _X = x.norm();
	if(Freq_1_m<_X && _X<Freq_1_M) return 1.;
	else return 0.;
}
template<unsigned int DIM, class T>
T func_twoHills(const Vector<DIM,T> x)
{
	double _X = x.norm();
	if((Freq_1_m<_X && _X<Freq_1_M) || (Freq_2_m<_X && _X<Freq_2_M)) return 1.;
	else return 0.;
}
template<unsigned int DIM, class T>
T func_threeHills(const Vector<DIM,T> x)
{
	double _X = x.norm();
	if((Freq_1_m<_X && _X<Freq_1_M) || (Freq_2_m<_X && _X<Freq_2_M) || (Freq_3_m<_X && _X<Freq_3_M)) return 1.;
	else return 0.;
}

int factorial(const int& n);

//Old definition
class UlichneyFunction
{
	private:
		double _deltaWidth;
		double _sqrtPi;
		double _scale;
		double _kFourier;
		double _aGauss;
		double _kDelta;

	public:
		UlichneyFunction(int n, int dim, double aGauss, double kDelta);
		double eval(double x) const;
		double getInf() const;
};

} /* function */
} /* stk */

#endif
