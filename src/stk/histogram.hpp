#ifndef __STK_HISTOGRAM__
#define __STK_HISTOGRAM__

#include <stk/vector.hpp>
#include <stk/array.hpp>

namespace stk
{

/**
 * @class Histogram
 * @todo Ajouter l'interpolation des valeurs
 */
template<int DIM, typename POS, typename VAL>
class Histogram : public Array<DIM, VAL>
{
	private:
		Vector<DIM, POS> m_spaceMin;
		Vector<DIM, POS> m_spaceMax;
		Vector<DIM, POS> m_spaceSize;
		
	
	public:
		Histogram();
		Histogram(const Histogram<DIM, POS, VAL>& histo);
		Histogram(const Vector<DIM, int>& i_size);
		
		void setBoundaries(
			const Vector<DIM, POS>& i_min,
			const Vector<DIM, POS>& i_max);
		void setBoundaries(const Histogram<DIM, POS, VAL>& histo);
		
		VAL& get(const Vector<DIM, POS>& i_pos);
		const VAL& get(const Vector<DIM, POS>& i_pos) const;
		VAL& getFromPos(const Vector<DIM, POS>& i_pos);
		const VAL& getFromPos(const Vector<DIM, POS>& i_pos) const;
		VAL getFromPosInterpolate(const Vector<DIM, POS>& i_pos);
		int getIndiceFromPos(const Vector<DIM, POS>& i_pos) const;
		Vector<DIM, POS> getPosFromIndice(int i_ind) const;
		
		const Vector<DIM, POS>& getMinPosition() const;
		const Vector<DIM, POS>& getMaxPosition() const;
		
		VAL integrate() const;
		
		void operator+=(const Histogram<DIM, POS, VAL>& i_histo);
		void operator-=(const Histogram<DIM, POS, VAL>& i_histo);
		void operator*=(double i_val);
		void operator/=(double i_val);
};

template<int DIM, typename POS, typename VAL>
Histogram<DIM, POS, VAL>::Histogram()
	: Array<DIM, VAL>()
{
	m_spaceMin = Vector<DIM, POS>(0);
	for(int i=0; i<DIM; i++) m_spaceMax[i] = this->getSize()[i]-1;
	m_spaceSize = m_spaceMax - m_spaceMin;
}

template<int DIM, typename POS, typename VAL>
Histogram<DIM, POS, VAL>::Histogram(const Histogram<DIM, POS, VAL>& histo)
	: Array<DIM, VAL>(histo)
{
	m_spaceMin = histo.getMinPosition();
	m_spaceMax = histo.getMaxPosition();
	m_spaceSize = m_spaceMax - m_spaceMin;
}

template<int DIM, typename POS, typename VAL>
Histogram<DIM, POS, VAL>::Histogram(const Vector<DIM, int>& i_size)
	: Array<DIM, VAL>(i_size)
{
	m_spaceMin = Vector<DIM, POS>(0);
	for(int i=0; i<DIM; i++) m_spaceMax[i] = this->getSize()[i]-1;
	m_spaceSize = m_spaceMax - m_spaceMin;
}
		
template<int DIM, typename POS, typename VAL>
void Histogram<DIM, POS, VAL>::setBoundaries(
	const Vector<DIM, POS>& i_min,
	const Vector<DIM, POS>& i_max)
{
	m_spaceMin = i_min;
	m_spaceMax = i_max;
	m_spaceSize = m_spaceMax - m_spaceMin;
}

template<int DIM, typename POS, typename VAL>
void Histogram<DIM, POS, VAL>::setBoundaries(const Histogram<DIM, POS, VAL>& histo)
{
	m_spaceMin = histo.getMinPosition();
	m_spaceMax = histo.getMaxPosition();
	m_spaceSize = m_spaceMax - m_spaceMin;
}

template<int DIM, typename POS, typename VAL>
int Histogram<DIM, POS, VAL>::getIndiceFromPos(const Vector<DIM, POS>& i_pos) const
{
	Vector<DIM, int> arrayPos;
	for(int i=0; i<DIM; i++)
	{
		arrayPos[i] = std::min((double)this->getSize()[i]-1, floor(this->getSize()[i]*(i_pos[i]-m_spaceMin[i])/m_spaceSize[i]));
	}
	
	return this->getIndiceFromCoord(arrayPos);
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> Histogram<DIM, POS, VAL>::getPosFromIndice(int i_ind) const
{
	Vector<DIM, POS> pos;
	
	int sz = this->getArraySize();
	int ind = i_ind;
	double ratio;
	
	for(int i=DIM-1; i>=0; i--)
	{
		sz /= this->getSize()[i];
		ratio = (double)(ind/sz)/(double)(this->getSize()[i]-1);
		pos[i] = m_spaceMin[i] + m_spaceSize[i]*ratio;
		ind = ind % sz;
	}
	
	return pos;
}

template<int DIM, typename POS, typename VAL>
VAL Histogram<DIM, POS, VAL>::getFromPosInterpolate(const Vector<DIM, POS>& i_pos)
{
	Vector<DIM, int> coord, coordInt, coordFrac, sign;
	for(int i=0; i<DIM; i++)
	{
		coord[i] = this->getSize()[i]*(i_pos[i]-m_spaceMin[i])/m_spaceSize[i];
		coordInt[i] = std::floor(coord[i]+0.5);
		coordFrac[i] = coord[i]-coordInt[i];
		if(coordFrac[i]!=0) sign[i] = std::floor(coordFrac[i]/std::abs(coordFrac[i]) + 0.5);
		else sign[i] = 0;
		coordFrac[i] = std::abs(coordFrac[i]);
	}

	VAL res = 0;
	int i, j;
	try
	{
		i=coordInt[0]; j=coordInt[1];
		if( 0<=i && i<this->getSize()[0] && 0<=j && j<this->getSize()[1] ) res += (1-coordFrac[0])*(1-coordFrac[1])*this->atCoord(stk::Vector<DIM, int>(i, j));
		i=coordInt[0]+sign[0]; j=coordInt[1];
		if( sign[0]!=0 && 0<=i && i<this->getSize()[0] && 0<=j && j<this->getSize()[1] ) res += (coordFrac[0])*(1-coordFrac[1])*this->atCoord(stk::Vector<DIM, int>(i, j));
		i=coordInt[0]; j=coordInt[1]+sign[1];
		if( sign[1]!=0 && 0<=i && i<this->getSize()[0] && 0<=j && j<this->getSize()[1] ) res += (1-coordFrac[0])*(coordFrac[1])*this->atCoord(stk::Vector<DIM, int>(i, j));
		i=coordInt[0]+sign[0]; j=coordInt[1]+sign[1];
		if( sign[0]!=0 && sign[1]!=0 && 0<=i && i<this->getSize()[0] && 0<=j && j<this->getSize()[1] ) res += (coordFrac[0])*(coordFrac[1])*this->atCoord(stk::Vector<DIM, int>(i, j));
	}
	catch(stk::exception::OutOfRange& e)
	{
		throw e
			.addVar("pos", i_pos)
			.addVar("posMin", m_spaceMin)
			.addVar("posMax", m_spaceMax)
			.addVar("size", this->getSize());
	}
	return res;
}

template<int DIM, typename POS, typename VAL>
const VAL& Histogram<DIM, POS, VAL>::getFromPos(const Vector<DIM, POS>& i_pos) const
{
	int ind = getIndiceFromPos(i_pos);
	
	try
	{
		return this->getFromIndice(ind);
	}
	catch(stk::exception::OutOfRange& e)
	{
		throw e
			.addVar("pos", i_pos)
			.addVar("posMin", m_spaceMin)
			.addVar("posMax", m_spaceMax)
			.addVar("size", this->getSize());
	}
}

template<int DIM, typename POS, typename VAL>
VAL& Histogram<DIM, POS, VAL>::getFromPos(const Vector<DIM, POS>& i_pos)
{
	int ind = getIndiceFromPos(i_pos);
	
	try
	{
		return this->getFromIndice(ind);
	}
	catch(stk::exception::OutOfRange& e)
	{
		throw e
			.addVar("pos", i_pos)
			.addVar("posMin", m_spaceMin)
			.addVar("posMax", m_spaceMax)
			.addVar("size", this->getSize());
	}
}

template<int DIM, typename POS, typename VAL>
VAL& Histogram<DIM, POS, VAL>::get(const Vector<DIM, POS>& i_pos)
{
	return this->getFromPos(i_pos);
}

template<int DIM, typename POS, typename VAL>
const VAL& Histogram<DIM, POS, VAL>::get(const Vector<DIM, POS>& i_pos) const
{
	return this->getFromPos(i_pos);
}
		
template<int DIM, typename POS, typename VAL>
const Vector<DIM, POS>& Histogram<DIM, POS, VAL>::getMinPosition() const
{
	return m_spaceMin;
}

template<int DIM, typename POS, typename VAL>
const Vector<DIM, POS>& Histogram<DIM, POS, VAL>::getMaxPosition() const
{
	return m_spaceMax;
}

template<int DIM, typename POS, typename VAL>
VAL Histogram<DIM, POS, VAL>::integrate() const
{
	return this->getSum() / m_spaceSize.volume();
}
		
template<int DIM, typename POS, typename VAL>
void Histogram<DIM, POS, VAL>::operator+=(const Histogram<DIM, POS, VAL>& i_histo)
{
	if(i_histo.getSize() != this->getSize())
	{
		throw exception::Message("Histogram::operator+= wrong size", STK_DBG_INFO);
	}
	
	for(int i=0; i<this->getArraySize(); i++)
	{
		this->getFromIndice(i) += i_histo.getFromIndice(i);
	}
}
		
template<int DIM, typename POS, typename VAL>
void Histogram<DIM, POS, VAL>::operator-=(const Histogram<DIM, POS, VAL>& i_histo)
{
	if(i_histo.getSize() != this->getSize())
	{
		throw exception::Message("Histogram::operator+= wrong size", STK_DBG_INFO);
	}
	
	for(int i=0; i<this->getArraySize(); i++)
	{
		this->getFromIndice(i) -= i_histo.getFromIndice(i);
	}
}
		
template<int DIM, typename POS, typename VAL>
void Histogram<DIM, POS, VAL>::operator*=(double i_val)
{
	for(int i=0; i<this->getArraySize(); i++)
	{
		this->getFromIndice(i) *= i_val;
	}
}
		
template<int DIM, typename POS, typename VAL>
void Histogram<DIM, POS, VAL>::operator/=(double i_val)
{
	for(int i=0; i<this->getArraySize(); i++)
	{
		this->getFromIndice(i) /= i_val;
	}
}

typedef Histogram<1, double, double> Histogram1dd;
typedef Histogram<1, double, Complexd> Histogram1dc;
typedef Histogram<1, double, int> Histogram1di;
typedef Histogram<2, double, double> Histogram2dd;
typedef Histogram<2, double, Complexd> Histogram2dc;
typedef Histogram<2, double, int> Histogram2di;

}

#endif
