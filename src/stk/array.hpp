#ifndef __STK_ARRAY__
#define __STK_ARRAY__

#include <vector>
#include <limits>

#include <stk/vector.hpp>
#include <stk/type.hpp>
#include <stk/exception.hpp>

namespace stk
{

/**
 * @class Array
 * Generic implementation of n-dimension array. All data are stored in
 * a 1D array. You can use 1D coordinate system, called "indice", or nD
 * coordinate system, called "coord".
 */
template<int DIM, typename VAL>
class Array
{
	private:
		Vector<DIM, int> m_size;
		int m_arraySize;
		VAL* m_data;
		
	
	public:
		Array();
		Array(const Array<DIM, VAL>& array);
		Array(const Vector<DIM, int>& i_size);
		~Array();

		void resize(const Vector<DIM, int>& i_size);
		
		const Vector<DIM, int>& getSize() const;
		int getArraySize() const;
		
		/**
		 * @deprecated
		 * @see atIndice
		 */
		VAL& getFromIndice(int i_ind);
		
		/**
		 * @deprecated
		 * @see atIndice
		 */
		const VAL& getFromIndice(int i_ind) const;
		
		/**
		 * @deprecated
		 * @see coord2indice
		 */
		int getIndiceFromCoord(const Vector<DIM, int>& i_pos) const;
		
		/**
		 * @deprecated
		 * @see indice2coord
		 */
		Vector<DIM, int> getCoordFromIndice(int i_ind) const;
		
		/**
		 * @deprecated
		 * @see atCoord
		 */
		VAL& getData(const Vector<DIM, int>& i_pos);
		
		/**
		 * @deprecated
		 * @see atCoord
		 */
		const VAL& getData(const Vector<DIM, int>& i_pos) const;
		
		Vector<DIM, int> indice2coord(int i_ind) const;
		int coord2indice(const Vector<DIM, int>& i_coord) const;
		VAL& atIndice(int i_ind);
		const VAL& atIndice(int i_ind) const;
		VAL& atCoord(const Vector<DIM, int>& i_coord);
		const VAL& atCoord(const Vector<DIM, int>& i_coord) const;
		
		VAL& operator[](int i);
		const VAL& operator[](int i) const;
		
		void fill(const VAL& i_value);
		VAL getSum() const;
		VAL getMean() const;
		VAL getMin() const;
		std::vector<VAL> getMin(unsigned int nbVal) const;
		VAL getMax() const;
		std::vector<VAL> getMax(unsigned int nbVal) const;
};

/* CONSTRUCTOR ********************************************************/

template<int DIM, typename VAL>
Array<DIM, VAL>::Array()
{
	m_size = Vector2i(1,1);
	m_arraySize = m_size.volume();
	
	m_data = new VAL[m_arraySize];
}

template<int DIM, typename VAL>
Array<DIM, VAL>::Array(const Array<DIM, VAL>& array)
{
	m_size = array.getSize();
	m_arraySize = m_size.volume();
	
	m_data = new VAL[m_arraySize];

	for(int i=0; i<m_arraySize; i++)
	{
		m_data[i] = array[i];
	}
}

template<int DIM, typename VAL>
Array<DIM, VAL>::Array(const Vector<DIM, int>& i_size)
{
	m_size = i_size;
	m_arraySize = m_size.volume();
	
	m_data = new VAL[m_arraySize];
}

template<int DIM, typename VAL>
Array<DIM, VAL>::~Array()
{
	delete[] m_data;
}

template<int DIM, typename VAL>
void Array<DIM, VAL>::resize(const Vector<DIM, int>& i_size)
{
	int arraySize = i_size.volume();
	VAL* data = new VAL[arraySize];

	delete[] m_data;
	m_data = data;
	m_arraySize = arraySize;
	m_size = i_size;
	
}

/* ACCESSOR ***********************************************************/

template<int DIM, typename VAL>
int Array<DIM, VAL>::coord2indice(const Vector<DIM, int>& i_pos) const
{
	int ind = i_pos[0];
	
	for(int i=0; i<DIM; i++)
	{
		if(i_pos[i] >= m_size[i] || i_pos[i] < 0) throw exception::OutOfRange(i_pos[i], 0, m_size[i], STK_DBG_INFO);
	}
	
	for(int i=1; i<DIM; i++)
	{
		ind += i_pos[i] * m_size[i-1];
	}
	
	return ind;
}

template<int DIM, typename VAL>
Vector<DIM, int> Array<DIM, VAL>::indice2coord(int i_ind) const
{
	Vector<DIM, int> pos;
	
	int sz = this->getArraySize();
	int ind = i_ind;
	double ratio;
	
	for(int i=DIM-1; i>=0; i--)
	{
		sz /= this->getSize()[i];
		pos[i] = ind/sz;
		ind = ind % sz;
	}
	
	return pos;
}

template<int DIM, typename VAL>
VAL& Array<DIM, VAL>::atCoord(const Vector<DIM, int>& i_coord)
{
	int ind = coord2indice(i_coord);
	
	if(ind < 0)
	{
		int dim = -ind;
		throw exception::OutOfRange(i_coord[dim], 0, m_size[dim], STK_DBG_INFO).addVar("coord", i_coord).addVar("dim", dim);
	}
	
	return m_data[ind];
}

template<int DIM, typename VAL>
const VAL& Array<DIM, VAL>::atCoord(const Vector<DIM, int>& i_coord) const
{
	int ind = coord2indice(i_coord);
	
	if(ind < 0)
	{
		int dim = -ind;
		throw exception::OutOfRange(i_coord[dim], 0, m_size[dim], STK_DBG_INFO).addVar("coord", i_coord).addVar("dim", dim);
	}
	
	return m_data[ind];
}

template<int DIM, typename VAL>
VAL& Array<DIM, VAL>::atIndice(int i_ind)
{
	if(i_ind < 0 || i_ind >= m_arraySize)
	{
		throw exception::OutOfRange(i_ind, 0, m_arraySize, STK_DBG_INFO);
	}
	return m_data[i_ind];
}

template<int DIM, typename VAL>
const VAL& Array<DIM, VAL>::atIndice(int i_ind) const
{
	if(i_ind < 0 || i_ind >= m_arraySize)
	{
		throw exception::OutOfRange(i_ind, 0, m_arraySize, STK_DBG_INFO);
	}
	return m_data[i_ind];
}

/* OTHER **************************************************************/

/* DEPRECATED *********************************************************/

template<int DIM, typename VAL>
int Array<DIM, VAL>::getIndiceFromCoord(const Vector<DIM, int>& i_pos) const
{
	return coord2indice(i_pos);
}

template<int DIM, typename VAL>
Vector<DIM, int> Array<DIM, VAL>::getCoordFromIndice(int i_ind) const
{
	return indice2coord(i_ind);
}

template<int DIM, typename VAL>
VAL& Array<DIM, VAL>::getData(const Vector<DIM, int>& i_pos)
{
	return atCoord(i_pos);
}

template<int DIM, typename VAL>
const VAL& Array<DIM, VAL>::getData(const Vector<DIM, int>& i_pos) const
{
	return atCoord(i_pos);
}

template<int DIM, typename VAL>
VAL& Array<DIM, VAL>::getFromIndice(int i_ind)
{
	return atIndice(i_ind);
}

template<int DIM, typename VAL>
const VAL& Array<DIM, VAL>::getFromIndice(int i_ind) const
{
	return atIndice(i_ind);
}

template<int DIM, typename VAL>
const Vector<DIM, int>& Array<DIM, VAL>::getSize() const
{
	return m_size;
}

template<int DIM, typename VAL>
int Array<DIM, VAL>::getArraySize() const
{
	return m_arraySize;
}

template<int DIM, typename VAL>
void Array<DIM, VAL>::fill(const VAL& i_value)
{
	for(int i=0; i<m_arraySize; i++) m_data[i] = i_value;
}

template<int DIM, typename VAL>
VAL& Array<DIM, VAL>::operator[](int i)
{
	return atIndice(i);
}

template<int DIM, typename VAL>
VAL Array<DIM, VAL>::getMin() const
{
	VAL min = this->atIndice(0);
	
	for(int i=1; i<this->getArraySize(); i++)
	{
		if(min > this->atIndice(i))
		{
			min = this->atIndice(i);
		}
	}
	
	return min;
}

template<int DIM, typename VAL>
std::vector<VAL> Array<DIM, VAL>::getMin(unsigned int nbVal) const
{
	std::vector<VAL> min(nbVal, std::numeric_limits<VAL>::max());

	for(unsigned int i=0; i<this->getArraySize(); i++)
	{
		for( int j=nbVal-1; j>=0; j--)
		{
			if(this->atIndice(i) >= min.at(j))
			{
				for( int k=nbVal-1; k>j+1; k--) min.at(k) = min.at(k-1);
				if(j+1<nbVal) min.at(j+1) = this->atIndice(i);
				break;
			}
			else if(j==0)
			{
				for( int k=nbVal-1; k>j; k--) min.at(k) = min.at(k-1);
				min.at(0) = this->atIndice(i);
				break;
			}
		}
	}
	return min;
}

template<int DIM, typename VAL>
VAL Array<DIM, VAL>::getMax() const
{
	VAL max = this->atIndice(0);

	for(int i=1; i<this->getArraySize(); i++)
	{
		if(max < this->atIndice(i))
		{
			max = this->atIndice(i);
		}
	}
	return max;
}

template<int DIM, typename VAL>
std::vector<VAL> Array<DIM, VAL>::getMax(unsigned int nbVal) const
{
	std::vector<VAL> max(nbVal, std::numeric_limits<VAL>::min());

	for(unsigned int i=0; i<this->getArraySize(); i++)
	{
		for( int j=nbVal-1; j>=0; j--)
		{
			if(this->atIndice(i) <= max.at(j))
			{
				for(unsigned int k=nbVal-1; k>j+1; k--) max.at(k) = max.at(k-1);
				if(j+1<nbVal) max.at(j+1) = this->atIndice(i);
				break;
			}
			else if(j==0)
			{
				for(unsigned int k=nbVal-1; k>j; k--) max.at(k) = max.at(k-1);
				max.at(0) = this->atIndice(i);
				break;
			}
		}
	}
	return max;
}

template<int DIM, typename VAL>
VAL Array<DIM, VAL>::getSum() const
{
	VAL sum = 0;
	
	for(int i=0; i<this->getArraySize(); i++)
	{
		sum += this->atIndice(i);
	}
	
	return sum;
}

template<int DIM, typename VAL>
VAL Array<DIM, VAL>::getMean() const
{
	return this->getSum()/this->getArraySize();
}

template<int DIM, typename VAL>
const VAL& Array<DIM, VAL>::operator[](int i) const
{
	return atIndice(i);
}

template<typename VAL>
std::ostream& operator<<(std::ostream &os, const Array<1, VAL>& i_array)
{
		for(int i=0; i<i_array.getSize()[0]; i++)
		{
			os << i_array.getData(Vector2i(i)) << std::endl;
		}
	
	return os;
}

template<typename VAL>
std::ostream& operator<<(std::ostream &os, const Array<2, VAL>& i_array)
{
	for(int j=0; j<i_array.getSize()[1]; j++)
	{
		for(int i=0; i<i_array.getSize()[0]; i++)
		{
			os << i_array.getData(Vector2i(i, j)) << "\t";
		}
		os << std::endl;
	}
	
	return os;
}

typedef Array<2, int> Array2i;
typedef Array<2, double> Array2d;
typedef Array<2, Complexd> Array2c;

}

#endif
