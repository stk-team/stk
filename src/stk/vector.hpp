#ifndef __STK_VECTOR__
#define __STK_VECTOR__

#include <iostream>
#include <cmath>
#include <cstdlib>

namespace stk
{

/**
 * @class Vector src/stk/vector.hpp
 * @brief Vector/Point representation of dimension DIM and with a precision type T
*/
template<int DIM, typename T>
class Vector
{
protected:
	T m_vals[DIM]; /**< \brief array of vector values*/

public:
/**
 * @brief Default constructor
 */
	Vector<DIM, T>()
	{
		for(unsigned int i=0; i<DIM; i++)
			m_vals[i] = 0;
	}

/**
 * @brief Constructor
 * @param i_vect vector to copy from
 */
	template<int DIM2>
	Vector<DIM, T>(const Vector<DIM2, T>& i_vect)
	{
		for(unsigned int i=0; i< (DIM < DIM2 ? DIM : DIM2) ; i++)
			m_vals[i] = i_vect[i];
		if(DIM > DIM2)
		{
			for(unsigned int i=DIM2; i<DIM ; i++)
				m_vals[i] = 0;
		}
	}

/**
 * @brief Constructor
 * @param i_vect vector of dimension DIM-1 to copy from
 * @param i_val val for the last dimension
 */
	Vector<DIM, T>(const Vector<DIM-1, T>& i_vect,const T& i_val)
	{
		for(unsigned int i=0; i<DIM-1; i++)
			m_vals[i] = i_vect[i];
		m_vals[DIM-1] = i_val;
	}

/**
 * @brief Constructor
 * @param i_vect vector of dimension DIM+1 to copy from
 * the copy erase the extra dimension value
 */
	Vector<DIM, T>(const Vector<DIM+1, T>& i_vect)
	{
		for(unsigned int i=0; i<DIM; i++)
			m_vals[i] = i_vect[i];
	}

/**
 * @brief Constructor
 * @param i_vect vector of dimension DIM+2 to copy from
 * the copy erase the extras dimensions values
 */
	Vector<DIM, T>(const Vector<DIM+2, T>& i_vect)
	{
		for(unsigned int i=0; i<DIM; i++)
			m_vals[i] = i_vect[i];
	}

/**
 * @brief Constructor
 * @param i_default_value default value for all component of the vector
 */
	Vector<DIM, T>(const T& i_default_value)
	{
		for(unsigned int i=0; i<DIM; i++)
			m_vals[i] = i_default_value;
	}

/**
 * @brief Constructor
 * @param i_v0 value for first vector component
 * @param i_v1 value for second vector component
 */
	Vector<DIM, T>(const T& i_v0, const T& i_v1)
	{
		if(DIM > 0) m_vals[0] = i_v0;
		if(DIM > 1) m_vals[1] = i_v1;
	}

/**
 * @brief Constructor
 * @param i_v0 value for first vector component
 * @param i_v1 value for second vector component
 * @param i_v2 value for third vector component
 */
	Vector<DIM, T>(const T& i_v0, const T& i_v1, const T& i_v2)
	{
		if(DIM > 0) m_vals[0] = i_v0;
		if(DIM > 1) m_vals[1] = i_v1;
		if(DIM > 2) m_vals[2] = i_v2;
	}

/**
 * @brief Constructor
 * @param i_v0 value for first vector component
 * @param i_v1 value for second vector component
 * @param i_v2 value for third vector component
 * @param i_v3 value for fourth vector component
 */
	Vector<DIM, T>(const T& i_v0, const T& i_v1, const T& i_v2, const T& i_v3)
	{
		if(DIM > 0) m_vals[0] = i_v0;
		if(DIM > 1) m_vals[1] = i_v1;
		if(DIM > 2) m_vals[2] = i_v2;
		if(DIM > 3) m_vals[3] = i_v3;
	}

/**
 * @brief Accessor to first vector component
 */
	T& x()
		{ return m_vals[0]; }

/**
 * @brief Accessor to second vector component
 */
	T& y()
		{ return m_vals[1]; }

/**
 * @brief Accessor to third vector component
 */
	T& z()
		{ return m_vals[2]; }

/**
 * @brief Accessor to first vector component (const return version)
 */
	const T& x() const
		{ return m_vals[0]; }

/**
 * @brief Accessor to second vector component (const return version)
 */
	const T& y() const
		{ return m_vals[1]; }

/**
 * @brief Accessor to third vector component (const return version)
 */
	const T& z() const
		{ return m_vals[2]; }

/**
 * @brief Mutator for first vector component
 * @param i_val value for first vector component
 */
	void x(const T& i_val)
		{ m_vals[0] = i_val; }

/**
 * @brief Mutator for second vector component
 * @param i_val value for second vector component
 */
	void y(const T& i_val)
		{ m_vals[1] = i_val; }

/**
 * @brief Mutator for third vector component
 * @param i_val value for third vector component
 */
	void z(const T& i_val)
		{ m_vals[2] = i_val; }

/**
 * @brief Mutator for n-th vector component
 * @param i_index array index to access component
 * @param i_val value for n-th vector component
 */
	void set(const unsigned int& i_index,const T& i_val)
		{ m_vals[i_index] = i_val; }

/**
 * @brief Accessor to n-th vector component
 * @param i_index array index to access component
 */
	T& get(const unsigned int& index)
		{ return m_vals[index]; }

/**
 * @brief Accessor to n-th vector component (const return version)
 * @param i_index array index to access component
 */
	const T& get(const unsigned int& index) const
		{ return m_vals[index]; }

/**
 * @brief Accessing operator [] to n-th vector component
 * @param i_index array index to access component
 */
	T& operator[](const unsigned int& i_index)
		{ return get(i_index); }

/**
 * @brief Accessing operator [] to n-th vector component (const return version)
 * @param i_index array index to access component
 */
	const T& operator[](const unsigned int& i_index) const
		{ return get(i_index); }

/**
 * @brief Get the maximum value into a vector
 */
	T getMax() const
	{
		T m = m_vals[0];
		for(int i=1; i<DIM; i++)
		{
			if(m_vals[i] > m) m = m_vals[i];
		}
		return m;
	}

/**
 * @brief Get the minimum value into a vector
 */
	T getMin() const
	{
		T m = m_vals[0];
		for(int i=1; i<DIM; i++)
		{
			if(m_vals[i] < m) m = m_vals[i];
		}
		return m;
	}


	Vector<DIM, T> operator+(const T& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] + b;
		return r;
	}

	Vector<DIM, T> operator-(const T& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] - b;
		return r;
	}

	Vector<DIM, T> operator+(const Vector<DIM, T>& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] + b[i];
		return r;
	}

	Vector<DIM, T> operator-(const Vector<DIM, T>& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] - b[i];
		return r;
	}

	Vector<DIM, T> operator*(const Vector<DIM, T>& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] * b[i];
		return r;
	}

	Vector<DIM, T> operator/(const Vector<DIM, T>& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] / b[i];
		return r;
	}

	Vector<DIM, T> operator*(const T& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] * b;
		return r;
	}

	Vector<DIM, T> operator/(const T& b) const
	{
		Vector<DIM, T> r;
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i] / b;
		return r;
	}

	bool operator<(const T& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] >= b)
				return false;
		return true;
	}

	bool operator>(const T& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] <= b)
				return false;
		return true;
	}

	bool operator<=(const T& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] > b)
				return false;
		return true;
	}

	bool operator>=(const T& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] < b)
				return false;
		return true;
	}

	bool operator<(const Vector<DIM,T>& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
		{
			if(m_vals[i] < b[i]) return true;
			else if(m_vals[i] > b[i]) return false;
		}
		return false;
	}

	bool operator>(const Vector<DIM,T>& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
		{
			if(m_vals[i] > b[i]) return true;
			else if(m_vals[i] < b[i]) return false;
		}
		return false;
	}

	bool operator<=(const Vector<DIM,T>& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] > b[i])
				return false;
		return true;
	}

	bool operator>=(const Vector<DIM,T>& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] < b[i])
				return false;
		return true;
	}

	bool operator==(const Vector<DIM,T>& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] != b[i])
				return false;
		return true;
	}

	bool operator!=(const Vector<DIM,T>& b) const
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] != b[i])
				return true;
		return false;
	}

	void operator+=(const Vector<DIM, T>& b)
		{ for(unsigned int i=0; i<DIM; i++) m_vals[i] += b[i]; }

	void operator-=(const Vector<DIM, T>& b)
		{ for(unsigned int i=0; i<DIM; i++) m_vals[i] -= b[i]; }

	void operator*=(const Vector<DIM, T>& b)
		{ for(unsigned int i=0; i<DIM; i++) m_vals[i] *= b[i]; }

	void operator/=(const Vector<DIM, T>& b)
		{ for(unsigned int i=0; i<DIM; i++) m_vals[i] /= b[i]; }

	void operator*=(const T& b)
		{ for(unsigned int i=0; i<DIM; i++) m_vals[i] *= b; }

	void operator/=(const T& b)
		{ for(unsigned int i=0; i<DIM; i++) m_vals[i] /= b; }

	bool operator!=(const Vector<DIM, T>& b)
	{
		for(unsigned int i=0; i<DIM; i++)
			if(m_vals[i] != b[i])
				return true;
		return false;
	}

/**
 * @brief compute the L-norm of the vector
 * @param i_space L-space dimension in which compute norm.\n
 *        0 : Infinity/maximum/Chebyshev norm\n
 *        Default = 2 (euclidian space)
 * @return norm of the vector
 */
	T norm(const unsigned int& i_space=2) const
	{
		T r = 0;
		if(i_space==0)
		{
			for(unsigned int i=0; i<DIM; i++)
				r = std::max(r, (T)std::abs(m_vals[i]));
		}
		else
		{
			for(unsigned int i=0; i<DIM; i++)
				r += pow(std::abs(m_vals[i]), i_space);
			r = pow(r, 1./i_space);
		}
		return r;
	}

/**
 * @brief normalize vector
 * @param i_space L-space dimension in which compute norm for normalization.\n Default = 2 (euclidian space)
 * @return normalized vector
 */
	Vector<DIM, T> normalize(const unsigned int& i_space=2) const
	{
		Vector<DIM, T> r;
		double n = norm(i_space);
		for(unsigned int i=0; i<DIM; i++) r[i] = m_vals[i]/n;
		return r;
	}

/**
 * @brief test if vector/point is inside an hypercube
 * @param i_min minimum bound of hypercube
 * @param i_max maximum bound of hypercube
 * @return True/False if inside/not inside the hypercube
 */
	bool insideHypercube(const Vector<DIM, T>& i_min, const Vector<DIM, T>& i_max) const
	{
		for(unsigned int i=0; i<DIM; i++)
		{
			if(m_vals[i] < i_min[i] || m_vals[i] > i_max[i])
				return false;
		}
		return true;
	}

/**
 * @brief compute distance from another vector/point
 * @param i_vect vector to witch compute the distance
 * @param i_space L-space dimension in which compute the distance. Default = 2 (euclidian space)
 * @return distance between vectors
 */
	double dist(const Vector<DIM, T>& i_vect, const unsigned int& i_space=2)
	{
		Vector<DIM, T> v = *this-i_vect;
		return v.norm(i_space);
	}

/**
 * @brief compute sum of all vector/point component
 * @return total of component
 */
	T total()
	{
		T r = 0;
		for(unsigned int i=0; i<DIM; i++)
			r += m_vals[i];
		return r;
	}

/**
 * @brief compute volume of hypercube from origin to the vector/point
 * @return volume of the hypercube
 */
	T volume() const
	{
		T res = 1;
		for(unsigned int i=0; i<DIM; i++)
			res *= m_vals[i];
		return res;
	}

	Vector<DIM, T>& operator=(const T& i_val)
	{
		m_vals[0] = i_val;
		return *this;
	}
};

/**
 * @brief compute the dot-product of to vector of equal dimension
 * @param i_V1 first vector
 * @param i_V2 second vector
 * @return scalar result of dot-product
 */
template<int DIM, typename T>
T dotProduct(const Vector<DIM, T>& i_V1, const Vector<DIM, T>& i_V2)
{
	T r = 0;
	for(unsigned int i=0; i<DIM; i++)
		r += i_V1[i] * i_V2[i];
	return r;
}

/**
 * @brief generate a random vector
 * @return the random vector
 * @todo how method reconize DIM and T?
 */
template<int DIM, typename T>
Vector<DIM, T> getRandVect()
{
	srand48(time(NULL));
	Vector<DIM, T> v;
	for(unsigned int i=0; i<DIM; i++)
		v[i] = drand48();
	return v;
}

template<class T>
Vector<3, T> convertVector(const Vector<4, T>& b)
{
	return Vector<3, T>(b[0], b[1], b[2])/b[3];
}

template<class T>
Vector<2, T> vectorFromEuler(double angle, double r)
{
	return Vector<2, T>(cos(angle)*r, sin(angle)*r);
}

template<int DIM, typename T>
std::ostream& operator<<(std::ostream &os, const Vector<DIM, T>& i_vect)
{
	if(DIM>0)
	{
		os << i_vect[0];
		for(unsigned int i=1; i<DIM; i++)
			os << '\t' << i_vect[i];
	}
	return os;
}

template<int DIM, typename T>
std::istream& operator>>(std::istream &os, Vector<DIM, T>& i_vect)
{
	for(unsigned int i=0; i<DIM; i++)
		os >> i_vect[i];
	return os;
}

typedef Vector<2, float> Vector2f;
typedef Vector<3, float> Vector3f;
typedef Vector<4, float> Vector4f;

typedef Vector<2, double> Vector2d;
typedef Vector<3, double> Vector3d;
typedef Vector<4, double> Vector4d;

typedef Vector<2, int> Vector2i;
typedef Vector<3, int> Vector3i;
typedef Vector<4, int> Vector4i;

typedef Vector<2, unsigned int> Vector2u;
typedef Vector<3, unsigned int> Vector3u;
typedef Vector<4, unsigned int> Vector4u;

typedef Vector<2, char> Vector2c;
typedef Vector<3, char> Vector3c;
typedef Vector<4, char> Vector4c;

}

#endif
