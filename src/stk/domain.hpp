#ifndef __STK_DOMAIN__
#define __STK_DOMAIN__

#include <sstream>
#include <stk/exception.hpp>
#include <stk/vector.hpp>

namespace stk
{

/* DOMAIN *************************************************************/

template<int DIM, typename POS>
class Domain
{
	protected:
		bool m_toroidal;
	
	public:
		Domain() : m_toroidal(false)
			{ }
		
		virtual Vector<DIM, POS> pos(const Vector<DIM, POS>& a) const
		{
			return a;
		}
		
		virtual Vector<DIM, POS> diff(const Vector<DIM, POS>& a, const Vector<DIM, POS>& b) const
		{
			return b-a;
		}
		
		virtual double dist(const Vector<DIM, POS>& a, const Vector<DIM, POS>& b) const
		{
			return diff(a, b).norm();
		}
		
		bool isToroidal() const
		{
			return m_toroidal;
		}
		
		bool toroidal(bool t)
		{
			m_toroidal = t;
		}
		
	//Abstract method
	public:
		virtual Domain<DIM, POS>* clone() const = 0;
		virtual Vector<DIM, POS> referencePoint() const = 0;
		virtual Vector<DIM, POS> centralPoint() const = 0;
		virtual double volume() const = 0;
		
		virtual std::vector< Vector<DIM, POS> > replicationVectors() const = 0;
		
		virtual Vector<DIM, POS> boundingBoxSize() const = 0;
		virtual Vector<DIM, POS> boundingBoxMin() const = 0;
		virtual Vector<DIM, POS> boundingBoxMax() const = 0;
		
		virtual POS distanceToBoundaries(const Vector<DIM, POS>& a) const = 0;
		virtual bool hitTest(const Vector<DIM, POS>& a) const = 0;
		
		virtual std::string toString() const = 0;
};

template<int DIM, typename POS>
std::ostream& operator<<(std::ostream &os, const Domain<DIM, POS>& i_pts)
{
	os << i_pts.toString();
	return os;
}

/* UNIT DOMAIN ********************************************************/

template<int DIM, typename POS>
class UnitDomain : public Domain<DIM, POS>
{
	public:
		UnitDomain()
			{ }
		
		virtual Domain<DIM, POS>* clone() const
		{
			return new UnitDomain<DIM, POS>(static_cast<const UnitDomain<DIM, POS>&>(*this));
		}
		
		virtual Vector<DIM, POS> pos(const Vector<DIM, POS>& a) const
		{
			if(this->isToroidal())
			{
				Vector<DIM, POS> r;
				for(int i=0; i<DIM; i++)
				{
					if(a[i] < 0.0)
					{
						r[i] = 1.0 + a[i] - (double)((int)a[i]);
						if(r[i] == 1.0) r[i] = 0.0;
					}
					else r[i] = a[i] - (double)((int)a[i]);
				}
				return r;
			}
			else return a;
		}
		
		virtual Vector<DIM, POS> diff(const Vector<DIM, POS>& a, const Vector<DIM, POS>& b) const
		{
			if(this->isToroidal())
			{
				Vector<DIM, POS> r = pos(b)-pos(a);
				
				for(int i=0; i<DIM; i++)
				{
					if(r[i] >= 0.5) r[i] -= 1.0;
					else if(r[i] < -0.5) r[i] += 1.0;
				}
				
				return r;
			}
			else return b-a;
		}
		
		virtual POS distanceToBoundaries(const Vector<DIM, POS>& a) const
		{
			if(this->isToroidal())
			{
				Vector<DIM, POS> r;
				
				for(int i=0; i<DIM; i++)
				{
					r[i] = std::min(std::max(a[i], (POS)0.0), (POS)1.0);
				}
				
				return (a-r).norm();
			}
			else return 0;
		}
		
		virtual Vector<DIM, POS> referencePoint() const
		{
			return Vector<DIM, POS>(0.0);
		}
		
		virtual double volume() const
		{
			return 1.0;
		}
		
		virtual Vector<DIM, POS> centralPoint() const
		{
			return Vector<DIM, POS>(0.5);
		}
		
		virtual std::vector< Vector<DIM, POS> > replicationVectors() const
		{
			std::vector< Vector<DIM, POS> > r;
			r.push_back(Vector<DIM, POS>(-1.0, 0.0));
			r.push_back(Vector<DIM, POS>(-1.0, 1.0));
			r.push_back(Vector<DIM, POS>(0.0, 1.0));
			r.push_back(Vector<DIM, POS>(1.0, 1.0));
			r.push_back(Vector<DIM, POS>(1.0, 0.0));
			r.push_back(Vector<DIM, POS>(1.0, -1.0));
			r.push_back(Vector<DIM, POS>(0.0, -1.0));
			r.push_back(Vector<DIM, POS>(-1.0, -1.0));
			return r;
		}
		
		virtual Vector<DIM, POS> boundingBoxSize() const
		{
			return Vector<DIM, POS>(1.0);
		}
		
		virtual Vector<DIM, POS> boundingBoxMin() const
		{
			return Vector<DIM, POS>(0.0);
		}
		
		virtual Vector<DIM, POS> boundingBoxMax() const
		{
			return Vector<DIM, POS>(1.0);
		}
		
		virtual bool hitTest(const Vector<DIM, POS>& a) const
		{
			for(int i=0; i<DIM; i++)
			{
				if(a[i] < 0.0) return false;
				if(a[i] >= 1.0) return false;
			}
			return true;
		}
		
		virtual std::string toString() const
		{
			std::stringstream out;
			out << "UnitDomain(" << (this->isToroidal() ? "toroidal" : "bounded") << ")";
			return out.str();
		}
};

/* RECTANGULAR DOMAIN *************************************************/

template<int DIM, typename POS>
class RectangularDomain : public Domain<DIM, POS>
{
	protected:
		Vector<DIM, POS> m_min;
		Vector<DIM, POS> m_max;
		Vector<DIM, POS> m_size;
	
	public:
		RectangularDomain(const Vector<DIM, POS>& i_min, const Vector<DIM, POS>& i_max)
		{
			m_min = i_min;
			m_max = i_max;
			m_size = m_max - m_min;
		}
		
		virtual Vector<DIM, POS> pos(const Vector<DIM, POS>& a) const
		{
			if(this->isToroidal())
			{
				Vector<DIM, POS> r;
				r = (a-this->m_min)/this->m_size;
				
				for(unsigned int i=0; i<DIM; i++)
				{			
					if(r[i] < 0.0)
					{
						r[i] = 1.0 + r[i] - (double)((int)r[i]);
						if(r[i] == 1.0) r[i] = 0.0;
					}
					else r[i] = r[i] - (double)((int)r[i]);	
				}
				
				r = r*this->m_size + this->m_min;
				return r;
			}
			else return a;
		}
		
		virtual Vector<DIM, POS> diff(const Vector<DIM, POS>& a, const Vector<DIM, POS>& b) const
		{
			if(this->isToroidal())
			{
				Vector<DIM, POS> r = pos(b)-pos(a);
				for(unsigned int i=0; i<DIM; i++)
				{
					if(r[i] >= this->m_size[i]/2.)
					{
						r[i] -= this->m_size[i];
					}
					else if(r[i] < -this->m_size[i]/2.)
					{
						r[i] += this->m_size[i];
					}
				}
				return r;
			}
			else return b-a;
		}
		
		virtual POS distanceToBoundaries(const Vector<DIM, POS>& a) const
		{
			if(this->isToroidal())
			{
				Vector<DIM, POS> r;
				
				for(int i=0; i<DIM; i++)
				{
					r[i] = std::min(std::max(a[i], m_min[i]), m_max[i]);
				}
				
				return (a-r).norm();
			}
			else return 0;
		}
		
		virtual Domain<DIM, POS>* clone() const
		{
			return new RectangularDomain<DIM, POS>(static_cast<const RectangularDomain<DIM, POS>&>(*this));
		}
		
		virtual Vector<DIM, POS> referencePoint() const
		{
			return Vector<DIM, POS>(m_min);
		}
		
		virtual double volume() const
		{
			return m_size.volume();
		}
		
		virtual Vector<DIM, POS> centralPoint() const
		{
			return m_min + m_size/2.0;
		}
		
		virtual std::vector< Vector<DIM, POS> > replicationVectors() const
		{
			std::vector< Vector<DIM, POS> > r;
			r.push_back(Vector<DIM, POS>(-1.0, 0.0)*m_size);
			r.push_back(Vector<DIM, POS>(-1.0, 1.0)*m_size);
			r.push_back(Vector<DIM, POS>(0.0, 1.0)*m_size);
			r.push_back(Vector<DIM, POS>(1.0, 1.0)*m_size);
			r.push_back(Vector<DIM, POS>(1.0, 0.0)*m_size);
			r.push_back(Vector<DIM, POS>(1.0, -1.0)*m_size);
			r.push_back(Vector<DIM, POS>(0.0, -1.0)*m_size);
			r.push_back(Vector<DIM, POS>(-1.0, -1.0)*m_size);
			return r;
		}
		
		virtual Vector<DIM, POS> boundingBoxSize() const
		{
			return Vector<DIM, POS>(m_size);
		}
		
		virtual Vector<DIM, POS> boundingBoxMin() const
		{
			return Vector<DIM, POS>(m_min);
		}
		
		virtual Vector<DIM, POS> boundingBoxMax() const
		{
			return Vector<DIM, POS>(m_max);
		}
		
		virtual bool hitTest(const Vector<DIM, POS>& a) const
		{
			for(int i=0; i<DIM; i++)
			{
				if(a[i] < m_min[i]) return false;
				if(a[i] >= m_max[i]) return false;
			}
			return true;
		}
		
		virtual std::string toString() const
		{
			std::stringstream out;
			out << "RectangularDomain(" << (this->isToroidal() ? "toroidal" : "bounded");
			out << ", ";
			
			for(int i=0; i<DIM; i++)
			{
				out << "[" << m_min[i] << ";"<< m_max[i] << "]";
				if(i<DIM-1) out << " x ";
			}
			
			out << ")";
			return out.str();
		}
};

/* BASE DOMAIN ******************************************************/

template<typename POS>
class BaseDomain : public Domain<2, POS>
{
	protected:
		Vector<2, POS> m_min;
		Vector<2, POS> m_v0;
		Vector<2, POS> m_v1;
		Vector<2, POS> m_center;
		Vector<2, POS> m_boxSize;
		Vector<2, POS> m_boxMin;
		Vector<2, POS> m_boxMax;
		double m_volume;
		
		Vector<2, POS> m_matrix0;
		Vector<2, POS> m_matrix1;
		Vector<2, POS> m_matrixInv0;
		Vector<2, POS> m_matrixInv1;
	
	public:
		BaseDomain(const Vector<2, POS>& i_min, const std::vector< Vector<2, POS> >& i_base)
		{
			m_min = i_min;
			m_v0 = i_base.at(0);
			m_v1 = i_base.at(1);
			
			Vector<2, POS> p0 = m_min + m_v0;
			Vector<2, POS> p1 = m_min + m_v1;
			Vector<2, POS> p2 = m_min + m_v0 + m_v1;
			
			for(int i=0; i<2; i++)
			{
				m_boxMin[i] = std::min(p0[i], p1[i]);
				m_boxMin[i] = std::min(m_boxMin[i], p2[i]);
				m_boxMin[i] = std::min(m_boxMin[i], m_min[i]);
				
				m_boxMax[i] = std::max(p0[i], p1[i]);
				m_boxMax[i] = std::max(m_boxMax[i], p2[i]);
				m_boxMax[i] = std::max(m_boxMax[i], m_min[i]);
			}
			m_boxSize = m_boxMax - m_boxMin;
			
			m_volume = m_v0[0]*m_v1[1] - m_v0[1]*m_v1[0];
			m_center = i_min+(m_v0 + m_v1)/2.0;
			
			m_matrixInv0[0] = m_v0[0];
			m_matrixInv0[1] = m_v1[0];
			m_matrixInv1[0] = m_v0[1];
			m_matrixInv1[1] = m_v1[1];
			double det = m_matrixInv0[0]*m_matrixInv1[1] - m_matrixInv0[1]*m_matrixInv1[0];
			if(std::abs(det) > 1e-30)
			{
				m_matrix0[0] = m_matrixInv1[1]/det;
				m_matrix0[1] = -m_matrixInv0[1]/det;
				m_matrix1[0] = -m_matrixInv1[0]/det;
				m_matrix1[1] = m_matrixInv0[0]/det;
			}
			else
			{
				std::stringstream out;
				out << "wrong base domain definition, ";
				out << "det = " << det << "  in matrix made from (" << m_v0[0] << ", " << m_v0[1] << ") and (" << m_v1[0] << ", " << m_v1[1] << ")";
				
				throw stk::exception::Message(out.str(), STK_DBG_INFO);
			}
		}
		
		virtual Domain<2, POS>* clone() const
		{
			return new BaseDomain<POS>(static_cast<const BaseDomain<POS>&>(*this));
		}
		
		virtual Vector<2, POS> pos(const Vector<2, POS>& a) const
		{
			if(this->isToroidal())
			{
				Vector<2, POS> r;
				Vector<2, POS> tmp;
				r = a-this->m_min;
				
				tmp[0] = (this->m_matrix0*r).total();
				tmp[1] = (this->m_matrix1*r).total();
				
				for(unsigned int i=0; i<2; i++)
				{
					if(tmp[i] < 0.0)
					{
						tmp[i] = 1.0 + tmp[i] - (double)((int)tmp[i]);
						if(tmp[i] == 1.0) tmp[i] = 0.0;
					}
					else tmp[i] = tmp[i] - (double)((int)tmp[i]);
				}
				
				r[0] = (this->m_matrixInv0*tmp).total();
				r[1] = (this->m_matrixInv1*tmp).total();
				r += this->m_min;
				return r;
			}
			else return a;
		}
		
		virtual Vector<2, POS> diff(const Vector<2, POS>& a, const Vector<2, POS>& b) const
		{
			if(this->isToroidal())
			{
				Vector<2, POS> r = pos(b)-pos(a);
				Vector<2, POS> tmp;
				
				tmp[0] = (this->m_matrix0*r).total();
				tmp[1] = (this->m_matrix1*r).total();
				
				tmp[0] += 0.5;
				tmp[1] += 0.5;
				for(unsigned int i=0; i<2; i++)
				{
					if(tmp[i] < 0.0)
					{
						tmp[i] = 1.0 + tmp[i] - (double)((int)tmp[i]);
						if(tmp[i] == 1.0) tmp[i] = 0.0;
					}
					else tmp[i] = tmp[i] - (double)((int)tmp[i]);
				}
				tmp[0] -= 0.5;
				tmp[1] -= 0.5;
				
				r[0] = (this->m_matrixInv0*tmp).total();
				r[1] = (this->m_matrixInv1*tmp).total();
				
				return r;
			}
			else return b-a;
		}
		
		virtual POS distanceToBoundaries(const Vector<2, POS>& a) const
		{
			if(this->isToroidal())
			{
				Vector<2, POS> r;
				Vector<2, POS> tmp;
				r = a-this->m_min;
				
				tmp[0] = (this->m_matrix0*r).total();
				tmp[1] = (this->m_matrix1*r).total();
				
				for(unsigned int i=0; i<2; i++)
				{
					tmp[i] = std::min(std::max(tmp[i], (POS)0), (POS)1);
				}
				
				r[0] = (this->m_matrixInv0*tmp).total();
				r[1] = (this->m_matrixInv1*tmp).total();
				r += this->m_min;
				
				return (a-r).norm();
			}
			else return 0;
		}
		
		virtual Vector<2, POS> referencePoint() const
		{
			return m_min;
		}
		
		virtual double volume() const
		{
			return m_volume;
		}
		
		virtual Vector<2, POS> centralPoint() const
		{
			return m_center;
		}
		
		virtual std::vector< Vector<2, POS> > replicationVectors() const
		{
			std::vector< Vector<2, POS> > r;
			r.push_back(m_v0*-1.0);
			r.push_back(m_v1 - m_v0);
			r.push_back(m_v1);
			r.push_back(m_v0 + m_v1);
			r.push_back(m_v0);
			r.push_back(m_v0 - m_v1);
			r.push_back(m_v1*-1.0);
			r.push_back((m_v0 + m_v1)*-1.0);
			return r;
		}
		
		virtual Vector<2, POS> boundingBoxSize() const
		{
			return m_boxSize;
		}
		
		virtual Vector<2, POS> boundingBoxMin() const
		{
			return m_boxMin;
		}
		
		virtual Vector<2, POS> boundingBoxMax() const
		{
			return m_boxMax;
		}
		
		virtual bool hitTest(const Vector<2, POS>& a) const
		{
			Vector<2, POS> r;
			Vector<2, POS> tmp;
			r = a-this->m_min;
			
			tmp[0] = (m_matrix0*r).total();
			tmp[1] = (m_matrix1*r).total();
			
			for(unsigned int i=0; i<2; i++)
			{
				if(tmp[i] < 0.0) return false;
				if(tmp[i] >= 1.0) return false;
			}
			
			return true;
		}
		
		Vector<2, POS> vector(int i) const
		{
			if(i == 0) return m_v0;
			else if(i == 1) return m_v1;
			else
			{
				std::stringstream out;
				out << "this domain contain only 2 vectors";
				
				throw stk::exception::Message(out.str(), STK_DBG_INFO);
			}
		}
		
		Vector<2, POS> transform(const Vector<2, POS>& a) const
		{
			Vector<2, POS> tmp;
			
			tmp[0] = (m_matrix0*a).total();
			tmp[1] = (m_matrix1*a).total();
			
			return tmp;
		}
		
		Vector<2, POS> inverseTransform(const Vector<2, POS>& a) const
		{
			Vector<2, POS> tmp;
			
			tmp[0] = (m_matrixInv0*a).total();
			tmp[1] = (m_matrixInv1*a).total();
			
			return tmp;
		}
		
		Vector<2, POS> posInBasis(const Vector<2, POS>& a) const
		{
			return transform(a);
		}
		
		virtual std::string toString() const
		{
			std::stringstream out;
			out << "BaseDomain(" << (this->isToroidal() ? "toroidal" : "bounded") << ")";
			return out.str();
		}
};

/* HEXAGONAL DOMAIN *************************************************/

template<typename POS>
class HexagonalDomain : public Domain<2, POS>
{
	protected:
		Vector<2, POS> m_u;
		Vector<2, POS> m_v;
		
		Vector<2, POS> m_x;
		Vector<2, POS> m_y;
		
		POS m_r;
		POS m_R;
		double m_area;
	
	public:
		HexagonalDomain()
		{
			m_u = Vector<2, POS>(1.0, 0.0);
			m_v = Vector<2, POS>(0.5, 0.866025404);
			
			m_x = Vector<2, POS>(1.0, 0.0);
			m_y = Vector<2, POS>(-0.577350269, 1.154700538);
			
			m_r = 0.5;
			m_R = m_r*2.0/std::sqrt(3.0);
			m_area = (m_r*m_r) * 3.0 * std::sqrt(3.0) / 2.0;
		}
		
		virtual Domain<2, POS>* clone() const
		{
			return new HexagonalDomain<POS>(static_cast<const HexagonalDomain<POS>&>(*this));
		}
		
		virtual Vector<2, POS> pos(const Vector<2, POS>& a) const
		{
			if(this->isToroidal())
			{
				Vector<2, POS> cUV;
				cUV[0] = std::floor(a[0]*m_x[0] + a[1]*m_y[0] + 0.5);
				cUV[1] = std::floor(a[0]*m_x[1] + a[1]*m_y[1] + 0.5);
				
				Vector<2, POS> cXY;
				cXY[0] = cUV[0]*m_u[0] + cUV[1]*m_v[0];
				cXY[1] = cUV[0]*m_u[1] + cUV[1]*m_v[1];
				
				Vector<2, POS> cDist = a - cXY;
				POS cNorm = cDist.norm();
				
				if(cNorm < 0.25) return cDist;
				
				Vector<2, POS> fDist = cDist;
				POS fNorm = cNorm;
				
				Vector<2, POS> nDist;
				if(cDist[0] < 0.0)
				{
					nDist = cDist + m_u;
					if(nDist.norm() < fNorm)
					{
						fDist = nDist;
						fNorm = fDist.norm();
					}
					
					nDist = cDist + m_v;
					if(nDist.norm() < fNorm)
					{
						fDist = nDist;
						fNorm = fDist.norm();
					}
				}
				else
				{
					nDist = cDist - m_u;
					if(nDist.norm() < fNorm)
					{
						fDist = nDist;
						fNorm = fDist.norm();
					}
					
					nDist = cDist - m_v;
					if(nDist.norm() < fNorm)
					{
						fDist = nDist;
						fNorm = fDist.norm();
					}
				}
				
				return fDist;
			}
			else return a;
		}
		
		virtual Vector<2, POS> diff(const Vector<2, POS>& a, const Vector<2, POS>& b) const
		{
			if(this->isToroidal())
			{
				return this->pos(b - a);
			}
			else return b-a;
		}
		
		virtual POS distanceToBoundaries(const Vector<2, POS>& a) const
		{
			return 0;
		}
		
		virtual Vector<2, POS> referencePoint() const
		{
			return Vector<2, POS>(0.0, 0.0);
		}
		
		virtual double volume() const
		{
			return m_area;
		}
		
		virtual Vector<2, POS> centralPoint() const
		{
			return Vector<2, POS>(0.0, 0.0);
		}
		
		virtual std::vector< Vector<2, POS> > replicationVectors() const
		{
			std::vector< Vector<2, POS> > r;
			r.push_back(m_u);
			r.push_back(m_v);
			r.push_back(m_v - m_u);
			r.push_back(m_u*-1.0);
			r.push_back(m_v*-1.0);
			r.push_back(m_u - m_v);
			return r;
		}
		
		virtual Vector<2, POS> boundingBoxSize() const
		{
			return Vector<2, POS>(2.0*m_r, 2.0*m_R);
		}
		
		virtual Vector<2, POS> boundingBoxMin() const
		{
			return Vector<2, POS>(-m_r, -m_R);
		}
		
		virtual Vector<2, POS> boundingBoxMax() const
		{
			return Vector<2, POS>(m_r, m_R);
		}
		
		virtual bool hitTest(const Vector<2, POS>& a) const
		{
			POS test;
			
			if(std::abs(dotProduct(m_u, a)) > m_r) return false;
			if(std::abs(dotProduct(m_v, a)) > m_r) return false;
			if(std::abs(dotProduct(m_v-m_u, a)) > m_r) return false;
			
			return true;
		}
		
		const POS& circumcircleRadius() const
		{
			return m_R;
		}
		
		const POS& incircleRadius() const
		{
			return m_r;
		}
		
		virtual std::string toString() const
		{
			std::stringstream out;
			out << "HexagonalDomain(" << (this->isToroidal() ? "toroidal" : "bounded") << ")";
			return out.str();
		}
};

}

#endif
