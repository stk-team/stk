#ifndef __STK_POINTSET__
#define __STK_POINTSET__

#include <vector>

#include <stk/point.hpp>
#include <stk/vector.hpp>
#include <stk/type.hpp>
#include <stk/histogram.hpp>
#include <stk/function/generic.hpp>
#include <stk/domain.hpp>

namespace stk
{

/**
 * @class PointSet
 * @brief Set of point represented by a std::vector
 */
template<int DIM, typename POS, typename VAL>
class PointSet
	: public std::vector< Point<DIM, POS, VAL> >
{
	private:
		//PointSet 2.0
		Domain<DIM, POS>* m_domain;
		bool m_sharedDomain;
		
		//PointSet 1.0
		VAL m_defaultVal;
	
	public:
		PointSet();
		PointSet(const PointSet<DIM, POS, VAL>& i_pts);
		
		void setParams(const PointSet<DIM, POS, VAL>& pts);
		
		/**
		 * Check if the pointset is toroidal
		 * @return Return true if toroidal, else false
		 * @see toroidal()
		 */
		bool isToroidal() const;
		
		void toroidal(bool t);

		/**
		 * Remove the toroidal periodic part from a vector
		 * @param i_vect vector to "simplify"
		 * @return Return this vector
		 */
		void minTorVect(Vector<DIM, POS>& i_vect) const;

		/**
		 * Remove the toroidal periodic part from a point
		 * @param i_pt point to "simplify"
		 * @see torCrop()
		 */
		void normalize(Vector<DIM, POS>& i_pt) const;
		Vector<DIM, POS> normalize(const Vector<DIM, POS>& i_pt) const;

		/**
		 * Remove the periodic part from all points of toroidal pointset
		 */
		void normalize();

		/**
		 * Compute the difference between two points with respect to toroidal metric if set
		 * @return vector
		 */
		Vector<DIM, POS> diff(const Vector<DIM, POS>& i_v1, const Vector<DIM, POS>& i_v2) const;

		/**
		 * Compute the distance (norm) between two points with respect to toroidal metric if set
		 * @return norm
		 */
		POS dist(const Vector<DIM, POS>& i_v1, const Vector<DIM, POS>& i_v2) const;
		POS dist(const Point<DIM, POS, VAL>& i_p1, const Point<DIM, POS, VAL>& i_p2) const;
		POS dist(unsigned int i_p1, unsigned int i_p2) const;
		POS dist(unsigned int i_p1, const Point<DIM, POS, VAL>& i_p2) const;
		POS dist(const Point<DIM, POS, VAL>& i_p1, unsigned int i_p2) const;
		POS dist(unsigned int i_p1, const Vector<DIM, POS>& i_p2) const;
		POS dist(const Vector<DIM, POS>& i_p1, unsigned int i_p2) const;

		/**
		 * Get the default value for this pointset (set to 1.0 by
		 * default)
		 * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
		 */
		VAL getDefaultVal() const;

		/**
		 * Get center position in space for all dimension
		 * @see getMinSpace
		 * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
		 */
		Vector<DIM, POS> centralPoint() const;

		/**
		 * Set maximum and minimum position in all dimension
		 * @param i_min A vector with minimum position for all dimension
		 * @param i_max A vector with maximum position for all dimension
		 * @see getMaxSpace
		 * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
		 */
		void setBoundaries(
				const Vector<DIM, POS>& i_min,
				const Vector<DIM, POS>& i_max);

		/**
		 * Get sum of all values
		 * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
		 */
		VAL integrate() const;

		//~ void applyFunction(const Function1d* func);
		//~ void applyInverseWindow(const Function1d* func);
		//~ void setMagnitude(const Function1d* func);
		//~ void diffMagnitude(const Function1d* func);

		void operator+=(const VAL& i_val);
		void operator-=(const VAL& i_val);
		void operator*=(const VAL& i_val);
		void operator/=(const VAL& i_val);
		
		void operator*=(const function::GenericFunction<DIM, POS, VAL>& FUNC);
		
		/* PointSet 2.0 ***********************************************/
		
		/* Statistics *************************************************/
	public:
		double meanDistance() const;
		
		/* Domain *****************************************************/
	private:
		void deleteDomain();
	
	public:
		bool isSharedDomain() const;
		const Domain<DIM, POS>* domain() const;
		Domain<DIM, POS>* domain();
		void sharedDomain(Domain<DIM, POS>* domain);
		void embedDomain(Domain<DIM, POS>* i_domain);
		void unitDomain();
		void unitToroidalDomain();
		void rectangularDomain(const Vector<DIM, POS>& i_min, const Vector<DIM, POS>& i_max);
		void rectangularToroidalDomain(const Vector<DIM, POS>& i_min, const Vector<DIM, POS>& i_max);
		void baseDomain(const Vector<DIM, POS>& i_min, const std::vector< Vector<DIM, POS> >& i_base);
		void baseToroidalDomain(const Vector<DIM, POS>& i_min, const std::vector< Vector<DIM, POS> >& i_base);
		void defaultDomain();
		
		bool hitDomainTest(const Vector<DIM, POS>& a) const;

		/**
		 * Get minimum position of the space bounding box
		 * @see getMaxSpace
		 * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
		 */
		Vector<DIM, POS> boundingBoxMin() const;

		/**
		 * Get maximum position of the space bounding box
		 * @see getMinSpace
		 * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
		 */
		Vector<DIM, POS> boundingBoxMax() const;

		/**
		 * Get size of the space bounding box
		 * @see getSpaceSize
		 * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
		 */
		Vector<DIM, POS> boundingBoxSize() const;
		
		/* Replication ************************************************/
	public:
		class replicationIterator
		{
			private:
				std::vector< Vector<DIM, POS> > m_shift;
				int m_shiftCur;
				const PointSet<DIM, POS, VAL>* m_pts;
				int m_ptsCur;
				POS m_distMax;
				bool m_distMaxEnabled;
				Vector<DIM, POS> m_pos;
				bool m_keepOriginal;
			
				void getVector()
				{
					while(m_shiftCur >= 0)
					{
						if(m_ptsCur >= m_pts->size())
						{
							m_shiftCur++;
							m_ptsCur = 0;
						}
						
						if(m_shiftCur >= m_shift.size())
						{
							m_ptsCur = -1; //End()
							m_shiftCur = -1; //End()
						}
						
						if(m_ptsCur >= 0)
						{
							m_pos = m_pts->at(m_ptsCur).pos() + m_shift.at(m_shiftCur);
							
							if(m_distMaxEnabled)
							{
								if(m_pts->domain()->distanceToBoundaries(m_pos) < m_distMax) break;
							}
							else break;
							
							m_ptsCur++;
						}
					}
				}
			
			public:
				replicationIterator()
				{
					m_keepOriginal = false;
					m_shiftCur = 0;
					m_ptsCur = 0;
					m_pts = NULL;
				}
				
				replicationIterator(const PointSet<DIM, POS, VAL>* i_pts)
				{
					m_shiftCur = -1;
					m_ptsCur = -1;
					m_pts = i_pts;
				}
				
				replicationIterator(const PointSet<DIM, POS, VAL>* i_pts, int i_cursor, bool keepOriginal)
				{
					m_shiftCur = 0;
					m_ptsCur = 0;
					m_pts = i_pts;
					m_keepOriginal = keepOriginal;
					m_distMaxEnabled = false;
					
					if(m_keepOriginal)
					{
						m_shift.push_back(Vector<DIM, POS>(0.0));
						
						std::vector< Vector<DIM, POS> > shift;
						shift = m_pts->domain()->replicationVectors();
						for(int i=0; i<shift.size(); i++)
						{
							m_shift.push_back(shift[i]);
						}
					}
					else
					{
						m_shift = m_pts->domain()->replicationVectors();
					}
					
					getVector();
				}
				
				replicationIterator(const PointSet<DIM, POS, VAL>* i_pts, int i_cursor, bool keepOriginal, POS distMax)
				{
					m_shiftCur = 0;
					m_ptsCur = 0;
					m_pts = i_pts;
					m_keepOriginal = keepOriginal;
					m_distMax = distMax;
					m_distMaxEnabled = true;
					
					if(m_keepOriginal)
					{
						m_shift.push_back(Vector<DIM, POS>(0.0));
						
						std::vector< Vector<DIM, POS> > shift;
						shift = m_pts->domain()->replicationVectors();
						for(int i=0; i<shift.size(); i++)
						{
							m_shift.push_back(shift[i]);
						}
					}
					else
					{
						m_shift = m_pts->domain()->replicationVectors();
					}
					
					getVector();
				}
			
				int id() const
				{
					return m_ptsCur;
				}
				
				Vector<DIM, POS> pos() const
				{
					return m_pos;
				}
			
				replicationIterator& operator++()
				{
					m_ptsCur++;
					getVector();
					
					return *this;
				}
			
				replicationIterator& operator++(int)
				{
					m_ptsCur++;
					getVector();
					return *this;
				}
				
				bool operator!=(const replicationIterator& i_iter) const
				{
					return (i_iter.m_shiftCur != m_shiftCur || i_iter.m_ptsCur != m_ptsCur || i_iter.m_pts != m_pts);
				}
				
				bool operator==(const replicationIterator& i_iter) const
				{
					return (i_iter.m_shiftCur != m_shiftCur && i_iter.m_ptsCur == m_ptsCur && i_iter.m_pts != m_pts);
				}
		};
	
	public:
		replicationIterator replicationBegin(bool keepOriginal) const;
		replicationIterator replicationBegin(bool keepOriginal,	const POS& distMax) const;
		replicationIterator replicationEnd() const;
};




/*
 * ********************************************************************
 * IMPLEMENTATION
 * ********************************************************************
 */

template<int DIM, typename POS, typename VAL>
PointSet<DIM, POS, VAL>::PointSet()
	: m_domain(NULL), m_sharedDomain(true)
{
	m_defaultVal = 1;
	
	defaultDomain();
}

template<int DIM, typename POS, typename VAL>
PointSet<DIM, POS, VAL>::PointSet(const PointSet<DIM, POS, VAL>& i_pts)
	: m_domain(NULL), m_sharedDomain(true), std::vector<Point<DIM, POS, VAL> >(i_pts)
{
	m_defaultVal = 1;
	
	//Copy domain
	if(i_pts.isSharedDomain())
	{
		sharedDomain(const_cast< Domain<DIM, POS>* >(i_pts.domain()));
	}
	else
	{
		embedDomain(i_pts.domain()->clone());
	}
}

template<int DIM, typename POS, typename VAL>
bool PointSet<DIM, POS, VAL>::isToroidal() const
{
	return m_domain->isToroidal();
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::toroidal(bool t)
{
	m_domain->toroidal(t);
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::normalize(const Vector<DIM, POS>& i_pt) const
{
	return m_domain->pos(i_pt);
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::normalize()
{
	for(unsigned int i=0; i<this->size(); i++)
		this->at(i).pos() = m_domain->pos(this->at(i).pos());
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::diff(const Vector<DIM, POS>& i_v1, const Vector<DIM, POS>& i_v2) const
{
	return m_domain->diff(i_v1, i_v2);
}

template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Vector<DIM, POS>& i_v1, const Vector<DIM, POS>& i_v2) const
{
	return m_domain->dist(i_v1, i_v2);
}

template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Point<DIM, POS, VAL>& i_p1, const Point<DIM, POS, VAL>& i_p2) const
{
	return m_domain->dist(i_p1.pos(), i_p2.pos());
}

template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(unsigned int i_p1, unsigned int i_p2) const
{
	return m_domain->dist(this->at(i_p1).pos(), this->at(i_p2).pos());
}

template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(unsigned int i_p1, const Point<DIM, POS, VAL>& i_p2) const
{
	return m_domain->dist(this->at(i_p1).pos(), i_p2.pos());
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Point<DIM, POS, VAL>& i_p1, unsigned int i_p2) const
{
	return m_domain->dist(i_p1.pos(), this->at(i_p2).pos());
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(unsigned int i_p1, const Vector<DIM, POS>& i_p2) const
{
	return m_domain->dist(this->at(i_p1).pos(), i_p2);
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Vector<DIM, POS>& i_p1, unsigned int i_p2) const
{
	return m_domain->dist(i_p1, this->at(i_p2).pos());
}

template<int DIM, typename POS, typename VAL>
VAL PointSet<DIM, POS, VAL>::getDefaultVal() const
{
	return m_defaultVal;
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::centralPoint() const
{
	return m_domain->centralPoint();
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::boundingBoxMin() const
{
	return m_domain->boundingBoxMin();
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::boundingBoxMax() const
{
	return m_domain->boundingBoxMax();
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::boundingBoxSize() const
{
	return m_domain->boundingBoxSize();
}

template<int DIM, typename POS, typename VAL>
bool PointSet<DIM, POS, VAL>::hitDomainTest(const Vector<DIM, POS>& a) const
{
	return m_domain->hitTest(a);
}

template<int DIM, typename POS, typename VAL>
VAL PointSet<DIM, POS, VAL>::integrate() const
{
	VAL res = 0;
	for(int i=0; i<this->size(); i++)
	{
		res += this->at(i).val();
	}
	return res;
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::operator+=(const VAL& i_val)
{
	for(int i=0; i<this->size(); i++)
	{
		this->at(i).val() += i_val;
	}
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::operator-=(const VAL& i_val)
{
	for(int i=0; i<this->size(); i++)
	{
		this->at(i).val() -= i_val;
	}
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::operator*=(const VAL& i_val)
{
	for(int i=0; i<this->size(); i++)
	{
		this->at(i).val() *= i_val;
	}
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::operator/=(const VAL& i_val)
{
	for(int i=0; i<this->size(); i++)
	{
		this->at(i).val() /= i_val;
	}
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::operator*=(const function::GenericFunction<DIM, POS, VAL>& func)
{
	for(int i=0; i<this->size(); i++)
	{
		this->at(i).val() *= func(this->at(i).pos());
	}
}

//~ template<int DIM, typename POS, typename VAL>
//~ void PointSet<DIM, POS, VAL>::applyWindow(const Function1d* func)
//~ {
	//~ Vector<N, P> center = (m_pMax + m_pMin)/2;
	//~ Vector<N, P> size = (m_pMax - m_pMin)/2;
	//~ double fact;
	//~ 
	//~ for(int i=0; i<this->size(); i++)
	//~ {
		//~ fact = func->eval(((this->at(i).pos - center)/size).norm());
		//~ this->at(i).val *= fact;
	//~ }
//~ }
//~ 
//~ template<int DIM, typename POS, typename VAL>
//~ void PointSet<DIM, POS, VAL>::applyInverseWindow(const Function1d* func)
//~ {
	//~ Vector<N, P> center = (m_pMax + m_pMin)/2;
	//~ Vector<N, P> size = (m_pMax - m_pMin)/2;
	//~ double fact;
	//~ 
	//~ for(int i=0; i<this->size(); i++)
	//~ {
		//~ fact = func->eval(((this->at(i).pos - center)/size).norm());
		//~ this->at(i).val /= fact;
	//~ }
//~ }
//~ 
//~ template<int DIM, typename POS, typename VAL>
//~ void PointSet<DIM, POS, VAL>::diffMagnitude(const Function1d* func)
//~ {
	//~ double fact;
	//~ for(int i=0; i<this->size(); i++)
	//~ {
		//~ fact = std::abs(this->at(i).val) - func->eval(this->at(i).pos.norm());
		//~ this->at(i).val = std::polar(fact,arg(this->at(i).val));
	//~ }
//~ }
//~ 
//~ template<int DIM, typename POS, typename VAL>
//~ void PointSet<DIM, POS, VAL>::setMagnitude(const Function1d* func)
//~ {
	//~ double fact;
	//~ for(int i=0; i<this->size(); i++)
	//~ {
		//~ fact = func->eval(this->at(i).pos.norm());
		//~ this->at(i).val = std::polar(fact,arg(this->at(i).val));
	//~ }
//~ }

/* PointSet 2.0 - Domain **********************************************/

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::deleteDomain()
{
	if(!isSharedDomain() && m_domain != NULL) delete m_domain;
}

template<int DIM, typename POS, typename VAL>
bool PointSet<DIM, POS, VAL>::isSharedDomain() const
{
	return m_sharedDomain;
}

template<int DIM, typename POS, typename VAL>
const Domain<DIM, POS>* PointSet<DIM, POS, VAL>::domain() const
{
	return m_domain;
}

template<int DIM, typename POS, typename VAL>
Domain<DIM, POS>* PointSet<DIM, POS, VAL>::domain()
{
	return m_domain;
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::sharedDomain(Domain<DIM, POS>* domain)
{
	deleteDomain();
	m_sharedDomain = true;
	m_domain = domain;
	normalize();
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::embedDomain(Domain<DIM, POS>* i_domain)
{
	deleteDomain();
	m_sharedDomain = false;
	m_domain = i_domain;
	normalize();
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::unitDomain()
{
	UnitDomain<DIM, POS>* domain = new UnitDomain<DIM, POS>();
	domain->toroidal(false);
	embedDomain(domain);
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::unitToroidalDomain()
{
	UnitDomain<DIM, POS>* domain = new UnitDomain<DIM, POS>();
	domain->toroidal(true);
	embedDomain(domain);
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::rectangularDomain(const Vector<DIM, POS>& i_min, const Vector<DIM, POS>& i_max)
{
	RectangularDomain<DIM, POS>* domain = new RectangularDomain<DIM, POS>(i_min, i_max);
	domain->toroidal(false);
	embedDomain(domain);
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::rectangularToroidalDomain(const Vector<DIM, POS>& i_min, const Vector<DIM, POS>& i_max)
{
	RectangularDomain<DIM, POS>* domain = new RectangularDomain<DIM, POS>(i_min, i_max);
	domain->toroidal(true);
	embedDomain(domain);
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::baseDomain(const Vector<DIM, POS>& i_min, const std::vector< Vector<DIM, POS> >& i_base)
{
	BaseDomain<POS>* domain = new BaseDomain<POS>(i_min, i_base);
	domain->toroidal(false);
	embedDomain(domain);
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::baseToroidalDomain(const Vector<DIM, POS>& i_min, const std::vector< Vector<DIM, POS> >& i_base)
{
	BaseDomain<POS>* domain = new BaseDomain<POS>(i_min, i_base);
	domain->toroidal(true);
	embedDomain(domain);
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::defaultDomain()
{
	unitToroidalDomain();
}

/* PointSet 2.0 - Replication *****************************************/

template<int DIM, typename POS, typename VAL>
typename PointSet<DIM, POS, VAL>::replicationIterator PointSet<DIM, POS, VAL>::replicationBegin(
	bool keepOriginal) const
{
	return PointSet<DIM, POS, VAL>::replicationIterator(this, 0, keepOriginal);
}

template<int DIM, typename POS, typename VAL>
typename PointSet<DIM, POS, VAL>::replicationIterator PointSet<DIM, POS, VAL>::replicationBegin(
	bool keepOriginal,
	const POS& distMax) const
{
	return PointSet<DIM, POS, VAL>::replicationIterator(this, 0, keepOriginal, meanDistance()*distMax);
}

template<int DIM, typename POS, typename VAL>
typename PointSet<DIM, POS, VAL>::replicationIterator PointSet<DIM, POS, VAL>::replicationEnd() const
{
	return typename PointSet<DIM, POS, VAL>::replicationIterator(this);
}

/* PointSet 2.0 - Statistics ******************************************/

template<int DIM, typename POS, typename PT>
double PointSet<DIM, POS, PT>::meanDistance() const
{
	if(DIM == 2)
	{
		return std::sqrt(m_domain->volume()/this->size());
	}
	else
	{
		return std::pow(m_domain->volume()/this->size(), 1.0/(double)DIM);
	}
}

template<int DIM, typename POS, typename VAL>
std::ostream& operator<<(std::ostream &os, const PointSet<DIM, POS, VAL>& i_pts)
{
	os << "PointSet" << DIM;
	
	if(typeid(int) == typeid(POS)) os << "i";
	else if(typeid(double) == typeid(POS)) os << "d";
	else if(typeid(float) == typeid(POS)) os << "f";
	else os << "?";
	
	if(typeid(int) == typeid(VAL)) os << "i";
	if(typeid(int) == typeid(VAL)) os << "ui";
	else if(typeid(double) == typeid(VAL)) os << "d";
	else if(typeid(float) == typeid(VAL)) os << "f";
	else if(typeid(std::complex<double>) == typeid(VAL)) os << "c";
	else if(typeid(std::complex<float>) == typeid(VAL)) os << "c";
	else os << "?";
	
	os << "(" << i_pts.size() << " pts, " << *(i_pts.domain()) << ")";
	
	return os;
}

/* Typedef ************************************************************/

typedef PointSet<2, double, Complexd > PointSet2dc;
typedef PointSet<2, double, double > PointSet2dd;
typedef PointSet<2, double, int > PointSet2di;
typedef PointSet<3, double, Complexd > PointSet3dc;
typedef PointSet<3, double, double > PointSet3dd;
typedef PointSet<3, double, int > PointSet3di;

}

#endif
