#ifndef __STK_POINTSET__
#define __STK_POINTSET__

#include <vector>
#include <stdexcept>

#include "point.hpp"
#include "vector.hpp"

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
        bool m_replication;
        bool m_toroidal;
        int m_repPtsSize;
        double m_repSize;
        VAL m_defaultVal;
        Vector<DIM, POS> m_pMin;
        Vector<DIM, POS> m_pMax;
        Vector<DIM, POS> m_pSize;

    public:
        PointSet();

        void setParams(const PointSet<DIM, POS, VAL>& pts);

        /**
         * Replicate a set of point in all direction. After the
         * replication, the number of point is replicate
         * (i_repSize*2+1)^DIM - 1 times
         * With no argument, choose a replication of wide 10 points around boundaries
         * @brief Replicate a set of point in all direction
         * @param i_repSize Number of replication in each direction
         * @see removeReplication
         * @see isReplicated
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        void replicate();
        void replicate(float i_repSize);

        /**
         * @return Return the number of point before the replication
         * @see replicate
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        int getNotReplicatedSize() const;

        /**
         * Remove points added by replicate and reset the replications
         * state
         * @see replicate
         * @see isReplicated
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        void removeReplication();

        /**
         * Check if the pointset is replicated
         * @return Return true if replicate() was calling before, else
         * false
         * @see replicate
         * @see removeReplication
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        bool isReplicated() const;

        /**
         * set if consider the pointset space as toroidal or not
         * @see isToroidal()
         */
        void toroidal(bool i_tor);
        /**
         * Check if the pointset is toroidal
         * @return Return true if toroidal, else false
         * @see toroidal()
         */
        bool isToroidal() const;

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
         * Get minimum position in space for all dimension
         * @see getMaxSpace
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        const Vector<DIM, POS>& getMinPosition() const;

        /**
         * Get maximum position in space for all dimension
         * @see getMinSpace
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        const Vector<DIM, POS>& getMaxPosition() const;

        /**
         * Get center position in space for all dimension
         * @see getMinSpace
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        Vector<DIM, POS> getMeanPosition() const;

        /**
         * Get space width for all dimension
         * @see getSpaceSize
         * @author Adrien PILLEBOUE <adrien.pilleboue@liris.cnrs.fr>
         */
        const Vector<DIM, POS>& getSpaceSize() const;

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

        void operator+=(const VAL& i_val);
        void operator-=(const VAL& i_val);
        void operator*=(const VAL& i_val);
        void operator/=(const VAL& i_val);
};




/*
 * ********************************************************************
 * IMPLEMENTATION
 * ********************************************************************
 */

template<int DIM, typename POS, typename VAL>
PointSet<DIM, POS, VAL>::PointSet()
{
    m_replication = false;
    m_toroidal = false;
    m_repPtsSize = 0;
    m_repSize = 0;
    m_defaultVal = 1;

    for(int i=0; i<DIM; i++)
    {
        m_pMin[i] = 0;
        m_pMax[i] = 1;
    }
    m_pSize = m_pMax - m_pMin;
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::setParams(const PointSet<DIM, POS, VAL>& pts)
{
    this->setBoundaries(pts.getMinPosition(), pts.getMaxPosition());
    this->m_defaultVal = pts.getDefaultVal();
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::replicate()
{
    double factor = 6.4*m_pSize.getMin()/std::sqrt(this->size());
    replicate(factor);
}

    template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::replicate(float i_repSize)
{
    if(m_replication)
    {
        if(m_repSize == i_repSize) return;
        else removeReplication();
    }

    m_repPtsSize = this->size();
    m_replication = true;
    m_toroidal = false;
    m_repSize = i_repSize;

    int ceilRepSize = ceil(i_repSize);
    Vector<DIM, POS> ind;

    for(int i=0; i<DIM; i++) ind[i] = -ceilRepSize;

    int d=0;
    float sum;
    while(d < DIM)
    {
        sum = 0;

        for(int i=0; i<DIM; i++) sum += std::abs(ind[i]);

        if(sum != 0.)
        {
            for(int i=0; i<m_repPtsSize; i++)
            {
                Point<DIM, POS, VAL> pt(this->at(i));
                pt.pos() += ind;

                if( (pt.pos()-0.5).norm(0)<(0.5+i_repSize) ) this->push_back(pt);
            }
        }

        ind[0] += m_pSize[0];

        d=0;

        while(ind[d] > ceilRepSize && d < DIM)
        {
            ind[d] = -ceilRepSize;
            d++;
            if(d < DIM) ind[d] += m_pSize[d];
        }
    }
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::removeReplication()
{
    this->erase(this->begin()+m_repPtsSize, this->end());
    m_replication = false;
}

template<int DIM, typename POS, typename VAL>
bool PointSet<DIM, POS, VAL>::isReplicated() const
{
    return m_replication;
}

    template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::toroidal(bool i_tor)
{
    m_toroidal = i_tor;
    if(isToroidal())
    {
        if(isReplicated())
            removeReplication();
        normalize();
    }
}

template<int DIM, typename POS, typename VAL>
bool PointSet<DIM, POS, VAL>::isToroidal() const
{
    return m_toroidal;
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::minTorVect(Vector<DIM, POS>& i_vect) const
{
    if(!isToroidal())
        throw std::runtime_error("this pointset is not toroidal!");

    for(unsigned int i=0; i<DIM; i++)
    {
        while(std::abs(i_vect[i]-m_pMin[i]) > m_pSize[i]/2.)
        {
            if(i_vect[i]-m_pMin[i] > m_pSize[i]/2.)
                i_vect[i] -= m_pSize[i];
            else
                i_vect[i] += m_pSize[i];
        }
    }
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::normalize(Vector<DIM, POS>& i_pt) const
{
    if(isToroidal())
    {
        for(unsigned int i=0; i<DIM; i++)
        {
            while(i_pt[i]>m_pMax[i])
            {
                i_pt[i] -= m_pSize[i];
            }
            while(i_pt[i]<m_pMin[i])
            {
                i_pt[i] += m_pSize[i];
            }
        }
    }
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::normalize(const Vector<DIM, POS>& i_pt) const
{
    if(isToroidal())
    {
        Vector<DIM, POS> res = i_pt;

        for(unsigned int i=0; i<DIM; i++)
        {
            while(res[i] > m_pMax[i])
            {
                res[i] -= m_pSize[i];
            }
            while(res[i] < m_pMin[i])
            {
                res[i] += m_pSize[i];
            }
        }

        return res;
    }
    else return i_pt;
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::normalize()
{
    for(unsigned int i=0; i<this->size(); i++)
        normalize(this->at(i).pos());
}

template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::diff(const Vector<DIM, POS>& i_v1, const Vector<DIM, POS>& i_v2) const
{
    Vector<DIM, POS> o_vect;
    if(isToroidal())
    {
        o_vect = normalize(i_v2)-normalize(i_v1);
        minTorVect(o_vect);
        return o_vect;
    }
    else return i_v2-i_v1;
}

template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Vector<DIM, POS>& i_v1, const Vector<DIM, POS>& i_v2) const
{
    return diff(i_v1, i_v2).norm();
}

template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Point<DIM, POS, VAL>& i_p1, const Point<DIM, POS, VAL>& i_p2) const
{
    return dist(i_p1.pos(), i_p2.pos());
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(unsigned int i_p1, unsigned int i_p2) const
{
    return dist(this->at(i_p1).pos(), this->at(i_p2).pos());
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(unsigned int i_p1, const Point<DIM, POS, VAL>& i_p2) const
{
    return dist(this->at(i_p1).pos(), i_p2.pos());
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Point<DIM, POS, VAL>& i_p1, unsigned int i_p2) const
{
    return dist(i_p1.pos(), this->at(i_p2).pos());
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(unsigned int i_p1, const Vector<DIM, POS>& i_p2) const
{
    return dist(this->at(i_p1).pos(), i_p2);
}
template<int DIM, typename POS, typename VAL>
POS PointSet<DIM, POS, VAL>::dist(const Vector<DIM, POS>& i_p1, unsigned int i_p2) const
{
    return dist(i_p1, this->at(i_p2).pos());
}

template<int DIM, typename POS, typename VAL>
VAL PointSet<DIM, POS, VAL>::getDefaultVal() const
{
    return m_defaultVal;
}


template<int DIM, typename POS, typename VAL>
Vector<DIM, POS> PointSet<DIM, POS, VAL>::getMeanPosition() const
{
    return (m_pMax + m_pMin)/2.0;
}

template<int DIM, typename POS, typename VAL>
const Vector<DIM, POS>& PointSet<DIM, POS, VAL>::getMinPosition() const
{
    return m_pMin;
}

template<int DIM, typename POS, typename VAL>
const Vector<DIM, POS>& PointSet<DIM, POS, VAL>::getMaxPosition() const
{
    return m_pMax;
}

template<int DIM, typename POS, typename VAL>
const Vector<DIM, POS>& PointSet<DIM, POS, VAL>::getSpaceSize() const
{
    return m_pSize;
}

template<int DIM, typename POS, typename VAL>
void PointSet<DIM, POS, VAL>::setBoundaries(
    const Vector<DIM, POS>& i_min,
    const Vector<DIM, POS>& i_max)
{
    m_pMin = i_min;
    m_pMax = i_max;
    m_pSize = m_pMax - m_pMin;
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

typedef PointSet<2, double, double > PointSet2dd;
typedef PointSet<2, double, int > PointSet2di;
typedef PointSet<3, double, double > PointSet3dd;
typedef PointSet<3, double, int > PointSet3di;

}

#endif
