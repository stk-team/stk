#ifndef __STK_OPTIMIZATION_LLOYD__
#define __STK_OPTIMIZATION_LLOYD__

#include <stk/pointset.hpp>
#include <stk/voronoi.hpp>
#include <stk/optimization/optimizer.hpp>

namespace stk
{

namespace optimization
{

template<typename VAL>
class Lloyd : public Optimizer<2, double, VAL>
{
	protected:
		double m_moveRatio;
		double m_weightEnabled;
	
	public:
		Lloyd(double mr = 1.0)
			: m_moveRatio(mr)
		{
			m_weightEnabled = false;
		}

		void enableWeight()
		{
			m_weightEnabled = true;
		}
		
		void optimize(PointSet<2, double, VAL>& _pts)
		{
			optimize(_pts, 1);
		}

		void optimize(PointSet<2, double, VAL>& _pts, const unsigned int& i_iter)
		{
			bool torState = _pts.isToroidal();
			_pts.toroidal(true);

			std::vector<stk::VoronoiCell<double, VAL> > voronoi;
			for(int j=0; j<i_iter; j++)
			{
				if(m_weightEnabled) stk::weightedVoronoi(_pts, voronoi);
				else stk::voronoi(_pts, voronoi);

				for(int i=0; i<_pts.size(); i++)
				{
					_pts[i].pos() += (voronoi[i].getCentroid() - _pts[i].pos()) * m_moveRatio;
				}
			}

			_pts.toroidal(torState);
			_pts.normalize();
		}
		
		void optimize(PointSet<2, double, VAL>& _pts, const std::vector<bool>& i_constFlag)
		{
			optimize(_pts, i_constFlag, 1);
		}

		void optimize(PointSet<2, double, VAL>& _pts, const std::vector<bool>& i_constFlag, const unsigned int& i_iter)
		{
			if(_pts.size() > i_constFlag.size())
			{
				throw stk::exception::Message("_pts.size() > i_constFlag.size()", STK_DBG_INFO);
			}
			bool torState = _pts.isToroidal();
			_pts.toroidal(true);
			
			std::vector<stk::VoronoiCell<double, VAL> > voronoi;
			for(int j=0; j<i_iter; j++)
			{
				if(m_weightEnabled) stk::weightedVoronoi(_pts, voronoi);
				else stk::voronoi(_pts, voronoi);

				for(int i=0; i<_pts.size(); i++)
				{
					if(!i_constFlag[i])
					{
						_pts[i].pos() += (voronoi[i].getCentroid() - _pts[i].pos()) * m_moveRatio;
					}
				}
			}

			_pts.toroidal(torState);
			_pts.normalize();
		}
};

}

}

#endif
