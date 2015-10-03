#ifndef __STK_FOURIER_TRANSFORM_CFT__
#define __STK_FOURIER_TRANSFORM_CFT__

#include <stk/fourier/transform-commons.hpp>

namespace stk
{

namespace fourier
{

/*
 ***********************************************************************
 * TRANSFORM CFT CORE
 ***********************************************************************
 */
 
template<int DIM, typename INPUT, typename POS, typename VAL>
void _transformCore(
	const INPUT& _input,
	const Vector<DIM, POS>& i_pos,
	VAL& o_val,
	Direction i_dir);

template<int DIM, typename POS, typename VALO, typename VALI>
void _transformCore(
	const PointSet<DIM, POS, VALI>& i_input,
	const Vector<DIM, POS>& i_oPos,
	std::complex<VALO>& o_oVal,
	Direction i_dir)
{
	double phFact = i_dir * 2.0 * M_PI;
	
	o_oVal = 0;
	for(unsigned int j=0; j<i_input.size(); j++)
	{
		const VALI& iVal = i_input[j].val();
		const Vector<DIM, POS>& iPos = i_input[j].pos();
		
		double ph = phFact * dotProduct(i_oPos, iPos);
		o_oVal += iVal * std::polar(1.0, ph);
	}
}

template<int DIM, typename POS, typename VALO, typename VALI>
void _transformCore(
	const PointSet<DIM, POS, VALI>& i_input,
	const Vector<DIM, POS>& i_oPos,
	VALO& o_oVal,
	Direction i_dir)
{
	double phFact = i_dir * 2.0 * M_PI;
	
	std::complex<VALO> res(0.0, 0.0);
	for(unsigned int j=0; j<i_input.size(); j++)
	{
		const VALI& iVal = i_input[j].val();
		const Vector<DIM, POS>& iPos = i_input[j].pos();
		
		double ph = phFact * dotProduct(i_oPos, iPos);
		res += iVal * std::polar(1.0, ph);
	}
	o_oVal = std::abs(res);
}

template<int DIM, typename POS, typename VALO, typename VALI>
void _transformCore(
	const Histogram<DIM, POS, VALI>& i_input,
	const Vector<DIM, POS>& i_oPos,
	std::complex<VALO>& o_oVal,
	Direction i_dir)
{
	double phFact = i_dir * 2.0 * M_PI;
	
	o_oVal = 0;
	for(unsigned int j=0; j<i_input.getArraySize(); j++)
	{
		const VALI& iVal = i_input.getFromIndice(j);
		Vector<DIM, POS> iPos = i_input.getPosFromIndice(j);
		
		double ph = phFact * dotProduct(i_oPos, iPos);
		o_oVal += iVal * std::polar(1.0, ph);
	}
}

template<int DIM, typename POS, typename VALO, typename VALI>
void _transformCore(
	const Histogram<DIM, POS, VALI>& i_input,
	const Vector<DIM, POS>& i_oPos,
	VALO& o_oVal,
	Direction i_dir)
{
	double phFact = i_dir * 2.0 * M_PI;
	
	std::complex<VALO> res(0.0, 0.0);
	for(unsigned int j=0; j<i_input.getArraySize(); j++)
	{
		const VALI& iVal = i_input.getFromIndice(j);
		Vector<DIM, POS> iPos = i_input.getPosFromIndice(j);
		
		double ph = phFact * dotProduct(i_oPos, iPos);
		res += iVal * std::polar(1.0, ph);
	}
	o_oVal = std::abs(res);
}

template<typename INPUT, int DIM, typename POS, typename VAL>
void _transformCft(
	const INPUT& i_input,
	PointSet<DIM, POS, VAL>& o_output,
	Direction i_dir)
{
	unsigned int oSz = o_output.size();
	
	for(unsigned int i=0; i<oSz; i++)
	{
		VAL& oVal = o_output[i].val();
		const Vector<DIM, POS>& oPos = o_output[i].pos();
		
		_transformCore(i_input, oPos, oVal, i_dir);
	}
}

template<typename INPUT, int DIM, typename POS, typename VAL>
void _transformCft(
	const INPUT& i_input,
	Histogram<DIM, POS, VAL>& o_output,
	Direction i_dir)
{	
	for(int i=0; i<o_output.getArraySize(); i++)
	{
		VAL& oVal = o_output.getFromIndice(i);
		const Vector<DIM, POS> oPos = o_output.getPosFromIndice(i);
		
		_transformCore(i_input, oPos, oVal, i_dir);
	}
}

}

}

#endif
