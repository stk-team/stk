#ifndef __STK_HUNGARIAN_ALGORITHM_RAW__
#define __STK_HUNGARIAN_ALGORITHM_RAW__

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

namespace stk
{

#define HUNGARIANALGO_EPS 1e-25

template<typename VAL>
class _HungarianAlgorithm
{
private:
	VAL* m_matrix;
	int* m_mask;
	std::vector<bool> m_mask_col;
	std::vector<bool> m_mask_row;

	int m_width;
	int m_height;
	int m_size;

	int m_curX;
	int m_curY;

	bool findUncoveredInMatrix(const VAL& item, int& i, int& j) const;
	int step1();
	int step2();
	int step3();
	int step4();
	int step5();

public:
	_HungarianAlgorithm(VAL* i_matrix, int i_width, int i_height);
	~_HungarianAlgorithm();
	void solve(std::map<int,int>& o_result);
};

template<typename VAL>
_HungarianAlgorithm<VAL>::_HungarianAlgorithm(VAL* i_matrix, int i_width, int i_height) :
	m_width(i_width),
	m_height(i_height),
	m_matrix(i_matrix),
	m_size(i_width*i_height)
{
	m_mask = new int[m_size];
	
	m_mask_col.resize(m_width);
	m_mask_row.resize(m_height);
	
	for(int i=0; i<m_width; i++) m_mask_col[i] = false;
	for(int i=0; i<m_height; i++) m_mask_row[i] = false;
	
	for(int i=0; i<m_size; i++)
	{
		m_matrix[i] = i_matrix[i];
	}
}

template<typename VAL>
_HungarianAlgorithm<VAL>::~_HungarianAlgorithm()
{
	delete[] m_mask;
}

template<typename VAL>
bool _HungarianAlgorithm<VAL>::findUncoveredInMatrix(const VAL& i_item, int& o_x, int& o_y) const
{
	for(int j=0; j<m_height; j++)
		if(!m_mask_row[j])
			for(int i=0; i<m_width; i++)
				if(!m_mask_col[i])
					if(m_matrix[j*m_width+i] == i_item)
					{
						o_x = i;
						o_y = j;
						return true;
					}

	return false;
}

template<typename VAL>
void _HungarianAlgorithm<VAL>::solve(std::map<int,int>& o_result)
{
	for(int i=0; i<m_size; i++)
		m_mask[i] = 0;

	int step = 1;
	int nbOps = 0;
	while(step > 0)
	{
		//std::cout << "===== step " << step << " =====" << std::endl;

		if(nbOps%5000 == 0 && nbOps != 0)
			std::cout << "still alive ! (nbOps = " << nbOps << ")" << std::endl;

		switch(step)
		{
			case 1:
				step = step1();
				break;
			case 2:
				step = step2();
				break;
			case 3:
				step = step3();
				break;
			case 4:
				step = step4();
				break;
			case 5:
				step = step5();
				break;
			default:
				step = 0;
		}

		nbOps++;

//		 for(int j=0; j<m_height; j++)
//		 {
//			 for(int i=0; i<m_width; i++)
//			 {
//				 std::cout << m_matrix[j*m_height+i] << "\t";
//			 }
//			 std::cout << std::endl;
//		 }
//		 for(int j=0; j<m_height; j++)
//		 {
//			 for(int i=0; i<m_width; i++)
//			 {
//				 std::cout << m_mask[j*m_height+i] << "\t";
//			 }
//			 std::cout << std::endl;
//		}
//		std::cout<<"row : ";
//		for(int j=0; j<m_height; j++)
//			std::cout<<m_mask_row[j]<<" ";
//		std::cout<<std::endl<<"col : ";
//		for(int j=0; j<m_width; j++)
//			std::cout<<m_mask_col[j]<<" ";
//		std::cout<<std::endl;
	}

	for(int j=0; j<m_height; j++)
		for(int i=0; i<m_width; i++)
			if(m_mask[j*m_height+i] == 1)
				o_result[i] = j;
}

template<typename VAL>
int _HungarianAlgorithm<VAL>::step1()
{
	for(int j=0; j<m_height; j++)
		for(int i=0; i<m_width; i++)
		{
			double& val = m_matrix[j*m_width+i];
			int& maskVal = m_mask[j*m_width+i];

			if(val == 0.0)
			{
				bool freeZero = true;
				for(int x=0; x<m_width; x++)
					if(m_mask[j*m_width+x] == 1)
					{
						freeZero = false;
						break;
					}

				if(freeZero)
					for(int x=0; x<m_width; x++)
						if(m_mask[x*m_width+i] == 1)
						{
							freeZero = false;
							break;
						}

				if(freeZero)
					m_mask[j*m_width+i] = 1;
			}
		}

	return 2;
}

template<typename VAL>
int _HungarianAlgorithm<VAL>::step2()
{
	int covercount = 0;
	for(int j=0; j<m_height; j++)
		for(int i=0; i<m_width; i++)
			if(m_mask[j*m_width+i] == 1)
			{
				covercount++;
				m_mask_col[i] = true;
			}

	if(m_width > m_height)
	{
		if(covercount >= m_width) return 0;
	}
	else
	{
		if(covercount >= m_height) return 0;
	}
	

	return 3;
}

template<typename VAL>
int _HungarianAlgorithm<VAL>::step3()
{
	if(findUncoveredInMatrix(0.0, m_curX, m_curY))
		m_mask[m_curY*m_width+m_curX] = 2;
	else
		return 5;

	for(int i=0; i<m_width; i++)
		if(m_mask[m_curY*m_width+i] == 1)
		{
			m_mask_row[m_curY] = true;
			m_mask_col[i] = false;
			return 3;
		}

	return 4;
}

template<typename VAL>
int _HungarianAlgorithm<VAL>::step4()
{
	std::vector<int> seq;
	seq.push_back(m_curY*m_width+m_curX);

	int i = m_curX;
	int j;
	bool found;

	do
	{
		found = false;
		for(j=0; j<m_height; j++)
			if(m_mask[j*m_width+i] == 1)
				if(std::find(seq.begin(), seq.end(), j*m_width+i) == seq.end())
				{
					seq.push_back(j*m_width+i);
					found = true;
					break;
				}

		if(found)
			for(i=0; i<m_width; i++)
				if(m_mask[j*m_width+i] == 2)
					if(std::find(seq.begin(), seq.end(), j*m_width+i) == seq.end())
					{
						seq.push_back(j*m_width+i);
						found = true;
						break;
					}
	}
	while(found);

	for(int k=0; k<seq.size(); k++)
		m_mask[seq[k]]--;

	for(int k=0; k<m_size; k++)
		if(m_mask[k] == 2) m_mask[k] = 0;

	for(int k=0; k<m_width; k++)
		m_mask_col[k] = false;
	for(int k=0; k<m_height; k++)
		m_mask_row[k] = false;

	return 2;
}

template<typename VAL>
int _HungarianAlgorithm<VAL>::step5()
{
	VAL h;

	for(int j=0; j<m_height; j++)
		if(!m_mask_row[j])
			for(int i=0; i<m_width; i++)
				if(!m_mask_col[i])
				{
					VAL& v = m_matrix[j*m_width+i];
					if(std::abs(h) < HUNGARIANALGO_EPS || h > v && std::abs(v) > HUNGARIANALGO_EPS)
						h = v;
				}

	for(int j=0; j<m_height; j++)
		if(m_mask_row[j])
			for(int i=0; i<m_width; i++)
				m_matrix[j*m_width+i] += h;

	for(int i=0; i<m_width; i++)
		if(!m_mask_col[i])
			for(int j=0; j<m_height; j++)
			{
				VAL& v = m_matrix[j*m_width+i];
				v -= h;
				if(std::abs(v) < HUNGARIANALGO_EPS)
					v = 0.0;
			}

	//~ std::cout << h << std::endl;

	return 3;
}

}

#endif
