#ifndef __STK_IO_HISTOGRAM__
#define __STK_IO_HISTOGRAM__

#include <fstream>
#include <sstream>

#include <stk/histogram.hpp>
#include <stk/vector.hpp>
#include <stk/exception.hpp>
#include <stk/io/commons.hpp>
#include <stk/io/array.hpp>

namespace stk
{

namespace io
{

template<typename POS, typename VAL>
void read(std::string i_filename, Histogram<1, POS, VAL>& o_histo)
{
	std::vector<VAL> data;

	//Open file
	std::ifstream file;
	file.open(i_filename.c_str());
	if(!file) throw exception::FileNotFound(i_filename, STK_DBG_INFO);
	
	//Read head line
	POS pos;
	POS pMin;
	POS pMax;
	VAL value;
	if(file >> pos >> value)
	{
		pMin = pos;
		pMax = pos;
		
		do
		{
			data.push_back(value);
			if(pMin > pos) pMin = pos;
			if(pMax < pos) pMax = pos;
		}
		while(file >> pos >> value);
	}

	o_histo.resize(data.size());
	for(int i=0; i<data.size(); i++)
	{
		o_histo[i] = data[i];
	}

	o_histo.setBoundaries(pMin, pMax);
}

template<int DIM, typename POS, typename VAL>
void write(
	std::string i_filename,
	const Histogram<DIM, POS, VAL>& i_array,
	FileType i_type = Dat);

template<int DIM, typename POS, typename VAL>
void writeDat(std::string i_filename, const Histogram<DIM, POS, VAL>& i_array)
{
	std::ofstream file;
	file.open(i_filename.c_str());
	if(!file) throw stk::exception::FileNotWritable(i_filename, STK_DBG_INFO);
	
	file.precision(10);
	for(int i=0; i < i_array.getArraySize(); i++)
	{
		file << i_array.getPosFromIndice(i) << '\t';
		file << i_array[i] << std::endl;
	}
	file.close();
}

}

}

#endif
