#ifndef __STK_IO_ARRAY__
#define __STK_IO_ARRAY__

#include <fstream>
#include <sstream>
#include <vector>

#include <stk/array.hpp>
#include <stk/vector.hpp>
#include <stk/exception.hpp>
#include <stk/io/commons.hpp>
#include <png.h>

namespace stk
{

namespace io
{

template<typename VAL>
void read(std::string i_filename, Array<1, VAL>& o_array)
{
	std::vector<VAL> data;

	//Open file
	std::ifstream file;
	file.open(i_filename.c_str());
	if(!file) throw exception::FileNotFound(i_filename, STK_DBG_INFO);
	
	//Read head line
	VAL value;
	while(file >> value)
	{
		data.push_back(value);
	}

	o_array.resize(data.size());
	for(int i=0; i<data.size(); i++)
	{
		o_array[i] = data[i];
	}
}

template<int DIM, typename VAL>
void write(
	std::string i_filename,
	const Array<DIM, VAL>& i_array,
	FileType i_type = Dat);

template<int DIM, typename VAL>
void writeDat(std::string i_filename, const Array<DIM, VAL>& i_array)
{
	std::ofstream file;
	file.open(i_filename.c_str());
	if(!file) throw stk::exception::FileNotWritable(i_filename, STK_DBG_INFO);
	
	file.precision(10);
	for(int i=0; i < i_array.getArraySize(); i++)
	{
		file << i_array.getData(i) << std::endl;
	}
	file.close();
}

template<typename VAL>
void writeDat(std::string i_filename, const Array<2, VAL>& i_array)
{
	std::ofstream file;
	file.open(i_filename.c_str());
	if(!file) throw stk::exception::FileNotWritable(i_filename, STK_DBG_INFO);
	
	file.precision(10);
	for(int j=0; j < i_array.getSize()[0]; j++)
	{
		for(int i=0; i < i_array.getSize()[1]; i++)
		{
			file << i_array.getData(Vector2i(i, j)) << '\t';
		}
		file << std::endl;
	}
	
	file.close();
}

void writePng(
	std::string i_filename,
	const Array<2, double>& i_array,
	double i_min = 0.0,
	double i_max = 1.0);

void draw(
	std::string i_filename,
	const Array<2, double>& i_array,
	double i_min = 0.0,
	double i_max = 1.0);

#ifdef CAIRO_ENABLED
void drawRegion(
	std::string i_filename,
	const Array2i& i_array,
	const Vector2i& i_size);
#endif

}

}

#endif
