#include <stk/exception.hpp>

#include <sstream>
#include <iostream>

#ifdef CAIRO_ENABLED
#include <cairo/cairo.h>
#endif

namespace stk
{

namespace exception
{

Message::Message(const stk::debug::Info& i_dbg) throw()
{
	m_dbg = i_dbg;
}

Message::Message(std::string i_msg, const stk::debug::Info& i_dbg) throw()
{
	m_dbg = i_dbg;
	setMessage(i_msg);
}

Message::~Message() throw()
{
	
}

const char* Message::what() const throw()
{
	
	//~ for(int i=0; i<m_vars.size(); i++)
	//~ {
		//~ out << "\t" << m_vars[i].first << " = " << m_vars[i].second << std::endl;
	//~ }
	
	
	return m_string.c_str();
}

void Message::setMessage(const std::string& i_msg)
{
	m_msg = i_msg;
	
	std::stringstream out;
	out << "\033[1m\033[31m" << m_msg << std::endl;
	out << "\033[0m\033[33m";
	out << "file : " << m_dbg.file << std::endl;
	out << "func : " << m_dbg.func << std::endl;
	out << "line : " << m_dbg.line << std::endl;
	out << "\033[0m";
	m_string = out.str();
}

FileNotFound::FileNotFound(std::string i_filename, const stk::debug::Info& i_dbg) throw()
	: Message(i_dbg)
{
	std::stringstream out;
	out << i_filename << " not found";
	setMessage(out.str());
}

FileNotWritable::FileNotWritable(std::string i_filename, const stk::debug::Info& i_dbg) throw()
	: Message(i_dbg)
{
	std::stringstream out;
	out << "can't write in " << i_filename;
	setMessage(out.str());
}

OutOfRange::OutOfRange(int i_pos, int i_min, int i_max, const stk::debug::Info& i_dbg) throw()
	: Message(i_dbg)
{
	std::stringstream out;
	out << "out of range : " << i_pos << " not in [" << i_min << ", " << i_max << "]";
	setMessage(out.str());
}

#ifdef CAIRO_ENABLED
CairoError::CairoError(const int& error, const stk::debug::Info& i_dbg) throw()
	: Message(i_dbg)
{
	std::stringstream out;
	out << "cairo error ";
	switch(error)
	{
		case CAIRO_STATUS_SUCCESS:
			out << "CAIRO_STATUS_SUCCESS";
			break;
		case CAIRO_STATUS_NO_MEMORY:
			out << "CAIRO_STATUS_NO_MEMORY";
			break;
		case CAIRO_STATUS_SURFACE_TYPE_MISMATCH:
			out << "CAIRO_STATUS_SURFACE_TYPE_MISMATCH";
			break;
		case CAIRO_STATUS_WRITE_ERROR:
			out << "CAIRO_STATUS_WRITE_ERROR";
			break;
		case CAIRO_STATUS_INVALID_SIZE:
			out << "CAIRO_STATUS_INVALID_SIZE";
			break;
		default:
			out << "#" << error;
	}
	
	setMessage(out.str());
}
#endif

}

}
