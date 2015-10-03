#ifndef __STK_EXCEPTION__
#define __STK_EXCEPTION__

#include <string>
#include <sstream>
#include <exception>
#include <vector>
#include <typeinfo>
#include <complex>

#include <stk/debug.hpp>

namespace stk
{

template<int DIM, typename POS> class Vector;

namespace exception
{

class Message : public std::exception
{
	private:
		std::string m_msg;
		debug::Info m_dbg;
		std::string m_string;
	
	protected:
		void setMessage(const std::string& i_msg);
	
	public:
		Message(const stk::debug::Info& dbg) throw();
		Message(std::string i_msg, const stk::debug::Info& dbg) throw();
		virtual ~Message() throw();
		virtual const char* what() const throw();
		
		template<typename T>
		Message& addVar(const std::string& a, const T& b)
		{
			std::stringstream out;
			out << "\033[36m" << a << " = " << b << "\033[0m" << std::endl;
			m_string += out.str();
			return *this;
		}
		
		template<int DIM, typename POS>
		Message& addVar(const std::string& a, const stk::Vector<DIM, POS>& b)
		{
			std::stringstream out;
			out << "\033[36m" << a << " = " << "Vector" << DIM << "(";
			for(int i=0; i<DIM; i++)
			{
				out << b[i];
				if(i < DIM-1) out << ", ";
			}
			out << ")" << "\033[0m" << std::endl;
			m_string += out.str();
			return *this;
		}
};

class FileNotFound : public Message
{
	public:
		FileNotFound(std::string filename, const stk::debug::Info& dbg) throw();
};

class FileNotWritable : public Message
{
	public:
		FileNotWritable(std::string filename, const stk::debug::Info& dbg) throw();
};

class OutOfRange : public Message
{
	public:
		OutOfRange(int i, int m, int M, const stk::debug::Info& dbg) throw();
};

#ifdef CAIRO_ENABLED
class CairoError : public Message
{
	public:
		CairoError(const int& error, const stk::debug::Info& dbg) throw();
};
#endif

}

}

#endif
