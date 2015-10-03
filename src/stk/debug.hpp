#ifndef __STK_DEBUG__
#define __STK_DEBUG__

#include <string>
#include <exception>

#define STK_DBG_INFO stk::debug::getInfo(__FILE__, __PRETTY_FUNCTION__, __LINE__)

namespace stk
{

namespace debug
{

class Info
{
	public:
		int line;
		std::string func;
		std::string file;
};

Info getInfo(const std::string& file, const std::string& func, int line);

}

}

#endif
