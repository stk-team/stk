#include <stk/debug.hpp>

#include <sstream>

namespace stk
{

namespace debug
{

Info getInfo(
	const std::string& file,
	const std::string& func,
	int line)
{
	Info info;
	info.line = line;
	info.func = func;
	info.file = file;
	
	return info;
}

}

}
