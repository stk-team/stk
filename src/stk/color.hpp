#ifndef __STK_COLOR__
#define __STK_COLOR__

#include <stk/vector.hpp>

namespace stk
{

namespace color
{

Vector3d rgb2hsv(const Vector3d& i_rgb);

Vector3d hsv2rgb(const Vector3d& i_hsv);

}

}

#endif
