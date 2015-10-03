#include <algorithm>

#include <stk/color.hpp>
#include <stk/exception.hpp>

namespace stk
{

namespace color
{

Vector3d rgb2hsv(const Vector3d& i_rgb)
{
	Vector3d hsv;
	Vector3d rgbNorm;
	double rgbMin = i_rgb.getMin();
	double rgbMax = i_rgb.getMax();
	double diff = rgbMax - rgbMin;

	//Value
	hsv[2] = rgbMax;
	if(hsv[2] == 0.0)
	{
		hsv[0] = 0.0;
		hsv[1] = 0.0;
		return hsv;
	}
	
	//Saturation
	hsv[1] = diff / rgbMax;
	if(hsv[1] == 0.0)
	{
		hsv[0] = 0.0;
		return hsv;
	}

	//Hue
	if(i_rgb[0] == rgbMax)		hsv[0] = i_rgb[1] - i_rgb[2];
	else if(i_rgb[1] == rgbMax)	hsv[0] = 2.0 + i_rgb[2] - i_rgb[1];
	else						hsv[0] = 4.0 + i_rgb[0] - i_rgb[1];

	hsv[0] /= diff;
	hsv[0] *= 60.0/360.0;
	if(hsv[0] < 0.0) hsv[0] += 1.0;
	
	return hsv;
}

Vector3d hsv2rgb(const Vector3d& i_hsv)
{
	if(i_hsv[1] == 0.0)
	{
		return Vector3d(i_hsv[2]);
	}

	double f, p, q, t;
	int i;

	i = floor(i_hsv[0]*360.0/60.0);
	f = i_hsv[0]*360.0/60.0 - (double) i;
	p = i_hsv[2] * (1.0 - i_hsv[1]);
	q = i_hsv[2] * (1.0 - i_hsv[1] * f);
	t = i_hsv[2] * (1.0 - i_hsv[1] * (1.0 - f));

	switch(i)
	{
		case 0:
			return Vector3d(i_hsv[2], t, p);
		case 1:
			return Vector3d(q, i_hsv[2], p);
		case 2:
			return Vector3d(p, i_hsv[2], t);
		case 3:
			return Vector3d(p, q, i_hsv[2]);
		case 4:
			return Vector3d(t, p, i_hsv[2]);
		default:
			return Vector3d(i_hsv[2], p, q);
	}
}

}

}
