/*
Please note that the first three functions are the same as in the other courseworks.
*/
#pragma once

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

template<typename T>
T cotan(T value)
{
	//Both of these are valid

	//return std::tan(PI / (T)2 - value); //(OR: 1 / tan(x))
	return 1.0 / std::tan(value);
}