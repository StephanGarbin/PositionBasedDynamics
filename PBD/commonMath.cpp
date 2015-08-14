#include "commonMath.h"

Eigen::Matrix3f kroneckerProduct(const Eigen::Vector3f& left, const Eigen::Vector3f& right)
{
	Eigen::Matrix3f result;

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			result(row, col) = left[row] * right[col];
		}
	}

	return result;
}