#include "commonMath.h"
#include <iostream>

int dsyevc3(double A[3][3], double w[3]);

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

float mcAuleyBrackets(float x)
{
	return 0.5f * (x + std::abs(x));
}

float tensorProductInner(const Eigen::Matrix3f& left, const Eigen::Matrix3f& right)
{
	return (left.transpose() * right).trace();
}

void eigenDecompositionCardano(Eigen::Matrix3f& A, Eigen::Matrix3f& eigenValues, Eigen::Matrix3f& eigenVectors)
{
	double m, c1, c0;

	// Determine coefficients of characteristic poynomial. We write
	//       | a   d   f  |
	//  A =  | d*  b   e  |
	//       | f*  e*  c  |
	double de = A(0, 1) * A(1, 2);                                    // d * e
	double dd = sqr(A(0, 1));                                         // d^2
	double ee = sqr(A(1, 2));                                         // e^2
	double ff = sqr(A(0, 2));                                         // f^2
	m = A(0, 0) + A(1, 1) + A(2, 2);
	c1 = (A(0, 0) * A(1, 1) + A(0, 0) * A(2, 2) + A(1, 1) * A(2, 2))        // a*b + a*c + b*c - d^2 - e^2 - f^2
		- (dd + ee + ff);
	c0 = A(2, 2) * dd + A(0, 0) * ee + A(1, 1) * ff - A(0, 0) * A(1, 1) * A(2, 2)
		- 2.0 * A(0, 2) * de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

	double p, sqrt_p, q, c, s, phi;
	p = sqr(m) - 3.0*c1;
	q = m*(p - (3.0 / 2.0)*c1) - (27.0 / 2.0)*c0;
	sqrt_p = sqrt(fabs(p));

	phi = 27.0 * (0.25*sqr(c1)*(p - c1) + c0*(q + 27.0 / 4.0*c0));
	phi = (1.0 / 3.0) * atan2(sqrt(fabs(phi)), q);

	c = sqrt_p*cos(phi);
	s = (1.0 / 1.73205080756887729352744634151)*sqrt_p*sin(phi);

	eigenValues(1, 1) = (1.0 / 3.0)*(m - c);
	eigenValues(2, 2) = eigenValues(1,1) + s;
	eigenValues(0, 0) = eigenValues(1,1) + c;
	eigenValues(1, 1) -= s;

	for (int e = 0; e < 2; ++e)
	{
		eigenVectors.col(e) = (A.col(0) - eigenValues(e, e) * Eigen::Vector3f(1.0f, 0.0f, 0.0f)).cross(A.col(1) - eigenValues(e, e) * Eigen::Vector3f(0.0f, 1.0f, 0.0f));
	}

	if (std::abs(eigenValues(0, 0) - eigenValues(1, 1)) == 0.0f)
	{
		eigenVectors.col(1) = eigenVectors.col(0).cross(A.col(0) - eigenValues(0, 0) * Eigen::Vector3f(1.0f, 0.0f, 0.0f));
	}

	eigenVectors.col(2) = eigenVectors.col(0).cross(eigenVectors.col(1));

	for (int i = 0; i < 3; ++i)
	{
		if (eigenVectors.col(i).squaredNorm() != 0)
		{
			eigenVectors.col(i) = eigenVectors.col(i).normalized();
		}
	}
}