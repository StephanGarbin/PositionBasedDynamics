#include "commonMath.h"
#include <iostream>

int dsyevj3(Eigen::Matrix3f& A, Eigen::Matrix3f& Q, Eigen::Matrix3f& w);

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
	sqrt_p = sqrt(abs(p));

	phi = 27.0 * (0.25*sqr(c1)*(p - c1) + c0*(q + 27.0 / 4.0*c0));
	phi = (1.0 / 3.0) * atan2(sqrt(abs(phi)), q);

	c = sqrt_p * cos(phi);
	s = 0.57735026918962576450914878050195745564760175127012687601860232648397767230293334569371539558574952522520871380513556767 * sqrt_p * sin(phi);

	Eigen::Matrix3d A_d;
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			A_d(row, col) = (double)A(row, col);
		}
	}
	Eigen::Matrix3d eigenValues_d;
	Eigen::Matrix3d eigenVectors_d;
	eigenValues_d.setZero();
	eigenValues_d(1, 1) = (1.0 / 3.0)*(m - c);
	eigenValues_d(2, 2) = eigenValues_d(1, 1) + s;
	eigenValues_d(0, 0) = eigenValues_d(1, 1) + c;
	eigenValues_d(1, 1) -= s;

	for (int e = 0; e < 2; ++e)
	{
		eigenVectors_d.col(e) = (A_d.col(0) - eigenValues_d(e, e) * Eigen::Vector3d(1.0, 0.0, 0.0)).cross(A_d.col(1) - eigenValues_d(e, e) * Eigen::Vector3d(0.0, 1.0, 0.0));
	}

	if (std::abs(eigenValues_d(0, 0) - eigenValues_d(1, 1)) == 0.0)
	{
		eigenVectors_d.col(1) = eigenVectors_d.col(0).cross(A_d.col(0) - eigenValues_d(0, 0) * Eigen::Vector3d(1.0, 0.0, 0.0));
	}

	eigenVectors_d.col(2) = eigenVectors_d.col(0).cross(eigenVectors_d.col(1));
	
	bool recompute = false;
	for (int i = 0; i < 3; ++i)
	{
		if (eigenVectors_d.col(i).squaredNorm() > 1e-10)
		{
			eigenVectors_d.col(i) = eigenVectors_d.col(i).normalized();
		}
		else
		{
			recompute = true;
		}
	}

	if (recompute)
	{
		//Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(A);
		//eigenValues = eigenSolver.pseudoEigenvalueMatrix(); //squared eigenvalues of F
		//eigenVectors = eigenSolver.pseudoEigenvectors(); //eigenvectors
		eigenValues.setZero();
		dsyevj3(A, eigenVectors, eigenValues);
	}
	else
	{
		for (int row = 0; row < 3; ++row)
		{
			for (int col = 0; col < 3; ++col)
			{
				eigenValues(row, col) = (float)eigenValues_d(row, col);
				eigenVectors(row, col) = (float)eigenVectors_d(row, col);
			}
		}
	}

}

void raySphereIntersect(const Eigen::Vector3f& sphereCentre, float sphereRadius,
	const Eigen::Vector3f& rayOrigin, const Eigen::Vector3f rayDirection, int& numIntersection, std::vector<Eigen::Vector3f>& intersectionPoints)
{
	float term = std::powf(rayDirection.dot(rayOrigin - sphereCentre), 2.0f) - std::sqrtf((rayOrigin - sphereCentre).squaredNorm()) + std::powf(sphereRadius, 2.0f);

	if (term < 0)
	{
		numIntersection = 0;
		intersectionPoints[0].setZero();
		intersectionPoints[1].setZero();
	}
	else if (term == 0.0f)
	{
		numIntersection = 1;
		float d = -rayDirection.dot(rayOrigin - sphereCentre);
		intersectionPoints[0] = rayOrigin + d * rayDirection;
		intersectionPoints[1].setZero();
	}
	else
	{
		numIntersection = 2;
		float d1 = -rayDirection.dot(rayOrigin - sphereCentre) + std::sqrtf(term);
		float d2 = -rayDirection.dot(rayOrigin - sphereCentre) - std::sqrtf(term);
		intersectionPoints[0] = rayOrigin + d1 * rayDirection;
		intersectionPoints[1] = rayOrigin + d2 * rayDirection;
	}

}

int dsyevj3(Eigen::Matrix3f& A, Eigen::Matrix3f& Q, Eigen::Matrix3f& w)
{
	const int n = 3;
	double sd, so;                  // Sums of diagonal resp. off-diagonal elements
	double s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
	double g, h, z, theta;          // More temporary storage
	double thresh;

	// Initialize Q to the identitity matrix
	Q.setIdentity();

	// Initialize w to diag(A)
	for (int i = 0; i < n; i++)
		w(i,i) = A(i,i);

	// Calculate SQR(tr(A))  
	sd = 0.0;
	for (int i = 0; i < n; i++)
		sd += fabs(w(i,i));
	sd = sqr(sd);

	// Main iteration loop
	for (int nIter = 0; nIter < 50; nIter++)
	{
		// Test for convergence 
		so = 0.0;
		for (int p = 0; p < n; p++)
		for (int q = p + 1; q < n; q++)
			so += fabs(A(p,q));
		if (so == 0.0)
			return 0;

		if (nIter < 4)
			thresh = 0.2 * so / sqr(n);
		else
			thresh = 0.0;

		// Do sweep
		for (int p = 0; p < n; p++)
		for (int q = p + 1; q < n; q++)
		{
			g = 100.0 * fabs(A(p,q));
			if (nIter > 4 && fabs(w(p,p)) + g == fabs(w(p,p))
				&& fabs(w(q,q)) + g == fabs(w(q,q)))
			{
				A(p,q) = 0.0;
			}
			else if (fabs(A(p,q)) > thresh)
			{
				// Calculate Jacobi transformation
				h = w(q,q) - w(p,p);
				if (fabs(h) + g == fabs(h))
				{
					t = A(p,q) / h;
				}
				else
				{
					theta = 0.5 * h / A(p,q);
					if (theta < 0.0)
						t = -1.0 / (sqrt(1.0 + sqr(theta)) - theta);
					else
						t = 1.0 / (sqrt(1.0 + sqr(theta)) + theta);
				}
				c = 1.0 / sqrt(1.0 + sqr(t));
				s = t * c;
				z = t * A(p,q);

				// Apply Jacobi transformation
				A(p,q) = 0.0;
				w(p) -= z;
				w(q) += z;
				for (int r = 0; r < p; r++)
				{
					t = A(r,p);
					A(r,p) = c*t - s*A(r,q);
					A(r,q) = s*t + c*A(r,q);
				}
				for (int r = p + 1; r < q; r++)
				{
					t = A(p,r);
					A(p,r) = c*t - s*A(r,q);
					A(r,q) = s*t + c*A(r,q);
				}
				for (int r = q + 1; r < n; r++)
				{
					t = A(p,r);
					A(p,r) = c*t - s*A(q,r);
					A(q,r) = s*t + c*A(q,r);
				}

				// Update eigenvectors
#ifndef EVALS_ONLY          
				for (int r = 0; r < n; r++)
				{
					t = Q(r, p);
					Q(r,p) = c*t - s*Q(r,q);
					Q(r,q) = s*t + c*Q(r,q);
				}
#endif
			}
		}
	}

	return -1;
}