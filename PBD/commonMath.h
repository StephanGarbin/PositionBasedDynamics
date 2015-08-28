#pragma once

#include <vector>
#include <Eigen\Dense>


Eigen::Matrix3f kroneckerProduct(const Eigen::Vector3f& left, const Eigen::Vector3f& right);

float mcAuleyBrackets(float x);

float tensorProductInner(const Eigen::Matrix3f& left, const Eigen::Matrix3f& right);

void eigenDecompositionCardano(Eigen::Matrix3f& A, Eigen::Matrix3f& eigenValues, Eigen::Matrix3f& eigenVectors);

void eigenDecompositionJacobi(Eigen::Matrix3f& A, Eigen::Matrix3f& eigenValues, Eigen::Matrix3f& eigenVectors);


template<typename T>
T sqr(T x)
{
	return x * x;
}


void raySphereIntersect(const Eigen::Vector3f& sphereCentre, float sphereRadius,
	const Eigen::Vector3f& rayOrigin, const Eigen::Vector3f rayDirection, int& numIntersection, std::vector<Eigen::Vector3f>& intersectionPoints);