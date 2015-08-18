#pragma once

#include <Eigen\Dense>


Eigen::Matrix3f kroneckerProduct(const Eigen::Vector3f& left, const Eigen::Vector3f& right);

float mcAuleyBrackets(float x);

float tensorProductInner(const Eigen::Matrix3f& left, const Eigen::Matrix3f& right);