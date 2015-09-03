#pragma once
#include <vector>
#include <string>
#include <Eigen\Dense>


void writeCustomTetAttributes(std::vector<float>& youngsModulus,
	std::vector<float>& anisotropyStrength, std::vector<Eigen::Vector3f>& anisotropyDirections,
	const std::string& fileName);


void readCustomTetAttributes(std::vector<float>& youngsModulus,
std::vector<float>& anisotropyStrength, std::vector<Eigen::Vector3f>& anisotropyDirections,
const std::string& fileName);



