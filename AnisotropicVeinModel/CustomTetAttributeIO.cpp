#include "CustomTetAttributeIO.h"

#include <fstream>
#include <iostream>

#include <boost\algorithm\string.hpp>
#include <boost\algorithm\string\trim_all.hpp>

void writeCustomTetAttributes(std::vector<float>& youngsModulus,
	std::vector<float>& anisotropyStrength, std::vector<Eigen::Vector3f>& anisotropyDirections,
	const std::string& fileName)
{
	std::ofstream file;
	file.open(fileName);

	for (int i = 0; i < youngsModulus.size(); ++i)
	{
		file << youngsModulus[i] << " " << anisotropyStrength[i] << " "
			<< anisotropyDirections[i][0] << " " << anisotropyDirections[i][1] << " " << anisotropyDirections[i][2] << std::endl;

	}
}


void readCustomTetAttributes(std::vector<float>& youngsModulus,
	std::vector<float>& anisotropyStrength, std::vector<Eigen::Vector3f>& anisotropyDirections,
	const std::string& fileName)
{
	std::ifstream file;
	file.open(fileName);

	std::string line;
	while (std::getline(file, line))
	{
		std::vector<std::string> floatVector;
		boost::split(floatVector, line, boost::is_any_of(" "));
		youngsModulus.push_back(std::stof(floatVector[0]));
		anisotropyStrength.push_back(std::stof(floatVector[1]));
		anisotropyDirections.push_back(Eigen::Vector3f(std::stof(floatVector[2]), std::stof(floatVector[3]), std::stof(floatVector[4])));
	}
}