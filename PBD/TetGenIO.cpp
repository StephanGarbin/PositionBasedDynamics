#include "TetGenIO.h"

#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

TetGenIO::TetGenIO()
{
}


TetGenIO::~TetGenIO()
{
}


bool
TetGenIO::readNodes(const std::string& fileName, std::vector<PBDParticle>& particles,
	float inverseMass, Eigen::Vector3f velocity)
{
	std::ifstream file;
	file.open(fileName);

	if (!file.is_open())
	{
		std::cout << "ERROR: Could not open file " << fileName << std::endl;
		return false;
	}
	else
	{
		std::cout << "READING FILE: " << fileName << std::endl;
	}

	std::string currentLine;

	Eigen::Vector3f position;

	//ignore first line
	std::getline(file, currentLine);

	while (std::getline(file, currentLine))
	{
		removeUnnecessarySpaces(currentLine);
		if (currentLine[0] == '#')
		{
			std::cout << "Ignored Comment: " << currentLine << std::endl;
			continue;
		}

		std::vector<std::string> inputs;

		boost::split(inputs, currentLine, boost::is_any_of(" "));

		position.x() = std::stod(inputs[1]);
		position.y() = std::stod(inputs[2]);
		position.z() = std::stod(inputs[3]);

		particles.emplace_back(position, velocity, inverseMass);
	}

	std::cout << "Read " << particles.size() << " nodes." << std::endl;
	return true;
}

bool
TetGenIO::readTetrahedra(const std::string& fileName, std::vector<PBDTetrahedra3d>& tetrahedra,
	const std::shared_ptr<std::vector<PBDParticle>>& particles)
{
	std::ifstream file;
	file.open(fileName);

	if (!file.is_open())
	{
		std::cout << "ERROR: Could not open file " << fileName << std::endl;
		return false;
	}
	else
	{
		std::cout << "READING FILE: " << fileName << std::endl;
	}

	std::string currentLine;


	//ignore first line
	std::getline(file, currentLine);

	while (std::getline(file, currentLine))
	{
		//std::cout << currentLine << std::endl;

		removeUnnecessarySpaces(currentLine);
		if (currentLine[0] == '#')
		{
			std::cout << "Ignored Comment: " << currentLine << std::endl;
			continue;
		}

		std::vector<int> vertexIndices(4);

		std::vector<std::string> inputs;
		boost::split(inputs, currentLine, boost::is_any_of(" "));

		//std::cout << inputs.size() << std::endl;

		//for (int r = 0; r < inputs.size(); ++r)
		//{
		//	std::cout << inputs[r] << std::endl;
		//}

		vertexIndices[0] = std::stoi(inputs[1]) - 1;
		vertexIndices[1] = std::stoi(inputs[4]) - 1;
		vertexIndices[2] = std::stoi(inputs[2]) - 1;
		vertexIndices[3] = std::stoi(inputs[3]) - 1;

		//std::cout << "Converted Indices" << std::endl;

		tetrahedra.emplace_back(std::move(vertexIndices), particles);
		//std::cout << "done: " << currentLine << std::endl;
	}

	std::cout << "Read " << tetrahedra.size() << " tets. " << std::endl;
	return true;
}	


void
TetGenIO::removeUnnecessarySpaces(std::string& input)
{
	bool change = true;
	std::string previous = input;

	while (change)
	{
		boost::replace_all(input, "  ", " ");

		if (previous == input)
		{
			change = false;
		}
		else
		{
			previous = input;
		}
	}
}