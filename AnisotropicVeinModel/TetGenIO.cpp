#include "TetGenIO.h"

#include <fstream>
#include <iostream>
#include <string>

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
		removeWhiteSpaceAtBeginning(currentLine);

		std::vector<std::string> inputs;

		int add = 0;

		boost::split(inputs, currentLine, boost::is_any_of(" "));
		//std::cout << inputs[1] << "; " << inputs[2] << "; " << inputs[3] << std::endl;

		position.x() = std::stod(inputs[1 + add]);
		position.y() = std::stod(inputs[2 + add]);
		position.z() = std::stod(inputs[3 + add]);

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

		removeWhiteSpaceAtBeginning(currentLine);

		std::vector<int> vertexIndices(5);

		std::vector<std::string> inputs;
		boost::split(inputs, currentLine, boost::is_any_of(" "));

		//std::cout << inputs.size() << std::endl;

		//for (int r = 0; r < inputs.size(); ++r)
		//{
		//	std::cout << inputs[r] << std::endl;
		//}

		int add = 0;

		//vertexIndices[0] = std::stoi(inputs[1 + add]);
		//vertexIndices[1] = std::stoi(inputs[4 + add]);
		//vertexIndices[2] = std::stoi(inputs[2 + add]);
		//vertexIndices[3] = std::stoi(inputs[3 + add]);
		vertexIndices[0] = std::stoi(inputs[1 + add]);
		vertexIndices[1] = std::stoi(inputs[4 + add]);
		vertexIndices[2] = std::stoi(inputs[2 + add]);
		vertexIndices[3] = std::stoi(inputs[3 + add]);

		//std::cout << inputs[1] << "; " << inputs[2] << "; " << inputs[3] << "; " << inputs[4] << std::endl;

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

void
TetGenIO::removeWhiteSpaceAtBeginning(std::string& input)
{
	while (input[0] == ' ')
	{
		input.erase(input.begin());
	}
}