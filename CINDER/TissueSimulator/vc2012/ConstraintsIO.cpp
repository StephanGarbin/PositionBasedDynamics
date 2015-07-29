#include "ConstraintsIO.h"

#include <fstream>
#include <iostream>
#include <string>

ConstraintsIO::ConstraintsIO()
{
}


ConstraintsIO::~ConstraintsIO()
{
}


bool
ConstraintsIO::readMayaVertexConstraints(std::vector<int>& vertexIds, const std::string& fileName)
{
	std::ifstream file;
	file.open(fileName);

	if (!file.is_open())
	{
		std::cout << "ERROR: Could not read constraints from file: " << fileName << std::endl;
		return false;
	}

	std::string currentLine;

	while (std::getline(file, currentLine))
	{
		for (int i = 0; i < currentLine.size(); ++i)
		{
			if (currentLine[i] == '[')
			{
				currentLine = currentLine.substr(i + 1, currentLine.size() - 1);
				break;
			}
		}

		for (int i = 0; i < currentLine.size(); ++i)
		{
			if (currentLine[i] == '[')
			{
				int endIdx = 1;
				while (currentLine[i + endIdx] != ']')
				{
					++endIdx;
				}

				vertexIds.push_back(std::atoi(currentLine.substr(i + 1, endIdx - 1).c_str()));

				i = i + endIdx;
			}
		}
	}

	std::cout << "Read " << vertexIds.size() << " constraints from disk." << std::endl;
	return true;
}
