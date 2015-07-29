#pragma once

#include <vector>

class ConstraintsIO
{
public:
	ConstraintsIO();
	~ConstraintsIO();

	static bool readMayaVertexConstraints(std::vector<int>& vertexIds, const std::string& fileName);
};

