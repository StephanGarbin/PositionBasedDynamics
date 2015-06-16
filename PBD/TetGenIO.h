#pragma once

#include <string>
#include <vector>

#include "PBDTetrahedra3d.h"
#include "PBDParticle.h"

#include <Eigen/Dense>

class TetGenIO
{
public:
	static bool readNodes(const std::string& fileName, std::vector<PBDParticle>& particles,
		float inverseMass, Eigen::Vector3f velocity);

	static bool readTetrahedra(const std::string& fileName, std::vector<PBDTetrahedra3d>& tetrahedra,
		const std::shared_ptr<std::vector<PBDParticle>>& particles);

private:

	static void removeUnnecessarySpaces(std::string& input);

	TetGenIO();
	~TetGenIO();
};

