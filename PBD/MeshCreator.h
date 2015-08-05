#pragma once

#include <vector>
#include <memory>

#include <Eigen\Dense>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"

//NOTE: The core code (i.e. what is found int the function 'generateTetBar') is DIRECTLY from Bender et al. (the 'Positions-Based Simulation of Continuous Materials' paper)
class MeshCreator
{
public:

	static void generateTetBar(std::shared_ptr<std::vector<PBDParticle>>& particles, std::vector<PBDTetrahedra3d>& tets,
		int width, int height, int depth);

	static void generateTetBarToFit(std::shared_ptr<std::vector<PBDParticle>>& particles, std::vector<PBDTetrahedra3d>& tets,
		int numDivsWidth, int numDivsHeight, int numDivsDepth,
		Eigen::Vector2f bottomLeft_above, Eigen::Vector2f topLeft_above, Eigen::Vector2f topRight_above,
		float thickness);

private:
	MeshCreator();
	~MeshCreator();
};

