#pragma once

#include <vector>
#include <memory>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"

//NOTE: This code is DIRECTLY from Bender et al. (the 'Positions-Based Simulation of Continuous Materials' paper)
class MeshCreator
{
public:

	static void generateTetBar(std::shared_ptr<std::vector<PBDParticle>>& particles, std::vector<PBDTetrahedra3d>& tets,
		int width, int height, int depth);

private:
	MeshCreator();
	~MeshCreator();
};

