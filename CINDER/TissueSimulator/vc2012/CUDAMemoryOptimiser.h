#pragma once

#include <vector>
#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"

class CUDAMemoryOptimiser
{
public:
	static void optimiseTetrahedralIndexingBasedOnNodeMemory(std::vector<PBDParticle>& nodes,
		std::vector<PBDTetrahedra3d>& tetrahedra);

	static void optimiseTetrahedralIndexingBasedOnMemoryChunks(std::vector<PBDParticle>& nodes,
		std::vector<PBDTetrahedra3d>& tetrahedra, int numThreadsPerBlock);

	static void sortNodesBasedOnConnectivity(std::vector<PBDParticle>& nodes,
		std::vector<PBDTetrahedra3d>& tetrahedra);

private:

	static int calculateMemoryAlignmentScore(const PBDTetrahedra3d& reference,
		const PBDTetrahedra3d& target);
	static int calculateMemoryAlignmentScore(const std::vector<PBDTetrahedra3d>& order);
	CUDAMemoryOptimiser();
	~CUDAMemoryOptimiser();
};

