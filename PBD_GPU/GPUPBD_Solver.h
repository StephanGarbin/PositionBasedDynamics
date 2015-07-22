#pragma once

#include <vector>

#include "PBDParticle.h";
#include "PBDTetrahedra3d.h";
#include "PBDSolverSettings.h"

class GPUPBD_Solver
{
public:
	GPUPBD_Solver();
	~GPUPBD_Solver();

	void setup(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles);

	void advanceSystem(std::shared_ptr<std::vector<PBDParticle>>& particles,
		const PBDSolverSettings& settings);

private:
	float m_isSetup;

	std::vector<float> m_inverseMasses;
	std::vector<int> m_indices;
	std::vector<float> m_undeformedVolumes;
	std::vector<float> m_referenceShapeMatrices;

	std::vector<float> m_positions;

	int CUDA_NUM_BLOCKS;
	int CUDA_NUM_THREADS_PER_BLOCK;
	int CUDA_NUM_PARTICLES;
	int CUDA_TRUE_NUM_CONSTRAINTS;

	void determineCUDALaunchParameters(int numParticles);

	void printError(const std::string& message);
};

