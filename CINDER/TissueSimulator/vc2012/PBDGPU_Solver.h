#pragma once

#include <vector>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "Parameters.h"

class PBDGPU_Solver
{
public:
	PBDGPU_Solver();
	~PBDGPU_Solver();

	void setup(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles);

	void advanceSystem(std::shared_ptr<std::vector<PBDParticle>>& particles,
		Parameters& settings);

private:

	void advanceVelocities(std::shared_ptr<std::vector<PBDParticle>>& particles,
		const Parameters& settings);

	void advancePositions(std::shared_ptr<std::vector<PBDParticle>>& particles,
		const Parameters& settings);

	void updateVelocities(std::shared_ptr<std::vector<PBDParticle>>& particles,
		const Parameters& settings);

	void setupCUDA();
	void deleteCUDA();

	float m_isSetup;

	std::vector<float> m_inverseMasses;
	std::vector<int> m_indices;
	std::vector<float> m_undeformedVolumes;
	std::vector<float> m_referenceShapeMatrices;
	std::vector<float> m_positions;
	std::vector<float> m_F;
	std::vector<float> m_U;
	std::vector<float> m_V;

	//CUDA device ptrs
	int* m_device_indices;
	float* m_device_positions;
	float* m_device_undeformedVolumes;
	float* m_device_referenceShapeMatrices;
	float* m_device_inverseMasses;
	float* m_device_F;
	float* m_device_U;
	float* m_device_V;

	int m_CUDA_NUM_BLOCKS;
	int m_CUDA_NUM_THREADS_PER_BLOCK;
	int m_CUDA_NUM_PARTICLES;
	int m_CUDA_TRUE_NUM_CONSTRAINTS;

	void determineCUDALaunchParameters(int numParticles);

	void printError(const std::string& message);
};

