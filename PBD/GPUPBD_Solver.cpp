#include "GPUPBD_Solver.h"

#include <iostream>
#include <memory>
#include <vector>

#include "CUDAPBD_SolverSettings.h"
#include "CUDA_WRAPPER.h"

GPUPBD_Solver::GPUPBD_Solver()
{
	m_isSetup = false;
}


GPUPBD_Solver::~GPUPBD_Solver()
{
	std::cout << "Destroying CUDA Solver..." << std::endl;
}

void
GPUPBD_Solver::determineCUDALaunchParameters(int numParticles)
{
	CUDA_TRUE_NUM_CONSTRAINTS = numParticles;

	CUDA_NUM_THREADS_PER_BLOCK = 64;

	CUDA_NUM_BLOCKS = (numParticles / CUDA_NUM_THREADS_PER_BLOCK) + 1;

	CUDA_NUM_PARTICLES = CUDA_NUM_BLOCKS * CUDA_NUM_THREADS_PER_BLOCK;

	std::cout << "Determined (from " << numParticles << " tendered tetrahedra):" << std::endl;
	std::cout << "	NUM_BLOCKS           : " << CUDA_NUM_BLOCKS << std::endl;
	std::cout << "	NUM_THREADS_PER_BLOCK: " << CUDA_NUM_THREADS_PER_BLOCK << std::endl;
	std::cout << "	This adds " << CUDA_NUM_PARTICLES - numParticles << " to the solver." << std::endl;
}

void
GPUPBD_Solver::setup(std::vector<PBDTetrahedra3d>& tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& particles)
{
	//0. Determine CUDA Launch parameters
	determineCUDALaunchParameters(tetrahedra.size());

	//1. Inverse masses
	for (int i = 0; i < particles->size(); ++i)
	{
		m_inverseMasses.push_back((*particles)[i].inverseMass());
	}


	//2. Indices
	for (int i = 0; i < tetrahedra.size(); ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			m_indices.push_back(tetrahedra[i].getVertexIndices()[j]);
		}
	}

	
	//3. Undeformed Volume
	for (int i = 0; i < tetrahedra.size(); ++i)
	{
		m_undeformedVolumes.push_back(tetrahedra[i].getUndeformedVolume());
	}

	//4. Reference Shape Matrices
	for (int i = 0; i < tetrahedra.size(); ++i)
	{
		for (int row = 0; row < 3; ++row)
		{
			for (int col = 0; col < 3; ++col)
			{
				m_referenceShapeMatrices.push_back(tetrahedra[i].getReferenceShapeMatrixInverseTranspose()(row, col));
			}
		}
	}

	m_positions.resize(m_inverseMasses.size());

	queryCUDADevices();

	std::cout << "CUDA Solver successfully initialised..." << std::endl;
	m_isSetup = true;
}

void
GPUPBD_Solver::advanceSystem(std::shared_ptr<std::vector<PBDParticle>>& particles,
	const PBDSolverSettings& settings)
{
	//check that the system is set up correctly
	if (!m_isSetup)
	{
		printError("Cannot advance system, setup() must be called first!");
	}

	//1. Get Positions
	for (int i = 0; i < particles->size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			m_positions[i * 3 + j] = (*particles)[i].position()[j];
		}
	}

	//2. Determine Settings
	CUDAPBD_SolverSettings solverSettings;
	solverSettings.lambda = settings.lambda;
	solverSettings.mu = settings.mu;
	solverSettings.numBlocks = CUDA_NUM_BLOCKS;
	solverSettings.numThreadsPerBlock = CUDA_NUM_THREADS_PER_BLOCK;
	solverSettings.numIterations = settings.numConstraintIts;
	solverSettings.trueNumberOfConstraints = CUDA_TRUE_NUM_CONSTRAINTS;

	//3. Advance System
	std::cout << "Solving System with CUDA..." << std::endl;
	CUDA_projectConstraints(m_indices, m_positions, m_inverseMasses,
		m_referenceShapeMatrices, m_undeformedVolumes, solverSettings);

	//4. Copy Positions back
	std::cout << "Copying Positions back to particles..." << std::endl;

	for (int i = 0; i < particles->size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			(*particles)[i].position()[j] = m_positions[i * 3 + j];
		}
	}

	std::cout << "...done" << std::endl;
}

void
GPUPBD_Solver::printError(const std::string& message)
{
	std::cerr << "GPUPBD Solver Error: " << message << std::endl;
}
