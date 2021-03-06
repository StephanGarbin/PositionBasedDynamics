#include "PBDGPU_Solver.h"

#include <iostream>
#include <memory>
#include <vector>

#include "Parameters.h"
#include "CUDA_WRAPPER.h"
#include "CUDA_GLOBALS.h"

PBDGPU_Solver::PBDGPU_Solver()
{
	m_isSetup = false;
}


PBDGPU_Solver::~PBDGPU_Solver()
{
	std::cout << "Destroying CUDA Solver..." << std::endl;
	deleteCUDA();
}

void
PBDGPU_Solver::determineCUDALaunchParameters(int numParticles)
{
	m_CUDA_TRUE_NUM_CONSTRAINTS = numParticles;

	m_CUDA_NUM_THREADS_PER_BLOCK = NUM_THREADS_PER_BLOCK;

	m_CUDA_NUM_BLOCKS = (numParticles / m_CUDA_NUM_THREADS_PER_BLOCK) + 1;

	m_CUDA_NUM_PARTICLES = m_CUDA_NUM_BLOCKS * m_CUDA_NUM_THREADS_PER_BLOCK;

	std::cout << "Determined (from " << numParticles << " tendered tetrahedra):" << std::endl;
	std::cout << "	NUM_BLOCKS           : " << m_CUDA_NUM_BLOCKS << std::endl;
	std::cout << "	NUM_THREADS_PER_BLOCK: " << m_CUDA_NUM_THREADS_PER_BLOCK << std::endl;
	std::cout << "	This adds " << m_CUDA_NUM_PARTICLES - numParticles << " to the solver." << std::endl;
}

void
PBDGPU_Solver::setup(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles)
{
	//0. Determine CUDA Launch parameters
	determineCUDALaunchParameters(tetrahedra.size());

	//1. Inverse masses

	//for (int i = 0; i < particles->size(); ++i)
	//{
	//	m_inverseMasses.push_back((*particles)[i].inverseMass());
	//}
	for (int j = 0; j < 4; ++j)
	{
		for (int i = 0; i < tetrahedra.size(); ++i)
		{
			m_inverseMasses.push_back((*particles)[tetrahedra[i].getVertexIndices()[j]].inverseMass());
		}
	}

	//2. Indices

	//for (int i = 0; i < tetrahedra.size(); ++i)
	//{
	//	for (int j = 0; j < 4; ++j)
	//	{
	//		m_indices.push_back(tetrahedra[i].getVertexIndices()[j]);
	//	}
	//}
	for (int j = 0; j < 4; ++j)
	{
		for (int i = 0; i < tetrahedra.size(); ++i)
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

	//for (int i = 0; i < tetrahedra.size(); ++i)
	//{
	//	for (int row = 0; row < 3; ++row)
	//	{
	//		for (int col = 0; col < 3; ++col)
	//		{
	//			m_referenceShapeMatrices.push_back(tetrahedra[i].getReferenceShapeMatrixInverseTranspose().transpose()(row, col));
	//		}
	//	}
	//}

	//for (int i = 0; i < tetrahedra.size(); ++i)
	//{
		for (int row = 0; row < 3; ++row)
		{
			for (int col = 0; col < 3; ++col)
			{
				for (int i = 0; i < tetrahedra.size(); ++i)
				{
					m_referenceShapeMatrices.push_back(tetrahedra[i].getReferenceShapeMatrixInverseTranspose().transpose()(row, col));
				}
			}
		}
	//}

	m_positions.resize(m_inverseMasses.size() * 3);

	//Resize arrays for diagonalisation
	m_F.resize(tetrahedra.size() * 3);
	m_U.resize(tetrahedra.size() * 9);
	m_V.resize(tetrahedra.size() * 9);

	//Resize arrays for anisotropy & viscoelasticity
	m_anisotropyDirection.resize(tetrahedra.size() * 3);
	m_viscosityMultiplier.resize(tetrahedra.size() * 9);

	queryCUDADevices();

	setupCUDA();

	std::cout << "CUDA Solver successfully initialised..." << std::endl;
	m_isSetup = true;
}

void
PBDGPU_Solver::advanceSystem(std::shared_ptr<std::vector<PBDParticle>>& particles,
Parameters& settings)
{
	//check that the system is set up correctly
	if (!m_isSetup)
	{
		printError("Cannot advance system, setup() must be called first!");
	}

	advanceVelocities(particles, settings);
	advancePositions(particles, settings);

	//1. Get Positions
	for (int i = 0; i < particles->size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			m_positions[i * 3 + j] = (*particles)[i].position()[j];
		}
	}
	//for (int j = 0; j < 3; ++j)
	//{
	//	for (int i = 0; i < particles->size(); ++i)
	//	{
	//		m_positions[i + j * particles->size()] = (*particles)[i].position()[j];
	//	}
	//}

	//2. Determine Settings
	settings.numBlocks = m_CUDA_NUM_BLOCKS;
	settings.numThreadsPerBlock = m_CUDA_NUM_THREADS_PER_BLOCK;
	settings.trueNumberOfConstraints = m_CUDA_TRUE_NUM_CONSTRAINTS;
	settings.numParticles = particles->size();

	//3. Advance System
	settings.calculateMu();
	settings.calculateLambda();

	settings.normaliseAnisotropyDirection();
	bool updateAnisotropy = m_anisotropyDirection[0] != settings.anisotropyDirection[0]
		|| m_anisotropyDirection[1] != settings.anisotropyDirection[1]
		|| m_anisotropyDirection[2] != settings.anisotropyDirection[2];

	if (updateAnisotropy)
	{
		for (int i = 0; i < m_anisotropyDirection.size() / 3; ++i)
		{
			m_anisotropyDirection[0 * m_anisotropyDirection.size() / 3 + i] = settings.anisotropyDirection[0];
			m_anisotropyDirection[1 * m_anisotropyDirection.size() / 3 + i] = settings.anisotropyDirection[1];
			m_anisotropyDirection[2 * m_anisotropyDirection.size() / 3 + i] = settings.anisotropyDirection[2];
		}
	}

	CUDA_updateBuffers(m_device_positions, m_positions, m_device_anisotropyDirection, m_anisotropyDirection, updateAnisotropy);

	CUDA_projectConstraints(m_device_indices, m_device_positions, m_device_inverseMasses,
		m_device_referenceShapeMatrices, m_device_undeformedVolumes,
		m_device_F, m_device_U, m_device_V,
		m_device_anisotropyDirection, m_device_viscosityMultiplier,
		settings);

	CUDA_getBuffers(m_device_positions, m_positions);

	//4. Copy Positions back
	for (int i = 0; i < particles->size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			(*particles)[i].position()[j] = m_positions[i * 3 + j];
		}
	}
	//for (int j = 0; j < 3; ++j)
	//{
	//	for (int i = 0; i < particles->size(); ++i)
	//	{
	//		(*particles)[i].position()[j] = m_positions[i + j * particles->size()];
	//	}
	//}

	//Update Velocities
	updateVelocities(particles, settings);

	//swap particles states
	for (auto& p : *particles)
	{
		p.swapStates();
	}
}

void
PBDGPU_Solver::printError(const std::string& message)
{
	std::cerr << "GPUPBD Solver Error: " << message << std::endl;
}

void
PBDGPU_Solver::advanceVelocities(std::shared_ptr<std::vector<PBDParticle>>& particles,
const Parameters& settings)
{
	for (auto& p : *particles)
	{
		float temp = settings.timeStep * p.inverseMass() * settings.gravity;
		p.velocity().x() = p.previousVelocity().x() + 0;
		p.velocity().y() = p.previousVelocity().y() + temp;
		p.velocity().z() = p.previousVelocity().z() + 0;
	}
}

void
PBDGPU_Solver::advancePositions(std::shared_ptr<std::vector<PBDParticle>>& particles,
const Parameters& settings)
{
	for (auto& p : *particles)
	{
		p.position() = p.previousPosition() + settings.timeStep * p.velocity();
		//std::cout << p.velocity().transpose() << std::endl;
	}
}

void
PBDGPU_Solver::updateVelocities(std::shared_ptr<std::vector<PBDParticle>>& particles,
const Parameters& settings)
{
	for (auto& p : *particles)
	{
		p.velocity() = (1.0 / settings.timeStep) * (p.position() - p.previousPosition());
	}
}


void
PBDGPU_Solver::setupCUDA()
{
	CUDA_allocateBuffers(&m_device_indices, &m_device_positions, &m_device_inverseMasses,
		&m_device_referenceShapeMatrices, &m_device_undeformedVolumes,
		&m_device_F, &m_device_U, &m_device_V,
		&m_device_anisotropyDirection, &m_device_viscosityMultiplier,
		m_indices, m_positions, m_inverseMasses, m_referenceShapeMatrices, m_undeformedVolumes,
		m_F, m_U, m_V,
		m_anisotropyDirection, m_viscosityMultiplier);
}

void
PBDGPU_Solver::deleteCUDA()
{
	CUDA_destroyBuffers(m_device_indices, m_device_positions, m_device_inverseMasses,
		m_device_referenceShapeMatrices, m_device_undeformedVolumes,
		m_device_F, m_device_U, m_device_V,
		m_device_anisotropyDirection, m_device_viscosityMultiplier);
}

