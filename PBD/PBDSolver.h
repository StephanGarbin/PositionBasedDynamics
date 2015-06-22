#pragma once

#include <vector>
#include <memory>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "PBDSolverSettings.h"

#include <boost/thread.hpp>

class PBDSolver
{
public:
	void advanceSystem(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
		std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences);
	PBDSolver();

	~PBDSolver();

private:

	void advanceVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void advancePositions(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraints(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraintsSOR(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
		std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences);

	void updateVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);


private:
	bool correctInversion(Eigen::Matrix3f& F, 
		Eigen::Matrix3f& FTransposeF,
		Eigen::Matrix3f& FInverseTranspose, Eigen::Matrix3f& PF,
		Eigen::Matrix3f& U, Eigen::Matrix3f& V,
		float I1, float I2, float logI3,
		const PBDSolverSettings& settings);
};


struct mutexStruct
{
	boost::mutex m_mutexSOR;
};

void projectConstraintsSOR_CORE(mutexStruct& sorMutex, std::vector<PBDTetrahedra3d>& tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
	std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences,
	int start, int end);

