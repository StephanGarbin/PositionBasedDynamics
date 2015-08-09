#pragma once

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "PBDSolverSettings.h"
#include "PBDProbabilisticConstraint.h"

#include <boost/thread.hpp>

class PBDSolver
{
public:
	void advanceSystem(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
		std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences,
		std::vector<PBDProbabilisticConstraint>& probabilisticConstraints);
	PBDSolver();

	~PBDSolver();

	void advanceVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void advancePositions(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraints(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraintsGeometricInversionHandling(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraintsOLD(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraintsVISCOELASTIC(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
		std::vector<PBDProbabilisticConstraint>& probabilisticConstraints);

	void projectConstraintsSOR(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
		std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences);

	void updateVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	float calculateTotalStrainEnergy(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings, int it,
		std::ofstream& file);

	void projectConstraintsDistance(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, int numIterations, float k);

	void projectConstraintsVolume(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, int numIterations, float k);

	void projectConstraintsNeoHookeanMixture(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
		std::vector<PBDProbabilisticConstraint>& probabilisticConstraints);

	void projectConstraintsMooneyRivlin(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	bool correctInversion(Eigen::Matrix3f& F, 
		Eigen::Matrix3f& FTransposeF,
		Eigen::Matrix3f& FInverseTranspose, Eigen::Matrix3f& PF,
		Eigen::Matrix3f& U, Eigen::Matrix3f& V,
		float I1, float I2, float logI3,
		float& strainEnergy, float volume,
		const PBDSolverSettings& settings);

	void computeGreenStrainAndPiolaStress(const Eigen::Matrix3f &F,
		const float restVolume,
		const float mu, const float lambda, Eigen::Matrix3f &epsilon, Eigen::Matrix3f &sigma, float &energy);

	inline void computeDeltaXPositionConstraint(float w1, float w2, float restDistance,
		const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, Eigen::Vector3f& temp, Eigen::Vector3f& deltaX);

	int m_currentFrame;
};

void
PBDSolver::computeDeltaXPositionConstraint(float w1, float w2, float restDistance,
const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, Eigen::Vector3f& temp, Eigen::Vector3f& deltaX)
{
	temp = x1 - x2;

	float squaredNorm = temp.squaredNorm();

	if (w1 + w2 == 0.0f || squaredNorm == 0.0f)
	{
		deltaX.setZero();
	}
	else
	{
		deltaX = (squaredNorm - restDistance) * (temp.normalized()) / (w1 + w2);
	}
}


struct mutexStruct
{
	boost::mutex m_mutexSOR;
};

void projectConstraintsSOR_CORE(mutexStruct& sorMutex, std::vector<PBDTetrahedra3d>& tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
	std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences,
	int start, int end);

void computeGreenStrainAndPiolaStressInversion(const Eigen::Matrix3f& F, const Eigen::Matrix3f& FTransposeF,
	Eigen::Matrix3f& U, Eigen::Matrix3f& V,
	const float restVolume,
	const float mu, const float lambda, Eigen::Matrix3f &epsilon, Eigen::Matrix3f &sigma, float &energy, int it);

