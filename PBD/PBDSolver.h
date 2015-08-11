#pragma once

#include <vector>
#include <memory>
#include <algorithm>

#include <iostream>
#include <fstream>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "Parameters.h"
#include "PBDSolverSettings.h"
#include "PBDProbabilisticConstraint.h"

#include <boost/thread.hpp>

class PBDSolver
{
public:
	void advanceSystem(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings,
		std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences,
		std::vector<PBDProbabilisticConstraint>& probabilisticConstraints);
	PBDSolver();

	~PBDSolver();

	void advanceVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings);

	void advancePositions(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings);

	void projectConstraints(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraintsGeometricInversionHandling(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraintsOLD(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	void projectConstraintsVISCOELASTIC(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings,
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

inline int smallestDistanceToOppositePlane(PBDTetrahedra3d& tet, float& resultDist, Eigen::Vector3f& normal)
{
	//1. 4-2, 4-3
	Eigen::Vector3f n = (tet.get_x(3).position() - tet.get_x(1).position()).cross(tet.get_x(3).position() - tet.get_x(2).position()).normalized();
	Eigen::Vector3f w = (tet.get_x(3).position() - tet.get_x(0).position());
	float dist1 = n.dot(w);


	n = (tet.get_x(3).position() - tet.get_x(0).position()).cross(tet.get_x(3).position() - tet.get_x(2).position()).normalized();
	w = (tet.get_x(3).position() - tet.get_x(1).position());
	float dist2 = n.dot(w);

	n = (tet.get_x(3).position() - tet.get_x(0).position()).cross(tet.get_x(3).position() - tet.get_x(1).position()).normalized();
	w = (tet.get_x(3).position() - tet.get_x(2).position());

	float dist3 = n.dot(w);

	n = (tet.get_x(1).position() - tet.get_x(0).position()).cross(tet.get_x(1).position() - tet.get_x(2).position()).normalized();
	w = (tet.get_x(3).position() - tet.get_x(1).position());

	float dist4 = n.dot(w);

	std::vector<float> distances;
	distances.push_back(dist1); distances.push_back(dist2); distances.push_back(dist3); distances.push_back(dist4);
	std::sort(distances.begin(), distances.end());

	//std::cout << dist1 << ", " << dist2 << ", " << dist3 << ", " << dist4 << std::endl;
	if (distances[0] == dist1)
	{
		resultDist = dist1;
		normal = (tet.get_x(3).position() - tet.get_x(1).position()).cross(tet.get_x(3).position() - tet.get_x(2).position()).normalized();
		return 0;
	}
	else if (distances[0] == dist2)
	{
		resultDist = dist2;
		normal = (tet.get_x(3).position() - tet.get_x(0).position()).cross(tet.get_x(3).position() - tet.get_x(2).position()).normalized();
		return 1;
	}
	else if (distances[0] == dist3)
	{
		resultDist = dist3;
		normal = (tet.get_x(3).position() - tet.get_x(0).position()).cross(tet.get_x(3).position() - tet.get_x(1).position()).normalized();
		return 2;
	}
	else if (distances[0] == dist4)
	{
		resultDist = dist4;
		normal = (tet.get_x(1).position() - tet.get_x(0).position()).cross(tet.get_x(1).position() - tet.get_x(2).position()).normalized();
		return 3;
	}
}