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
#include "CollisionMesh.h"

#include <boost/thread.hpp>

class PBDSolver
{
public:
	void advanceSystem(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings,
		std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences,
		std::vector<PBDProbabilisticConstraint>& probabilisticConstraints,
		std::vector<CollisionMesh>& collisionGeometry);
	PBDSolver();

	~PBDSolver();

	void advanceVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings);

	void advancePositions(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings);

	void projectConstraintsVISCOELASTIC(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings,
		std::vector<PBDProbabilisticConstraint>& probabilisticConstraints,
		std::vector<CollisionMesh>& collisionGeometry);

	void updateVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings);

	float calculateTotalStrainEnergy(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings, int it,
		std::ofstream& file);

	void projectConstraintsDistance(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, int numIterations, float k);

	void projectConstraintsVolume(std::vector<PBDTetrahedra3d>& tetrahedra,
		std::shared_ptr<std::vector<PBDParticle>>& particles, int numIterations, float k);

	bool correctInversion(Eigen::Matrix3f& F, 
		Eigen::Matrix3f& FTransposeF,
		Eigen::Matrix3f& FInverseTranspose, Eigen::Matrix3f& PF,
		Eigen::Matrix3f& U, Eigen::Matrix3f& V,
		float I1, float I2, float logI3,
		float& strainEnergy, float volume,
		const PBDSolverSettings& settings);

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

inline int smallestDistanceToOppositePlane(PBDTetrahedra3d& tet, float& resultDist, Eigen::Vector3f& normal)
{
	std::vector<Eigen::Vector3f> normals;

	//1. 4-2, 4-3
	normals.push_back((tet.get_x(3).position() - tet.get_x(2).position()).cross(tet.get_x(3).position() - tet.get_x(1).position()).normalized());
	Eigen::Vector3f w = (tet.get_x(3).position() - tet.get_x(0).position());
	float dist0 = normals[0].dot(w);


	normals.push_back((tet.get_x(3).position() - tet.get_x(2).position()).cross(tet.get_x(3).position() - tet.get_x(0).position()).normalized());
	w = (tet.get_x(3).position() - tet.get_x(1).position());
	float dist1 = normals[1].dot(w);

	normals.push_back((tet.get_x(3).position() - tet.get_x(0).position()).cross(tet.get_x(3).position() - tet.get_x(1).position()).normalized());
	w = (tet.get_x(3).position() - tet.get_x(2).position());

	float dist2 = normals[2].dot(w);

	normals.push_back((tet.get_x(1).position() - tet.get_x(0).position()).cross(tet.get_x(1).position() - tet.get_x(2).position()).normalized());
	w = (tet.get_x(3).position() - tet.get_x(1).position());

	float dist3 = normals[3].dot(w);

	std::vector<float> distances;
	if (tet.get_x(0).inverseMass() != 0.0f)
	{
		distances.push_back(dist0);
	}
	if (tet.get_x(1).inverseMass() != 0.0f)
	{
		distances.push_back(dist1);
	}
	if (tet.get_x(2).inverseMass() != 0.0f)
	{
		distances.push_back(dist2);
	}
	if (tet.get_x(3).inverseMass() != 0.0f)
	{
		distances.push_back(dist3);
	}

	std::sort(distances.begin(), distances.end());

	//std::cout << dist0 << ", " << dist1 << ", " << dist2 << ", " << dist3 << std::endl;
	if (distances[0] == dist0)
	{
		resultDist = dist0;
		normal = normals[0];
		return 0;
	}
	else if (distances[0] == dist1)
	{
		resultDist = dist1;
		normal = normals[1];
		return 1;
	}
	else if (distances[0] == dist2)
	{
		resultDist = dist2;
		normal = normals[2];
		return 2;
	}
	else if (distances[0] == dist3)
	{
		resultDist = dist3;
		normal = normals[3];
		return 3;
	}
}