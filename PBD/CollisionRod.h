#pragma once

#include <vector>
#include <string>
#include <memory>

#include <Eigen\Dense>

#include "PBDParticle.h"
#include "AbcReaderTransform.h"

class CollisionRod
{
public:
	CollisionRod();
	~CollisionRod();

	void readFromAbc(const std::string& fileName, const std::vector<std::string>& topBottomTransformNames);

	void resolveParticleCollisions(std::vector<PBDParticle>& particles, int systemFrame, float timeStep,
		int numSpheres, float sphereRadius);

	Eigen::Vector3f& getCollisionMeshTranslation()
	{
		return m_translation;
	}

	void glRender(int systemFrame, float timeStep, int numSpheres, float sphereRadius);

private:

	void computeDeltaXPositionConstraint(float w1, float w2, float restDistance,
		const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, Eigen::Vector3f& temp, Eigen::Vector3f& deltaX);

	std::shared_ptr<AbcReaderTransform> m_reader;
	Eigen::Vector3f m_translation;

	std::vector<Eigen::Vector3f> m_collisionSphereCentres;
};
