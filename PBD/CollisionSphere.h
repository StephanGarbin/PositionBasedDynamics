#pragma once

#include <vector>
#include <string>
#include <memory>

#include <Eigen\Dense>

#include "PBDParticle.h"
#include "AbcReaderTransform.h"

class CollisionSphere
{
public:
	CollisionSphere();
	~CollisionSphere();

	void readFromAbc(const std::string& fileName, const std::string& transformName);

	void resolveParticleCollisions(std::vector<PBDParticle>& particles, int systemFrame, float timeStep,
		float sphereRadius);

	void calculateNewSphereCentre(int systemFrame, float timeStep);

	void resolveParticleCollisions_SAFE(std::vector<PBDParticle>& particles, int systemFrame, float timeStep,
		float sphereRadius,
		int start, int end);

	void checkForSinglePointIntersection_SAFE(const Eigen::Vector3f& point, float& penetrationDistance, bool& penetrates, float sphereRadius);

	Eigen::Vector3f& getCollisionMeshTranslation()
	{
		return m_translation;
	}

	void glRender(int systemFrame, float timeStep, float sphereRadius);

	int& getFrameLimit(){ return m_frameLimit; }
private:

	void computeDeltaXPositionConstraint(float w1, float w2, float restDistance,
		const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, Eigen::Vector3f& temp, Eigen::Vector3f& deltaX);

	std::shared_ptr<AbcReaderTransform> m_reader;
	Eigen::Vector3f m_translation;

	Eigen::Vector3f m_collisionSphereCentre;

	int m_frameLimit;

	int m_lastProcessedFrame;

	Eigen::Vector3f m_previousCollisionSphereCentre;
};

