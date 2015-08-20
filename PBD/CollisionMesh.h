#pragma once

#include <vector>
#include <string>
#include <memory>

#include <Eigen\Dense>

#include "PBDParticle.h"
#include "AbcReader.h"

class CollisionMesh
{
public:
	CollisionMesh();
	~CollisionMesh();

	void readFromAbc(const std::string& fileName);

	void resolveParticleCollisions(std::vector<PBDParticle>& particles);

	Eigen::Vector3f& getCollisionMeshTranslation()
	{
		return m_translation;
	}

private:
	std::shared_ptr<AbcReader> m_reader;
	Eigen::Vector3f m_translation;
};

