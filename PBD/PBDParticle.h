#pragma once

#include <vector>
#include <Eigen/Dense>

class PBDParticle
{
public:
	PBDParticle();
	PBDParticle(const Eigen::Vector3f& position, const Eigen::Vector3f& velocity, float inverseMass);
	~PBDParticle();

	Eigen::Vector3f& position() { return m_position; }
	Eigen::Vector3f& velocity() { return m_velocity; }

	Eigen::Vector3f& previousPosition() { return m_previousPosition; }
	Eigen::Vector3f& previousVelocity() { return m_previousVelocity; }

	Eigen::Vector3f& pastPosition() { return m_pastPosition; }
	Eigen::Vector3f& pastVelocity() { return m_pastVelocity; }

	float& inverseMass() { return m_inverseMass;  }

	void swapStates();

	int getNumContainingTetrahedra() { return m_tetrahedra.size(); }

	std::vector<int>& getContainingTetIdxs(){ return m_tetrahedra; }

private:
	std::vector<int> m_tetrahedra;
	float m_inverseMass;

	Eigen::Vector3f m_pastPosition;
	Eigen::Vector3f m_pastVelocity;

	Eigen::Vector3f m_previousPosition;
	Eigen::Vector3f m_previousVelocity;

	Eigen::Vector3f m_position;
	Eigen::Vector3f m_velocity;
};

