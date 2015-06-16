#pragma once

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

	float& inverseMass() { return m_inverseMass;  }

	void swapStates();

private:

	float m_inverseMass;

	Eigen::Vector3f m_previousPosition;
	Eigen::Vector3f m_previousVelocity;

	Eigen::Vector3f m_position;
	Eigen::Vector3f m_velocity;
};

