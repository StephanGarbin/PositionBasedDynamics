#include "PBDParticle.h"

PBDParticle::PBDParticle()
{
	m_inverseMass = 0.0;
}

PBDParticle::PBDParticle(const Eigen::Vector3f& position, const Eigen::Vector3f& velocity, float inverseMass)
{
	m_position = position;
	m_previousPosition = position;
	m_pastPosition = position;
	m_velocity = velocity;
	m_previousVelocity = velocity;
	m_pastVelocity = velocity;
	m_inverseMass = inverseMass;
}


PBDParticle::~PBDParticle()
{
}

void
PBDParticle::swapStates()
{
	m_pastPosition = m_previousPosition;
	m_pastVelocity = m_previousVelocity;
	m_previousPosition = m_position;
	m_previousVelocity = m_velocity;
}
