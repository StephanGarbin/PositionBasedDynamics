#pragma once
#include <vector>
#include <memory>
#include <Eigen\Dense>

#include "PBDParticle.h"


class PBDProbabilisticConstraint
{
public:
	PBDProbabilisticConstraint();
	~PBDProbabilisticConstraint();

	void updateCovariance(const Eigen::Matrix3f& newCovariance)
	{
		m_covariance = newCovariance;
	}

	Eigen::Matrix3f& getCovariance()
	{
		return m_covariance;
	}

	void project(std::vector<PBDParticle>& particles);

	void initialise(std::vector<PBDParticle>& particles, float radius);

	Eigen::Vector3f& getConstraintPosition()
	{
		return m_constraintPosition;
	}

	float getInitialRadius()
	{
		return m_initialRadius;
	}

private:

	Eigen::Vector3f m_constraintPosition;

	Eigen::Matrix3f m_covariance;

	std::vector<int> m_particleInfluences;
	std::vector<float> m_initialDistances;

	float m_initialRadius;
};

