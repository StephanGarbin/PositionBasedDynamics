#pragma once
#include <vector>
#include <Eigen\Dense>


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

	void project();

	void initialise(float radius);

private:
	Eigen::Matrix3f m_covariance;

	std::vector<int> m_particleInfluences;
	std::vector<float> m_initialDistances;
};

