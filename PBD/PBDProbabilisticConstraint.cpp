#include "PBDProbabilisticConstraint.h"
#include "PBDSolver.h"

#include <iostream>


PBDProbabilisticConstraint::PBDProbabilisticConstraint()
{
}


PBDProbabilisticConstraint::~PBDProbabilisticConstraint()
{
}

void
PBDProbabilisticConstraint::project(std::vector<PBDParticle>& particles)
{
	//correct particle positions (assuming that our current point is fixed
	Eigen::Vector3f temp;
	Eigen::Vector3f deltaX;

	PBDSolver solver;

	float w1 = 0.0;
	Eigen::Vector3f p1 = m_constraintPosition;

	float w2;
	Eigen::Vector3f p2;
	for (int i = 0; i < m_particleInfluences.size(); ++i)
	{
		w2 = particles[m_particleInfluences[i]].inverseMass();
		p2 = particles[m_particleInfluences[i]].position();

		if ((p1 - p2).squaredNorm() <= m_initialDistances[i])
		{
			continue;
		}

		solver.computeDeltaXPositionConstraint(w1, w2, m_initialDistances[i], p1, p2, temp, deltaX);

		particles[m_particleInfluences[i]].position() += deltaX * w2;
	}
}


void
PBDProbabilisticConstraint::initialise(std::vector<PBDParticle>& particles, float radius)
{
	m_initialRadius = radius;

	std::cout << "Initialising Prob Constraint" << std::endl;
	std::cout << m_constraintPosition << std::endl;
	for (int p = 0; p < particles.size(); ++p)
	{
		float distance = (particles[p].position() - m_constraintPosition).squaredNorm();
		//std::cout << distance << std::endl;
		if (distance <= radius)
		{
			m_particleInfluences.push_back(p);
			m_initialDistances.push_back(distance);
		}
	}

	std::cout << "Added " << m_particleInfluences.size() << " particles to probabilistic constraint." << std::endl;
}