#include "PBDSolver.h"

#include <iostream>
#include <fstream>

#include <tbb\parallel_for.h>
#include <tbb\mutex.h>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "EnergyConstraint.h"


PBDSolver::PBDSolver()
{
	m_currentFrame = 0;
}


PBDSolver::~PBDSolver()
{
}

void
PBDSolver::advanceSystem(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings,
std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences,
std::vector<PBDProbabilisticConstraint>& probabilisticConstraints,
std::vector<CollisionMesh>& collisionGeometry)
{
	//Advance Velocities
	advanceVelocities(tetrahedra, particles, settings);

	//Advance Positions
	advancePositions(tetrahedra, particles, settings);

	if (!settings.disableConstraintProjection)
	{
		//Project Constraints
		//if (!settings.useSOR)
		//{
		//	projectConstraints(tetrahedra, particles, settings);
		//}
		//else
		//{
		//	projectConstraintsSOR(tetrahedra, particles, settings, temporaryPositions, numConstraintInfluences);
		//}
		//projectConstraintsOLD(tetrahedra, particles, settings);
		//projectConstraintsNeoHookeanMixture(tetrahedra, particles, settings, probabilisticConstraints);
		//projectConstraintsMooneyRivlin(tetrahedra, particles, settings);
		//projectConstraintsDistance(tetrahedra, particles, settings.numConstraintIts, settings.youngsModulus);
		//projectConstraintsGeometricInversionHandling(tetrahedra, particles, settings);
		//projectConstraintsVolume(tetrahedra, particles, settings.numConstraintIts, settings.youngsModulus);
		projectConstraintsVISCOELASTIC(tetrahedra, particles, settings, probabilisticConstraints, collisionGeometry);
	}

	//Update Velocities
	updateVelocities(tetrahedra, particles, settings);

	//swap particles states
	for (auto& p : *particles)
	{
		//std::cout << p.previousPosition() - p.position() << std::endl;
		p.swapStates();
	}

	++m_currentFrame;
}


void
PBDSolver::advanceVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings)
{
	for (auto& p : *particles)
	{
		//gravity
		float temp = settings.deltaT * p.inverseMass() * settings.gravity;
		p.velocity().x() = p.previousVelocity().x() + 0;
		p.velocity().y() = p.previousVelocity().y() + temp;
		p.velocity().z() = p.previousVelocity().z() + 0;

		//other external forces
		p.velocity() += settings.deltaT * settings.forceMultiplicationFactor * settings.externalForce;
	}
}

void
PBDSolver::advancePositions(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings)
{
	for (auto& p : *particles)
	{
		p.position() = p.previousPosition() + settings.deltaT * p.velocity() * p.inverseMass();
		//std::cout << p.velocity().transpose() << std::endl;
	}

	//if (settings.currentFrame == 1)
	//{
	//	for (auto& p : *particles)
	//	{
	//		Eigen::Vector3f temp(0.0f, 1.5f, 0.0f);
	//		//p.position() += settings.deltaT * settings.positionDeltaMultiplicationFactor * settings.externalPositionDelta;
	//		//p.position() += settings.getCurrentTime() * settings.positionDeltaMultiplicationFactor * settings.externalPositionDelta;
	//		p.position() = temp;
	//		//std::cout << p.velocity().transpose() << std::endl;
	//	}
	//}
}

void
PBDSolver::updateVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	for (auto& p : *particles)
	{
		p.velocity() = (1.0 / settings.deltaT) * (p.position() - p.previousPosition());
	}
}

float
PBDSolver::calculateTotalStrainEnergy(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings, int it,
std::ofstream& file)
{
	Eigen::Matrix3f F;
	Eigen::Matrix3f FInverseTranspose;
	Eigen::Matrix3f FTransposeF;

	Eigen::Matrix3f PF;
	Eigen::Matrix3f gradientTemp;
	Eigen::MatrixXf gradient; gradient.resize(3, 4);

	Eigen::Matrix3f U;
	Eigen::Matrix3f V;

	Eigen::Matrix3f sigma;
	Eigen::Matrix3f epsilon;

	Eigen::Vector3f deltaX;

	float strainEnergyTotal = 0.0;

	for (int t = 0; t < tetrahedra.size(); ++t)
	{
		float lagrangeM;
		float strainEnergy;

		//Compute volume
		float Volume = tetrahedra[t].getUndeformedVolume();

		//Get deformation gradient
		F = tetrahedra[t].getDeformationGradient();

		if (F.isIdentity())
		{
			continue;
		}

		FInverseTranspose = F.inverse().transpose();
		FTransposeF = F.transpose() * F;

		//Compute Isotropic Invariants
		float I1 = (FTransposeF).trace();
		float I3 = (FTransposeF).determinant();

		float logI3 = log(I3);

		//computeGreenStrainAndPiolaStressInversion(F, FTransposeF, U, V, Volume, settings.mu, settings.lambda, PF, epsilon, strainEnergy, 100);

		strainEnergyTotal += strainEnergy;
	}

	if (settings.printStrainEnergy)
	{
		if (it < 10)
		{
			std::cout << "Strain Energy Before [  " << it << "]: "
				<< strainEnergyTotal << std::endl;
		}
		else if (it < 100)
		{
			std::cout << "Strain Energy Before [ " << it << "]: "
				<< strainEnergyTotal << std::endl;
		}
		else if (it < 1000)
		{
			std::cout << "Strain Energy Before [" << it << "]: "
				<< strainEnergyTotal << std::endl;
		}
	}

	if (settings.printStrainEnergyToFile)
	{
		file << strainEnergyTotal << std::endl;
	}

	return strainEnergyTotal;
}

bool
PBDSolver::correctInversion(Eigen::Matrix3f& F,
Eigen::Matrix3f& FTransposeF,
Eigen::Matrix3f& FInverseTranspose, Eigen::Matrix3f& PF,
Eigen::Matrix3f& U, Eigen::Matrix3f& V,
float I1, float I2, float logI3,
float& strainEnergy, float volume,
const PBDSolverSettings& settings)
{
	//Compute Eigendecomposition of F
	Eigen::Vector3f S;

	Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);
	U = eigenSolver.pseudoEigenvalueMatrix();
	S[0] = U(0, 0);
	S[1] = U(1, 1);
	S[2] = U(2, 2);
	V = eigenSolver.pseudoEigenvectors();

	//Make sure all eigenvalues are > 0
	for (int i = 0; i < 3; ++i)
	{
		if (S(i, i) < 0.0f)
		{
			S(i, i) = 0.0f;
		}
	}

	//if det V is smaller than 0
	if (V.determinant() < 0.0)
	{
		int pos = 0;
		float smallestValue = 111111111111111111111.0f;
		for (int i = 0; i < 3; ++i)
		{
			if (S(i, i) < smallestValue)
			{
				pos = i;
				smallestValue = S(i, i);
			}
		}

		V.col(pos) = -V.col(pos);
	}

	Eigen::Matrix3f Fhat;
	Fhat.setZero();
	int numEntriesBelowThreshold = 0;

	int pos = 0;
	for (int i = 0; i < 3; ++i)
	{
		Fhat(i, i) = std::sqrtf(S(i, i));
		if (Fhat(i, i) < 1.0e-4f)
		{
			pos = i;
			++numEntriesBelowThreshold;
		}
	}

	if (numEntriesBelowThreshold == 0)
	{
		//DIRECTLY FROM BENDER ET AL
		Eigen::Vector3f hatFInv(1.0f / Fhat(0, 0), 1.0f / Fhat(1, 1), 1.0f / Fhat(2, 2));
		U = F * V;
		for (unsigned char l = 0; l < 3; l++)
		{
			for (unsigned char m = 0; m < 3; m++)
			{
				U(m, l) *= hatFInv[l];
			}
		}
	}
	else if (numEntriesBelowThreshold > 1)
	{
		U.setIdentity();
	}
	else
	{
		//DIRECTLY FROM BENDER ET AL
		U = F * V;
		for (unsigned char l = 0; l < 3; l++)
		{
			if (l != pos)
			{
				for (unsigned char m = 0; m < 3; m++)
				{
					U(m, l) *= 1.0f / Fhat(l, l);
				}
			}
		}

		Eigen::Vector3f v[2];
		unsigned char index = 0;
		for (unsigned char l = 0; l < 3; l++)
		{
			if (l != pos)
			{
				v[index++] = Eigen::Vector3f(U(0, l), U(1, l), U(2, l));
			}
		}
		Eigen::Vector3f vec = v[0].cross(v[1]);
		vec.normalize();
		U(0, pos) = vec[0];
		U(1, pos) = vec[1];
		U(2, pos) = vec[2];
	}

	if (U.determinant() < 0.0)
	{
		int pos = 2;
		float smallestValue = 11111111111;
		for (int i = 0; i < 3; ++i)
		{
			pos = (S(i, i) < smallestValue) ? i : pos;
			smallestValue = (S(i, i) < smallestValue) ? S(i, i) : smallestValue;
		}

		Fhat(pos, pos) *= -1.0;
		U.col(pos) *= -1.0;
	}

	const float minSVal = 0.577f;

	for (unsigned char i = 0; i < 3; i++)
	{
		if (Fhat(i, i) < minSVal)
		{
			Fhat(i, i) = minSVal;
		}
	}

	//FInverseTranspose = Fhat.inverse().transpose();

	//I1 = (Fhat).trace();
	//float I3 = (Fhat).determinant();

	//logI3 = log(I3);

	////Compute Stress tensor
	//PF = settings.mu * Fhat - settings.mu * FInverseTranspose
	//	+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

	////Compute Strain Energy density field
	//strainEnergy = volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

	// epsilon for hatF
	Eigen::Vector3f epsilonHatF(0.5f*(Fhat(0, 0) * Fhat(0, 0) - 1.0f), 0.5f*(Fhat(1, 1) * Fhat(1, 1) - 1.0f), 0.5f*(Fhat(2, 2) * Fhat(2, 2) - 1.0f));

	const float trace = epsilonHatF[0] + epsilonHatF[1] + epsilonHatF[2];
	const float ltrace = settings.lambda * trace;
	Eigen::Vector3f sigmaVec = epsilonHatF * 2.0f * settings.mu;
	sigmaVec[0] += ltrace;
	sigmaVec[1] += ltrace;
	sigmaVec[2] += ltrace;
	sigmaVec[0] = Fhat(0, 0) * sigmaVec[0];
	sigmaVec[1] = Fhat(1, 1) * sigmaVec[1];
	sigmaVec[2] = Fhat(2, 2) * sigmaVec[2];

	Eigen::Matrix3f sigmaDiag, epsDiag;

	sigmaDiag.row(0) = Eigen::Vector3f(sigmaVec[0], 0.0f, 0.0f);
	sigmaDiag.row(1) = Eigen::Vector3f(0.0f, sigmaVec[1], 0.0f);
	sigmaDiag.row(2) = Eigen::Vector3f(0.0f, 0.0f, sigmaVec[2]);

	epsDiag.row(0) = Eigen::Vector3f(epsilonHatF[0], 0.0f, 0.0f);
	epsDiag.row(1) = Eigen::Vector3f(0.0f, epsilonHatF[1], 0.0f);
	epsDiag.row(2) = Eigen::Vector3f(0.0f, 0.0f, epsilonHatF[2]);

	Eigen::Matrix3f  epsilon = U * epsDiag * V.transpose();
	PF = U * sigmaDiag * V.transpose();

	float psi = 0.0f;
	for (unsigned char j = 0; j < 3; j++)
	for (unsigned char k = 0; k < 3; k++)
		psi += epsilon(j, k) * epsilon(j, k);
	psi = settings.mu*psi + 0.5f*settings.lambda * trace*trace;
	strainEnergy = volume*psi;

	return true;
}

void
PBDSolver::projectConstraintsDistance(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, int numIterations, float k)
{
	float w1;
	float w2;
	float restDistance;

	Eigen::Vector3f x1;
	Eigen::Vector3f x2;

	Eigen::Vector3f deltaX;

	Eigen::Vector3f temp;

	float stiffnessMultiplier;


	if (k == 1.0f)
	{
		stiffnessMultiplier = 1.0f;
	}
	else
	{
		if (k > 1.0f)
		{
			stiffnessMultiplier = 1.0f;
		}
		else
		{
			stiffnessMultiplier = 1.0f - std::pow((1.0f - k), 1.0f / (float)numIterations);
		}
	}


	for (int it = 0; it < numIterations; ++it)
	{
		for (int t = 0; t < tetrahedra.size(); ++t)
		{
			w1 = tetrahedra[t].get_x(0).inverseMass();
			x1 = tetrahedra[t].get_x(0).position();

			w2 = tetrahedra[t].get_x(2).inverseMass();
			x2 = tetrahedra[t].get_x(2).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(0),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(0).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(2).position() += deltaX * w2 * stiffnessMultiplier;

			w2 = tetrahedra[t].get_x(3).inverseMass();
			x2 = tetrahedra[t].get_x(3).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(1),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(0).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(3).position() += deltaX * w2 * stiffnessMultiplier;

			w2 = tetrahedra[t].get_x(1).inverseMass();
			x2 = tetrahedra[t].get_x(1).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(2),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(0).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(1).position() += deltaX * w2 * stiffnessMultiplier;

			//---------------------------------------------------

			w1 = tetrahedra[t].get_x(2).inverseMass();
			x1 = tetrahedra[t].get_x(2).position();

			w2 = tetrahedra[t].get_x(3).inverseMass();
			x2 = tetrahedra[t].get_x(3).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(3),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(2).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(3).position() += deltaX * w2 * stiffnessMultiplier;


			w2 = tetrahedra[t].get_x(1).inverseMass();
			x2 = tetrahedra[t].get_x(1).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(4),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(2).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(1).position() += deltaX * w2 * stiffnessMultiplier;

			//---------------------------------------------------

			w1 = tetrahedra[t].get_x(3).inverseMass();
			x1 = tetrahedra[t].get_x(3).position();

			w2 = tetrahedra[t].get_x(1).inverseMass();
			x2 = tetrahedra[t].get_x(1).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(5),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(3).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(1).position() += deltaX * w2 * stiffnessMultiplier;
		}
	}

}


bool solveVolumeConstraint(
	const Eigen::Vector3f &p0, float invMass0,
	const Eigen::Vector3f &p1, float invMass1,
	const Eigen::Vector3f &p2, float invMass2,
	const Eigen::Vector3f &p3, float invMass3,
	const float restVolume,
	const float negVolumeStiffness,
	const float posVolumeStiffness,
	Eigen::Vector3f &corr0, Eigen::Vector3f &corr1, Eigen::Vector3f &corr2, Eigen::Vector3f &corr3);

void
PBDSolver::projectConstraintsVolume(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, int numIterations, float k)
{
	float w1;
	float w2;
	float restDistance;

	Eigen::Vector3f x1;
	Eigen::Vector3f x2;

	Eigen::Vector3f deltaX;

	Eigen::Vector3f temp;
	float stiffnessMultiplier;

	if (k == 1.0f)
	{
		stiffnessMultiplier = 1.0f;
	}
	else
	{
		if (k > 1.0f)
		{
			stiffnessMultiplier = 1.0f;
		}
		else
		{
			stiffnessMultiplier = 1.0f - std::pow((1.0f - k), 1.0f / (float)numIterations);
		}
	}

	std::vector<int> idxs = { 0, 1, 2, 3 };

	for (int it = 0; it < numIterations; ++it)
	{
		for (int t = 0; t < tetrahedra.size(); ++t)
		{

			w1 = tetrahedra[t].get_x(0).inverseMass();
			x1 = tetrahedra[t].get_x(0).position();

			w2 = tetrahedra[t].get_x(2).inverseMass();
			x2 = tetrahedra[t].get_x(2).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(0),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(0).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(2).position() += deltaX * w2 * stiffnessMultiplier;

			w2 = tetrahedra[t].get_x(3).inverseMass();
			x2 = tetrahedra[t].get_x(3).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(1),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(0).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(3).position() += deltaX * w2 * stiffnessMultiplier;

			w2 = tetrahedra[t].get_x(1).inverseMass();
			x2 = tetrahedra[t].get_x(1).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(2),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(0).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(1).position() += deltaX * w2 * stiffnessMultiplier;

			//---------------------------------------------------

			w1 = tetrahedra[t].get_x(2).inverseMass();
			x1 = tetrahedra[t].get_x(2).position();

			w2 = tetrahedra[t].get_x(3).inverseMass();
			x2 = tetrahedra[t].get_x(3).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(3),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(2).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(3).position() += deltaX * w2 * stiffnessMultiplier;


			w2 = tetrahedra[t].get_x(1).inverseMass();
			x2 = tetrahedra[t].get_x(1).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(4),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(2).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(1).position() += deltaX * w2 * stiffnessMultiplier;

			//---------------------------------------------------

			w1 = tetrahedra[t].get_x(3).inverseMass();
			x1 = tetrahedra[t].get_x(3).position();

			w2 = tetrahedra[t].get_x(1).inverseMass();
			x2 = tetrahedra[t].get_x(1).position();
			computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(5),
				x1, x2, temp, deltaX);
			tetrahedra[t].get_x(3).position() += -deltaX * w1 * stiffnessMultiplier;
			tetrahedra[t].get_x(1).position() += deltaX * w2 * stiffnessMultiplier;
		}

		for (int t = 0; t < tetrahedra.size(); ++t)
		{
			//if (std::abs(tetrahedra[t].getVolume() - tetrahedra[t].getUndeformedVolume()) < 1e-12f)
			//{
			//	continue;
			//}

			Eigen::Vector3f	 corr0;
			Eigen::Vector3f	 corr1;
			Eigen::Vector3f	 corr2;
			Eigen::Vector3f	 corr3;

			Eigen::Vector3f p0 = tetrahedra[t].get_x(0).position();
			Eigen::Vector3f p1 = tetrahedra[t].get_x(1).position();
			Eigen::Vector3f p2 = tetrahedra[t].get_x(2).position();
			Eigen::Vector3f p3 = tetrahedra[t].get_x(3).position();

			float m1 = tetrahedra[t].get_x(0).inverseMass();
			float m2 = tetrahedra[t].get_x(1).inverseMass();
			float m3 = tetrahedra[t].get_x(2).inverseMass();
			float m4 = tetrahedra[t].get_x(3).inverseMass();

			//m1 = 1.0f;
			//m2 = 1.0f;
			//m3 = 1.0f;
			//m4 = 1.0f;

			if (solveVolumeConstraint(p0, m1, p1, m2, p2, m3, p3, m4, tetrahedra[t].getUndeformedVolumeAlternative(),
				1.0f,
				1.0f,
				corr0, corr1, corr2, corr3))
			{
				if (m1 != 0.0f)
				{
					tetrahedra[t].get_x(0).position() += corr0;
					//std::cout << corr0 << std::endl;
				}
				if (m2 != 0.0f)
				{
					tetrahedra[t].get_x(1).position() += corr1;
				}
				if (m3 != 0.0f)
				{
					tetrahedra[t].get_x(2).position() += corr2;
				}
				if (m4 != 0.0f)
				{
					tetrahedra[t].get_x(3).position() += corr3;
				}
			}
		}
	}
}


bool solveVolumeConstraint(
	const Eigen::Vector3f &p0, float invMass0,
	const Eigen::Vector3f &p1, float invMass1,
	const Eigen::Vector3f &p2, float invMass2,
	const Eigen::Vector3f &p3, float invMass3,
	const float restVolume,
	const float negVolumeStiffness,
	const float posVolumeStiffness,
	Eigen::Vector3f &corr0, Eigen::Vector3f &corr1, Eigen::Vector3f &corr2, Eigen::Vector3f &corr3)
{
	float volume = 1.0f / 6.0f * (p1 - p0).cross(p2 - p0).dot(p3 - p0);

	//if (volume == restVolume)
	//{
	//	return false;
	//}


	//std::cout << "Points: " << std::endl;
	//std::cout << p0 << std::endl;
	//std::cout << p1 << std::endl;
	//std::cout << p2 << std::endl;
	//std::cout << p3 << std::endl;
	//std::cout << "Inverse Masses: " << invMass0 << ", " << invMass1 << ", " << invMass2 << ", " << invMass3 << std::endl;

	Eigen::Vector3f d1 = p1 - p0;
	Eigen::Vector3f d2 = p2 - p0;
	Eigen::Vector3f d3 = p3 - p0;

	corr0.setZero(); corr1.setZero(); corr2.setZero(); corr3.setZero();

	if (posVolumeStiffness == 0.0f && volume > 0.0f)
		return false;

	if (negVolumeStiffness == 0.0f && volume < 0.0f)
		return false;


	Eigen::Vector3f grad0 = (p1 - p2).cross(p3 - p2);
	Eigen::Vector3f grad1 = (p2 - p0).cross(p3 - p0);
	Eigen::Vector3f grad2 = (p0 - p1).cross(p3 - p1);
	Eigen::Vector3f grad3 = (p1 - p0).cross(p2 - p0);

	float lambda =
		invMass0 * grad0.squaredNorm() +
		invMass1 * grad1.squaredNorm() +
		invMass2 * grad2.squaredNorm() +
		invMass3 * grad3.squaredNorm();

	if (fabs(lambda) < 1.0e-9)
		return false;

	//if (volume < 0.0f)
	//	lambda = negVolumeStiffness * (volume - restVolume) / lambda;
	//else
	//	lambda = posVolumeStiffness * (volume - restVolume) / lambda;

	lambda = 1.0f * (volume - restVolume) / lambda;

	//std::cout << "Volume: " << volume << "; Undeformed Volume: " << restVolume << std::endl;

	//corr0 = -lambda * invMass0 * grad0;
	//corr1 = -lambda * invMass1 * grad1;
	//corr2 = -lambda * invMass2 * grad2;
	//corr3 = -lambda * invMass3 * grad3;

	//std::cout << "Correction Terms: " << std::endl;
	//std::cout << corr0 << std::endl;
	//std::cout << corr1 << std::endl;
	//std::cout << corr2 << std::endl;
	//std::cout << corr3 << std::endl;
	//std::cout << "___________________________________________" << std::endl;

	return true;
}


void
PBDSolver::projectConstraintsVISCOELASTIC(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, PBDSolverSettings& settings,
std::vector<PBDProbabilisticConstraint>& probabilisticConstraints,
std::vector<CollisionMesh>& collisionGeometry)
{
	Eigen::Vector3f a(0.0f, 1.0f, 0.0f);

	float w1;
	float w2;
	float restDistance;

	Eigen::Vector3f x1;
	Eigen::Vector3f x2;

	Eigen::Vector3f temp;

	Eigen::Matrix3f TEMP_1;

	Eigen::Matrix3f F_orig;
	Eigen::Matrix3f F;
	Eigen::Matrix3f FInverseTranspose;
	Eigen::Matrix3f FTransposeF;

	Eigen::Matrix3f PF;
	Eigen::Matrix3f PF_vol;
	Eigen::Matrix3f gradientTemp;
	Eigen::MatrixXf gradient; gradient.resize(3, 4);

	Eigen::Matrix3f U;
	Eigen::Matrix3f U_orig;
	Eigen::Matrix3f V;
	bool isInverted;

	Eigen::Matrix3f S;
	Eigen::Matrix3f Fhat;

	Eigen::Vector3f deltaX;

	//Mooney-Rivlin
	Eigen::Matrix3f C;


	std::ofstream strainEnergyfile;

	float averageDeltaXLengthAccumulator = 0.0f;
	float averageDeltaXLengthAccumulatorCounter = 0.0f;

	if (settings.printStrainEnergyToFile)
	{
		std::stringstream ss;
		ss << "C:/Users/Stephan/Documents/MATLAB/dissertation/pbd/strainEnergyDebug/strainEnergy_" << m_currentFrame << ".txt";
		strainEnergyfile.open(ss.str());
		ss.clear();
	}

	if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
	{
		calculateTotalStrainEnergy(tetrahedra, particles, settings, -1, strainEnergyfile);
	}
	bool inversionHandled = false;

	isInverted = false;

	for (int it = 0; it < settings.numConstraintIts; ++it)
	{
		for (int t = 0; t < tetrahedra.size(); ++t)
		{
			float lagrangeM;
			float strainEnergy;
			float Volume;

			//Get deformation gradient
			F_orig = tetrahedra[t].getDeformationGradient();

			FTransposeF = F_orig.transpose() * F_orig;

			//Deal with inverted elements / degenerate elements
			//if (std::fabs(F.determinant()) < 1.0e-4f)
			if (true)
			{
				Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);
				S = eigenSolver.pseudoEigenvalueMatrix(); //squared eigenvalues of F
				V = eigenSolver.pseudoEigenvectors(); //eigenvectors
				//eigenDecompositionCardano(FTransposeF, S, V);

				for (int i = 0; i < 3; ++i)
				{
					if (S(i, i) < 0.0f)
					{
						S(i, i) = 0.0f;
					}
				}

				//remove reflection from V if necessary
				if (V.determinant() < 0.0f)
				{
					float minElementValue = 10000000000000;
					int minElementIdx = 111;
					for (int i = 0; i < 3; ++i)
					{
						if (S(i, i) < minElementValue)
						{
							minElementValue = V(i, i);
							minElementIdx = i;
						}
					}
					V.col(minElementIdx) *= -1.0f;
					//V.col(0) *= -1.0f;
				}

				//determine entries of F
				F.setZero();
				for (int i = 0; i < 3; ++i)
				{
					F(i, i) = std::sqrtf(S(i, i));
				}

				U = F_orig * V * F.inverse();
				U_orig = U;
				int numEntriesNearZero = 0;
				int idx = 0;
				for (int i = 0; i < 3; ++i)
				{
					if (std::fabs(F(i, i)) < 1.0e-4f)
					{
						++numEntriesNearZero;
						idx = i;
					}
				}

				if (numEntriesNearZero > 0)
				{
					if (numEntriesNearZero == 1)
					{
						if (idx == 0)
						{
							U.col(0) = U.col(1).cross(U.col(2)).normalized();
						}
						else if (idx == 1)
						{
							U.col(1) = U.col(0).cross(U.col(2)).normalized();
						}
						else
						{
							U.col(2) = U.col(0).cross(U.col(1)).normalized();
						}
					}
					else
					{
						U.setIdentity();
					}
				}

				//remove reflection from U if necessary
				if (U.determinant() < 0.0f)
				{
					float minElementValue = 1000000000000000;
					int minElementIdx = 111;
					for (int i = 0; i < 3; ++i)
					{
						if (F(i, i) < minElementValue)
						{
							minElementValue = F(i, i);
							minElementIdx = i;
						}
					}

					F(minElementIdx, minElementIdx) *= -1.0f;
					U.col(minElementIdx) *= -1.0f;
				}

				const float minXVal = 0.577f;

				for (unsigned char j = 0; j < 3; j++)
				{
					if (F(j, j) < minXVal)
						F(j, j) = minXVal;
				}

				const float maxXVal = 10000.0f;

				for (unsigned char j = 0; j < 3; j++)
				{
					if (F(j, j) > maxXVal)
						F(j, j) = maxXVal;
				}
			}

			FInverseTranspose = F.inverse().transpose();
			FTransposeF = F.transpose() * F;

			switch (settings.materialModel)
			{
			case PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN:
			{
				Volume = tetrahedra[t].getUndeformedVolume();

				//Compute Isotropic Invariants
				float I1 = (FTransposeF).trace();
				float I3 = (FTransposeF).determinant();

				float logI3 = log(I3);

				/*PF = settings.mu * F - settings.mu * FInverseTranspose;
				PF_vol =  ((settings.lambda * logI3) / 2.0) * FInverseTranspose;*/
				PF = settings.mu * F - settings.mu * FInverseTranspose
				 + ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

				//Compute Strain Energy density field
				strainEnergy = Volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

			}
			break;
			case PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER:
			{
				Volume = tetrahedra[t].getUndeformedVolume();

				//Compute Isotropic Invariants
				float I1 = (FTransposeF).trace();
				float I3 = (FTransposeF).determinant();

				float logI3 = log(I3);

				Eigen::Vector3f rotated_a = V.transpose() * settings.MR_a;

				//PF = settings.mu * F - settings.mu * FInverseTranspose;
				//PF_vol = ((settings.lambda * logI3) / 2.0) * FInverseTranspose;
				PF = settings.mu * F - settings.mu * FInverseTranspose
				+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;


				//Compute Strain Energy density field
				strainEnergy = Volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

				//STRETCH ('pseudo-invariant' of C)
				float lambda = std::sqrtf((rotated_a.transpose() * (1.0f * FTransposeF)).dot(rotated_a));

				strainEnergy += (settings.anisotropyParameter / 2.0f) * std::pow(lambda - 1.0f, 2.0f);

				PF += F * std::pow(F.determinant(), -2.0 / 3.0)
					* (settings.anisotropyParameter * (lambda - 1.0f)
					* (kroneckerProduct(rotated_a, rotated_a) + (lambda / 3.0f) * FTransposeF.inverse()));
			}
			break;
			case PBDSolverSettings::CONSTITUTIVE_MODEL::RUBIN_BODNER:
			{
				std::cout << "Rubin-Bodner Implementation Removed Temporarily!" << std::endl;
			}
			break;
			default:
				std::cout << "ERROR: No Constitutive Model Selected!" << std::endl;
				break;
			}
			
			//VISCOELASTICITY -----------------------------------------------------------------------------------------------------------
			////if (settings.alpha != 0.0f && settings.rho != 0.0f)
			//{
			//	//FInverseTranspose = F.inverse();

			//	PF *= F.inverse();
			//	//PF_vol *= FInverseTranspose;

			//	/*FInverseTranspose = F_orig.inverse().transpose();

			//	PF *= FInverseTranspose.transpose();
			//	*/
			//	Eigen::Matrix3f vMult;
			//	
			//	if (settings.useFullPronySeries)
			//	{
			//		Eigen::Matrix3f temp;
			//		temp.setZero();

			//		for (int pComponent = 0; pComponent < settings.fullAlpha.size(); ++pComponent)
			//		{
			//			temp += (2.0f * settings.deltaT * settings.fullAlpha[pComponent] * PF
			//				+ settings.fullRho[pComponent] * tetrahedra[t].getFullUpsilon(pComponent)) / (settings.deltaT + settings.fullRho[pComponent]);

			//			tetrahedra[t].getFullUpsilon(pComponent) = (2.0f * settings.deltaT * settings.fullAlpha[pComponent] * PF
			//				+ settings.fullRho[pComponent] * tetrahedra[t].getFullUpsilon(pComponent)) / (settings.deltaT + settings.fullRho[pComponent]);
			//		}

			//		vMult = temp;
			//	}
			//	else
			//	{
			//		vMult = tetrahedra[t].getUpsilon() * V.transpose().inverse();

			//		vMult = (2.0f * settings.deltaT * settings.alpha * PF + settings.rho * vMult) / (settings.deltaT + settings.rho);

			//		tetrahedra[t].getUpsilon() = U * vMult * V.transpose();
			//	}

			//	PF = (PF) * 2.0f - vMult * 1.0f;

			//	/*PF = F_orig * PF;*/

			//	if (settings.trackS)
			//	{
			//		settings.tracker.S.push_back(PF);
			//	}

			//	PF = F * PF;
			//}

			PF = U * PF * V.transpose();
			/*PF = U * PF * V.transpose();*/

			//PF GRADIENT ---------------------------------------------------------------------------------------------------------------

			gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
			gradient.col(0) = gradientTemp.col(0);
			gradient.col(1) = gradientTemp.col(1);
			gradient.col(2) = gradientTemp.col(2);
			gradient.col(3) = -gradientTemp.rowwise().sum();


			//PBD MAIN ROUTINE-----------------------------------------------------------------------------------------------------------

			float denominator = 0.0;

			for (int cI = 0; cI < 4; ++cI)
			{
				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					denominator += tetrahedra[t].get_x(cI).inverseMass()
						* gradient.col(cI).lpNorm<2>();
				}
			}

			//prevent division by zero if there is no deformation
			if (denominator < 1e-15)
			{
				continue;
			}

			lagrangeM = -(strainEnergy / denominator);


			if (std::isnan(lagrangeM) || std::isinf(lagrangeM))
			{
				//std::cout << "U orig: " << std::endl;
				//std::cout << U_orig << std::endl;
				//std::cout << "U: " << std::endl;
				//std::cout << U << std::endl;
				//std::cout << "V: " << std::endl;
				//std::cout << V << std::endl;
				//std::cout << "F: " << std::endl;
				//std::cout << F << std::endl;
				//std::cout << "S: " << std::endl;
				//std::cout << S << std::endl;
				//std::cout << "det: " << F.determinant() << std::endl;
				//std::cout << "vol: " << Volume << std::endl;
				//std::cout << "PF: " << std::endl;
				//std::cout << PF << std::endl;
				//std::cout << "Gradient: " << std::endl;
				//std::cout << gradient << std::endl;
				//std::cout << "LaM: " << lagrangeM << std::endl;
				//std::cout << "Strain Energy: " << strainEnergy << std::endl;
				//std::cout << "Log(I3): " << log(FTransposeF.determinant()) << std::endl;

				//std::cout << "---------------------------------" << std::endl;
				continue;
			}
			else
			{
				for (int cI = 0; cI < 4; ++cI)
				{
					if (tetrahedra[t].get_x(cI).inverseMass() != 0)
					{
						if (!settings.disablePositionCorrection)
						{
							deltaX = (tetrahedra[t].get_x(cI).inverseMass()
								* lagrangeM) * gradient.col(cI);

							tetrahedra[t].get_x(cI).position() += deltaX;

							if (settings.trackAverageDeltaXLength)
							{
								float value = deltaX.squaredNorm();
								if (!std::isnan(value) && !std::isinf(value))
								{
									averageDeltaXLengthAccumulator += value;
									averageDeltaXLengthAccumulatorCounter += 1.0f;
								}
							}
						}
					}
				}
			}
		}

		if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
		{
			calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
		}

		if (settings.enableGroundPlaneCollision)
		{
			for (int p = 0; p < particles->size(); ++p)
			{
				if ((*particles)[p].position()[1] < 0.0f)
				{
					(*particles)[p].position()[1] = 0.0f;
				}
			}
		}

		//COLLISION HANDLING
		for (int c = 0; c < collisionGeometry.size(); ++c)
		{
			collisionGeometry[c].resolveParticleCollisions(*particles);
		}

	}
	//------------------------------------

	//if (settings.trackS)
	//{
	//	settings.tracker.S.push_back(PF * F_orig.inverse());
	//}
	if (settings.trackF)
	{
		settings.tracker.F.push_back(F);
	}
	if (settings.trackPF)
	{
		settings.tracker.PF.push_back(PF);
	}

	if (settings.trackSpecificPosition)
	{
		settings.tracker.specificPosition.push_back((*particles)[settings.trackSpecificPositionIdx].position());
	}

	if (settings.trackAverageDeltaXLength)
	{
		if (averageDeltaXLengthAccumulatorCounter != 0.0f)
		{
			settings.tracker.averageDeltaXLength.push_back(averageDeltaXLengthAccumulator / averageDeltaXLengthAccumulatorCounter);
		}
		else
		{
			settings.tracker.averageDeltaXLength.push_back(0.0f);
		}
		
		
	}

	//for (int pC = 0; pC < probabilisticConstraints.size(); ++pC)
	//{
	//	probabilisticConstraints[pC].project(*particles);
	//}

	if (settings.printStrainEnergyToFile)
	{
		strainEnergyfile.close();
	}

	if (inversionHandled)
	{
		std::cout << "Inversion handled successfully!" << std::endl;
	}
}
