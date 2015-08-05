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
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences)
{
	//Advance Velocities
	advanceVelocities(tetrahedra, particles, settings);

	//Advance Positions
	advancePositions(tetrahedra, particles, settings);

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
	//projectConstraintsNeoHookeanMixture(tetrahedra, particles, settings);
	//projectConstraintsMooneyRivlin(tetrahedra, particles, settings);
	//projectConstraintsDistance(tetrahedra, particles, settings.numConstraintIts, settings.youngsModulus);
	//projectConstraintsGeometricInversionHandling(tetrahedra, particles, settings);
	projectConstraintsVolume(tetrahedra, particles, settings.numConstraintIts, settings.youngsModulus);

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
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	for (auto& p : *particles)
	{
		float temp = settings.deltaT * p.inverseMass() * settings.gravity;
		p.velocity().x() = p.previousVelocity().x() + 0;
		p.velocity().y() = p.previousVelocity().y() + temp;
		p.velocity().z() = p.previousVelocity().z() + 0;
	}
}

void
PBDSolver::advancePositions(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	for (auto& p : *particles)
	{
		p.position() = p.previousPosition() + settings.deltaT * p.velocity() * p.inverseMass();
		//std::cout << p.velocity().transpose() << std::endl;
	}
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

	for (int t = 0; t < settings.numTetrahedra; ++t)
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

		computeGreenStrainAndPiolaStressInversion(F, FTransposeF, U, V, Volume, settings.mu, settings.lambda, PF, epsilon, strainEnergy, 100);

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

void
PBDSolver::projectConstraints(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
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

	Eigen::Vector3f deltaXGeometric;

	std::ofstream strainEnergyfile;

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

	float local_lambda = settings.lambda;
	float local_mu = settings.mu;
	float currentYM = settings.youngsModulus;

	bool increasingYM = false;

	//if (settings.youngsModulus > 2.0f)
	//{
	//	increasingYM = true;
	//	currentYM = 2.0f;
	//}

	//int numIncreaseIts = settings.youngsModulus / 2;

	//for (int increaseIt = 1; increaseIt <= numIncreaseIts; ++increaseIt)
	{
		//currentYM = (settings.youngsModulus / numIncreaseIts) * increaseIt;

		for (int it = 0; it < settings.numConstraintIts; ++it)
		{
			for (int t = 0; t < settings.numTetrahedra; ++t)
			//for (int t = settings.numTetrahedra - 1; t >= 0; --t)
			{
				for (int innerT = 0; innerT < settings.numTetrahedraIterations; ++innerT)
				{

					local_lambda = settings.calculateLambda(currentYM, settings.poissonRatio);
					local_mu = settings.calculateMu(currentYM, settings.poissonRatio);

					float lagrangeM;
					float strainEnergy;
					bool repeatIteration = false;

					//Compute volume
					float Volume = std::abs(tetrahedra[t].getUndeformedVolume());

					//Get deformation gradient
					F = tetrahedra[t].getDeformationGradient();

					FInverseTranspose = F.inverse().transpose();
					FTransposeF = F.transpose() * F;

					//Compute Isotropic Invariants
					float I1 = (FTransposeF).trace();
					float I3 = (FTransposeF).determinant();

					float logI3 = log(I3);

					if (F.isIdentity())
					{
						continue;
					}

					//if (!correctInversion(F, FTransposeF, FInverseTranspose, PF, U, V, I1, 0.0, logI3, strainEnergy, Volume, settings))
					//{
					//	continue;
					//}
					//correctInversion(F, FTransposeF, FInverseTranspose, PF, U, V, I1, 0.0, logI3, strainEnergy, Volume, settings);
					//if (F.determinant() < 1.0e-04f || tetrahedra[t].getVolume() < 0.00001)
					//if (F.determinant() < 0.0f)
					{
						computeGreenStrainAndPiolaStressInversion(F, FTransposeF, U, V, Volume, local_mu, local_lambda, PF, epsilon, strainEnergy, it);
					}
					//else
					{
						//computeGreenStrainAndPiolaStress(F, Volume, settings.mu, settings.lambda, PF, epsilon, strainEnergy);
					}

					//std::cout << "Current Volume: " << tetrahedra[t].getVolume();
					//if (tetrahedra[t].getVolume() < 0.00001)
					//{
					//	std::cout << "Degenerate/Inverted Tetrahedron at " << t << "; V =  " << Volume << std::endl;
					//}

					//if (strainEnergy > 0.01)
					//{
					//	std::cout << strainEnergy << std::endl;
					//}

					if (settings.correctStrongForcesWithSubteps)
					{
						if (std::abs(strainEnergy) > 1.0e-4f)
						{
							strainEnergy /= strainEnergy / 1.0e-4f;
							repeatIteration = true;
						}
					}

					gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
					gradient.col(0) = gradientTemp.col(0);
					gradient.col(1) = gradientTemp.col(1);
					gradient.col(2) = gradientTemp.col(2);
					gradient.col(3) = -gradientTemp.rowwise().sum();

					//std::cout << "Strain Energy: " << strainEnergy << std::endl;

					//Compute Lagrange Multiplier
					float denominator = 0.0;

					for (int cI = 0; cI < 4; ++cI)
					{
						if (tetrahedra[t].get_x(cI).inverseMass() != 0)
						{
							denominator += tetrahedra[t].get_x(cI).inverseMass()
								* gradient.col(cI).squaredNorm();
						}
					}


					//if (denominator < 1.0e-9f)
					//	continue;

					//if (std::fabs(denominator) < 1.0e-6f)
					//{
					//	continue;
					//}
					//else
					{
						lagrangeM = -(strainEnergy / denominator);
					}

					for (int cI = 0; cI < 4; ++cI)
					{
						//float w1;
						//float w2;
						//Eigen::Vector3f x1;
						//Eigen::Vector3f x2;

						//if (settings.useGeometricConstraintLimits)
						//{
						//	switch (cI)
						//	{
						//	case 0:
						//		break;
						//	case 1:
						//		break;
						//	case 2:
						//		break;
						//	case 3:
						//		break;
						//	default:

						//	}
						//}

						if (tetrahedra[t].get_x(cI).inverseMass() != 0)
						{
							deltaX = (tetrahedra[t].get_x(cI).inverseMass()
								* lagrangeM) * gradient.col(cI);

							tetrahedra[t].get_x(cI).position() += deltaX;
						}
					}

					if (repeatIteration)
					{
						//++t;
						--t;
					}
				}
			}

			if (increasingYM)
			{
				if (currentYM >= settings.youngsModulus)
				{
					currentYM = settings.youngsModulus;
					increasingYM = false;
				}
				else
				{
					currentYM += 0.15;
				}
			}

			if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
			{
				calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
			}
		}
	}

	if (settings.printStrainEnergyToFile)
	{
		strainEnergyfile.close();
	}
}

void
PBDSolver::projectConstraintsOLD(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	Eigen::Matrix3f F;
	Eigen::Matrix3f FInverseTranspose;
	Eigen::Matrix3f FTransposeF;

	Eigen::Matrix3f PF;
	Eigen::Matrix3f gradientTemp;
	Eigen::MatrixXf gradient; gradient.resize(3, 4);

	Eigen::Matrix3f U;
	Eigen::Matrix3f V;
	bool isInverted;

	Eigen::Matrix3f S;
	Eigen::Matrix3f Fhat;

	Eigen::Vector3f deltaX;

	std::ofstream strainEnergyfile;

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

	for (int it = 0; it < settings.numConstraintIts; ++it)
	{
		std::cout << "IT: " << it << std::endl;

		for (int t = 0; t < settings.numTetrahedra; ++t)
		{
			float lagrangeM;

			//Get deformation gradient
			F = tetrahedra[t].getDeformationGradient();

			//std::cout << t << ": " << std::endl;
			//std::cout << tetrahedra[t].getDeformedShapeMatrix() << std::endl;
			//std::cout << F << std::endl;
			//std::cout << tetrahedra[t].getReferenceShapeMatrixInverseTranspose().transpose() << std::endl;
			//if (t == 47)
			//{
			//	//std::cout << tetrahedra[t].getVertexIndices()[0] << ", ";
			//	//std::cout << tetrahedra[t].getVertexIndices()[1] << ", ";
			//	//std::cout << tetrahedra[t].getVertexIndices()[2] << ", ";
			//	//std::cout << tetrahedra[t].getVertexIndices()[3] << std::endl;
			//	//std::cout << (*particles)[tetrahedra[t].getVertexIndices()[0]].position() << std::endl;
			//	//std::cout << (*particles)[tetrahedra[t].getVertexIndices()[1]].position() << std::endl;
			//	//std::cout << (*particles)[tetrahedra[t].getVertexIndices()[2]].position() << std::endl;
			//	//std::cout << (*particles)[tetrahedra[t].getVertexIndices()[3]].position() << std::endl;
			//	printf("F [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", t,
			//		F(0, 0), F(0, 1), F(0, 2),
			//		F(1, 0), F(1, 1), F(1, 2),
			//		F(2, 0), F(2, 1), F(2, 2));
			//	FTransposeF = tetrahedra[t].getReferenceShapeMatrixInverseTranspose().transpose();
			//	printf("RefInverse [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", t,
			//		FTransposeF(0, 0), FTransposeF(0, 1), FTransposeF(0, 2),
			//		FTransposeF(1, 0), FTransposeF(1, 1), FTransposeF(1, 2),
			//		FTransposeF(2, 0), FTransposeF(2, 1), FTransposeF(2, 2));
			//	FTransposeF = tetrahedra[t].getDeformedShapeMatrix();
			//	printf("Deformed [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", t,
			//		FTransposeF(0, 0), FTransposeF(0, 1), FTransposeF(0, 2),
			//		FTransposeF(1, 0), FTransposeF(1, 1), FTransposeF(1, 2),
			//		FTransposeF(2, 0), FTransposeF(2, 1), FTransposeF(2, 2));
			//}
			if (F.isIdentity())
			{
				continue;
			}

			isInverted = F.determinant() < 0.0;

			//check for inversion
			if (isInverted)
			{
				//inversionHandled = true;
				////1. Compute Eigendecomposition
				//Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);
				//S = eigenSolver.pseudoEigenvalueMatrix();
				//V = eigenSolver.pseudoEigenvectors();
				//for (int i = 0; i < 3; ++i)
				//{
				//	if (S(i, i) < 0.0f)
				//	{
				//		S(i, i) = 0.0f;
				//	}
				//}
				////2.  Detect if V is a reflection .
				////    Make a rotation out of it by multiplying one column with -1.
				//const float detV = V.determinant();
				//if (detV < 0.0)
				//{
				//	float minLambda = FLT_MAX;
				//	unsigned char pos = 0;
				//	for (unsigned char l = 0; l < 3; l++)
				//	{
				//		if (S(l,l) < minLambda)
				//		{
				//			pos = l;
				//			minLambda = S(l, l);
				//		}
				//	}
				//	V(0, pos) = -V(0, pos);
				//	V(1, pos) = -V(1, pos);
				//	V(2, pos) = -V(2, pos);
				//}
				////3. Compute Fhat
				//Fhat.setZero();
				//Fhat(0, 0) = sqrtf(S(0, 0));
				//Fhat(1, 1) = sqrtf(S(1, 1));
				//Fhat(2, 2) = sqrtf(S(2, 2));
				////4. Compute U
				//U = F * V * Fhat.inverse();

				////CORRECT U
				////
				//// Check for values of hatF near zero
				////
				//unsigned char chk = 0;
				//unsigned char pos = 0;
				//for (unsigned char l = 0; l < 3; l++)
				//{
				//	if (fabs(Fhat(l, l)) < 1.0e-4f)
				//	{
				//		pos = l;
				//		chk++;
				//	}
				//}

				//if (chk > 0)
				//{
				//	if (chk > 1)
				//	{
				//		U.setIdentity();
				//	}
				//	else
				//	{
				//		U = F * V;
				//		for (unsigned char l = 0; l < 3; l++)
				//		{
				//			if (l != pos)
				//			{
				//				for (unsigned char m = 0; m < 3; m++)
				//				{
				//					U(m, l) *= 1.0f / Fhat(l, l);
				//				}
				//			}
				//		}

				//		Eigen::Vector3f v[2];
				//		unsigned char index = 0;
				//		for (unsigned char l = 0; l < 3; l++)
				//		{
				//			if (l != pos)
				//			{
				//				v[index++] = Eigen::Vector3f(U(0, l), U(1, l), U(2, l));
				//			}
				//		}
				//		Eigen::Vector3f vec = v[0].cross(v[1]);
				//		vec.normalize();
				//		U(0, pos) = vec[0];
				//		U(1, pos) = vec[1];
				//		U(2, pos) = vec[2];
				//	}
				//}
				//else
				//{
				//	Eigen::Vector3f hatFInv(1.0f / Fhat(0, 0), 1.0f / Fhat(1, 1), 1.0f / Fhat(2, 2));
				//	U = F * V;
				//	for (unsigned char l = 0; l < 3; l++)
				//	{
				//		for (unsigned char m = 0; m < 3; m++)
				//		{
				//			U(m, l) *= hatFInv[l];
				//		}
				//	}
				//}

				////5. Check if U is also a rotation and correct
				//if (U.determinant() < 0)
				//{
				//	//find minimal element of U
				//	int minElementFhat = 0;
				//	float minElementFhatValue = Fhat(0, 0);
				//	for (int e = 0; e < 3; ++e)
				//	{
				//		if (Fhat(e, e) < minElementFhatValue)
				//		{
				//			minElementFhat = e;
				//			minElementFhatValue = Fhat(e, e);
				//		}
				//	}

				//	Fhat(minElementFhat, minElementFhat) *= -1.0;
				//	U.col(minElementFhat) *= -1.0;
				//}

			}

			FInverseTranspose = F.inverse().transpose();
			FTransposeF = F.transpose() * F;

			//Compute Isotropic Invariants
			float I1 = (FTransposeF).trace();
			float I3 = (FTransposeF).determinant();

			float logI3 = log(I3);

			//Compute Stress tensor
			if (isInverted)
			{
				//PF = settings.mu * Fhat - settings.mu * Fhat.inverse().transpose()
				//	+ ((settings.lambda * logI3) / 2.0) * Fhat.inverse().transpose();

				//PF = U * PF * V;
				//I1 = (Fhat.transpose() * Fhat).trace();
				//I3 = (Fhat.transpose() * Fhat).determinant();
				//logI3 = log(I3);
			}
			else
			{
				PF = settings.mu * F - settings.mu * FInverseTranspose
					+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;
			}

			//Compute volume
			float Volume = tetrahedra[t].getUndeformedVolume();


			//std::cout << "Current Volume: " << tetrahedra[t].getVolume();
			if (tetrahedra[t].getVolume() < 0.00001)
			{
				std::cout << "Degenerate/Inverted Tetrahedron at " << t << "; V =  " << Volume << std::endl;
			}

			gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
			gradient.col(0) = gradientTemp.col(0);
			gradient.col(1) = gradientTemp.col(1);
			gradient.col(2) = gradientTemp.col(2);
			gradient.col(3) = -gradientTemp.rowwise().sum();

			//Compute Strain Energy density field
			float strainEnergy = Volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

			//if (strainEnergy < -1e-5f)
			//{
			//	strainEnergy = -1e-5f;
			//}
			//std::cout << "Strain Energy: " << strainEnergy << std::endl;
			//Compute Lagrange Multiplier

			float denominator = 0.0;

			for (int cI = 0; cI < 4; ++cI)
			{
				if (t == 47)
				{
					std::cout << "[" << cI << "]: " << tetrahedra[t].get_x(cI).inverseMass() << "; "
						<< gradient.col(cI).lpNorm<2>() << std::endl;
				}

				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					denominator += tetrahedra[t].get_x(cI).inverseMass()
						* gradient.col(cI).lpNorm<2>();

					//if (std::fabs(denominator) > 1e-06)
					//{
					//	std::cout << "Condition met!" << std::endl;
					//}

				}
			}

			//if (std::fabs(denominator) < 1e-06)
			//{
			//	//std::cout << "Skipping!" << std::endl;
			//	continue;
			//}
			//else
			{
				lagrangeM = -(strainEnergy / denominator);
			}
			//if (t == 47)
			//{
			//	printf("I1 = %8.16f \n",I1);
			//	printf("I3 = %8.16f \n",I3);
			//	printf("lagrangeMultiplier = %8.16f \n", lagrangeM);
			//	printf("strainEnergy = %8.16f \n", strainEnergy);
			//	printf("denominator = %8.16f \n", denominator);
			//	printf("PF [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n", 47,
			//		PF(0, 0), PF(0, 1), PF(0, 2),
			//		PF(1, 0), PF(1, 1), PF(1, 2),
			//		PF(2, 0), PF(2, 1), PF(2, 2));
			//	printf("Gradient [%d]: \n %4.8f, %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f, %4.8f \n", 47,
			//		gradient(0, 0), gradient(0, 1), gradient(0, 2), gradient(0, 3),
			//		gradient(1, 0), gradient(1, 1), gradient(1, 2), gradient(1, 3),
			//		gradient(2, 0), gradient(2, 1), gradient(2, 2), gradient(2, 3));
			//}
			if (std::isnan(lagrangeM))
			{
				//if (isInverted)
				//{
				//	std::cout << "NAN! ";
				//}
				//else
				//{
				//	std::cout << "NON-INV NAN! ";
				//}
				//
				//std::cout << "Deformation Gradient" << std::endl;
				//std::cout << F << std::endl;
				//std::cout << "Inverse of deformation gradient:" << std::endl;
				//std::cout << F.inverse().transpose() << std::endl;
				//std::cout << "Stress Tensor" << std::endl;
				//std::cout << PF << std::endl;
				//std::cout << "Tensor Gradient " << std::endl;
				//std::cout << gradient << std::endl;
				//std::cout << "Strain Energy: " << strainEnergy << std::endl;
				//std::cout << "Denominator: " << denominator << std::endl;
				//std::cout << "Lagrange Multiplier: " << lagrangeM << std::endl;
				////std::cout << "Inverse Mass: " << tetrahedra[t].get_x(c).inverseMass() << std::endl;
				//std::cout << "Undeformed Volume: " << V << std::endl;
				//
				//std::cout << "STEPS: " << std::endl;
				//std::cout << (settings.mu * F) << std::endl;
				//std::cout << settings.mu * F.inverse().transpose() << std::endl;
				//std::cout << log(I3) << std::endl;
				//std::cout << F.inverse().transpose() << std::endl;

				lagrangeM = 0.0;
			}

			for (int cI = 0; cI <4; ++cI)
			{
				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					deltaX = (tetrahedra[t].get_x(cI).inverseMass()
						* lagrangeM) * gradient.col(cI);

					tetrahedra[t].get_x(cI).position() += deltaX;
					if (t == 47)
					{
						std::cout << "[ " << cI << "] : " << std::endl;
						std::cout << deltaX << std::endl;
					}
				}
			}

		}


		if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
		{
			calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
		}

	}

	if (settings.printStrainEnergyToFile)
	{
		strainEnergyfile.close();
	}

	if (inversionHandled)
	{
		std::cout << "Inversion handled successfully!" << std::endl;
	}
}

void
PBDSolver::projectConstraintsSOR(std::vector<PBDTetrahedra3d>& tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
	std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences)
{
		Eigen::Matrix3f F;
		Eigen::Matrix3f FInverseTranspose;
		Eigen::Matrix3f FTransposeF;

		Eigen::Matrix3f PF;
		Eigen::Matrix3f gradientTemp;
		Eigen::MatrixXf gradient; gradient.resize(3, 4);

		Eigen::Matrix3f U;
		Eigen::Matrix3f V;
		bool isInverted;

		Eigen::Vector3f deltaX;

		std::ofstream strainEnergyfile;

		if (settings.printStrainEnergyToFile)
		{
			std::stringstream ss;
			ss << "C:/Users/Stephan/Documents/MATLAB/dissertation/pbd/strainEnergyDebug/strainEnergySOR_" << m_currentFrame << ".txt";
			strainEnergyfile.open(ss.str());
			ss.clear();
		}

		if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
		{
			calculateTotalStrainEnergy(tetrahedra, particles, settings, -1, strainEnergyfile);
		}

		for (int it = 0; it < settings.numConstraintIts; ++it)
		{
			//reset the accmulator arrays
			for (int t = 0; t < particles->size(); ++t)
			{
				temporaryPositions[t].setZero();
				numConstraintInfluences[t] = 0;
			}

			boost::thread_group threads;

			int numThreads = 4;

			int stepSize = settings.numTetrahedra / numThreads;

			mutexStruct mutexInstance;
			//std::cout << settings.numTetrahedra << std::endl;
			for (int t = 0; t < settings.numTetrahedra; t += stepSize)
			{
				threads.create_thread(boost::bind(projectConstraintsSOR_CORE, boost::ref(mutexInstance), boost::ref(tetrahedra), boost::ref(particles), boost::ref(settings),
					boost::ref(temporaryPositions), boost::ref(numConstraintInfluences),
					0, settings.numTetrahedra - 1));
				//if (t < settings.numTetrahedra - stepSize - 1)
				//{
				//	threads.create_thread(boost::bind(projectConstraintsSOR_CORE, boost::ref(mutexInstance), boost::ref(tetrahedra), boost::ref(particles), boost::ref(settings),
				//		boost::ref(temporaryPositions), boost::ref(numConstraintInfluences),
				//		t, t + stepSize - 1));
				//	//std::cout << t << "; " << t + stepSize - 1 << std::endl;
				//}
				//else
				//{
				//	threads.create_thread(boost::bind(projectConstraintsSOR_CORE, boost::ref(mutexInstance), boost::ref(tetrahedra), boost::ref(particles), boost::ref(settings),
				//		boost::ref(temporaryPositions), boost::ref(numConstraintInfluences),
				//		t, settings.numTetrahedra - 1));
				//	//std::cout << t << "; " << settings.numTetrahedra - 1 << std::endl;
				//	continue;
				//}
			};
			
			threads.join_all();

			for (int t = 0; t < settings.numTetrahedra; ++t)
			{
				for (int cI = 0; cI < 4; ++cI)
				{

					if (tetrahedra[t].get_x(cI).inverseMass() != 0)
					{
						if (numConstraintInfluences[tetrahedra[t].getVertexIndices()[cI]] != 0)
						{
							//tetrahedra[t].get_x(cI).position() += (temporaryPositions[tetrahedra[t].getVertexIndices()[cI]] * settings.w) / (numConstraintInfluences[tetrahedra[t].getVertexIndices()[cI]] / 4.0f); 
							tetrahedra[t].get_x(cI).position() += (temporaryPositions[tetrahedra[t].getVertexIndices()[cI]]);
							//std::cout << temporaryPositions[tetrahedra[t].getVertexIndices()[cI]] << std::endl;
						}
					}
				}
			}

			if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
			{
				calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
			}
		}

		if (settings.printStrainEnergyToFile)
		{
			strainEnergyfile.close();
		}
	}


void
projectConstraintsSOR_CORE(mutexStruct& sorMutex, std::vector<PBDTetrahedra3d>& tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings,
	std::vector<Eigen::Vector3f>& temporaryPositions, std::vector<int>& numConstraintInfluences,
	int start, int end)
{
	Eigen::Matrix3f F;
	Eigen::Matrix3f FInverseTranspose;
	Eigen::Matrix3f FTransposeF;

	Eigen::Matrix3f PF;
	Eigen::Matrix3f gradientTemp;
	Eigen::MatrixXf gradient; gradient.resize(3, 4);

	Eigen::Matrix3f U;
	Eigen::Matrix3f V;
	Eigen::Matrix3f epsilon;
	bool isInverted;

	Eigen::Vector3f deltaX;

	for (size_t t = start; t != end; ++t)
	{
		float lagrangeM;

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

		//Compute volume
		float Volume = tetrahedra[t].getUndeformedVolume();

		float strainEnergy = 0;


		if (F.isIdentity())
		{
			continue;
		}


		//Compute Stress tensor
		//PF = settings.mu * F - settings.mu * FInverseTranspose
		//	+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;
		computeGreenStrainAndPiolaStressInversion(F, FTransposeF, U, V, Volume, settings.mu, settings.lambda, PF, epsilon, strainEnergy, -1);

		


		//std::cout << "Current Volume: " << tetrahedra[t].getVolume();
		if (tetrahedra[t].getVolume() < 0.00001)
		{
			std::cout << "Degenerate/Inverted Tetrahedron at " << t << "; V =  " << Volume << std::endl;
		}

		gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
		gradient.col(0) = gradientTemp.col(0);
		gradient.col(1) = gradientTemp.col(1);
		gradient.col(2) = gradientTemp.col(2);
		gradient.col(3) = -gradientTemp.rowwise().sum();

		//Compute Strain Energy density field
		//float strainEnergy = Volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

		//Compute Lagrange Multiplier
		float denominator = 0.0;

		for (int cI = 0; cI < 4; ++cI)
		{
			if (tetrahedra[t].get_x(cI).inverseMass() != 0)
			{
				denominator += tetrahedra[t].get_x(cI).inverseMass()
					* gradient.col(cI).lpNorm<2>();
			}
		}

		{
			lagrangeM = -(strainEnergy / denominator);
		}

		if (std::isnan(lagrangeM))
		{
			//std::cout << "NAN!" << std::endl;
			lagrangeM = 0.0;
		}

		
		boost::mutex::scoped_lock lock(sorMutex.m_mutexSOR);
		for (int cI = 0; cI < 4; ++cI)
		{
			if (tetrahedra[t].get_x(cI).inverseMass() != 0)
			{
				deltaX = (tetrahedra[t].get_x(cI).inverseMass()
					* lagrangeM) * gradient.col(cI);

				temporaryPositions[tetrahedra[t].getVertexIndices()[cI]] += deltaX;
				numConstraintInfluences[tetrahedra[t].getVertexIndices()[cI]] += 1;
			}
		}
	}

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
computeGreenStrainAndPiolaStressInversion(const Eigen::Matrix3f& F, const Eigen::Matrix3f& FTransposeF,
Eigen::Matrix3f& U, Eigen::Matrix3f& V,
const float restVolume,
const float mu, const float lambda, Eigen::Matrix3f &epsilon, Eigen::Matrix3f &sigma, float &energy, int it)
{
	//Compute Eigendecomposition
	Eigen::Vector3f S;
	
	Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);
	epsilon = eigenSolver.pseudoEigenvalueMatrix();
	S[0] = epsilon(0, 0);
	S[1] = epsilon(1, 1);
	S[2] = epsilon(2, 2);
	V = eigenSolver.pseudoEigenvectors();

	for (int i = 0; i < 3; ++i)
	{
		if (S[i] < 0.0f)
		{
			S[i] = 0.0f;
		}
	}

	// Detect if V is a reflection .
	// Make a rotation out of it by multiplying one column with -1.
	const float detV = V.determinant();
	if (detV < 0.0)
	{
		float minLambda = FLT_MAX;
		unsigned char pos = 0;
		for (unsigned char l = 0; l < 3; l++)
		{
			if (S[l] < minLambda)
			{
				pos = l;
				minLambda = S[l];
			}
		}
		V(0, pos) = -V(0, pos);
		V(1, pos) = -V(1, pos);
		V(2, pos) = -V(2, pos);
	}

	Eigen::Vector3f hatF;
	hatF[0] = sqrtf(S[0]);
	hatF[1] = sqrtf(S[1]);
	hatF[2] = sqrtf(S[2]);

	Eigen::Matrix3f VT;
	VT = V.transpose();

	//
	// Check for values of hatF near zero
	//
	unsigned char chk = 0;
	unsigned char pos = 0;
	for (unsigned char l = 0; l < 3; l++)
	{
		if (fabs(hatF[l]) < 1.0e-4f)
		{
			pos = l;
			chk++;
		}
	}

	if (chk > 0)
	{
		if (chk > 1)
		{
			U.setIdentity();
		}
		else
		{
			U = F * V;
			for (unsigned char l = 0; l < 3; l++)
			{
				if (l != pos)
				{
					for (unsigned char m = 0; m < 3; m++)
					{
						U(m, l) *= 1.0f / hatF[l];
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
	}
	else
	{
		Eigen::Vector3f hatFInv(1.0f / hatF[0], 1.0f / hatF[1], 1.0f / hatF[2]);
		U = F * V;
		for (unsigned char l = 0; l < 3; l++)
		{
			for (unsigned char m = 0; m < 3; m++)
			{
				U(m, l) *= hatFInv[l];
			}
		}
	}

	const float detU = U.determinant();

	// U is a reflection => tet is inverted
	if (detU < 0.0f)
	{
		//std::cout << "Inverted tet!\n";
		float minLambda = FLT_MAX;
		unsigned char pos = 0;
		for (unsigned char l = 0; l < 3; l++)
		{
			if (hatF[l] < minLambda)
			{
				pos = l;
				minLambda = hatF[l];
			}
		}
		hatF[pos] = -hatF[pos];
		U(0, pos) = -U(0, pos);
		U(1, pos) = -U(1, pos);
		U(2, pos) = -U(2, pos);
	}

	// Clamp small singular values
	//const float minXVal = 0.577f;
	const float minXVal = 0.177f;
	
	// epsilon for hatF
	Eigen::Vector3f epsilonHatF(0.5f*(hatF[0] * hatF[0] - 1.0f), 0.5f*(hatF[1] * hatF[1] - 1.0f), 0.5f*(hatF[2] * hatF[2] - 1.0f));

	const float trace = epsilonHatF[0] + epsilonHatF[1] + epsilonHatF[2];
	const float ltrace = lambda*trace;
	Eigen::Vector3f sigmaVec = epsilonHatF * 2.0f*mu;
	sigmaVec[0] += ltrace;
	sigmaVec[1] += ltrace;
	sigmaVec[2] += ltrace;
	sigmaVec[0] = hatF[0] * sigmaVec[0];
	sigmaVec[1] = hatF[1] * sigmaVec[1];
	sigmaVec[2] = hatF[2] * sigmaVec[2];

	Eigen::Matrix3f sigmaDiag, epsDiag;

	sigmaDiag.row(0) = Eigen::Vector3f(sigmaVec[0], 0.0f, 0.0f);
	sigmaDiag.row(1) = Eigen::Vector3f(0.0f, sigmaVec[1], 0.0f);
	sigmaDiag.row(2) = Eigen::Vector3f(0.0f, 0.0f, sigmaVec[2]);

	epsDiag.row(0) = Eigen::Vector3f(epsilonHatF[0], 0.0f, 0.0f);
	epsDiag.row(1) = Eigen::Vector3f(0.0f, epsilonHatF[1], 0.0f);
	epsDiag.row(2) = Eigen::Vector3f(0.0f, 0.0f, epsilonHatF[2]);

	epsilon = U*epsDiag*VT;
	sigma = U*sigmaDiag*VT;

	float psi = 0.0f;
	for (unsigned char j = 0; j < 3; j++)
	{
		for (unsigned char k = 0; k < 3; k++)
		{
			psi += epsilon(j, k) * epsilon(j, k);
		}
	}

	psi = mu*psi + 0.5f*lambda * trace*trace;

	//if (psi < 0.0001f)
	//std::cout << psi << std::endl;

	//if (psi > 1.0e+6f)
	//{
	//	std::cout << "Mu: " << mu << "; Psi: " << psi << "; lambda: " << lambda << "; trace^2 " << trace << std::endl;
	//	std::cout << hatF << std::endl;
	//	std::cout << "F" << std::endl;
	//	std::cout << F << std::endl;
	//	std::cout << "epsilon" << std::endl;
	//	std::cout << epsilon << std::endl;
	//	psi = 1e-6f;
	//}
	energy = restVolume*psi;

	//if (it < 5)
	//{
	//	if (energy < -1.0e-6f)
	//	{
	//		energy = -1.0e-6f;
	//	}

	//	if (energy > 1.0e-6f)
	//	{
	//		energy = 1.0e-6f;
	//	}
	//}
}

void
PBDSolver::computeGreenStrainAndPiolaStress(const Eigen::Matrix3f &F,
	const float restVolume,
	const float mu, const float lambda, Eigen::Matrix3f &epsilon, Eigen::Matrix3f &sigma, float &energy)
{
	// epsilon = 1/2 F^T F - I

	epsilon(0, 0) = 0.5f*(F(0, 0) * F(0, 0) + F(1, 0) * F(1, 0) + F(2, 0) * F(2, 0) - 1.0f);		// xx
	epsilon(1, 1) = 0.5f*(F(0, 1) * F(0, 1) + F(1, 1) * F(1, 1) + F(2, 1) * F(2, 1) - 1.0f);		// yy
	epsilon(2, 2) = 0.5f*(F(0, 2) * F(0, 2) + F(1, 2) * F(1, 2) + F(2, 2) * F(2, 2) - 1.0f);		// zz
	epsilon(0, 1) = 0.5f*(F(0, 0) * F(0, 1) + F(1, 0) * F(1, 1) + F(2, 0) * F(2, 1));			// xy
	epsilon(0, 2) = 0.5f*(F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2));			// xz
	epsilon(1, 2) = 0.5f*(F(0, 1) * F(0, 2) + F(1, 1) * F(1, 2) + F(2, 1) * F(2, 2));			// yz
	epsilon(1, 0) = epsilon(0, 1);
	epsilon(2, 0) = epsilon(0, 2);
	epsilon(2, 1) = epsilon(1, 2);

	// P(F) = F(2 mu E + lambda tr(E)I) => E = green strain
	const float trace = epsilon(0, 0) + epsilon(1, 1) + epsilon(2, 2);
	const float ltrace = lambda*trace;
	sigma = epsilon * 2.0f*mu;
	sigma(0, 0) += ltrace;
	sigma(1, 1) += ltrace;
	sigma(2, 2) += ltrace;
	sigma = F * sigma;

	float psi = 0.0;
	for (unsigned char j = 0; j < 3; j++)
	for (unsigned char k = 0; k < 3; k++)
		psi += epsilon(j, k) * epsilon(j, k);
	psi = mu*psi + 0.5f*lambda * trace*trace;
	energy = restVolume * psi;
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


void
PBDSolver::projectConstraintsNeoHookeanMixture(std::vector<PBDTetrahedra3d>& tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	float w1;
	float w2;
	float restDistance;

	Eigen::Vector3f x1;
	Eigen::Vector3f x2;

	Eigen::Vector3f temp;


	Eigen::Matrix3f F;
	Eigen::Matrix3f FInverseTranspose;
	Eigen::Matrix3f FTransposeF;

	Eigen::Matrix3f PF;
	Eigen::Matrix3f gradientTemp;
	Eigen::MatrixXf gradient; gradient.resize(3, 4);

	Eigen::Matrix3f U;
	Eigen::Matrix3f V;
	bool isInverted;

	Eigen::Matrix3f S;
	Eigen::Matrix3f Fhat;

	Eigen::Vector3f deltaX;

	std::ofstream strainEnergyfile;

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

	//projectConstraintsDistance(tetrahedra, particles, 10, 1);

	for (int it = 0; it < settings.numConstraintIts; ++it)
	{
		//if (it == 10)
		//{
		//	projectConstraintsDistance(tetrahedra, particles, 1);
		//}

		for (int t = 0; t < settings.numTetrahedra; ++t)
		{
			float lagrangeM;

			//Get deformation gradient
			F = tetrahedra[t].getDeformationGradient();

			if (F.isIdentity())
			{
				continue;
			}

			isInverted = F.determinant() < 0.0;

			//check for inversion
			//if (isInverted)
			//{
			//	inversionHandled = true;
			//	//1. Compute Eigendecomposition
			//	Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);
			//	S = eigenSolver.pseudoEigenvalueMatrix();
			//	V = eigenSolver.pseudoEigenvectors();
			//	for (int i = 0; i < 3; ++i)
			//	{
			//		if (S(i, i) < 0.0f)
			//		{
			//			S(i, i) = 0.0f;
			//		}
			//	}
			//	//2.  Detect if V is a reflection .
			//	//    Make a rotation out of it by multiplying one column with -1.
			//	const float detV = V.determinant();
			//	if (detV < 0.0)
			//	{
			//		float minLambda = FLT_MAX;
			//		unsigned char pos = 0;
			//		for (unsigned char l = 0; l < 3; l++)
			//		{
			//			if (S(l, l) < minLambda)
			//			{
			//				pos = l;
			//				minLambda = S(l, l);
			//			}
			//		}
			//		V(0, pos) = -V(0, pos);
			//		V(1, pos) = -V(1, pos);
			//		V(2, pos) = -V(2, pos);
			//	}
			//	//3. Compute Fhat
			//	Fhat.setZero();
			//	Fhat(0, 0) = sqrtf(S(0, 0));
			//	Fhat(1, 1) = sqrtf(S(1, 1));
			//	Fhat(2, 2) = sqrtf(S(2, 2));
			//	//4. Compute U
			//	U = F * V * Fhat.inverse();
			//	//CORRECT U
			//	//
			//	// Check for values of hatF near zero
			//	//
			//	unsigned char chk = 0;
			//	unsigned char pos = 0;
			//	for (unsigned char l = 0; l < 3; l++)
			//	{
			//		if (fabs(Fhat(l, l)) < 1.0e-4f)
			//		{
			//			pos = l;
			//			chk++;
			//		}
			//	}
			//	if (chk > 0)
			//	{
			//		if (chk > 1)
			//		{
			//			U.setIdentity();
			//		}
			//		else
			//		{
			//			U = F * V;
			//			for (unsigned char l = 0; l < 3; l++)
			//			{
			//				if (l != pos)
			//				{
			//					for (unsigned char m = 0; m < 3; m++)
			//					{
			//						U(m, l) *= 1.0f / Fhat(l, l);
			//					}
			//				}
			//			}
			//			Eigen::Vector3f v[2];
			//			unsigned char index = 0;
			//			for (unsigned char l = 0; l < 3; l++)
			//			{
			//				if (l != pos)
			//				{
			//					v[index++] = Eigen::Vector3f(U(0, l), U(1, l), U(2, l));
			//				}
			//			}
			//			Eigen::Vector3f vec = v[0].cross(v[1]);
			//			vec.normalize();
			//			U(0, pos) = vec[0];
			//			U(1, pos) = vec[1];
			//			U(2, pos) = vec[2];
			//		}
			//	}
			//	else
			//	{
			//		Eigen::Vector3f hatFInv(1.0f / Fhat(0, 0), 1.0f / Fhat(1, 1), 1.0f / Fhat(2, 2));
			//		U = F * V;
			//		for (unsigned char l = 0; l < 3; l++)
			//		{
			//			for (unsigned char m = 0; m < 3; m++)
			//			{
			//				U(m, l) *= hatFInv[l];
			//			}
			//		}
			//	}
			//	//5. Check if U is also a rotation and correct
			//	if (U.determinant() < 0)
			//	{
			//		//find minimal element of U
			//		int minElementFhat = 0;
			//		float minElementFhatValue = Fhat(0, 0);
			//		for (int e = 0; e < 3; ++e)
			//		{
			//			if (Fhat(e, e) < minElementFhatValue)
			//			{
			//				minElementFhat = e;
			//				minElementFhatValue = Fhat(e, e);
			//			}
			//		}
			//		Fhat(minElementFhat, minElementFhat) *= -1.0;
			//		U.col(minElementFhat) *= -1.0;
			//	}
			//}

			FInverseTranspose = F.inverse().transpose();
			FTransposeF = F.transpose() * F;

			//Compute Isotropic Invariants
			float I1 = (FTransposeF).trace();
			float I3 = (FTransposeF).determinant();

			float logI3 = log(I3);

			//Compute Stress tensor
			//if (isInverted)
			//{
			//	PF = settings.mu * Fhat - settings.mu * Fhat.inverse().transpose()
			//		+ ((settings.lambda * logI3) / 2.0) * Fhat.inverse().transpose();
			//	PF = U * PF * V;
			//	I1 = (Fhat.transpose() * Fhat).trace();
			//	I3 = (Fhat.transpose() * Fhat).determinant();
			//	logI3 = log(I3);
			//}
			//else
			{
				PF = settings.mu * F - settings.mu * FInverseTranspose
					+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;
			}

			//Compute volume
			float Volume = tetrahedra[t].getUndeformedVolume();


			//std::cout << "Current Volume: " << tetrahedra[t].getVolume();
			if (tetrahedra[t].getVolume() < 0.00001)
			{
				std::cout << "Degenerate/Inverted Tetrahedron at " << t << "; V =  " << Volume << std::endl;
			}

			gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
			gradient.col(0) = gradientTemp.col(0);
			gradient.col(1) = gradientTemp.col(1);
			gradient.col(2) = gradientTemp.col(2);
			gradient.col(3) = -gradientTemp.rowwise().sum();

			//Compute Strain Energy density field
			float strainEnergy = Volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

			//if (strainEnergy < -1e-5f)
			//{
			//	strainEnergy = -1e-5f;
			//}

			//std::cout << "Strain Energy: " << strainEnergy << std::endl;

			//Compute Lagrange Multiplier

			float denominator = 0.0;

			for (int cI = 0; cI < 4; ++cI)
			{
				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					denominator += tetrahedra[t].get_x(cI).inverseMass()
						* gradient.col(cI).lpNorm<2>();

					//if (std::fabs(denominator) > 1e-06)
					//{
					//	std::cout << "Condition met!" << std::endl;
					//}

				}
			}

			//if (std::fabs(denominator) < 1e-06)
			//{
			//	//std::cout << "Skipping!" << std::endl;
			//	continue;
			//}
			//else
			{
				lagrangeM = -(strainEnergy / denominator);
			}

			//if (std::isnan(lagrangeM))
			//{
			//	//if (isInverted)
			//	//{
			//	//	std::cout << "NAN! ";
			//	//}
			//	//else
			//	//{
			//	//	std::cout << "NON-INV NAN! ";
			//	//}
			//	//
			//	//std::cout << "Deformation Gradient" << std::endl;
			//	//std::cout << F << std::endl;
			//	//std::cout << "Inverse of deformation gradient:" << std::endl;
			//	//std::cout << F.inverse().transpose() << std::endl;
			//	//std::cout << "Stress Tensor" << std::endl;
			//	//std::cout << PF << std::endl;
			//	//std::cout << "Tensor Gradient " << std::endl;
			//	//std::cout << gradient << std::endl;
			//	//std::cout << "Strain Energy: " << strainEnergy << std::endl;
			//	//std::cout << "Denominator: " << denominator << std::endl;
			//	//std::cout << "Lagrange Multiplier: " << lagrangeM << std::endl;
			//	////std::cout << "Inverse Mass: " << tetrahedra[t].get_x(c).inverseMass() << std::endl;
			//	//std::cout << "Undeformed Volume: " << V << std::endl;
			//	//
			//	//std::cout << "STEPS: " << std::endl;

			//	//std::cout << (settings.mu * F) << std::endl;
			//	//std::cout << settings.mu * F.inverse().transpose() << std::endl;
			//	//std::cout << log(I3) << std::endl;
			//	//std::cout << F.inverse().transpose() << std::endl;

			//	//lagrangeM = 0.0;
			//}

			if (std::isnan(lagrangeM) || std::isinf(lagrangeM) || isInverted || std::abs(lagrangeM) > 0.5f)
			{
				w1 = tetrahedra[t].get_x(0).inverseMass();
				x1 = tetrahedra[t].get_x(0).position();

				w2 = tetrahedra[t].get_x(2).inverseMass();
				x2 = tetrahedra[t].get_x(2).position();
				computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(0),
					x1, x2, temp, deltaX);
				tetrahedra[t].get_x(0).position() += -deltaX * w1;
				tetrahedra[t].get_x(2).position() += deltaX * w2;

				w2 = tetrahedra[t].get_x(3).inverseMass();
				x2 = tetrahedra[t].get_x(3).position();
				computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(1),
					x1, x2, temp, deltaX);
				tetrahedra[t].get_x(0).position() += -deltaX * w1;
				tetrahedra[t].get_x(3).position() += deltaX * w2;

				w2 = tetrahedra[t].get_x(1).inverseMass();
				x2 = tetrahedra[t].get_x(1).position();
				computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(2),
					x1, x2, temp, deltaX);
				tetrahedra[t].get_x(0).position() += -deltaX * w1;
				tetrahedra[t].get_x(1).position() += deltaX * w2;

				//---------------------------------------------------

				w1 = tetrahedra[t].get_x(2).inverseMass();
				x1 = tetrahedra[t].get_x(2).position();

				w2 = tetrahedra[t].get_x(3).inverseMass();
				x2 = tetrahedra[t].get_x(3).position();
				computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(3),
					x1, x2, temp, deltaX);
				tetrahedra[t].get_x(2).position() += -deltaX * w1;
				tetrahedra[t].get_x(3).position() += deltaX * w2;


				w2 = tetrahedra[t].get_x(1).inverseMass();
				x2 = tetrahedra[t].get_x(1).position();
				computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(4),
					x1, x2, temp, deltaX);
				tetrahedra[t].get_x(2).position() += -deltaX * w1;
				tetrahedra[t].get_x(1).position() += deltaX * w2;

				//---------------------------------------------------

				w1 = tetrahedra[t].get_x(3).inverseMass();
				x1 = tetrahedra[t].get_x(3).position();

				w2 = tetrahedra[t].get_x(1).inverseMass();
				x2 = tetrahedra[t].get_x(1).position();
				computeDeltaXPositionConstraint(w1, w2, tetrahedra[t].getUndeformedSideLength(5),
					x1, x2, temp, deltaX);
				tetrahedra[t].get_x(3).position() += -deltaX * w1;
				tetrahedra[t].get_x(1).position() += deltaX * w2;
			}
			else
			{
				for (int cI = 0; cI < 4; ++cI)
				{
					if (tetrahedra[t].get_x(cI).inverseMass() != 0)
					{
						deltaX = (tetrahedra[t].get_x(cI).inverseMass()
							* lagrangeM) * gradient.col(cI);

						tetrahedra[t].get_x(cI).position() += deltaX;

						//std::cout << "[ " << cI << "] : " << std::endl;
						//std::cout << deltaX << std::endl;
					}
				}
			}

		}


		if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
		{
			calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
		}

	}

	if (settings.printStrainEnergyToFile)
	{
		strainEnergyfile.close();
	}

	if (inversionHandled)
	{
		std::cout << "Inversion handled successfully!" << std::endl;
	}
}

void
PBDSolver::projectConstraintsMooneyRivlin(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	Eigen::Matrix3f F;
	Eigen::Matrix3f FInverseTranspose;
	Eigen::Matrix3f FTransposeF;

	Eigen::Matrix3f PF;
	Eigen::Matrix3f gradientTemp;
	Eigen::MatrixXf gradient; gradient.resize(3, 4);

	Eigen::Matrix3f U;
	Eigen::Matrix3f V;
	bool isInverted;

	Eigen::Matrix3f S;
	Eigen::Matrix3f Fhat;

	Eigen::Vector3f deltaX;

	std::ofstream strainEnergyfile;

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

	for (int it = 0; it < settings.numConstraintIts; ++it)
	{
		for (int t = 0; t < settings.numTetrahedra; ++t)
		{
			float lagrangeM;

			//Get deformation gradient
			F = tetrahedra[t].getDeformationGradient();

			if (F.isIdentity())
			{
				continue;
			}

			isInverted = F.determinant() < 0.0;

			//check for inversion
			if (isInverted)
			{
				inversionHandled = true;
				//1. Compute Eigendecomposition
				Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);
				S = eigenSolver.pseudoEigenvalueMatrix();
				V = eigenSolver.pseudoEigenvectors();
				for (int i = 0; i < 3; ++i)
				{
					if (S(i, i) < 0.0f)
					{
						S(i, i) = 0.0f;
					}
				}

				//2.  Detect if V is a reflection .
				//    Make a rotation out of it by multiplying one column with -1.
				const float detV = V.determinant();
				if (detV < 0.0)
				{
					float minLambda = FLT_MAX;
					unsigned char pos = 0;
					for (unsigned char l = 0; l < 3; l++)
					{
						if (S(l, l) < minLambda)
						{
							pos = l;
							minLambda = S(l, l);
						}
					}
					V(0, pos) = -V(0, pos);
					V(1, pos) = -V(1, pos);
					V(2, pos) = -V(2, pos);
				}

				//3. Compute Fhat
				Fhat.setZero();
				Fhat(0, 0) = sqrtf(S(0, 0));
				Fhat(1, 1) = sqrtf(S(1, 1));
				Fhat(2, 2) = sqrtf(S(2, 2));

				//4. Compute U
				U = F * V * Fhat.inverse();

				//CORRECT U
				//
				// Check for values of hatF near zero
				//
				unsigned char chk = 0;
				unsigned char pos = 0;
				for (unsigned char l = 0; l < 3; l++)
				{
					if (fabs(Fhat(l, l)) < 1.0e-4f)
					{
						pos = l;
						chk++;
					}
				}

				if (chk > 0)
				{
					if (chk > 1)
					{
						U.setIdentity();
					}
					else
					{
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
				}
				else
				{
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

				//5. Check if U is also a rotation and correct
				if (U.determinant() < 0)
				{
					//find minimal element of U
					int minElementFhat = 0;
					float minElementFhatValue = Fhat(0, 0);
					for (int e = 0; e < 3; ++e)
					{
						if (Fhat(e, e) < minElementFhatValue)
						{
							minElementFhat = e;
							minElementFhatValue = Fhat(e, e);
						}
					}

					Fhat(minElementFhat, minElementFhat) *= -1.0;
					U.col(minElementFhat) *= -1.0;
				}

			}

			FInverseTranspose = F.inverse().transpose();
			FTransposeF = F.transpose() * F;

			//Compute Isotropic Invariants
			float I1 = (FTransposeF).trace();
			float I2 = 0.5f * (std::pow(I1, 2.0) - (FTransposeF * FTransposeF).trace());
			float I3 = (FTransposeF).determinant();

			float logI3 = log(I3);

			//Compute Stress tensor

			PF = settings.mu * F - settings.mu * FInverseTranspose
				+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

			//Compute volume
			float Volume = tetrahedra[t].getUndeformedVolume();


			//std::cout << "Current Volume: " << tetrahedra[t].getVolume();
			if (tetrahedra[t].getVolume() < 0.00001)
			{
				std::cout << "Degenerate/Inverted Tetrahedron at " << t << "; V =  " << Volume << std::endl;
			}

			gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
			gradient.col(0) = gradientTemp.col(0);
			gradient.col(1) = gradientTemp.col(1);
			gradient.col(2) = gradientTemp.col(2);
			gradient.col(3) = -gradientTemp.rowwise().sum();

			//Compute Strain Energy density field
			float strainEnergy = Volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

			//if (strainEnergy < -1e-5f)
			//{
			//	strainEnergy = -1e-5f;
			//}

			//std::cout << "Strain Energy: " << strainEnergy << std::endl;

			//Compute Lagrange Multiplier

			float denominator = 0.0;

			for (int cI = 0; cI < 4; ++cI)
			{
				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					denominator += tetrahedra[t].get_x(cI).inverseMass()
						* gradient.col(cI).lpNorm<2>();

					//if (std::fabs(denominator) > 1e-06)
					//{
					//	std::cout << "Condition met!" << std::endl;
					//}

				}
			}

			//if (std::fabs(denominator) < 1e-06)
			//{
			//	//std::cout << "Skipping!" << std::endl;
			//	continue;
			//}
			//else
			{
				lagrangeM = -(strainEnergy / denominator);
			}

			if (std::isnan(lagrangeM))
			{
				//if (isInverted)
				//{
				//	std::cout << "NAN! ";
				//}
				//else
				//{
				//	std::cout << "NON-INV NAN! ";
				//}
				//
				//std::cout << "Deformation Gradient" << std::endl;
				//std::cout << F << std::endl;
				//std::cout << "Inverse of deformation gradient:" << std::endl;
				//std::cout << F.inverse().transpose() << std::endl;
				//std::cout << "Stress Tensor" << std::endl;
				//std::cout << PF << std::endl;
				//std::cout << "Tensor Gradient " << std::endl;
				//std::cout << gradient << std::endl;
				//std::cout << "Strain Energy: " << strainEnergy << std::endl;
				//std::cout << "Denominator: " << denominator << std::endl;
				//std::cout << "Lagrange Multiplier: " << lagrangeM << std::endl;
				////std::cout << "Inverse Mass: " << tetrahedra[t].get_x(c).inverseMass() << std::endl;
				//std::cout << "Undeformed Volume: " << V << std::endl;
				//
				//std::cout << "STEPS: " << std::endl;

				//std::cout << (settings.mu * F) << std::endl;
				//std::cout << settings.mu * F.inverse().transpose() << std::endl;
				//std::cout << log(I3) << std::endl;
				//std::cout << F.inverse().transpose() << std::endl;

				lagrangeM = 0.0;
			}

			for (int cI = 0; cI <4; ++cI)
			{
				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					deltaX = (tetrahedra[t].get_x(cI).inverseMass()
						* lagrangeM) * gradient.col(cI);

					tetrahedra[t].get_x(cI).position() += deltaX;

					//std::cout << "[ " << cI << "] : " << std::endl;
					//std::cout << deltaX << std::endl;
				}
			}

		}


		if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
		{
			calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
		}

	}

	if (settings.printStrainEnergyToFile)
	{
		strainEnergyfile.close();
	}

	if (inversionHandled)
	{
		std::cout << "Inversion handled successfully!" << std::endl;
	}
}


void
PBDSolver::projectConstraintsGeometricInversionHandling(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	float w1;
	float w2;
	float restDistance;

	Eigen::Vector3f x1;
	Eigen::Vector3f x2;

	Eigen::Vector3f temp;


	Eigen::Matrix3f F;
	Eigen::Matrix3f FInverseTranspose;
	Eigen::Matrix3f FTransposeF;

	Eigen::Matrix3f PF;
	Eigen::Matrix3f gradientTemp;
	Eigen::MatrixXf gradient; gradient.resize(3, 4);

	Eigen::Matrix3f U;
	Eigen::Matrix3f V;
	bool isInverted;

	Eigen::Matrix3f S;
	Eigen::Matrix3f Fhat;

	Eigen::Vector3f deltaX;

	std::ofstream strainEnergyfile;

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

	//projectConstraintsDistance(tetrahedra, particles, 10, 1);

	for (int it = 0; it < settings.numConstraintIts; ++it)
	{
		//if (it == 10)
		//{
		//	projectConstraintsDistance(tetrahedra, particles, 1);
		//}

		for (int t = 0; t < settings.numTetrahedra; ++t)
		{
			float lagrangeM;

			//Get deformation gradient
			F = tetrahedra[t].getDeformationGradient();

			if (F.isIdentity())
			{
				continue;
			}

			isInverted = F.determinant() < 0.0;

			FInverseTranspose = F.inverse().transpose();
			FTransposeF = F.transpose() * F;

			//Compute Isotropic Invariants
			float I1 = (FTransposeF).trace();
			float I3 = (FTransposeF).determinant();

			float logI3 = log(I3);

			//Compute Stress tensor
			
			PF = settings.mu * F - settings.mu * FInverseTranspose
				+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

			//Compute volume
			float Volume = tetrahedra[t].getUndeformedVolume();


			//std::cout << "Current Volume: " << tetrahedra[t].getVolume();
			//if (tetrahedra[t].getVolume() < 0.00001)
			//{
			//	std::cout << "Degenerate/Inverted Tetrahedron at " << t << "; V =  " << Volume << std::endl;
			//}

			gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
			gradient.col(0) = gradientTemp.col(0);
			gradient.col(1) = gradientTemp.col(1);
			gradient.col(2) = gradientTemp.col(2);
			gradient.col(3) = -gradientTemp.rowwise().sum();

			//Compute Strain Energy density field
			float strainEnergy = Volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

			//Compute Lagrange Multiplier

			float denominator = 0.0;

			for (int cI = 0; cI < 4; ++cI)
			{
				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					denominator += tetrahedra[t].get_x(cI).inverseMass()
						* gradient.col(cI).lpNorm<2>();

					//if (std::fabs(denominator) > 1e-06)
					//{
					//	std::cout << "Condition met!" << std::endl;
					//}
				}
			}

			//if (std::fabs(denominator) < 1e-06)
			//{
			//	//std::cout << "Skipping!" << std::endl;
			//	continue;
			//}
			//else
			{
				lagrangeM = -(strainEnergy / denominator);
			}

			//if (std::isnan(lagrangeM))
			//{
			//	//if (isInverted)
			//	//{
			//	//	std::cout << "NAN! ";
			//	//}
			//	//else
			//	//{
			//	//	std::cout << "NON-INV NAN! ";
			//	//}
			//	//
			//	//std::cout << "Deformation Gradient" << std::endl;
			//	//std::cout << F << std::endl;
			//	//std::cout << "Inverse of deformation gradient:" << std::endl;
			//	//std::cout << F.inverse().transpose() << std::endl;
			//	//std::cout << "Stress Tensor" << std::endl;
			//	//std::cout << PF << std::endl;
			//	//std::cout << "Tensor Gradient " << std::endl;
			//	//std::cout << gradient << std::endl;
			//	//std::cout << "Strain Energy: " << strainEnergy << std::endl;
			//	//std::cout << "Denominator: " << denominator << std::endl;
			//	//std::cout << "Lagrange Multiplier: " << lagrangeM << std::endl;
			//	////std::cout << "Inverse Mass: " << tetrahedra[t].get_x(c).inverseMass() << std::endl;
			//	//std::cout << "Undeformed Volume: " << V << std::endl;
			//	//
			//	//std::cout << "STEPS: " << std::endl;

			//	//std::cout << (settings.mu * F) << std::endl;
			//	//std::cout << settings.mu * F.inverse().transpose() << std::endl;
			//	//std::cout << log(I3) << std::endl;
			//	//std::cout << F.inverse().transpose() << std::endl;

			//	//lagrangeM = 0.0;
			//}

			//if (true)
			//if (std::isnan(lagrangeM) || std::isinf(lagrangeM) || isInverted || std::abs(lagrangeM) > 0.01f || std::abs(F.determinant()) < 1e-06f)
			{
				Eigen::Vector3f p0 = tetrahedra[t].get_x(0).position();
				Eigen::Vector3f p1 = tetrahedra[t].get_x(1).position();
				Eigen::Vector3f p2 = tetrahedra[t].get_x(2).position();
				Eigen::Vector3f p3 = tetrahedra[t].get_x(3).position();
				gradient.col(0) = (p1 - p2).cross(p3 - p2);
				gradient.col(1) = (p2 - p0).cross(p3 - p0);
				gradient.col(2) = (p0 - p1).cross(p3 - p1);
				gradient.col(3) = (p1 - p0).cross(p2 - p0);

				//gradient.col(0) = (tetrahedra[t].get_x(1).position() - tetrahedra[t].get_x(2).position()).cross(tetrahedra[t].get_x(3).position() - tetrahedra[t].get_x(2).position());
				//gradient.col(1) = (tetrahedra[t].get_x(2).position() - tetrahedra[t].get_x(0).position()).cross(tetrahedra[t].get_x(3).position() - tetrahedra[t].get_x(0).position());
				//gradient.col(2) = (tetrahedra[t].get_x(0).position() - tetrahedra[t].get_x(1).position()).cross(tetrahedra[t].get_x(3).position() - tetrahedra[t].get_x(1).position());
				//gradient.col(3) = (tetrahedra[t].get_x(1).position() - tetrahedra[t].get_x(0).position()).cross(tetrahedra[t].get_x(2).position() - tetrahedra[t].get_x(0).position());
				lagrangeM = 0.0f;
				for (int cI = 0; cI < 4; ++cI)
				{
					lagrangeM += gradient.col(cI).squaredNorm() * tetrahedra[t].get_x(cI).inverseMass();
				}

				if (std::abs(lagrangeM) < 1e-9f)
				{
					continue;
				}


				//float volume = 1.0f / 6.0f * (p1 - p0).cross(p2 - p0).dot(p3 - p0);

				lagrangeM = (tetrahedra[t].getVolume() - tetrahedra[t].getUndeformedVolume()) / lagrangeM;

				for (int cI = 0; cI < 4; ++cI)
				{
					deltaX = -lagrangeM * tetrahedra[t].get_x(cI).inverseMass() * gradient.col(cI);

					tetrahedra[t].get_x(cI).position() += deltaX;
				}

			}
			//else
			//{
			//	for (int cI = 0; cI < 4; ++cI)
			//	{
			//		if (tetrahedra[t].get_x(cI).inverseMass() != 0)
			//		{
			//			deltaX = (tetrahedra[t].get_x(cI).inverseMass()
			//				* lagrangeM) * gradient.col(cI);

			//			tetrahedra[t].get_x(cI).position() += deltaX;

			//			if (std::isnan(deltaX[0]) || std::isinf(deltaX[0])
			//				|| std::isnan(deltaX[1]) || std::isinf(deltaX[1])
			//				|| std::isnan(deltaX[2]) || std::isinf(deltaX[2]))
			//			{
			//					std::cout << "Deformation Gradient" << std::endl;
			//					std::cout << F << std::endl;
			//					std::cout << "Inverse of deformation gradient:" << std::endl;
			//					std::cout << F.inverse().transpose() << std::endl;
			//					std::cout << "Stress Tensor" << std::endl;
			//					std::cout << PF << std::endl;
			//					std::cout << "Tensor Gradient " << std::endl;
			//					std::cout << gradient << std::endl;
			//					std::cout << "Strain Energy: " << strainEnergy << std::endl;
			//					std::cout << "Denominator: " << denominator << std::endl;
			//					std::cout << "Lagrange Multiplier: " << lagrangeM << std::endl;
			//					std::cout << "Inverse Mass: " << tetrahedra[t].get_x(cI).inverseMass() << std::endl;
			//					std::cout << "Undeformed Volume: " << V << std::endl;
			//					
			//					std::cout << "STEPS: " << std::endl;

			//					std::cout << (settings.mu * F) << std::endl;
			//					std::cout << settings.mu * F.inverse().transpose() << std::endl;
			//					std::cout << log(I3) << std::endl;
			//					std::cout << F.inverse().transpose() << std::endl;

			//			}

			//			//std::cout << "[ " << cI << "] : " << std::endl;
			//			//std::cout << deltaX << std::endl;
			//		}
			//	}
			//}

		}


		if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
		{
			calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
		}

	}

	if (settings.printStrainEnergyToFile)
	{
		strainEnergyfile.close();
	}

	if (inversionHandled)
	{
		std::cout << "Inversion handled successfully!" << std::endl;
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

			if(solveVolumeConstraint(p0, m1, p1, m2, p2, m3, p3, m4, tetrahedra[t].getUndeformedVolumeAlternative(),
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