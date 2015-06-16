#include "PBDSolver.h"

#include <iostream>

#include "EnergyConstraint.h"


PBDSolver::PBDSolver()
{
}


PBDSolver::~PBDSolver()
{
}

void
PBDSolver::advanceSystem(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	//Advance Velocities
	advanceVelocities(tetrahedra, particles, settings);

	//Advance Positions
	advancePositions(tetrahedra, particles, settings);

	//Project Constraints
	projectConstraints(tetrahedra, particles, settings);

	//Update Velocities
	updateVelocities(tetrahedra, particles, settings);

	//swap particles states
	for (auto& p : *particles)
	{
		p.swapStates();
	}
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
		p.position() = p.previousPosition() + settings.deltaT * p.velocity();
	}
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
	bool isInverted;

	Eigen::Vector3f deltaX;

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
				Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

				U = svd.matrixU();
				V = svd.matrixV();

				F = svd.singularValues().asDiagonal().toDenseMatrix();

				F(2, 2) *= -1;
				V.col(2) *= -1;

				std::cout << "Handling Inversion! " << std::endl;
			}



			FInverseTranspose = F.inverse().transpose();
			FTransposeF = F.transpose() * F;

			//Compute Isotropic Invariants
			float I1 = (FTransposeF).trace();
			float I3 = (FTransposeF).determinant();

			float logI3 = log(I3);

			//Compute Stress tensor
			PF = settings.mu * F - settings.mu * FInverseTranspose
				+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;


			if (isInverted)
			{
				PF = U * PF * V;
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
				lagrangeM = - (strainEnergy / denominator);
			}

			if (std::isnan(lagrangeM))
			{
				std::cout << "NAN!" << std::endl;
				//std::cout << "Deformation Gradient" << std::endl;
				//std::cout << F << std::endl;
				//std::cout << "Inverse of deformation gradient:" << std::endl;
				//std::cout << F.inverse().transpose() << std::endl;
				//std::cout << "Stress Tensor" << std::endl;
				//std::cout << PF << std::endl;
				//std::cout << "Tensor Gradient " << std::endl;
				//std::cout << gradient << std::endl;
				//std::cout << "Strain Energy: " << strainEnergy << std::endl;
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