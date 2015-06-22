#include "PBDSolver.h"

#include <iostream>

#include <tbb\parallel_for.h>
#include <tbb\mutex.h>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "EnergyConstraint.h"


PBDSolver::PBDSolver()
{
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
	projectConstraints(tetrahedra, particles, settings);
	//projectConstraintsSOR(tetrahedra, particles, settings, temporaryPositions, numConstraintInfluences);

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
PBDSolver::updateVelocities(std::vector<PBDTetrahedra3d>& tetrahedra,
std::shared_ptr<std::vector<PBDParticle>>& particles, const PBDSolverSettings& settings)
{
	for (auto& p : *particles)
	{
		p.velocity() = (1.0 / settings.deltaT) * (p.position() - p.previousPosition());
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

	Eigen::Vector3f deltaX;
	
	float lagrangeM;

	for (int it = 0; it < settings.numConstraintIts; ++it)
	{
		for (int t = 0; t < settings.numTetrahedra; ++t)
		{
			FInverseTranspose = F.inverse().transpose();
			FTransposeF = F.transpose() * F;

			//Compute Isotropic Invariants
			float I1 = (FTransposeF).trace();
			float I3 = (FTransposeF).determinant();

			float logI3 = log(I3);

			//Get deformation gradient
			F = tetrahedra[t].getDeformationGradient();

			if (F.isIdentity())
			{
				continue;
			}

			if (!correctInversion(F, FInverseTranspose, FTransposeF, PF, U, V, I1, 0.0, logI3, settings))
			{
				continue;
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
				//std::cout << "NAN!" << std::endl;
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

				continue;
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

		for (int it = 0; it < settings.numConstraintIts; ++it)
		{
			//reset the accmulator arrays
			for (int t = 0; t < particles->size(); ++t)
			{
				temporaryPositions[t].setZero();
				numConstraintInfluences[t] = 0;
			}

			boost::thread_group threads;

			int numThreads = 1;

			int stepSize = settings.numTetrahedra / numThreads;

			mutexStruct mutexInstance;
			//std::cout << settings.numTetrahedra << std::endl;
			for (int t = 0; t < settings.numTetrahedra; t += stepSize)
			{
				if (t < settings.numTetrahedra - stepSize - 1)
				{
					threads.create_thread(boost::bind(projectConstraintsSOR_CORE, boost::ref(mutexInstance), boost::ref(tetrahedra), boost::ref(particles), boost::ref(settings),
						boost::ref(temporaryPositions), boost::ref(numConstraintInfluences),
						t, t + stepSize - 1));
					//std::cout << t << "; " << t + stepSize - 1 << std::endl;
				}
				else
				{
					threads.create_thread(boost::bind(projectConstraintsSOR_CORE, boost::ref(mutexInstance), boost::ref(tetrahedra), boost::ref(particles), boost::ref(settings),
						boost::ref(temporaryPositions), boost::ref(numConstraintInfluences),
						t, settings.numTetrahedra - 1));
					//std::cout << t << "; " << settings.numTetrahedra - 1 << std::endl;
					continue;
				}
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
							tetrahedra[t].get_x(cI).position() += temporaryPositions[tetrahedra[t].getVertexIndices()[cI]] / numConstraintInfluences[tetrahedra[t].getVertexIndices()[cI]];
							//std::cout << temporaryPositions[tetrahedra[t].getVertexIndices()[cI]] << std::endl;
						}
					}
				}
			}
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

		isInverted = F.determinant() < 0.0;

		if (F.isZero())
		{
			continue;
			std::cout << "Zero F" << std::endl;
		}

		//check for inversion
		if (isInverted)
		{
			std::cout << "Inverted" << std::endl;
			//std::cout << F << std::endl;

			Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);

			U = svd.matrixU();
			V = svd.matrixV();

			F = svd.singularValues().asDiagonal().toDenseMatrix();

			/*std::cout << "F before: " << std::endl << F << std::endl;*/
			//std::cout << U << std::endl;
			//std::cout << V << std::endl;
			F(2, 2) *= -1.0;
			//std::cout << "U before: " << std::endl << U << std::endl;
			V.col(2) *= -1.0;
			//std::cout << "U after: " << std::endl << U << std::endl;

			//std::cout << "det U: " << U.determinant() << std::endl;
			//std::cout << "det V: " << V.determinant() << std::endl;

			//if (V.determinant() < 0.0)
			//{
			//	V.col(0) *= -1;
			//}

			//F = U * F * V;
			continue;
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

		//std::cout << PF << std::endl;

		if (isInverted)
		{
			for (int row = 0; row < 3; ++row)
			{
				for (int col = 0; col < 3; ++col)
				{
					//if (std::abs(PF(row, col)) > 0.0001)
					//{
					//	if (PF(row, col) > 0)
					//	{
					//		PF(row, col) = 0.0001;
					//	}
					//	else
					//	{
					//		PF(row, col) = -0.0001;
					//	}
					//}
					if (PF(row, col) < -0.0001)
					{
						PF(row, col) = -0.0001;
					}
				}
			}
			PF = V * PF * U;
			//std::cout << PF << std::endl;
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
const PBDSolverSettings& settings)
{
	float determinantF = F.determinant();

	//NOT INVERTED
	if (determinantF > 0)
	{
		//Compute Stress tensor
		PF = settings.mu * F - settings.mu * FInverseTranspose
			+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

		return true;
	}

	//INVERTED
	Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);

	Eigen::Matrix3f S = eigenSolver.eigenvalues().asDiagonal().toDenseMatrix();
	V = eigenSolver.eigenvectors();

	if (V.determinant() < 0.0)
	{
		int pos = 2;
		float smallestValue = 11111111111;
		for (int i = 0; i < 3; ++i)
		{
			pos = (S(i, i) < smallestValue) ? i : pos;
			smallestValue = (S(i, i) < smallestValue) ? S(i,i) : smallestValue;
		}

		V.col(pos) *= -1.0;
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
					U(m, l) *= 1.0f / Fhat(l,l);
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

	PF = U * Fhat * V.transpose();


	return true;
}
