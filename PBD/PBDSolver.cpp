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

	Eigen::Matrix3f sigma;
	Eigen::Matrix3f epsilon;

	Eigen::Vector3f deltaX;

	for (int it = 0; it < settings.numConstraintIts; ++it)
	{
		for (int t = 0; t < settings.numTetrahedra; ++t)
		{
			float lagrangeM;
			float strainEnergy;

			//Compute volume
			float Volume = tetrahedra[t].getUndeformedVolume();

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

			//if (F.determinant() < 1.0e-04f || tetrahedra[t].getVolume() < 0.00001)
			if (F.determinant() < 0.0f)
			{
				computeGreenStrainAndPiolaStressInversion(F, Volume, settings.mu, settings.lambda, sigma, epsilon, strainEnergy);
				//computeGradCGreen(Volume, tetrahedra[t].getReferenceShapeMatrixInverseTranspose(), sigma, gradient);
			}
			else
			{
				computeGreenStrainAndPiolaStress(F, Volume, settings.mu, settings.lambda, sigma, epsilon, strainEnergy);
				//computeGradCGreen(Volume, tetrahedra[t].getReferenceShapeMatrixInverseTranspose(), sigma, gradient);
			}

			float strainEnergyClamp = 0.00001;

			if (strainEnergy < -strainEnergyClamp)
			{
				strainEnergy = -strainEnergyClamp;
			}

			if (strainEnergy > strainEnergyClamp)
			{
				strainEnergy = strainEnergyClamp;
			}


			//std::cout << "Current Volume: " << tetrahedra[t].getVolume();
			if (tetrahedra[t].getVolume() < 0.00001)
			{
				std::cout << "Degenerate/Inverted Tetrahedron at " << t << "; V =  " << Volume << std::endl;
			}

			gradientTemp = Volume * sigma * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
			gradient.col(0) = gradientTemp.col(0);
			gradient.col(1) = gradientTemp.col(1);
			gradient.col(2) = gradientTemp.col(2);
			gradient.col(3) = -gradientTemp.rowwise().sum();

			//gradientTemp = Volume * PF * tetrahedra[t].getReferenceShapeMatrixInverseTranspose();
			//gradient.col(0) = gradientTemp.col(0);
			//gradient.col(1) = gradientTemp.col(1);
			//gradient.col(2) = gradientTemp.col(2);
			//gradient.col(3) = -gradientTemp.rowwise().sum();

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

			//if (std::fabs(denominator) < 1.0e-6f)
			//{
			//	continue;
			//}
			//else
			{
				lagrangeM = - (strainEnergy / denominator);
			}



			for (int cI = 0; cI < 4; ++cI)
			{
				if (tetrahedra[t].get_x(cI).inverseMass() != 0)
				{
					deltaX = (tetrahedra[t].get_x(cI).inverseMass()
						* lagrangeM) * gradient.col(cI);

					tetrahedra[t].get_x(cI).position() += deltaX;
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
float& strainEnergy, float volume,
const PBDSolverSettings& settings)
{
	float determinantF = F.determinant();

	//No deformation
	if (F.isIdentity())
	{
		return false;
	}

	//NOT INVERTED
	if (determinantF > 0)
	{
		//Compute Stress tensor
		PF = settings.mu * F - settings.mu * FInverseTranspose
			+ ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

		//Compute Strain Energy density field
		strainEnergy = volume * (0.5 * settings.mu * (I1 - logI3 - 3.0) + (settings.lambda / 8.0) * std::pow(logI3, 2.0));

		return true;
	}

	//std::cout << "Reflection" << std::endl;

	//INVERTED
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);

	Eigen::Matrix3f S = eigenSolver.eigenvalues().asDiagonal().toDenseMatrix();
	V = eigenSolver.eigenvectors();

	for (int i = 0; i < 3; ++i)
	{
		if (S(i, i) < 0.0f)
		{
			S(i, i) = 0.0f;
		}
	}


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
	Eigen::Matrix3f  sigma = U * sigmaDiag * V.transpose();

	float psi = 0.0f;
	for (unsigned char j = 0; j < 3; j++)
	for (unsigned char k = 0; k < 3; k++)
		psi += epsilon(j, k) * epsilon(j, k);
	psi = settings.mu*psi + 0.5f*settings.lambda * trace*trace;
	strainEnergy = volume*psi;

	return true;
}

void
PBDSolver::computeGreenStrainAndPiolaStressInversion(const Eigen::Matrix3f& F,
	const float restVolume,
	const float mu, const float lambda, Eigen::Matrix3f &epsilon, Eigen::Matrix3f &sigma, float &energy)
{

	Eigen::Matrix3f FT_F;
	FT_F = F.transpose() * F;

	// Inversion handling
	Eigen::Matrix3f V;
	Eigen::Vector3f S;

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenSolver(FT_F);

	S = eigenSolver.eigenvalues();
	V = eigenSolver.eigenvectors();

	if (S[0] < 0.0f) S[0] = 0.0f;		// safety for sqrt
	if (S[1] < 0.0f) S[1] = 0.0f;
	if (S[2] < 0.0f) S[2] = 0.0f;

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

	Eigen::Matrix3f U;
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
	const float minXVal = 0.577f;

	for (unsigned char j = 0; j < 3; j++)
	{
		if (hatF[j] < minXVal)
			hatF[j] = minXVal;
	}

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
	for (unsigned char k = 0; k < 3; k++)
		psi += epsilon(j, k) * epsilon(j, k);
	psi = mu*psi + 0.5f*lambda * trace*trace;
	energy = restVolume*psi;
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
PBDSolver::computeGradCGreen(float restVolume, const Eigen::Matrix3f &invRestMat, const Eigen::Matrix3f &sigma, Eigen::MatrixXf& J)
{
	Eigen::Matrix3f H;
	Eigen::Matrix3f T;
	T = invRestMat.transpose();
	H = sigma * T * restVolume;

	J(0,0) = H(0, 0);
	J(1,0) = H(0, 1);
	J(2,0) = H(0, 2);
	J(0,1) = H(1, 0);
	J(1,1) = H(1, 1);
	J(2,1) = H(1, 2);
	J(0,2) = H(2, 0);
	J(1,2) = H(2, 1);
	J(2,2) = H(2, 2);

	J.col(3) = -J.col(0) - J.col(1) - J.col(2);
}
