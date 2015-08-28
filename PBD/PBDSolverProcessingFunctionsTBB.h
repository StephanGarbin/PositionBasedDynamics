#pragma once

#include <vector>
#include <string>

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "CollisionMesh.h"
#include "CollisionRod.h"
#include "PBDProbabilisticConstraint.h"
#include "PBDSolverSettings.h"


#include <tbb\parallel_for.h>
#include <tbb\mutex.h>
#include <tbb\queuing_mutex.h>
#include <tbb\blocked_range.h>


//typedef tbb::queuing_mutex currentMutex_t;

struct PBDSolverTBB
{
	PBDSolverTBB(std::vector<PBDTetrahedra3d>& in_tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& in_particles, PBDSolverSettings& in_settings,
	std::vector<PBDProbabilisticConstraint>& in_probabilisticConstraints,
	std::vector<CollisionMesh>& in_collisionGeometry,
	std::vector<CollisionRod>& in_collisionGeometry2,
	std::vector<CollisionSphere>& in_collisionGeometry3,
	tbb::queuing_mutex& in_mutex) : tetrahedra(in_tetrahedra), particles(in_particles),
	settings(in_settings), probabilisticConstraints(in_probabilisticConstraints),
	collisionGeometry(in_collisionGeometry), collisionGeometry2(in_collisionGeometry2), collisionGeometry3(in_collisionGeometry3), mutex(in_mutex)
	{
		//nothing else to do
	}

	std::vector<PBDTetrahedra3d>& tetrahedra;
	std::shared_ptr<std::vector<PBDParticle>>& particles;
	PBDSolverSettings& settings;
	std::vector<PBDProbabilisticConstraint>& probabilisticConstraints;
	std::vector<CollisionMesh>& collisionGeometry;
	std::vector<CollisionRod>& collisionGeometry2;
	std::vector<CollisionSphere>& collisionGeometry3;

	tbb::queuing_mutex& mutex;

	void operator()(const tbb::blocked_range<size_t>& r) const
	{
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

		//for (int it = 0; it < settings.numConstraintIts; ++it)
		{
			for (size_t t = r.begin(); t != r.end(); ++t)
			{
				float lagrangeM;
				float strainEnergy;
				float Volume;

				//Get deformation gradient
				F_orig = tetrahedra[t].getDeformationGradient();

				FTransposeF = F_orig.transpose() * F_orig;

				if (!settings.disableInversionHandling)
				{
					//Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(FTransposeF);
					//S = eigenSolver.pseudoEigenvalueMatrix(); //squared eigenvalues of F
					//V = eigenSolver.pseudoEigenvectors(); //eigenvectors
					eigenDecompositionCardano(FTransposeF, S, V);

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

					const float maxXVal = 500.0f;

					for (unsigned char j = 0; j < 3; j++)
					{
						if (std::abs(F(j, j)) > maxXVal)
						{
							if (F(j, j) < 0.0f)
							{
								F(j, j) = -maxXVal;
							}
							else
							{
								F(j, j) = maxXVal;
							}
						}
					}
				}

				if (settings.disableInversionHandling)
				{
					F = tetrahedra[t].getDeformationGradient();
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
					//PF = settings.mu * F - settings.mu * FInverseTranspose
					//   + ((settings.lambda * logI3) / 2.0) * FInverseTranspose;
					PF = settings.mu * F - settings.mu * FInverseTranspose;

					PF_vol = ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

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

					 if (settings.disableInversionHandling)
					 {
						 rotated_a = settings.MR_a;
					 }

					 //PF = settings.mu * F - settings.mu * FInverseTranspose;
					 //PF_vol = ((settings.lambda * logI3) / 2.0) * FInverseTranspose;
					 PF = settings.mu * F - settings.mu * FInverseTranspose;

					 PF_vol = ((settings.lambda * logI3) / 2.0) * FInverseTranspose;

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
				//if (settings.alpha != 0.0f && settings.rho != 0.0f)
				{
					//FInverseTranspose = F.inverse();
					/*PF = U * PF * V.transpose();*/

					//F_orig = tetrahedra[t].getDeformationGradient();

					PF = F.inverse() * PF;
					PF_vol = F.inverse() * PF_vol;

					//PF *= F.inverse();
					//PF_vol *= F.inverse();

					//PF_vol *= FInverseTranspose;

					/*FInverseTranspose = F_orig.inverse().transpose();

					PF *= FInverseTranspose.transpose();
					*/
					Eigen::Matrix3f vMult;

					if (settings.useFullPronySeries)
					{
						Eigen::Matrix3f temp;
						temp.setZero();

						for (int pComponent = 0; pComponent < settings.fullAlpha.size(); ++pComponent)
						{
							temp += (2.0f * settings.deltaT * settings.fullAlpha[pComponent] * PF
								+ settings.fullRho[pComponent] * tetrahedra[t].getFullUpsilon(pComponent)) / (settings.deltaT + settings.fullRho[pComponent]);

							tetrahedra[t].getFullUpsilon(pComponent) = (2.0f * settings.deltaT * settings.fullAlpha[pComponent] * PF
								+ settings.fullRho[pComponent] * tetrahedra[t].getFullUpsilon(pComponent)) / (settings.deltaT + settings.fullRho[pComponent]);
						}

						vMult = temp;
					}
					else
					{
						vMult = tetrahedra[t].getUpsilon();

						vMult = (2.0f * settings.deltaT * settings.alpha * PF + settings.rho * vMult) / (settings.deltaT + settings.rho);

						//if (it == settings.numConstraintIts - 1)
						{
							tetrahedra[t].getUpsilon() = vMult;
						}
					}

					PF = 2.0f * PF_vol + 2.0f * PF - vMult * 1.0f;

					//if (settings.trackS)
					//{
					//	settings.tracker.S.push_back(PF);
					//}

					//PF = F_orig * PF;
					//PF *= F;
					PF = F * PF;
				}

				//PF = U * PF * V.transpose();

				if (!settings.disableInversionHandling)
				{
					PF = U * PF * V.transpose();
				}

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
				if (denominator < 1e-20)
				{
					continue;
				}

				lagrangeM = -(strainEnergy / denominator);


				if (std::isnan(lagrangeM) || std::isinf(lagrangeM))
				{
					continue;
				}
				else
				{
					//std::cout << strainEnergy << ",";

					//Acquire Lock
					tbb::queuing_mutex::scoped_lock lock(mutex);

					for (int cI = 0; cI < 4; ++cI)
					{
						if (tetrahedra[t].get_x(cI).inverseMass() != 0)
						{
							if (!settings.disablePositionCorrection)
							{
								deltaX = (tetrahedra[t].get_x(cI).inverseMass()
									* lagrangeM) * gradient.col(cI);

								Eigen::Vector3f proposedEndpoint = tetrahedra[t].get_x(cI).position() + deltaX;
								float penetrationDistance;
								bool penetrates;
								for (int cS = 0; cS < collisionGeometry3.size(); ++cS)
								{
									collisionGeometry3[cS].checkForSinglePointIntersection_SAFE(proposedEndpoint, penetrationDistance, penetrates,
										settings.collisionSpheresRadius[cS]);

									if (penetrates)
									{
										Eigen::Vector3f temp = penetrationDistance * (tetrahedra[t].get_x(cI).position() - proposedEndpoint).normalized();
										if (std::sqrtf(temp.squaredNorm()) > std::sqrtf((tetrahedra[t].get_x(cI).position() - proposedEndpoint).squaredNorm()))
										{
											proposedEndpoint = tetrahedra[t].get_x(cI).position();
										}
										else
										{
											proposedEndpoint += penetrationDistance * (tetrahedra[t].get_x(cI).position() - proposedEndpoint).normalized();
										}
									}
								}
								tetrahedra[t].get_x(cI).position() = proposedEndpoint;

								//tetrahedra[t].get_x(cI).position() += deltaX;

								//if (settings.trackAverageDeltaXLength)
								//{
								//	float value = deltaX.squaredNorm();
								//	if (!std::isnan(value) && !std::isinf(value))
								//	{
								//		averageDeltaXLengthAccumulator += value;
								//		averageDeltaXLengthAccumulatorCounter += 1.0f;
								//	}
								//}
							}
						}
					}
				}
			}

			if (settings.printStrainEnergy || settings.printStrainEnergyToFile)
			{
				//calculateTotalStrainEnergy(tetrahedra, particles, settings, it, strainEnergyfile);
			}

			//if (settings.trackAverageDeltaXLength)
			//{
			//	if (averageDeltaXLengthAccumulatorCounter != 0.0f)
			//	{
			//		settings.tracker.averageDeltaXLength.push_back(averageDeltaXLengthAccumulator / averageDeltaXLengthAccumulatorCounter);
			//	}
			//	else
			//	{
			//		settings.tracker.averageDeltaXLength.push_back(0.0f);
			//	}
			//}

			//for (int pC = 0; pC < probabilisticConstraints.size(); ++pC)
			//{
			//	probabilisticConstraints[pC].project(*particles);
			//}

			//if (settings.printStrainEnergyToFile)
			//{
			//	strainEnergyfile.close();
			//}

			//if (inversionHandled)
			//{
			//	std::cout << "Inversion handled successfully!" << std::endl;
			//}
		}
	}
};
