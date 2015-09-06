#pragma once

#include <iostream>
#include <vector>

#include "PBDSolverTracker.h"
#include "commonMath.h"

struct PBDSolverSettings
{
	bool disableConstraintProjection;
	bool disablePositionCorrection;

	float deltaT;

	int currentFrame;

	//EXTERNAL GLOBAL FORCES
	float gravity;
	Eigen::Vector3f externalForce;

	float forceMultiplicationFactor;

	//EXTERNAL POSITION
	Eigen::Vector3f externalPositionDelta;
	float positionDeltaMultiplicationFactor;

	int numConstraintIts;

	int numTetrahedraIterations;

	//Lame coefficients
	float mu;
	float lambda;

	float youngsModulus;
	float poissonRatio;

	float minYoungsModulus;

	//Viscoelasticity
	float alpha;
	float rho;

	bool useFullPronySeries;
	std::vector<float> fullAlpha;
	std::vector<float> fullRho;

	//For SOR
	float w;

	float inverseMass;

	bool useSOR;


	bool useGeometricConstraintLimits;
	bool correctStrongForcesWithSubteps;

	//debug print
	bool printStrainEnergy;
	bool printStrainEnergyToFile;
	bool printStressComponentsToFile;

	//TRACKING

	//Tracker for values in the solver
	PBDSolverTracker tracker;

	bool trackS;
	bool trackF;
	bool trackPF;
	bool trackAverageDeltaXLength;
	bool trackSpecificPosition;
	int trackSpecificPositionIdx;

	//Collision
	bool enableGroundPlaneCollision;
	float groundplaneHeight;

	//Collision Rods
	std::vector<int> collisionSpheresNum;
	std::vector<float> collisionSpheresRadius;

	bool useMultiThreadedSolver;

	bool useSecondOrderUpdates;

	enum CONSTITUTIVE_MODEL
	{
		NEO_HOOKEAN,
		NEO_HOOKEAN_FIBER,
		RUBIN_BODNER
	};

	CONSTITUTIVE_MODEL materialModel;

	float MR_J;
	float MR_A;
	float MR_B;
	float MR_K;
	float MR_T;

	float MR_f_active;
	float MR_f_passive;
	float MR_alpha;

	float anisotropyParameter;
	Eigen::Vector3f MR_a;
	Eigen::Matrix3f MR_A0;


	bool disableInversionHandling;

	bool usePerTetMaterialAttributes;

	void initialise()
	{
		
		materialModel = NEO_HOOKEAN;
		externalForce.setZero();
		forceMultiplicationFactor = 0.0f;
		externalPositionDelta.setZero();
		positionDeltaMultiplicationFactor = 0.0f;
		anisotropyParameter = 0.0f;

		disableConstraintProjection = false;
		disablePositionCorrection = false;
		useFullPronySeries = false;

		printStrainEnergy = false;
		printStrainEnergyToFile = false;
		printStressComponentsToFile = false;

		trackS = false;
		trackF = false;
		trackPF = false;
		trackSpecificPosition = false;
		trackAverageDeltaXLength = false;
		enableGroundPlaneCollision = false;
		groundplaneHeight = 0.0f;

		disableInversionHandling = false;
		useMultiThreadedSolver = true;
		useSecondOrderUpdates = false;
		usePerTetMaterialAttributes = false;
		minYoungsModulus = 0.0f;
	}

	void print()
	{
		std::cout << "deltaT: " << deltaT << std::endl;

		std::cout << "gravity: " << gravity << std::endl;

		std::cout << "num Constraint its: " << numConstraintIts << std::endl;;

		std::cout << "mu: " << mu << std::endl;
		std::cout << "lambda: " << lambda << std::endl;

		std::cout << "Young's Modulus: " << youngsModulus << std::endl;
		std::cout << "Poisson Ratio: " << poissonRatio << std::endl;
	}

	float calculateLambda(float youngsModulus, float poissonRatio) const
	{
		return (youngsModulus * poissonRatio) / ((1.0f + poissonRatio) * (1.0f - 2.0f * poissonRatio));
	}

	float calculateMu(float youngsModulus, float poissonRatio) const
	{
		return youngsModulus / (2.0f * (1.0f + poissonRatio));
	}

	void calculateLambda()
	{
		lambda =  (youngsModulus * poissonRatio) / ((1.0f + poissonRatio) * (1.0f - 2.0f * poissonRatio));
	}

	void calculateMu()
	{
		mu = youngsModulus / (2.0f * (1.0f + poissonRatio));
	}

	void calculateFiberStructureTensor()
	{
		MR_a.normalize();

		MR_A0 = kroneckerProduct(MR_a, MR_a);
	}

	float getCurrentTime() const
	{
		return deltaT * (float)currentFrame;
	}
};