#pragma once

#include <iostream>

#include "PBDSolverTracker.h"

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

	//Viscoelasticity
	float alpha;
	float rho;

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

	void initialise()
	{
		externalForce.setZero();
		forceMultiplicationFactor = 0.0f;
		externalPositionDelta.setZero();
		positionDeltaMultiplicationFactor = 0.0f;

		disableConstraintProjection = false;
		disablePositionCorrection = false;

		printStrainEnergy = false;
		printStrainEnergyToFile = false;
		printStressComponentsToFile = false;
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

	float getCurrentTime() const
	{
		return deltaT * (float)currentFrame;
	}
};