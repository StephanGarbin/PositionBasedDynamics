#pragma once

#include <cinder\Vector.h>
#include <vector>
#include <string>

struct Parameters
{
	//Material Model
	enum CONSTITUTIVE_MODEL
	{
		NEO_HOOKEAN,
		NEO_HOOKEAN_FIBER,
		NEO_HOOKEAN_FIBER_VISCOELASTIC
	};

	int materialModel;

	std::vector<std::string> materialModelNames;

	//Anisotropy
	cinder::Vec3f anisotropyDirection;
	float anisotropyStrength;
	void normaliseAnisotropyDirection()
	{
		anisotropyDirection.safeNormalize();
	}

	//Viscoelasticity
	float alpha;
	float rho;

	//Elasticity
	float youngsModulus;
	float poissonRatio;
	float inverseMass;
	float timeStep;
	int numConstraintIterations;
	int numGPUBlockIterations;
	float lambda;
	float mu;
	float gravity;

	std::string tetGenNodeFile;
	std::string tetGenElementFile;
	std::string positionConstraintFile;
	std::string meshOutputFile;

	//for GPU SOLVER
	int numThreadsPerBlock;
	int numBlocks;
	int trueNumberOfConstraints;
	int numParticles;

	void initialiseToDefaults()
	{
		materialModel = 0;
		alpha = 0.0f;
		rho = 0.0f;
		anisotropyStrength = 0.0f;
		anisotropyDirection.zero();

		youngsModulus = 20.0f;
		poissonRatio = 0.3;
		inverseMass = 1.0f;
		timeStep = 0.005f;
		numConstraintIterations = 5;
		numGPUBlockIterations = 1;

		gravity = -9.8f;

		tetGenNodeFile = "barout.node";
		tetGenElementFile = "barout.ele";
		positionConstraintFile = "barLowVertexConstraints.txt";
		meshOutputFile = "GPUSolverResult.abc";

		materialModelNames.push_back("Compr.Neo Hookean");
		materialModelNames.push_back("Compr. Ani. Neo-Hookean");
		materialModelNames.push_back("Compr. Ani. Visc. Neo-Hookean");
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
		lambda = (youngsModulus * poissonRatio) / ((1.0f + poissonRatio) * (1.0f - 2.0f * poissonRatio));
	}

	void calculateMu()
	{
		mu = youngsModulus / (2.0f * (1.0f + poissonRatio));
	}
};