#pragma once

struct Parameters
{
	float youngsModulus;
	float poissonRatio;
	float inverseMass;
	float timeStep;
	int numConstraintIterations;
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

	void initialiseToDefaults()
	{
		youngsModulus = 1.0f;
		poissonRatio = 0.3;
		inverseMass = 1.0f;
		timeStep = 0.005f;
		numConstraintIterations = 300;

		gravity = -9.8f;

		tetGenNodeFile = "barout.node";
		tetGenElementFile = "barout.ele";
		positionConstraintFile = "barLowVertexConstraints.txt";
		meshOutputFile = "GPUSolverResult.abc";
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