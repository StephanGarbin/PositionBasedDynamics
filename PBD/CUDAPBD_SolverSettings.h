#pragma once

struct CUDAPBD_SolverSettings
{
	int numIterations;
	int numThreadsPerBlock;
	int numBlocks;

	int trueNumberOfConstraints;

	int numConstraints;

	float lambda;
	float mu;
};