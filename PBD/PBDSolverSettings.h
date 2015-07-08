#pragma once


struct PBDSolverSettings
{
	float deltaT;

	float gravity;

	int numConstraintIts;

	int numTetrahedra;

	//Lame coefficients
	float mu;
	float lambda;

	float youngsModulus;
	float poissonRatio;

	//For SOR
	float w;

	bool useSOR;


	bool useGeometricConstraintLimits;
	bool correctStrongForcesWithSubteps;

	//debug print
	bool printStrainEnergy;
	bool printStrainEnergyToFile;

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
		return (youngsModulus * poissonRatio) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio));
	}

	float calculateMu(float youngsModulus, float poissonRatio) const
	{
		return youngsModulus / (2.0 * (1.0 + poissonRatio));
	}
};