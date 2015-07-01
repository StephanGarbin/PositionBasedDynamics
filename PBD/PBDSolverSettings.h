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
};