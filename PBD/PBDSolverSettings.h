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
};