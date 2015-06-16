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
};