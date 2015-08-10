#pragma once

#include "PBDSolverSettings.h"

struct Parameters
{
	//Global Settings
	bool disableSolver;
	int TEST_IDX;
	int TEST_VERSION;

	int frame2ApplyInitialDeformation;

	//---------------------------------

	int numMillisecondsToWaitBetweenFrames;

	int timingPrintInterval;
	int maxFrames;

	int globalHeight;
	int globalWidth;

	float executionTimeSum;

	//Camera
	float baryCentre[3];
	float radius;
	float rotation[3];
	float zoom;

	//Solver Type
	bool useFEMSolver;

	//Debug IO
	bool writeToAlembic;
	bool printStrainEnergyToFile;

	//Inversion Handling Test
	bool testingInversionHandling ;
	int dimToCollapse;

	//Mesh generation
	bool generateMeshInsteadOfDoingIO;
	bool generateMeshFromTrackingData;
	bool useTrackingConstraints;

	//Constraint IO
	bool readVertexConstraintData;

	bool applyInitialDeformationToMesh;

	//Settings for the solver
	PBDSolverSettings solverSettings;

	void initialiseToDefaults()
	{
		disableSolver = false;
		applyInitialDeformationToMesh = false;

		numMillisecondsToWaitBetweenFrames = 0;
		executionTimeSum = 0.0f;


		timingPrintInterval = 100;
		solverSettings.currentFrame = 1;
		maxFrames = 1000;

		useFEMSolver = false;
		writeToAlembic = true;
		printStrainEnergyToFile = false;

		testingInversionHandling = false;
		dimToCollapse = 1;

		generateMeshInsteadOfDoingIO = true;
		generateMeshFromTrackingData = false;

		readVertexConstraintData = false;
	}

	void initialiseCamera()
	{
		//above
		rotation[0] = 90.0f;
		rotation[1] = 360.0f;
		rotation[2] = 120.0f;
		zoom = 0.0f;

		//small generate bar side
		rotation[0] = 115.0f;
		rotation[1] = 265.0f;
		rotation[2] = 129.0f;
		zoom = 0.076f;
	}

	int getCurrentFrame()
	{
		return solverSettings.currentFrame;
	}

	void increaseCurrentFrame()
	{
		++solverSettings.currentFrame;
	}
};