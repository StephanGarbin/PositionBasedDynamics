#pragma once


struct Parameters
{
	int numMillisecondsToWaitBetweenFrames;

	int timingPrintInterval;
	int currentFrame;
	int maxFrames;

	int globalHeight;
	int globalWidth;

	float executionTimeSum;

	float baryCentre[3];
	float radius;
	float rotation[3];
	float zoom;

	bool useFEMSolver;
	bool writeToAlembic;
	bool printStrainEnergyToFile;

	bool testingInversionHandling ;
	int dimToCollapse;

	bool generateMeshInsteadOfDoingIO;

	bool generateMeshFromTrackingData;

	bool readVertexConstraintData;

	bool disableSolver;

	void initialiseToDefaults()
	{
		disableSolver = false;

		numMillisecondsToWaitBetweenFrames = 0;
		executionTimeSum = 0.0f;


		timingPrintInterval = 100;
		currentFrame = 1;
		maxFrames = 100000;

		useFEMSolver = false;
		writeToAlembic = true;
		printStrainEnergyToFile = false;

		testingInversionHandling = false;
		dimToCollapse = 1;

		generateMeshInsteadOfDoingIO = true;
		generateMeshFromTrackingData = true;

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
		//rotation[0] = 115.0f;
		//rotation[1] = 265.0f;
		//rotation[2] = 129.0f;
		//zoom = 0.076f;
	}
};