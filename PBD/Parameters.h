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

	void initialiseToDefaults()
	{
		numMillisecondsToWaitBetweenFrames = 0;
		executionTimeSum = 0.0f;


		timingPrintInterval = 100;
		currentFrame = 1;
		maxFrames = 10000;

		useFEMSolver = false;
		writeToAlembic = true;
		printStrainEnergyToFile = false;

		testingInversionHandling = false;
		dimToCollapse = 1;

		generateMeshInsteadOfDoingIO = true;
	}
};