#pragma once

#include "PBDSolverSettings.h"

struct Parameters
{
	//Global Settings
	bool disableSolver;
	int TEST_IDX;
	int TEST_VERSION;

	//Convenience functions
	int frame2ApplyInitialDeformation;
	int frame2DisApplyInitialDeformation;

	int continuousDeformationRelaxationFrame;
	float continuousDeformationStrainIncreaseFactor;

	bool invertSingleElementAtStart;
	float invertSingleElementAtStartAmount;

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
	bool collapseMeshAtStart;
	int dimToCollapse;

	//Mesh generation
	bool generateMeshInsteadOfDoingIO;
	bool generateMeshFromTrackingData;
	bool useTrackingConstraints;

	//Constraint IO
	bool readVertexConstraintData;

	bool applyInitialDeformationToMesh;
	bool applyContinuousDeformationToMesh;

	//Settings for the solver
	PBDSolverSettings solverSettings;

	//Fiber Mesh
	bool createFiberMesh;

	//Image IO
	bool doImageIO;
	std::string imageFileName;

	//Collision Geometry;
	bool readCollisionGeometry;
	std::vector<std::string> collisionGeometryFiles;

	bool translateCollisionGeometry;
	Eigen::Vector3f collisionGeometryTranslationAmount;
	int collisionGeometryTranslateUntilFrame;


	bool applyPressure;
	Eigen::Vector3f pressureForce;
	Eigen::Vector3f pressureCentre;
	float pressureRadius;
	int pressureMaxPositionIdx;
	int pressureStartFrame;
	int pressureEndFrame;


	bool renderCollisionGoemetry;

	bool overrideLookAt;
	Eigen::Vector3f lookatEye;

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
		createFiberMesh = false;

		collapseMeshAtStart = false;
		dimToCollapse = 1;

		generateMeshInsteadOfDoingIO = true;
		generateMeshFromTrackingData = false;

		readVertexConstraintData = false;

		doImageIO = false;
		readCollisionGeometry = false;
		translateCollisionGeometry = false;
		applyPressure = false;

		renderCollisionGoemetry = false;

		overrideLookAt = false;
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