#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

#include "Parameters.h"
#include "IOParameters.h"
#include "PBDSolverSettings.h"

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "CustomTetAttributeIO.h"

#include "CollisionRod.h"
#include "MovingHardConstraints.h"


bool parseTestParams(char* argument,
	int& testIdx, int& testVersion)
{
	std::string arg(argument);

	int firstUnderScore = arg.find_first_of("_");
	int secondUnderScore = arg.find_last_of("_");

	if (firstUnderScore == arg.size() || secondUnderScore == arg.size())
	{
		std::cout << "ERROR: Test is specified in an invalid way!. (Valid Example: [ TEST_1_0 ])" << std::endl;
		return false;
	}

	testIdx = std::stoi(arg.substr(firstUnderScore + 1, secondUnderScore - firstUnderScore));

	testVersion = std::stoi(arg.substr(secondUnderScore + 1, arg.size() - firstUnderScore));

	std::cout << "Running Test IDX [ " << testIdx << " ] - VERSION [ " << testVersion << "] " << std::endl;

	return true;
}

bool parseGeneralParams(Parameters& params, const int argc, char* argv[])
{
	//ELASTICITY
	params.solverSettings.youngsModulus = std::stof(argv[2]);
	params.solverSettings.poissonRatio = std::stof(argv[3]);

	//VISCO_ELASTICITY
	params.solverSettings.alpha = std::stof(argv[4]);
	params.solverSettings.rho = std::stof(argv[5]);

	//GENERAL
	params.solverSettings.inverseMass = std::stof(argv[5]);

	//SOLVER
	params.solverSettings.numConstraintIts = std::stoi(argv[6]);
	params.solverSettings.deltaT = std::stof(argv[7]);
	params.useFEMSolver = (std::string(argv[8]) == "USE_FEM");
	
	//IO
	params.writeToAlembic = (std::string(argv[9]) == "SAVE_MESH");
}

void initTest_0(Parameters& params, IOParameters& paramsIO);
void initTest_1(Parameters& params, IOParameters& paramsIO);
void initTest_2(Parameters& params, IOParameters& paramsIO);
void initTest_3(Parameters& params, IOParameters& paramsIO);

void initTest_4(Parameters& params, IOParameters& paramsIO);
void initTest_5(Parameters& params, IOParameters& paramsIO);

void initTest_6(Parameters& params, IOParameters& paramsIO);

void initTest_7(Parameters& params, IOParameters& paramsIO);

void initTest_8(Parameters& params, IOParameters& paramsIO);

//Average DeltaX Length Measure
void initTest_9(Parameters& params, IOParameters& paramsIO);

//Suction Test
void initTest_10(Parameters& params, IOParameters& paramsIO);

void initTest_11(Parameters& params, IOParameters& paramsIO);

//Sphere Indentation Test
void initTest_12(Parameters& params, IOParameters& paramsIO);

//Resolution
void initTest_13(Parameters& params, IOParameters& paramsIO);

//NAZIM PHANTOM
void initTest_14(Parameters& params, IOParameters& paramsIO);

//Torus Collision Test
void initTest_19(Parameters& params, IOParameters& paramsIO);

//INVERSION 2
void initTest_20(Parameters& params, IOParameters& paramsIO);

//Anisotropic
void initTest_21(Parameters& params, IOParameters& paramsIO);

bool parseTerminalParameters(const int argc, char* argv[],
	Parameters& params, IOParameters& paramsIO)
{
	//Check that we have enough arguments
	if (argc < 2)
	{
		std::cout << "Specify TEST_<IDX>_<VERSION> (where IDX is a placeholder for the test to run," << std::endl;
		std::cout << " and VERSION a placeholder for the sub-category of that test)." << std::endl;
		std::cout << "Additionally, optionally provide the following: " << std::endl;
		std::cout << "ELASTICITY--------------------" << std::endl;
		std::cout << "	- Youngs Modulus" << std::endl;
		std::cout << "	- Poisson Ratio" << std::endl;
		std::cout << "VISCO-ELASTICITY--------------------" << std::endl;
		std::cout << "	- Alpha" << std::endl;
		std::cout << "	- Rho" << std::endl;
		std::cout << "GENERAL--------------------" << std::endl;
		std::cout << "	- Inverse Mass" << std::endl;
		std::cout << "SOLVER--------------------" << std::endl;
		std::cout << "	- Num Constraint Its" << std::endl;
		std::cout << "	- Time Step Size" << std::endl;
		std::cout << "	- USE_FEM" << std::endl;
		std::cout << "IO--------------------" << std::endl;
		std::cout << "	- SAVE_MESH" << std::endl;
		return false;
	}

	//Defaults
	params.initialiseToDefaults();
	paramsIO.initialiseToDefaults();


	int testIDX;
	int testVERSION;
	//Check which test to run
	if (!parseTestParams(argv[1], testIDX, testVERSION))
	{
		return false;
	}

	params.TEST_IDX = testIDX;
	params.TEST_VERSION = testVERSION;

	//Initialise Test Parameters
	switch (testIDX)
	{
	case 0:
		initTest_0(params, paramsIO);
		break;
	case 1:
		initTest_1(params, paramsIO);
		break;
	case 2:
		initTest_2(params, paramsIO);
		break;
	case 3:
		initTest_3(params, paramsIO);
		break;
	case 4:
		initTest_4(params, paramsIO);
		break;
	case 5:
		initTest_5(params, paramsIO);
		break;
	case 6:
		initTest_6(params, paramsIO);
		break;
	case 7:
		initTest_7(params, paramsIO);
		break;
	case 8:
		initTest_8(params, paramsIO);
		break;
	case 9:
		initTest_9(params, paramsIO);
		break;
	case 10:
		initTest_10(params, paramsIO);
		break;
	case 11:
		initTest_11(params, paramsIO);
		break;
	case 12:
		initTest_12(params, paramsIO);
		break;
	case 13:
		initTest_13(params, paramsIO);
		break;
	case 14:
		initTest_14(params, paramsIO);
		break;
	case 19:
		initTest_19(params, paramsIO);
		break;
	case 20:
		initTest_20(params, paramsIO);
		break;
	case 21:
		initTest_21(params, paramsIO);
		break;
	default:
		break;
	}

	//Check if we have to also parse some other parameters
	if (argc > 2)
	{
		parseGeneralParams(params, argc, argv);
	}

	params.solverSettings.tracker.generateFileNames(params.TEST_IDX, params.TEST_VERSION);

	return true;
}



bool doIO(Parameters& params, IOParameters& paramsIO, std::vector<int>& vertexConstraintIndices,
	std::vector<PBDTetrahedra3d>& tetrahedra,
	std::shared_ptr<std::vector<PBDParticle>>& particles,
	std::vector<std::vector<Eigen::Vector2f>>& trackingData,
	std::vector<CollisionRod>& collisionRodGeometry,
	std::vector<CollisionSphere>& collisionSphereGeometry,
	std::vector<MovingHardConstraints>& movingHardConstraints)
{
	if (params.TEST_IDX == 0)
	{
		MeshCreator::generateTetBar(particles, tetrahedra, 10, 4, 4);
		return true;
	}
	else if (params.TEST_IDX == 1 || params.TEST_IDX == 2 || params.TEST_IDX == 3 || params.TEST_IDX == 4)
	{
		MeshCreator::generateSingleTet(particles, tetrahedra, 0, 0, 0);

		if (params.TEST_IDX == 3)
		{
			tetrahedra[0].setFullUpsilonCount(params.solverSettings.fullAlpha.size());
		}

		return true;
	}
	else if (params.TEST_IDX == 5)
	{
		MeshCreator::generateTetBar(particles, tetrahedra, 5, 3, 3);
		return true;
	}
	else if (params.TEST_IDX == 6)
	{
		MeshCreator::generateTetBar(particles, tetrahedra, 10, 6, 6);
		return true;
	}
	else if (params.TEST_IDX == 7 || params.TEST_IDX == 8)
	{
		MeshCreator::generateTetBar(particles, tetrahedra, 5, 3, 3);
		return true;
	}
	else if (params.TEST_IDX == 9)
	{
		MeshCreator::generateTetBar(particles, tetrahedra, 10, 4, 4);
		return true;
	}
	else if (params.TEST_IDX == 10)
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.setZero();
		TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);
	}
	else if (params.TEST_IDX == 11)
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.setZero();
		TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);

		if (params.TEST_VERSION != 5)
		{
			collisionRodGeometry.resize(1);
			std::vector<std::string> rodTransformNames;
			rodTransformNames.push_back("top1_bakedToWorld");
			rodTransformNames.push_back("bottom_bakedToWorld");

			collisionRodGeometry[0].readFromAbc("rod_liver1.abc", rodTransformNames);
		}

		if (paramsIO.constraintFile != "DUMMY")
		{
			ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, paramsIO.constraintFile);
			for (int i = 0; i < vertexConstraintIndices.size(); ++i)
			{
				(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
			}
		}
		else
		{
			std::cout << "WARNING: Constraint file ignored!" << std::endl;
		}
	}
	else if (params.TEST_IDX == 12)
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.setZero();
		TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);

		//pin ground vertices
		for (int p = 0; p < particles->size(); ++p)
		{
			if((*particles)[p].position()[1] <= 110.273f)
			{
				(*particles)[p].inverseMass() = 0.0;
			}
		}

		//for (int p = 0; p < particles->size(); ++p)
		//{
		//	std::cout << p << ": " << std::endl;
		//	std::cout << (*particles)[p].position() << std::endl;
		//}

		if (paramsIO.constraintFile != "DUMMY")
		{
			ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, paramsIO.constraintFile);
			for (int i = 0; i < vertexConstraintIndices.size(); ++i)
			{
				(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
			}
		}
		else
		{
			std::cout << "WARNING: Constraint file ignored!" << std::endl;
		}

		collisionSphereGeometry.resize(1);
		collisionSphereGeometry[0].readFromAbc("phantomDeformationSphere.abc","deforminationSphere");
		collisionSphereGeometry[0].getFrameLimit() = 28;
		return true;
		
	}
	else if (params.TEST_IDX == 13)
	{
		if (params.TEST_VERSION == 0)
		{
			MeshCreator::generateTetBar(particles, tetrahedra, 10, 4, 4);
		}
		else if (params.TEST_VERSION == 1)
		{
			MeshCreator::generateTetBar(particles, tetrahedra, 10, 8, 8);
		}
		else if (params.TEST_VERSION == 2)
		{
			MeshCreator::generateTetBar(particles, tetrahedra, 20, 8, 8);
		}
		else if (params.TEST_VERSION == 3)
		{
			MeshCreator::generateTetBar(particles, tetrahedra, 20, 16, 16);
		}
		else if (params.TEST_VERSION == 4)
		{
			MeshCreator::generateTetBar(particles, tetrahedra, 20, 32, 32);
		}
		else if (params.TEST_VERSION == 5)
		{
			MeshCreator::generateTetBar(particles, tetrahedra, 12, 11, 12);
		}
		return true;
	}
	else if (params.TEST_IDX == 14)
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.setZero();
		TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);

		ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, paramsIO.constraintFile);
		for (int i = 0; i < vertexConstraintIndices.size(); ++i)
		{
			(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
		}

		movingHardConstraints.resize(1);
		movingHardConstraints[0].readConstraintFiles(paramsIO.movingHardConstraintsIndexFiles);
		movingHardConstraints[0].initialisePositionMasses(*particles);
		movingHardConstraints[0].readFromAbc(paramsIO.movingHardConstraintsAbcArchive,
			paramsIO.movingHardConstraintsLocatorNames);
	}
	else if (params.TEST_IDX == 19)
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.setZero();
		TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);

		if (paramsIO.constraintFile != "DUMMY")
		{
			ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, paramsIO.constraintFile);
			for (int i = 0; i < vertexConstraintIndices.size(); ++i)
			{
				(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
			}
		}
		else
		{
			std::cout << "WARNING: Constraint file ignored!" << std::endl;
		}

		collisionSphereGeometry.resize(5);
		collisionSphereGeometry[0].readFromAbc("zeroZeroCollisionSphere.abc", "collisionSphere");
		collisionSphereGeometry[0].getFrameLimit() = 99;
		collisionSphereGeometry[1].readFromAbc("zeroZeroCollisionSphere1.abc", "collisionSphere1");
		collisionSphereGeometry[1].getFrameLimit() = 99;
		collisionSphereGeometry[2].readFromAbc("zeroZeroCollisionSphere2.abc", "collisionSphere2");
		collisionSphereGeometry[2].getFrameLimit() = 99;
		collisionSphereGeometry[3].readFromAbc("zeroZeroCollisionSphere3.abc", "collisionSphere3");
		collisionSphereGeometry[3].getFrameLimit() = 99;
		collisionSphereGeometry[4].readFromAbc("zeroZeroCollisionSphere4.abc", "collisionSphere4");
		collisionSphereGeometry[4].getFrameLimit() = 99;
		return true;

	}
	else if (params.TEST_IDX == 20)
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.setZero();
		TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);
	}
	else if (params.TEST_IDX == 21)
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.setZero();
		TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);

		ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, paramsIO.constraintFile);
		for (int i = 0; i < vertexConstraintIndices.size(); ++i)
		{
			(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
		}

		//Now read custom tet attributes
		std::vector<float> cYoungsModulus;
		std::vector<float> cAnisotropyStrength;
		std::vector<Eigen::Vector3f> cAnisotropyDirection;
		readCustomTetAttributes(cYoungsModulus, cAnisotropyStrength, cAnisotropyDirection, paramsIO.customTetAttributeFile);
		for (int t = 0; t < tetrahedra.size(); ++t)
		{
			if (params.TEST_VERSION > 0)
			{
				tetrahedra[t].getPerTetYoungsModulus() = cYoungsModulus[t];
			}
			if (params.TEST_VERSION > 1)
			{
				tetrahedra[t].getPerTetAnisotropyStrength() = cAnisotropyStrength[t];
				tetrahedra[t].getPerTetAnisotropyDirection() = cAnisotropyDirection[t];
			}
		}

		return true;
	}
	else
	{
		Eigen::Vector3f initialVelocity;
		initialVelocity.x() = 0; initialVelocity.y() = 0; initialVelocity.z() = 0;

		// HARD VERTEX CONSTRAINTS
		if (params.readVertexConstraintData)
		{
			ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, paramsIO.constraintFile);
			for (int i = 0; i < vertexConstraintIndices.size(); ++i)
			{
				(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
			}
		}

		// TRACKING DATA
		trackingData.resize(paramsIO.trackerFiles.size());
		for (int i = 0; i < paramsIO.trackerFiles.size(); ++i)
		{
			TrackerIO::readTrackerAnimationNuke(paramsIO.trackerFiles[i], trackingData[i]);
		}

		if (!params.generateMeshInsteadOfDoingIO)
		{
			// MESH
			TetGenIO::readNodes(paramsIO.nodeFile, *particles, params.solverSettings.inverseMass, initialVelocity);
			TetGenIO::readTetrahedra(paramsIO.elementFile, tetrahedra, particles);
		}
		else
		{
			//GENERATE the mesh data
			if (params.generateMeshFromTrackingData)
			{
				MeshCreator::generateTetBarToFit(particles, tetrahedra, 13, 7, 3,
					trackingData[0][0], trackingData[1][0], trackingData[2][0], 20.0f);
			}
			else
			{
				MeshCreator::generateTetBar(particles, tetrahedra, 10, 6, 6);
			}
		}

		return true;
	}
}

//VISCOELASTICITY TESTS
void initTest_0(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 100000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.disableSolver = false;

	params.solverSettings.poissonRatio = 0.4f;
	params.solverSettings.youngsModulus = 10.0f;
	params.solverSettings.numConstraintIts = 10;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.81f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 0.1f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 1.0f);
	//params.solverSettings.anisotropyParameter = 1.0f;

	params.solverSettings.materialModel = PBDSolverSettings::NEO_HOOKEAN_FIBER;

	params.solverSettings.disableInversionHandling = false;
	params.zoom = 0.328f;

	params.solverSettings.useSecondOrderUpdates = false;
	params.solverSettings.useMultiThreadedSolver = true;

	params.solverSettings.enableGroundPlaneCollision = true;
	params.solverSettings.groundplaneHeight = -1.0f;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;

	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.25f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 1.0f;
	}
}

void initTest_1(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 1000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	params.zoom = 1.318f;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 1.0f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = 0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;


	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	params.applyInitialDeformationToMesh = true;
	params.frame2ApplyInitialDeformation = 100;
	params.frame2DisApplyInitialDeformation = 500;

	params.solverSettings.numConstraintIts = 1;
	params.numMillisecondsToWaitBetweenFrames = 0;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.trackS = true;
	params.solverSettings.trackF = true;
	params.solverSettings.trackPF = true;

	params.solverSettings.forceMultiplicationFactor = 1.0f;
	params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.rho = 0.0f;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;

	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.25f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 1.0f;
	}
}

void initTest_2(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 1500;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	params.zoom = 1.318f;

	params.solverSettings.poissonRatio = 0.499f;
	params.solverSettings.youngsModulus = 200.001f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = 0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;


	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	params.applyInitialDeformationToMesh = false;
	params.applyContinuousDeformationToMesh = true;
	params.continuousDeformationRelaxationFrame = 500;
	params.continuousDeformationStrainIncreaseFactor = 0.001f;
	params.frame2ApplyInitialDeformation = -1;

	params.solverSettings.numConstraintIts = 1;
	params.numMillisecondsToWaitBetweenFrames = 0;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.trackS = true;
	params.solverSettings.trackF = true;
	params.solverSettings.trackPF = true;

	params.solverSettings.forceMultiplicationFactor = 0.0f;

	params.solverSettings.rho = 0.99f;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;

	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.25f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 1.0f;
	}
}

void initTest_3(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 1000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	params.zoom = 1.318f;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 1.0f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = 0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;


	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	params.applyInitialDeformationToMesh = true;
	params.frame2ApplyInitialDeformation = 100;
	params.frame2DisApplyInitialDeformation = 500;

	params.solverSettings.numConstraintIts = 1;
	params.numMillisecondsToWaitBetweenFrames = 0;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.trackS = true;
	params.solverSettings.trackF = true;
	params.solverSettings.trackPF = true;

	params.solverSettings.forceMultiplicationFactor = 1.0f;
	params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.useFullPronySeries = true;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.fullAlpha = { 0.60f, 0.2f, 0.2f };
		params.solverSettings.fullRho = { 0.019f, 0.9f, 0.001f };
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.fullAlpha = { 0.30f, 0.3f, 0.4f };
		params.solverSettings.fullRho = { 0.001f, 0.00001f, 0.9001f };
	}
	//else if (params.TEST_VERSION == 2)
	//{
	//	params.solverSettings.alpha = 0.5f;
	//}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

//INVERSION / DEGENERATE ELEMENT HANDLING
void initTest_4(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 1000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	params.zoom = 1.318f;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 1.0f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = 0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;


	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	//params.applyInitialDeformationToMesh = true;
	//params.frame2ApplyInitialDeformation = 100;
	//params.frame2DisApplyInitialDeformation = 500;

	params.solverSettings.numConstraintIts = 1;
	params.numMillisecondsToWaitBetweenFrames = 200;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	//params.solverSettings.trackS = true;
	//params.solverSettings.trackF = true;
	//params.solverSettings.trackPF = true;

	//params.solverSettings.forceMultiplicationFactor = 1.0f;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.useFullPronySeries = false;
	params.solverSettings.alpha = 0.0f;
	params.solverSettings.rho = 0.0f;

	params.invertSingleElementAtStart = true;

	if (params.TEST_VERSION == 0)
	{
		params.invertSingleElementAtStartAmount = 2.000001f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.invertSingleElementAtStartAmount = 2.000001f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.invertSingleElementAtStartAmount = 2.000001f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.numMillisecondsToWaitBetweenFrames = 500;
	}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

void initTest_5(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 1000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	//params.zoom = 1.318f;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 10.0f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = 0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.zoom = 1.5f;

	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	//params.applyInitialDeformationToMesh = true;
	//params.frame2ApplyInitialDeformation = 100;
	//params.frame2DisApplyInitialDeformation = 500;

	params.solverSettings.numConstraintIts = 1;
	params.numMillisecondsToWaitBetweenFrames = 50;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	//params.solverSettings.trackS = true;
	//params.solverSettings.trackF = true;
	//params.solverSettings.trackPF = true;

	//params.solverSettings.forceMultiplicationFactor = 1.0f;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.useFullPronySeries = false;
	params.solverSettings.alpha = 0.0f;
	params.solverSettings.rho = 0.0f;

	params.doImageIO = true;

	//params.invertSingleElementAtStart = true;

	if (params.TEST_VERSION == 0)
	{
		params.collapseMeshAtStart = true;
		params.dimToCollapse = 0;
		params.imageFileName = "inversionLowRes_dim0";

	}
	else if (params.TEST_VERSION == 1)
	{
		params.collapseMeshAtStart = true;
		params.dimToCollapse = 1;
		params.imageFileName = "inversionLowRes_dim1";
	}
	else if (params.TEST_VERSION == 2)
	{
		params.collapseMeshAtStart = true;
		params.dimToCollapse = 2;
		params.imageFileName = "inversionLowRes_dim2";
	}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

void initTest_6(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 100000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 10.0f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.81f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.zoom = 0.5f;


	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	//params.applyInitialDeformationToMesh = true;
	//params.frame2ApplyInitialDeformation = 100;
	//params.frame2DisApplyInitialDeformation = 500;

	params.solverSettings.numConstraintIts = 25;
	params.numMillisecondsToWaitBetweenFrames = 0;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	//params.solverSettings.trackS = true;
	//params.solverSettings.trackF = true;
	//params.solverSettings.trackPF = true;

	//params.solverSettings.forceMultiplicationFactor = 1.0f;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.useFullPronySeries = false;
	params.solverSettings.alpha = 0.0f;
	params.solverSettings.rho = 0.0f;

	//params.invertSingleElementAtStart = true;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;

		//params.solverSettings.disablePositionCorrection = true;
		params.solverSettings.MR_f_active = 1.0f;
		params.solverSettings.MR_f_passive = 1.0f;
		params.solverSettings.MR_alpha = 0.0f;
		params.solverSettings.MR_A = 3.0f;
		params.solverSettings.MR_B = 1.0f;
		params.solverSettings.MR_J = 1.0f;
		params.solverSettings.MR_K = 6.0f;
		params.solverSettings.MR_T = 1.0f;
		params.solverSettings.MR_a = Eigen::Vector3f(1.0f, 0.0f, 0.0f).normalized();

		params.solverSettings.anisotropyParameter = 1.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;

		//params.solverSettings.disablePositionCorrection = true;
		params.solverSettings.MR_f_active = 1.0f;
		params.solverSettings.MR_f_passive = 1.0f;
		params.solverSettings.MR_alpha = 0.0f;
		params.solverSettings.MR_A = 3.0f;
		params.solverSettings.MR_B = 1.0f;
		params.solverSettings.MR_J = 1.0f;
		params.solverSettings.MR_K = 6.0f;
		params.solverSettings.MR_T = 1.0f;
		params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 0.0f).normalized();

		params.solverSettings.MR_A0.setZero();
		params.solverSettings.MR_A0(1, 1) = 1.0f;

		params.solverSettings.anisotropyParameter = 1.0f;
	}
	//else if (params.TEST_VERSION == 2)
	//{
	//	params.solverSettings.alpha = 0.5f;
	//}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

void initTest_7(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 100000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 10.0f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.81f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.zoom = 0.5f;


	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	//params.applyInitialDeformationToMesh = true;
	//params.frame2ApplyInitialDeformation = 100;
	//params.frame2DisApplyInitialDeformation = 500;

	params.solverSettings.numConstraintIts = 25;
	params.numMillisecondsToWaitBetweenFrames = 0;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	//params.solverSettings.trackS = true;
	//params.solverSettings.trackF = true;
	//params.solverSettings.trackPF = true;

	//params.solverSettings.forceMultiplicationFactor = 1.0f;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.useFullPronySeries = false;
	params.solverSettings.alpha = 0.0f;
	params.solverSettings.rho = 0.0f;

	//params.invertSingleElementAtStart = true;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;

		//params.solverSettings.disablePositionCorrection = true;
		params.solverSettings.MR_f_active = 1.0f;
		params.solverSettings.MR_f_passive = 1.0f;
		params.solverSettings.MR_alpha = 0.0f;
		params.solverSettings.MR_A = 3.0f;
		params.solverSettings.MR_B = 1.0f;
		params.solverSettings.MR_J = 1.0f;
		params.solverSettings.MR_K = 6.0f;
		params.solverSettings.MR_T = 1.0f;
		params.solverSettings.MR_a = Eigen::Vector3f(1.0f, 0.0f, 0.0f).normalized();

		params.solverSettings.anisotropyParameter = 1.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;

		//params.solverSettings.disablePositionCorrection = true;
		params.solverSettings.MR_f_active = 1.0f;
		params.solverSettings.MR_f_passive = 1.0f;
		params.solverSettings.MR_alpha = 0.0f;
		params.solverSettings.MR_A = 3.0f;
		params.solverSettings.MR_B = 1.0f;
		params.solverSettings.MR_J = 1.0f;
		params.solverSettings.MR_K = 6.0f;
		params.solverSettings.MR_T = 1.0f;
		params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 0.0f).normalized();

		params.solverSettings.MR_A0.setZero();
		params.solverSettings.MR_A0(1, 1) = 1.0f;

		params.solverSettings.anisotropyParameter = 1.0f;
	}
	//else if (params.TEST_VERSION == 2)
	//{
	//	params.solverSettings.alpha = 0.5f;
	//}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

void initTest_8(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 100000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 10.0f;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.81f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.zoom = 0.5f;


	params.disableSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.solverSettings.disablePositionCorrection = false;
	//params.applyInitialDeformationToMesh = true;
	//params.frame2ApplyInitialDeformation = 100;
	//params.frame2DisApplyInitialDeformation = 500;

	params.solverSettings.numConstraintIts = 25;
	params.numMillisecondsToWaitBetweenFrames = 0;
	//params.solverSettings.forceMultiplicationFactor = 10000.2;
	//params.solverSettings.externalForce.x() = 1.0f;
	//params.solverSettings.trackS = true;
	//params.solverSettings.trackF = true;
	//params.solverSettings.trackPF = true;

	//params.solverSettings.forceMultiplicationFactor = 1.0f;
	//params.solverSettings.externalForce.x() = 1.0f;
	params.solverSettings.useFullPronySeries = false;
	params.solverSettings.alpha = 0.0f;
	params.solverSettings.rho = 0.0f;

	//params.invertSingleElementAtStart = true;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;

		//params.solverSettings.disablePositionCorrection = true;
		params.solverSettings.MR_f_active = 1.0f;
		params.solverSettings.MR_f_passive = 1.0f;
		params.solverSettings.MR_alpha = 0.0f;
		params.solverSettings.MR_A = 3.0f;
		params.solverSettings.MR_B = 1.0f;
		params.solverSettings.MR_J = 1.0f;
		params.solverSettings.MR_K = 6.0f;
		params.solverSettings.MR_T = 1.0f;
		params.solverSettings.MR_a = Eigen::Vector3f(1.0f, 0.0f, 0.0f).normalized();

		params.solverSettings.anisotropyParameter = 1.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.numConstraintIts = 25;
		params.maxFrames = 10000;
		params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::RUBIN_BODNER;

		//params.solverSettings.disablePositionCorrection = true;

	}
	//else if (params.TEST_VERSION == 2)
	//{
	//	params.solverSettings.alpha = 0.5f;
	//}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

void initTest_9(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 50000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.disableSolver = false;

	params.solverSettings.poissonRatio = 0.3f;
	params.solverSettings.youngsModulus = 5.0f;
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.81f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 1.0f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.solverSettings.trackAverageDeltaXLength = true;

	params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;

	params.zoom = 0.5f;
	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 0.0f);

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;

	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.25f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 1.0f;
	}
}

void initTest_10(Parameters& params, IOParameters& paramsIO)
{
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.disableSolver = false;

	params.solverSettings.poissonRatio = 0.45f;
	params.solverSettings.youngsModulus = 550.0f;
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = 0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 1.0f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	paramsIO.nodeFile = "tissueCube.1.node";
	paramsIO.elementFile = "tissueCube.1.ele";

	params.solverSettings.enableGroundPlaneCollision = true;
	params.readCollisionGeometry = true;
	params.collisionGeometryFiles.push_back("suctionCollisionTestMesh.abc");

	params.translateCollisionGeometry = true;
	params.collisionGeometryTranslationAmount = Eigen::Vector3f(0.0f, -0.001f, 0.0f);
	params.collisionGeometryTranslateUntilFrame = 40;

	params.applyPressure = true;
	params.pressureStartFrame = 60;
	params.pressureEndFrame = 360;
	params.pressureCentre = Eigen::Vector3f(0.995f, 0.508f, -0.993f);
	params.pressureRadius = 0.5f;
	params.pressureMaxPositionIdx = 1729;
	params.pressureForce = Eigen::Vector3f(0.0f, 3.0f, 0.0f);

	params.writeToAlembic = true;
	params.maxFrames = 700;

	params.solverSettings.rho = 0.5f;

	params.solverSettings.trackSpecificPosition = true;
	params.solverSettings.trackSpecificPositionIdx = 462;

	params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN;
	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 0.0f);

	params.solverSettings.useMultiThreadedSolver = true;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.25f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 0.99f;
	}
}

void initTest_11(Parameters& params, IOParameters& paramsIO)
{
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.renderCollisionGoemetry = true;

	params.solverSettings.poissonRatio = 0.40f;
	params.solverSettings.youngsModulus = 0.01f;
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 1.0f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	paramsIO.nodeFile = "liverPig1.node";
	paramsIO.elementFile = "liverPig1.ele";
	//paramsIO.constraintFile = "DUMMY";
	paramsIO.constraintFile = "pigLiver1AbsoluteConstraints.txt";

	//params.solverSettings.enableGroundPlaneCollision = true;
	//params.solverSettings.groundplaneHeight = -50.0f;

	//params.writeToAlembic = true;
	params.maxFrames = 1000;

	params.solverSettings.rho = 0.5f;

	params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;
	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
	params.solverSettings.anisotropyParameter = 0.0f;


	//rod
	params.solverSettings.collisionSpheresNum.push_back(4);
	params.solverSettings.collisionSpheresRadius.push_back(2.0f);

	params.writeToAlembic = false;
	params.disableSolver = false;
	params.solverSettings.useMultiThreadedSolver = true;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.25f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 0.99f;
	}
}

void initTest_12(Parameters& params, IOParameters& paramsIO)
{
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.renderCollisionGoemetry = true;

	params.solverSettings.poissonRatio = 0.40f;
	params.solverSettings.youngsModulus = 0.01f; //0.02f
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -0.00f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 0.58f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	paramsIO.nodeFile = "LiverInitialLowResolution_lowest.1.node";
	paramsIO.elementFile = "LiverInitialLowResolution_lowest.1.ele";
	paramsIO.constraintFile = "phantomIdentationGood_Lowest.txt";
	//paramsIO.nodeFile = "LiverInitialLowResolution_00625.1.node";
	//paramsIO.elementFile = "LiverInitialLowResolution_00625.1.ele";
	//paramsIO.constraintFile = "DUMMY";
	//paramsIO.constraintFile = "phantomIdentationGood_00625.txt";
	//paramsIO.constraintFile = "phantomIdentationLow_00625.txt";

	//params.solverSettings.enableGroundPlaneCollision = true;
	//params.solverSettings.groundplaneHeight = -50.0f;

	params.writeToAlembic = true;
	params.maxFrames = 6000;

	params.solverSettings.rho = 0.05f;

	params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;
	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
	params.solverSettings.anisotropyParameter = 0.0f;
	params.solverSettings.disableConstraintProjection = false;

	//rod
	params.solverSettings.collisionSpheresRadius.push_back(41.129f);

	params.writeToAlembic = true;
	params.disableSolver = false;
	params.solverSettings.useMultiThreadedSolver = true;

	//params.solverSettings.useSecondOrderUpdates = true;

	params.solverSettings.enableGroundPlaneCollision = false;
	params.solverSettings.groundplaneHeight = 103.41f;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.0f;
		paramsIO.nodeFile = "LiverInitialLowResolution_00625.1.node";
		paramsIO.elementFile = "LiverInitialLowResolution_00625.1.ele";
		paramsIO.constraintFile = "phantomIdentationLow_00625.txt";
	}
	else if (params.TEST_VERSION == 2)
	{
		paramsIO.nodeFile = "LiverInitialLowResolution_lowest2.1.node";
		paramsIO.elementFile = "LiverInitialLowResolution_lowest2.1.ele";
		paramsIO.constraintFile = "phantomIdentationGood_Lowest2.txt";
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 3)
	{
		paramsIO.nodeFile = "LiverInitialLowResolution_lowest3.1.node";
		paramsIO.elementFile = "LiverInitialLowResolution_lowest3.1.ele";
		paramsIO.constraintFile = "phantomIdentationGood_Lowest3.txt";
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 0.99f;
	}
	else if (params.TEST_VERSION == 5)
	{
		std::cout << "WARNING: This will produce an undeformed configuration for comparison!" << std::endl;
		params.solverSettings.disableConstraintProjection = true;
	}
}

void initTest_13(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 100000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.disableSolver = false;

	params.solverSettings.poissonRatio = 0.4f;
	params.solverSettings.youngsModulus = 10.0f;
	params.solverSettings.numConstraintIts = 10;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.81f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 0.1f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 1.0f);
	//params.solverSettings.anisotropyParameter = 1.0f;

	params.solverSettings.materialModel = PBDSolverSettings::NEO_HOOKEAN_FIBER;

	params.solverSettings.disableInversionHandling = false;
	params.zoom = 0.328f;

	params.solverSettings.useSecondOrderUpdates = false;
	params.solverSettings.useMultiThreadedSolver = true;

	params.solverSettings.enableGroundPlaneCollision = true;
	params.solverSettings.groundplaneHeight = -1.0f;
	params.solverSettings.alpha = 0.0f;
	//if (params.TEST_VERSION == 0)
	//{
	//	params.solverSettings.alpha = 0.0f;

	//}
	//else if (params.TEST_VERSION == 1)
	//{
	//	params.solverSettings.alpha = 0.25f;
	//}
	//else if (params.TEST_VERSION == 2)
	//{
	//	params.solverSettings.alpha = 0.5f;
	//}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

void initTest_14(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 1000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.disableSolver = false;

	params.solverSettings.poissonRatio = 0.4f;
	params.solverSettings.youngsModulus = 0.05f;
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -0.0f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 0.1f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 1.0f);
	//params.solverSettings.anisotropyParameter = 1.0f;

	params.solverSettings.materialModel = PBDSolverSettings::NEO_HOOKEAN_FIBER;

	params.solverSettings.disableInversionHandling = false;
	params.zoom = 0.328f;

	params.solverSettings.useSecondOrderUpdates = false;
	params.solverSettings.useMultiThreadedSolver = true;

	paramsIO.nodeFile = "parenchyma_reduce_1.node";
	paramsIO.elementFile = "parenchyma_reduce_1.ele";
	paramsIO.constraintFile = "config1Hard.txt";

	paramsIO.movingHardConstraintsAbcArchive = "config1_locatorAnim.abc";
	paramsIO.movingHardConstraintsIndexFiles.push_back("config1Pulley.txt");
	paramsIO.movingHardConstraintsLocatorNames.push_back("locatorPulley");
	paramsIO.movingHardConstraintsIndexFiles.push_back("config1Table.txt");
	paramsIO.movingHardConstraintsLocatorNames.push_back("tablePulley");

	params.overrideLookAt = true;
	params.lookatEye = Eigen::Vector3f(177.564f, 237.164f, -303.59f);

	params.writeToAlembic = true;
	params.solverSettings.disableConstraintProjection = false;
	//if (params.TEST_VERSION == 0)
	//{
	//	params.solverSettings.alpha = 0.0f;

	//}
	//else if (params.TEST_VERSION == 1)
	//{
	//	params.solverSettings.alpha = 0.25f;
	//}
	//else if (params.TEST_VERSION == 2)
	//{
	//	params.solverSettings.alpha = 0.5f;
	//}
	//else if (params.TEST_VERSION == 3)
	//{
	//	params.solverSettings.alpha = 0.75f;
	//}
	//else if (params.TEST_VERSION == 4)
	//{
	//	params.solverSettings.alpha = 1.0f;
	//}
}

void initTest_19(Parameters& params, IOParameters& paramsIO)
{
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.solverSettings.disableConstraintProjection = false;
	params.renderCollisionGoemetry = true;

	params.solverSettings.poissonRatio = 0.40f;
	params.solverSettings.youngsModulus = 5.1f; //0.02f
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.81f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 0.58f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.zoom = 1.43f;

	paramsIO.nodeFile = "torus.node";
	paramsIO.elementFile = "torus.ele";


	params.renderCollisionGoemetry = true;

	paramsIO.constraintFile = "DUMMY";

	params.writeToAlembic = true;
	params.maxFrames = 3000;

	params.solverSettings.rho = 0.58f;

	params.solverSettings.materialModel = PBDSolverSettings::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER;
	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
	params.solverSettings.anisotropyParameter = 0.0f;
	params.solverSettings.disableConstraintProjection = false;

	//rod
	params.solverSettings.collisionSpheresRadius.push_back(1.538f);
	params.solverSettings.collisionSpheresRadius.push_back(1.538f);
	params.solverSettings.collisionSpheresRadius.push_back(1.538f);
	params.solverSettings.collisionSpheresRadius.push_back(1.538f);
	params.solverSettings.collisionSpheresRadius.push_back(1.538f);

	params.writeToAlembic = true;
	params.disableSolver = false;
	params.solverSettings.useMultiThreadedSolver = true;

	params.solverSettings.enableGroundPlaneCollision = true;
	params.solverSettings.groundplaneHeight = 0.0f;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.99f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.0f;
		params.solverSettings.youngsModulus = 10.1f; //0.02f
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.0f;
		params.solverSettings.youngsModulus = 15.5f; //0.02f
	}
	else if (params.TEST_VERSION == 4)
	{

	}

}


void initTest_20(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 5000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.disableSolver = false;

	params.solverSettings.poissonRatio = 0.4f;
	params.solverSettings.youngsModulus = 0.05f;
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.8f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 0.1f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 1.0f);
	//params.solverSettings.anisotropyParameter = 1.0f;

	params.solverSettings.materialModel = PBDSolverSettings::NEO_HOOKEAN_FIBER;

	params.solverSettings.disableInversionHandling = false;
	params.zoom = 0.328f;

	params.solverSettings.useSecondOrderUpdates = false;
	params.solverSettings.useMultiThreadedSolver = true;

	paramsIO.nodeFile = "parenchyma_reduce_1.node";
	paramsIO.elementFile = "parenchyma_reduce_1.ele";

	params.solverSettings.enableGroundPlaneCollision = true;
	params.solverSettings.groundplaneHeight = 90.0f;

	params.overrideLookAt = true;
	params.lookatEye = Eigen::Vector3f(177.564f, 237.164f, -303.59f);

	params.writeToAlembic = true;
	params.solverSettings.disableConstraintProjection = false;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.alpha = 0.25f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 1.0f;
	}
}

void initTest_21(Parameters& params, IOParameters& paramsIO)
{
	params.maxFrames = 5000;
	params.writeToAlembic = false;
	params.useTrackingConstraints = false;
	params.readVertexConstraintData = false;
	params.useFEMSolver = false;
	params.disableSolver = false;

	params.solverSettings.poissonRatio = 0.4f;
	params.solverSettings.youngsModulus = 0.01f;
	params.solverSettings.numConstraintIts = 5;
	params.solverSettings.deltaT = 0.005f;
	params.solverSettings.inverseMass = 1.0f;
	params.solverSettings.printStrainEnergy = false;
	params.solverSettings.printStrainEnergyToFile = false;
	params.solverSettings.gravity = -9.8f;
	params.solverSettings.externalForce.setZero();
	params.solverSettings.forceMultiplicationFactor = 0.0f;
	params.solverSettings.rho = 0.1f;
	params.solverSettings.numTetrahedraIterations = 0;
	params.solverSettings.correctStrongForcesWithSubteps = false;
	params.solverSettings.useGeometricConstraintLimits = false;

	params.solverSettings.MR_a = Eigen::Vector3f(0.0f, 1.0f, 1.0f);
	//params.solverSettings.anisotropyParameter = 1.0f;

	params.solverSettings.materialModel = PBDSolverSettings::NEO_HOOKEAN_FIBER;

	params.solverSettings.disableInversionHandling = false;
	params.zoom = 0.328f;

	params.solverSettings.useSecondOrderUpdates = false;
	params.solverSettings.useMultiThreadedSolver = true;

	paramsIO.nodeFile = "ani/cube2.1.node";
	paramsIO.elementFile = "ani/cube2.1.ele";
	paramsIO.customTetAttributeFile = "ani/cube2CustomAttrs.txt";
	paramsIO.constraintFile = "ani/cube2HardConstraints.txt";

	params.solverSettings.enableGroundPlaneCollision = true;
	params.solverSettings.groundplaneHeight = 00.0f;

	params.overrideLookAt = true;
	params.lookatEye = Eigen::Vector3f(9.875f, 11.122f, -5.686f);

	params.writeToAlembic = true;
	params.solverSettings.disableConstraintProjection = false;

	if (params.TEST_VERSION == 0)
	{
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 1)
	{
		params.solverSettings.usePerTetMaterialAttributes = true;
		params.solverSettings.alpha = 0.0f;
	}
	else if (params.TEST_VERSION == 2)
	{
		params.solverSettings.usePerTetMaterialAttributes = true;
		params.solverSettings.alpha = 0.5f;
	}
	else if (params.TEST_VERSION == 3)
	{
		params.solverSettings.alpha = 0.75f;
	}
	else if (params.TEST_VERSION == 4)
	{
		params.solverSettings.alpha = 1.0f;
	}
}

