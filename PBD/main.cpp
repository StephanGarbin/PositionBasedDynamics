#include <iostream>
#include <chrono>
#include <thread>

#include <tbb\tick_count.h>
#include <sstream>

#include "TetGenIO.h"
#include "ConstraintsIO.h"
#include "VegaIO.h"
#include "TrackerIO.h"

#include "AbcWriter.h"
#include "SurfaceMeshHandler.h"

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "PBDSolver.h"
#include "PBDSolverSettings.h"
#include "PBDProbabilisticConstraint.h"

#include "GLUTHelper.h"
#include "AntTweakBar.h"

#include "FEMSimulator.h"
#include "MeshCreator.h"

#include "Parameters.h"
#include "IOParameters.h"

std::vector<PBDTetrahedra3d> tetrahedra;
std::shared_ptr<std::vector<PBDParticle>> particles = std::make_shared<std::vector<PBDParticle>>();
std::vector<Eigen::Vector3f> currentPositions;
std::vector<Eigen::Vector3f> initialPositions;
std::vector<int> numConstraintInfluences;
std::vector<PBDProbabilisticConstraint> probabilisticConstraints;
std::vector<std::vector<Eigen::Vector2f>> trackingData;

PBDSolver solver;
FEMSimulator FEMsolver;

std::shared_ptr<SurfaceMeshHandler> smHandler;

PBDSolverSettings settings;
Parameters parameters;
IOParameters ioParameters;

void mainLoop();


void applyFEMDisplacementsToParticles()
{
	std::vector<double>& displacements = FEMsolver.getCurrentDisplacements();

	for (int i = 0; i < displacements.size(); i += 3)
	{
		(*particles)[i / 3].position()[0] = initialPositions[i / 3][0] + displacements[i];
		(*particles)[i / 3].position()[1] = initialPositions[i / 3][1] + displacements[i + 1];
		(*particles)[i / 3].position()[2] = initialPositions[i / 3][2] + displacements[i + 2];

		//std::cout << (*particles)[i / 3].position() << std::endl;
		displacements[i] = 0;
		displacements[i + 1] = 0;
		displacements[i + 2] = 0;
	}
}

void setInitialPositionsFromParticles()
{
	initialPositions.resize(particles->size());
	for (int i = 0; i < particles->size(); ++i)
	{
		initialPositions[i] = (*particles)[i].position();
	}
}

void getCurrentPositionFromParticles()
{
	if (currentPositions.size() == 0)
	{
		currentPositions.resize(particles->size());
	}

	for (int i = 0; i < (*particles).size(); ++i)
	{
		currentPositions[i] = (*particles)[i].position();
	}
}

void updateProbabilisticConstraints()
{
	Eigen::Vector3f currentConstraintPosition = 
		TrackerIO::getInterpolatedConstraintPosition(trackingData[3], 1.0f / 24.0f, settings.deltaT, parameters.currentFrame * settings.deltaT);

	float scaleFactor = 0.01;

	//1. subtract displacement
	float widthSubtraction = trackingData[0][0].x();
	float heightSubtraction = trackingData[0][0].y();

	std::cout << widthSubtraction << ", " << heightSubtraction << std::endl;

	currentConstraintPosition[0] -= widthSubtraction;
	currentConstraintPosition[2] -= heightSubtraction;

	//2. Rescale
	currentConstraintPosition *= scaleFactor;

	probabilisticConstraints[0].getConstraintPosition() = currentConstraintPosition;
}

void createProabilisticConstraints()
{
	probabilisticConstraints.resize(1);
	updateProbabilisticConstraints();
	probabilisticConstraints[0].initialise(*particles, 1000.1);
}

void setCamera()
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity;

	gluPerspective( /* field of view in degree */ 90.0,
		/* aspect ratio */ 1.0,
		/* Z near */ 1.0, /* Z far */ 500000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity;

	gluLookAt(-(parameters.baryCentre[0] + parameters.radius * parameters.zoom),
		parameters.baryCentre[1] + parameters.radius * parameters.zoom,
		-(parameters.baryCentre[2] + parameters.radius * parameters.zoom),  /* eye is at (0,0,5) */
		parameters.baryCentre[0], parameters.baryCentre[1], parameters.baryCentre[2],      /* center is at (0,0,0) */
		0.0, 1.0, 0.0);      /* up is in positive Y direction */

	glTranslated(parameters.baryCentre[0], parameters.baryCentre[1], parameters.baryCentre[2]);
	glRotated(parameters.rotation[0], 1.0, 0.0, 0.0);
	glRotated(parameters.rotation[1], 0.0, 1.0, 0.0);
	glRotated(parameters.rotation[2], 0.0, 0.0, 1.0);
	glTranslated(-parameters.baryCentre[0], -parameters.baryCentre[1], -parameters.baryCentre[2]);
}

void determineLookAt()
{
	//1. Compute Barycentre
	Eigen::Vector3f baryCentreTemp;
	baryCentreTemp.setZero();

	for (int i = 0; i < particles->size(); ++i)
	{
		baryCentreTemp += (*particles)[i].position();
	}

	baryCentreTemp /= (float)particles->size();
	parameters.baryCentre[0] = baryCentreTemp[0];
	parameters.baryCentre[1] = baryCentreTemp[1];
	parameters.baryCentre[2] = baryCentreTemp[2];

	//2. Compute Radius
	float radiusTemp = 0;
	for (int i = 0; i < particles->size(); ++i)
	{
		Eigen::Vector3f distanceVector = (*particles)[i].position() - baryCentreTemp;
		if (distanceVector.squaredNorm() > radiusTemp)
		{
			radiusTemp = distanceVector.squaredNorm();
		}
	}

	parameters.radius = radiusTemp;

	std::cout << "Barycentre: " << std::endl;
	std::cout << baryCentreTemp << std::endl;
	std::cout << "Radius: " << std::endl;
	std::cout << parameters.radius << std::endl << std::endl;
}

void idleLoopGlut(void)
{
	mainLoop();
}

void mainLoopGlut(void)
{
	mainLoop();
}

void mainLoop()
{
	settings.calculateLambda();
	settings.calculateMu();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);

	setCamera();

	//Render Tets
	for (int t = 0; t < tetrahedra.size(); ++t)
	{
		tetrahedra[t].glRender(0.5, 0.5, 0.5);
	}

	glEnable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset(-1, -1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//Render Tets
	for (int t = 0; t < tetrahedra.size(); ++t)
	{
		tetrahedra[t].glRender(1.0, 1.0, 1.0);
	}
	glDisable(GL_POLYGON_OFFSET_LINE);

	//Advance Solver
	tbb::tick_count start = tbb::tick_count::now();
	if (!parameters.disableSolver)
	{
		if (!parameters.useFEMSolver)
		{
			updateProbabilisticConstraints();

			solver.advanceSystem(tetrahedra, particles, settings, currentPositions, numConstraintInfluences,
				probabilisticConstraints);
		}
		else
		{
			FEMsolver.doTimeStep(true);
			applyFEMDisplacementsToParticles();
		}
	}
	tbb::tick_count end = tbb::tick_count::now();
	parameters.executionTimeSum += (end - start).seconds();
	if (parameters.currentFrame % parameters.timingPrintInterval == 0)
	{
		std::cout << "Average simulation Time: " << parameters.executionTimeSum / parameters.currentFrame << "s." << std::endl;
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	TwDraw();
	glutSwapBuffers();
	if (!parameters.disableSolver)
	{
		++parameters.currentFrame;
	}

	if (parameters.writeToAlembic)
	{
		getCurrentPositionFromParticles();
		smHandler->setSample(currentPositions);
	}

	if (parameters.maxFrames <= parameters.currentFrame)
	{
		std::cout << "Leaving Glut Main Loop..." << std::endl;
		glutLeaveMainLoop();
	}

	Sleep(parameters.numMillisecondsToWaitBetweenFrames);
	//std::cout << "Current Frame: " << currentFrame << std::endl;
}

// Callback function called by GLUT when window size changes
void Reshape(int width, int height)
{
	// Set OpenGL viewport and camera
	glViewport(0, 0, width, height);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//gluPerspective(40, (double)width / height, 1, 10);
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	//gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
	//glTranslatef(0, 0.6f, -1);

	// Send the new window size to AntTweakBar
	TwWindowSize(width, height);
}

// Function called at exit
void TerminateAll(void)
{
	TwTerminate();
}

void collapseMesh()
{
	float minDimValue = -100000000;

	for (int p = 0; p < particles->size(); ++p)
	{
		if ((*particles)[p].position()[parameters.dimToCollapse] > minDimValue)
		{
			minDimValue = (*particles)[p].position()[parameters.dimToCollapse];
		}
	}

	for (int p = 0; p < particles->size(); ++p)
	{
		if ((*particles)[p].inverseMass() != 0.0f)
		{
			(*particles)[p].position()[parameters.dimToCollapse] = minDimValue;
			(*particles)[p].previousPosition()[parameters.dimToCollapse] = minDimValue;
		}
	}

	getCurrentPositionFromParticles();
}

void doIO(float inverseMass, std::vector<int>& vertexConstraintIndices)
{
	Eigen::Vector3f initialVelocity;
	initialVelocity.x() = 0; initialVelocity.y() = 0; initialVelocity.z() = 0;

	// HARD VERTEX CONSTRAINTS
	if (parameters.readVertexConstraintData)
	{
		ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, ioParameters.constraintFile);
		for (int i = 0; i < vertexConstraintIndices.size(); ++i)
		{
			(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
		}
	}

	// TRACKING DATA
	trackingData.resize(ioParameters.trackerFiles.size());
	for (int i = 0; i < ioParameters.trackerFiles.size(); ++i)
	{
		TrackerIO::readTrackerAnimationNuke(ioParameters.trackerFiles[i], trackingData[i]);
	}

	if (!parameters.generateMeshInsteadOfDoingIO)
	{
		// MESH
		TetGenIO::readNodes(ioParameters.nodeFile, *particles, inverseMass, initialVelocity);
		TetGenIO::readTetrahedra(ioParameters.elementFile, tetrahedra, particles);
	}
	else
	{
		//GENERATE the mesh data
		if (parameters.generateMeshFromTrackingData)
		{
			MeshCreator::generateTetBarToFit(particles, tetrahedra, 10, 7, 3,
				trackingData[0][0], trackingData[1][0], trackingData[2][0], 20.0f);
		}
		else
		{
			MeshCreator::generateTetBar(particles, tetrahedra, 10, 6, 6);
		}
	}
}

int main(int argc, char* argv[])
{
	float youngsModulus;

	float poissonRatio;

	int numConstraintIts;

	float invM;

	float timeStep;

	if (argc == 1)
	{
		std::cout << "Please provide the following: " << std::endl;
		std::cout << "	- Youngs Modulus" << std::endl;
		std::cout << "	- Poisson Ratio" << std::endl;
		std::cout << "	- Inverse Mass" << std::endl;
		std::cout << "	- Num Constraint Its" << std::endl;
		std::cout << "	- Time Step Size" << std::endl;
		std::cout << "	- USE_FEM" << std::endl;
		std::cout << "	- SAVE_MESH" << std::endl;
		return 0;
	}

	if (argc > 1)
	{
		//Young's modulus
		youngsModulus = std::stof(std::string(argv[1]));
		//Poisson ratio
		poissonRatio = std::stof(std::string(argv[2]));
		//Inverse Mass
		invM = std::stof(std::string(argv[3]));

		numConstraintIts = std::stoi(std::string(argv[4]));

		timeStep = std::stof(std::string(argv[5]));

		parameters.useFEMSolver = std::string(argv[6]) == "USE_FEM";

		if (argc > 7)
		{
			if (std::string(argv[7]) == "SAVE_MESH")
			{
				parameters.writeToAlembic = true;
				std::cout << "Saving output mesh..." << std::endl;
			}
		}
	}
	else
	{
		std::cout << "ERROR! Not enough paramters found..." << std::endl;
		return 0;
	}

	if (parameters.useFEMSolver)
	{
		smHandler = std::make_shared<SurfaceMeshHandler>("WRITE_TETS", "deformedMeshFEM.abc");
	}
	else
	{
		smHandler = std::make_shared<SurfaceMeshHandler>("WRITE_TETS", "deformedMesh.abc");
	}

	settings.youngsModulus = youngsModulus;
	settings.poissonRatio = poissonRatio;
	settings.deltaT = timeStep;
	//settings.gravity = -9.81f;
	settings.gravity = 0.0f;
	settings.numConstraintIts = numConstraintIts;
	settings.w = 1.0;
	settings.printStrainEnergy = false;
	settings.printStrainEnergyToFile = parameters.printStrainEnergyToFile;
	settings.correctStrongForcesWithSubteps = false;
	settings.numTetrahedraIterations = 1;
	settings.calculateLambda();
	settings.calculateMu();
	settings.print();

	parameters.initialiseToDefaults();
	ioParameters.initialiseToDefaults();

	std::vector<int> vertexConstraintIndices;
	doIO(invM, vertexConstraintIndices);

	settings.numTetrahedra = tetrahedra.size();
	std::cout << "Num Tets: " << tetrahedra.size() << "; Num Nodes: " << particles->size() << std::endl;

	numConstraintInfluences.resize(particles->size());
	currentPositions.resize(particles->size());

	createProabilisticConstraints();

	GLUTSettings glutSettings;
	glutSettings.height = 500;
	glutSettings.width = 500;
	glutSettings.windowName = "PBD FEM";
	glutSettings.GLVersionMajor = 3;
	glutSettings.GLVersionMinor = 0;
	glutSettings.positionX = 100;
	glutSettings.positionY = 100;

	GLUTHelper helper;
	helper.initWindow(argc, argv, glutSettings);
	helper.setIdleFunc(idleLoopGlut);

	determineLookAt();
	parameters.initialiseCamera();

	if (parameters.useFEMSolver)
	{
		std::cout << "Writing .veg file..." << std::endl;
		VegaIO::writeVegFile(tetrahedra, particles, "FEMMesh.veg");

		std::cout << "Setting initial node positions..." << std::endl;
		setInitialPositionsFromParticles();

		std::cout << "Initialising FEM solver..." << std::endl;
		FEMsolver.initSolver("FEMMesh.veg", vertexConstraintIndices, settings.youngsModulus, settings.poissonRatio,
			settings.deltaT);

		std::cout << "Entering simulation loop..." << std::endl;
	}

	if (parameters.writeToAlembic)
	{
		smHandler->initTopology(*particles, tetrahedra);
		std::cout << "Initialised Topology for Alembic Output!" << std::endl;
	}

	if (parameters.testingInversionHandling)
	{
		collapseMesh();
		std::cout << "Collapsed mesh to test inversion Handling!" << std::endl;
	}

	//TweakBar Interface
	TwInit(TW_OPENGL, NULL);
	TwWindowSize(glutSettings.height, glutSettings.width);
	TwBar* solverSettings;
	solverSettings = TwNewBar("Solver Settings");

	TwDefine(" GLOBAL help='FEM based PBD Solver Demo.' ");
	TwAddVarRW(solverSettings, "stepSize", TW_TYPE_FLOAT, &settings.deltaT,
		" label='Step Size' min=0.0001 max=10 step=0.001 keyIncr=s keyDecr=S help='Internal Solver Step Size (0.005 is stable)' ");

	TwAddVarRW(solverSettings, "constraintIts", TW_TYPE_INT32, &settings.numConstraintIts,
		" label='Constraint Iterations' min=1 max=100 step=1 keyIncr=s keyDecr=S help='Internal Solver Constraint Iterations (5 is stable)' ");


	//TwBar* materialSettings;
	//materialSettings = TwNewBar("Material Settings");

	TwAddVarRW(solverSettings, "YoungsModulus", TW_TYPE_FLOAT, &settings.youngsModulus,
		" label='Youngs Modulus' min=0.0 max=1000.0 step=0.01 keyIncr=s keyDecr=S help='Stiffness' ");

	TwAddVarRW(solverSettings, "PoissonRatio", TW_TYPE_FLOAT, &settings.poissonRatio,
		" label='Poisson Ratio' min=0.0 max=0.5 step=0.01 keyIncr=s keyDecr=S help='Poisson Ratio' ");

	TwAddVarRW(solverSettings, "rotationX", TW_TYPE_FLOAT, &parameters.rotation[0],
		" label='Cam Rotation X' min=0.0 max=360.0 step=1 keyIncr=s keyDecr=S help='Rotation about X' ");
	TwAddVarRW(solverSettings, "rotationY", TW_TYPE_FLOAT, &parameters.rotation[1],
		" label='Cam Rotation Y' min=0.0 max=360.0 step=1 keyIncr=s keyDecr=S help='Rotation about X' ");
	TwAddVarRW(solverSettings, "rotationZ", TW_TYPE_FLOAT, &parameters.rotation[2],
		" label='Cam Rotation Z' min=0.0 max=360.0 step=1 keyIncr=s keyDecr=S help='Rotation about X' ");

	TwAddVarRW(solverSettings, "zoom", TW_TYPE_FLOAT, &parameters.zoom,
		" label='Cam Zoom' min=0.0 max=100 step=0.001 keyIncr=s keyDecr=S help='Zoom' ");


	// Set GLUT callbacks
	glutReshapeFunc(Reshape);
	atexit(TerminateAll);  // Called after glutMainLoop ends

	// Set GLUT event callbacks
	// - Directly redirect GLUT mouse button events to AntTweakBar
	glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
	// - Directly redirect GLUT mouse motion events to AntTweakBar
	glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
	glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// - Directly redirect GLUT key events to AntTweakBar
	glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
	// - Directly redirect GLUT special key events to AntTweakBar
	glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
	// - Send 'glutGetModifers' function pointer to AntTweakBar;
	//   required because the GLUT key event functions do not report key modifiers states.
	TwGLUTModifiersFunc(glutGetModifiers);

	helper.enterDisplayLoop(mainLoopGlut);

	return 0;
}