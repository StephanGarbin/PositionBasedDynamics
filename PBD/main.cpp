#define _USE_MATH_DEFINES

#include <iostream>
#include <chrono>
#include <thread>

#include <tbb\tick_count.h>
#include <sstream>

#include "TetGenIO.h"
#include "ConstraintsIO.h"
#include "VegaIO.h"
#include "TrackerIO.h"
#include "cImageIO.h"

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

#include "AppHelper.h"
#include "FiberMesh.h"
#include "CollisionMesh.h"
#include "CollisionRod.h"
#include "CollisionSphere.h"

std::vector<PBDTetrahedra3d> tetrahedra;
std::shared_ptr<std::vector<PBDParticle>> particles = std::make_shared<std::vector<PBDParticle>>();
std::vector<Eigen::Vector3f> currentPositions;
std::vector<Eigen::Vector3f> initialPositions;
std::vector<int> numConstraintInfluences;
std::vector<PBDProbabilisticConstraint> probabilisticConstraints;
std::vector<std::vector<Eigen::Vector2f>> trackingData;
std::shared_ptr<FiberMesh> fiberMesh;
std::vector<CollisionMesh> collisionGeometry;
std::vector<CollisionRod> collisionGeometry2;
std::vector<CollisionSphere> collisionGeometry3;

PBDSolver solver;
FEMSimulator FEMsolver;

std::shared_ptr<SurfaceMeshHandler> smHandler;

Parameters parameters;
IOParameters ioParameters;

void mainLoop();
void applyInitialDeformationToMesh();
void disapplyInitialDeformationToMesh();

void applyContinuousDeformationToMesh();

void collapseMesh();

void createFiberMesh();

void invertSingleElementAtStart();

void applyCollisitionGeometryTranslation();

void applyPressure();

std::string generateFileName(const std::string& base, const std::string& extension, int idx, int version)
{
	std::stringstream ss;
	ss << base << "_" << idx << "_" << version << "." << extension;
	std::string result = ss.str();
	ss.clear();

	return result;
}

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
	std::vector<int> activeConstraintIndices = { 0, 1 };
	std::vector<int> trackingDataIndices = { 3, 2 };

	for (int i = 0; i < activeConstraintIndices.size(); ++i)
	{
		Eigen::Vector3f currentConstraintPosition =
			TrackerIO::getInterpolatedConstraintPosition(trackingData[trackingDataIndices[i]], 1.0f / 24.0f,
			parameters.solverSettings.deltaT, parameters.getCurrentFrame() * parameters.solverSettings.deltaT);

		float scaleFactor = 0.01;

		//1. subtract displacement
		float widthSubtraction = trackingData[0][0].x();
		float heightSubtraction = trackingData[0][0].y();

		currentConstraintPosition[0] -= widthSubtraction;
		currentConstraintPosition[2] -= heightSubtraction;

		//2. Rescale
		currentConstraintPosition *= scaleFactor;

		probabilisticConstraints[activeConstraintIndices[i]].getConstraintPosition() = currentConstraintPosition;
	}
}

void createProabilisticConstraints()
{
	probabilisticConstraints.resize(2);
	updateProbabilisticConstraints();
	probabilisticConstraints[0].initialise(*particles, 0.1f);
	probabilisticConstraints[1].initialise(*particles, 0.2f);
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
	//Apply initial deformation if necessary
	if (parameters.applyInitialDeformationToMesh)
	{
		if (parameters.getCurrentFrame() == parameters.frame2ApplyInitialDeformation)
		{
			applyInitialDeformationToMesh();
		}

		if (parameters.getCurrentFrame() == parameters.frame2DisApplyInitialDeformation)
		{
			disapplyInitialDeformationToMesh();
		}
	}

	if (parameters.applyContinuousDeformationToMesh)
	{
		applyContinuousDeformationToMesh();
	}

	if (parameters.invertSingleElementAtStart)
	{
		invertSingleElementAtStart();
	}

	if (parameters.collapseMeshAtStart)
	{
		collapseMesh();
	}

	if (parameters.createFiberMesh)
	{
		createFiberMesh();
	}

	if (parameters.translateCollisionGeometry)
	{
		applyCollisitionGeometryTranslation();
	}

	if (parameters.applyPressure)
	{
		applyPressure();
	}

	parameters.solverSettings.calculateLambda();
	parameters.solverSettings.calculateMu();
	parameters.solverSettings.calculateFiberStructureTensor();

	//GLenum err = glGetError();
	//if (err != GL_NO_ERROR)
	//{
	//	std::cerr << "GL Error (Before glPolygonMode): " << err << std::endl;
	//}

	//glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);

	setCamera();

	//Render Tets
	//for (int t = 0; t < tetrahedra.size(); ++t)
	//{
	//	tetrahedra[t].glRender(0.5, 0.5, 0.5);
	//}

	//GLenum err = glGetError();
	//if (err != GL_NO_ERROR)
	//{
	//	std::cerr << "GL Error (Before glPolygonMode): " << err << std::endl;
	//}

	//glEnable(GL_POLYGON_OFFSET_LINE);
	//glPolygonOffset(-1, -1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//Render Tets
	for (int t = 0; t < tetrahedra.size(); ++t)
	{
		tetrahedra[t].glRender(1.0, 1.0, 1.0);
	}
	//glDisable(GL_POLYGON_OFFSET_LINE);

	if (parameters.renderCollisionGoemetry)
	{
		for (int i = 0; i < collisionGeometry2.size(); ++i)
		{
			collisionGeometry2[i].glRender(parameters.getCurrentFrame(), parameters.solverSettings.deltaT,
				parameters.solverSettings.collisionSpheresNum[i], parameters.solverSettings.collisionSpheresRadius[i]);
		}

		for (int i = 0; i < collisionGeometry3.size(); ++i)
		{
			collisionGeometry3[i].glRender(parameters.getCurrentFrame(), parameters.solverSettings.deltaT,
				parameters.solverSettings.collisionSpheresRadius[i]);
		}
	}

	for (int i = 0; i < probabilisticConstraints.size(); ++i)
	{
		glPushMatrix();
		glTranslated(probabilisticConstraints[i].getConstraintPosition().x(),
			probabilisticConstraints[i].getConstraintPosition().y(),
			probabilisticConstraints[i].getConstraintPosition().z());

		glutWireSphere(probabilisticConstraints[i].getInitialRadius(), 8, 8);
		glPopMatrix();
	}

	//Advance Solver
	tbb::tick_count start = tbb::tick_count::now();
	if (!parameters.disableSolver)
	{
		if (!parameters.useFEMSolver)
		{
			if (parameters.useTrackingConstraints)
			{
				updateProbabilisticConstraints();
			}
			solver.advanceSystem(tetrahedra, particles, parameters.solverSettings, currentPositions, numConstraintInfluences,
				probabilisticConstraints, collisionGeometry, collisionGeometry2, collisionGeometry3);
		}
		else
		{
			FEMsolver.doTimeStep(true);
			applyFEMDisplacementsToParticles();
		}
	}
	tbb::tick_count end = tbb::tick_count::now();
	parameters.executionTimeSum += (end - start).seconds();
	if (parameters.getCurrentFrame() % parameters.timingPrintInterval == 0)
	{
		std::cout << "Average simulation Time: " << parameters.executionTimeSum / parameters.getCurrentFrame() << "s."
			<< "FRAME: [ " << parameters.getCurrentFrame() << " ]" << std::endl;
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	if (parameters.doImageIO)
	{
		saveImage(700, 700, parameters.getCurrentFrame(), parameters.imageFileName);
	}

	TwDraw();

	glutSwapBuffers();
	if (!parameters.disableSolver)
	{
		parameters.increaseCurrentFrame();
	}

	if (parameters.writeToAlembic)
	{
		getCurrentPositionFromParticles();
		smHandler->setSample(currentPositions);
	}

	if (parameters.maxFrames <= parameters.getCurrentFrame())
	{
		//Write all debug information
		parameters.solverSettings.tracker.writeAll();

		std::cout << "Leaving Glut Main Loop..." << std::endl;
		glutLeaveMainLoop();
	}

	/*if (parameters.doImageIO)
	{
		saveImage(700, 700, parameters.getCurrentFrame(), parameters.imageFileName);
	}*/

	Sleep(parameters.numMillisecondsToWaitBetweenFrames);
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
	if (parameters.getCurrentFrame() != 1)
	{
		return;
	}

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

void applyInitialDeformationToMesh()
{
	(*particles)[1].position().y() += 0.5f;
	(*particles)[1].previousPosition().y() += 0.5f;
}

void disapplyInitialDeformationToMesh()
{
	(*particles)[1].position().y() -= 0.5f;
	(*particles)[1].previousPosition().y() -= 0.5f;
}

void applyContinuousDeformationToMesh()
{
	if (parameters.getCurrentFrame() < parameters.continuousDeformationRelaxationFrame)
	{
		//float increaseFactor = parameters.continuousDeformationStrainIncreaseFactor;
		//(*particles)[1].position().y() += parameters.continuousDeformationStrainIncreaseFactor;
		(*particles)[1].previousPosition().z() -= parameters.continuousDeformationStrainIncreaseFactor * (parameters.continuousDeformationRelaxationFrame - parameters.getCurrentFrame());
	}
	else if (parameters.getCurrentFrame() < parameters.continuousDeformationRelaxationFrame * 2)
	{
		//(*particles)[1].position().y() -= parameters.continuousDeformationStrainIncreaseFactor;
		(*particles)[1].previousPosition().z() += parameters.continuousDeformationStrainIncreaseFactor * ((parameters.getCurrentFrame() - parameters.continuousDeformationRelaxationFrame));
	}
	//else
	//{
	//	(*particles)[1].inverseMass() = 1.0f;
	//}
}

void invertSingleElementAtStart()
{
	if (parameters.getCurrentFrame() == 1)
	{
		switch (parameters.TEST_VERSION)
		{
		case 0:
			(*particles)[1].inverseMass() = 1.0f;
			(*particles)[1].position().y() -= parameters.invertSingleElementAtStartAmount;
			(*particles)[1].previousPosition().y() -= parameters.invertSingleElementAtStartAmount;
			break;
		case 1:
			(*particles)[2].inverseMass() = 1.0f;
			(*particles)[2].position().z() -= parameters.invertSingleElementAtStartAmount;
			(*particles)[2].previousPosition().z() -= parameters.invertSingleElementAtStartAmount;
			break;
		case 2:
			(*particles)[3].inverseMass() = 1.0f;
			(*particles)[3].position().x() += parameters.invertSingleElementAtStartAmount;
			(*particles)[3].previousPosition().x() += parameters.invertSingleElementAtStartAmount;
			break;
		case 3:
			(*particles)[0].position().x() = 0.3f;
			(*particles)[0].position().y() = 0.3f;
			(*particles)[0].position().z() = 0.6f;

			(*particles)[1].position().x() = 0.3f;
			(*particles)[1].position().y() = 0.0f;
			(*particles)[1].position().z() = 0.599135f;
			
			(*particles)[2].position().x() = 0.0f;
			(*particles)[2].position().y() = 0.0f;
			(*particles)[2].position().z() = 0.0f;
			
			(*particles)[3].position().x() = 0.3f;
			(*particles)[3].position().y() = 0.0f;
			(*particles)[3].position().z() = 0.6f;

			(*particles)[0].previousPosition() = (*particles)[0].position();
			(*particles)[1].previousPosition() = (*particles)[1].position();
			(*particles)[2].previousPosition() = (*particles)[2].position();
			(*particles)[3].previousPosition() = (*particles)[3].position();

			(*particles)[0].inverseMass() = 1.0f;
			(*particles)[1].inverseMass() = 1.0f;
			(*particles)[2].inverseMass() = 0.0f;
			(*particles)[3].inverseMass() = 1.0f;
		}
	}
}

void createFiberMesh()
{
	if (parameters.getCurrentFrame() != 1)
	{
		return;
	}

	fiberMesh = std::make_shared <FiberMesh>(particles, &tetrahedra);

	//	Eigen::Vector3f origin(0.0f, 0.0f, 0.0f);
	//	Eigen::Vector3f dimension(5.0f, 3.0f, 3.0f);
	//	Eigen::Vector3f rotation(0.0f, 0.0f, M_PI / 2.0f);
	//
	//	fiberMesh->generateFibersToFillCube(origin, rotation, dimension,
	//		2, 2, 0.0f, false);
	//}
}

void applyCollisitionGeometryTranslation()
{
	if (parameters.getCurrentFrame() < parameters.collisionGeometryTranslateUntilFrame)
	{
		collisionGeometry[0].getCollisionMeshTranslation() += parameters.collisionGeometryTranslationAmount;
	}
}

void applyPressure()
{
	if (parameters.getCurrentFrame() > parameters.pressureStartFrame
		&& parameters.getCurrentFrame() < parameters.pressureEndFrame)
	{
		for (int i = 0; i < particles->size(); ++i)
		{
			float dist = ((*particles)[i].position() - parameters.pressureCentre).squaredNorm();
			if (dist < parameters.pressureRadius)
			{
				(*particles)[i].previousVelocity() += parameters.pressureForce * parameters.solverSettings.deltaT;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	if (!parseTerminalParameters(argc, argv, parameters, ioParameters))
	{
		return 0;
	}

	std::cout << "Parameters setup completed..." << std::endl;

	//We potentially need these for the FEM solver
	std::vector<int> vertexConstraintIndices;

	if (!doIO(parameters, ioParameters, vertexConstraintIndices,
		tetrahedra, particles, trackingData, collisionGeometry2, collisionGeometry3))
	{
		return 0;
	}

	std::cout << "IO completed..." << std::endl;

	std::cout << "MESH COMPLEXITY: " << std::endl;
	std::cout << "Num Tets: " << tetrahedra.size() << "; Num Nodes: " << particles->size() << std::endl;
	std::cout << "----------------------------------------------" << std::endl;

	numConstraintInfluences.resize(particles->size());
	currentPositions.resize(particles->size());
	
	if (parameters.useTrackingConstraints)
	{
		createProabilisticConstraints();
	}


	GLUTSettings glutSettings;
	glutSettings.height = 600;
	glutSettings.width = 600;
	glutSettings.windowName = "PBD FEM";
	glutSettings.GLVersionMajor = 3;
	glutSettings.GLVersionMinor = 0;
	glutSettings.positionX = 100;
	glutSettings.positionY = 100;
	GLUTHelper helper;
	helper.initWindow(argc, argv, glutSettings);
	helper.setIdleFunc(idleLoopGlut);
	determineLookAt();
	//parameters.initialiseCamera();
	if (parameters.useFEMSolver)
	{
		std::cout << "FEM SOLVER INIT:" << std::endl;
		std::cout << "Writing .veg file..." << std::endl;
		VegaIO::writeVegFile(tetrahedra, particles, "FEMMesh.veg");

		std::cout << "Setting initial node positions..." << std::endl;
		setInitialPositionsFromParticles();

		std::cout << "Initialising FEM solver..." << std::endl;
		FEMsolver.initSolver("FEMMesh.veg", vertexConstraintIndices,
			parameters.solverSettings.youngsModulus, parameters.solverSettings.poissonRatio,
			parameters.solverSettings.deltaT);
		std::cout << "----------------------------------------------" << std::endl;
	}

	if (parameters.writeToAlembic)
	{
		if (parameters.useFEMSolver)
		{
			smHandler = std::make_shared<SurfaceMeshHandler>("WRITE_TETS", generateFileName("deformedMeshFEM", "abc", parameters.TEST_IDX, parameters.TEST_VERSION));
		}
		else
		{
			smHandler = std::make_shared<SurfaceMeshHandler>("WRITE_TETS", generateFileName("deformedMesh", "abc", parameters.TEST_IDX, parameters.TEST_VERSION));
		}
		smHandler->initTopology(*particles, tetrahedra);
		std::cout << "Initialised Topology for Alembic Output!" << std::endl;
	}

	if (parameters.readCollisionGeometry)
	{
		collisionGeometry.resize(parameters.collisionGeometryFiles.size());
		for (int c = 0; c < parameters.collisionGeometryFiles.size(); ++c)
		{
			collisionGeometry[c].readFromAbc(parameters.collisionGeometryFiles[c]);
		}

		std::cout << "Read collision geometry files!" << std::endl;
	}

	//TweakBar Interface
	TwInit(TW_OPENGL, NULL);
	TwWindowSize(glutSettings.height, glutSettings.width);
	TwBar* solverSettings;
	solverSettings = TwNewBar("Solver Settings");

	TwDefine(" GLOBAL help='FEM based PBD Solver Demo.' ");
	TwAddVarRW(solverSettings, "stepSize", TW_TYPE_FLOAT, &parameters.solverSettings.deltaT,
		" label='Step Size' min=0.0001 max=10 step=0.001 keyIncr=s keyDecr=S help='Internal Solver Step Size (0.005 is stable)' ");

	TwAddVarRW(solverSettings, "constraintIts", TW_TYPE_INT32, &parameters.solverSettings.numConstraintIts,
		" label='Constraint Iterations' min=1 max=100 step=1 keyIncr=s keyDecr=S help='Internal Solver Constraint Iterations (5 is stable)' ");

	TwAddSeparator(solverSettings, NULL, NULL);
	TwAddButton(solverSettings, "comment0", NULL, NULL, " label='Elasticity' ");
	TwAddVarRW(solverSettings, "YoungsModulus", TW_TYPE_FLOAT, &parameters.solverSettings.youngsModulus,
		" label='Youngs Modulus' min=0.0 max=1000.0 step=0.01 keyIncr=s keyDecr=S help='Stiffness' ");

	TwAddVarRW(solverSettings, "PoissonRatio", TW_TYPE_FLOAT, &parameters.solverSettings.poissonRatio,
		" label='Poisson Ratio' min=0.0 max=0.5 step=0.01 keyIncr=s keyDecr=S help='Poisson Ratio' ");

	TwAddSeparator(solverSettings, NULL, NULL);
	TwAddButton(solverSettings, "comment1", NULL, NULL, " label='Viscoelasticity' ");
	TwAddVarRW(solverSettings, "Alpha", TW_TYPE_FLOAT, &parameters.solverSettings.alpha,
		" label='Alpha' min=0.0 max=1.0 step=0.01 keyIncr=s keyDecr=S help='Alpha' ");
	TwAddVarRW(solverSettings, "Rho", TW_TYPE_FLOAT, &parameters.solverSettings.rho,
		" label='Rho' min=0.0 max=100 step=0.01 keyIncr=s keyDecr=S help='Rho' ");

	TwAddSeparator(solverSettings, NULL, NULL);
	TwAddButton(solverSettings, "comment2", NULL, NULL, " label='Anisotropy' ");
	TwAddVarRW(solverSettings, "AnisotropyStrength", TW_TYPE_FLOAT, &parameters.solverSettings.anisotropyParameter,
		" label='Strength' min=0.0 max=10000 step=0.01 keyIncr=s keyDecr=S help='Strength of Anisotropic Material Component' ");

	TwAddVarRW(solverSettings, "AnistropyDir", TW_TYPE_DIR3F, &parameters.solverSettings.MR_a[0],
		" label='Dir' help='Direction of Material Anisotropy' ");

	TwAddSeparator(solverSettings, NULL, NULL);

	TwAddVarRW(solverSettings, "Gravity", TW_TYPE_FLOAT, &parameters.solverSettings.gravity,
		" label='Gravity' min=-100.0 max=100 step=0.01 keyIncr=s keyDecr=S help='Gravity' ");
		
	TwAddSeparator(solverSettings, NULL, NULL);

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