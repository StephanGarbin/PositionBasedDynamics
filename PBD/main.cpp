#include <iostream>
#include <chrono>
#include <thread>

#include <tbb\tick_count.h>
#include <sstream>

#include "TetGenIO.h"
#include "ConstraintsIO.h"
#include "VegaIO.h"

#include "AbcWriter.h"
#include "SurfaceMeshHandler.h"

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "PBDSolver.h"
#include "PBDSolverSettings.h"

#include "GLUTHelper.h"
#include "AntTweakBar.h"

#include "lodepng.h"

#include "FEMSimulator.h"

std::vector<PBDTetrahedra3d> tetrahedra;
std::shared_ptr<std::vector<PBDParticle>> particles = std::make_shared<std::vector<PBDParticle>>();
std::vector<Eigen::Vector3f> currentPositions;
std::vector<Eigen::Vector3f> initialPositions;

std::vector<int> numConstraintInfluences;
PBDSolver solver;
FEMSimulator FEMsolver;

std::shared_ptr<SurfaceMeshHandler> smHandler;

PBDSolverSettings settings;

int numMilliseconds = 1000;

double sumExecutionTime;
int timingPrintInterval = 100;
int currentFrame = 1;
int maxFrames = 12000;

int globalHeight;
int globalWidth;

float baryCentre[3];
float radius;
float rotation[3];
float zoom;

void mainLoop();

float youngsModulus;
float poissonRatio;

float lambda;
float mu;

bool useFEMSolver = false;
bool writeToAlembic = true;
bool printStrainEnergyToFile = false;

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

void saveFrameBufferAsPng()
{
	unsigned width = globalWidth;
	unsigned height = globalHeight;
	

	std::stringstream ss;
	ss << "images/solverOutput_" << currentFrame << ".png";
	
	std::vector<double> RawImage;
	RawImage.resize(globalWidth * globalHeight * 4);

	glReadPixels(0, 0, width, height, GL_RGBA, GL_DOUBLE, &RawImage[0]);


	std::vector<unsigned char> image;
	image.resize(globalWidth * globalHeight * 4);
	for (int i = 0; i < globalHeight * globalWidth * 4; ++i)
	{
		image[i] = (unsigned char)(RawImage[i] * 255.0f);
	}

	std::vector<unsigned char> png;

	unsigned error = lodepng::encode(png, image, width, height);
	if (!error) lodepng::save_file(png, ss.str().c_str());

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;

	ss.clear();
}

void calculateLambdaAndMu()
{
	mu = youngsModulus / (2.0 * (1.0 + poissonRatio));
	lambda = (youngsModulus * poissonRatio) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio));

	settings.lambda = lambda;
	settings.mu = mu;

	settings.youngsModulus = youngsModulus;
	settings.poissonRatio = poissonRatio;
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

	gluLookAt(-(baryCentre[0] + radius * zoom), baryCentre[1] + radius * zoom, -(baryCentre[2] + radius * zoom),  /* eye is at (0,0,5) */
		baryCentre[0], baryCentre[1], baryCentre[2],      /* center is at (0,0,0) */
		0.0, 1.0, 0.0);      /* up is in positive Y direction */

	glTranslated(baryCentre[0], baryCentre[1], baryCentre[2]);
	glRotated(rotation[0], 1.0, 0.0, 0.0);
	glRotated(rotation[1], 0.0, 1.0, 0.0);
	glRotated(rotation[2], 0.0, 0.0, 1.0);
	glTranslated(-baryCentre[0], -baryCentre[1], -baryCentre[2]);
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
	baryCentre[0] = baryCentreTemp[0];
	baryCentre[1] = baryCentreTemp[1];
	baryCentre[2] = baryCentreTemp[2];

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

	radius = radiusTemp;

	//std::cout << "Barycentre: " << std::endl;
	//std::cout << baryCentreTemp << std::endl;
	//std::cout << "Radius: " << std::endl;
	//std::cout << radius << std::endl << std::endl;
}

void LookAtMesh()
{
	//gluLookAt(5, 5, 5, 0.0, 0.0, 0.0, 0, 1, 0);
	//glTranslatef(0, 0.6f, -1);

	//gluLookAt(baryCentre[0] + radius, baryCentre[1] + radius, baryCentre[2] - radius,
	//	baryCentre[0], baryCentre[1], baryCentre[2], 0.0, 1.0, 0.0);

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
	calculateLambdaAndMu();

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
	if (!useFEMSolver)
	{
		solver.advanceSystem(tetrahedra, particles, settings, currentPositions, numConstraintInfluences);
	}
	else
	{
		FEMsolver.doTimeStep(true);
		applyFEMDisplacementsToParticles();
	}
	tbb::tick_count end = tbb::tick_count::now();
	sumExecutionTime += (end - start).seconds();
	if (currentFrame % timingPrintInterval == 0)
	{
		std::cout << "Average simulation Time: " << sumExecutionTime / currentFrame << "s." << std::endl;
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	TwDraw();
	glutSwapBuffers();
	++currentFrame;

	if (writeToAlembic)
	{
		getCurrentPositionFromParticles();
		smHandler->setSample(currentPositions);
	}

	if (maxFrames <= currentFrame)
	{
		std::cout << "Leaving Glut Main Loop..." << std::endl;
		glutLeaveMainLoop();
	}
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
	globalHeight = height;
	globalWidth = width;
}

// Function called at exit
void TerminateAll(void)
{
	TwTerminate();
}


int main(int argc, char* argv[])
{
	//Young's modulus
	float k;

	//Poisson ratio
	float v;

	//constraint iterations
	int numConstraintIts;

	//Inverse Mass
	float invM;

	bool useFEM;
	bool useSOR = false;

	float timeStep;

	if (argc > 1)
	{
		//Young's modulus
		k = std::stof(std::string(argv[1]));

		//Poisson ratio
		v = std::stof(std::string(argv[2]));

		numConstraintIts = std::stoi(std::string(argv[3]));

		invM = std::stof(std::string(argv[4]));

		timeStep = std::stof(std::string(argv[5]));

		useFEM = std::string(argv[6]) == "USE_FEM";

		if (argc > 7)
		{
			if (std::string(argv[7]) == "SAVE_MESH")
			{
				writeToAlembic = true;
				std::cout << "Saving output mesh..." << std::endl;
			}
		}

		youngsModulus = k;
		poissonRatio = v;
	}
	else
	{
		k = 1.0;
		v = 0.4333;
		numConstraintIts = 5;
		invM = 1;
		useSOR = false;
		timeStep = 0.005;
	}

	if (useSOR)
	{
		std::cout << "Using SOR Solver..." << std::endl;
	}

	if (useFEM)
	{
		useFEMSolver = true;
		smHandler = std::make_shared<SurfaceMeshHandler>("WRITE_TETS", "deformedMeshFEM.abc");
	}
	else
	{
		smHandler = std::make_shared<SurfaceMeshHandler>("WRITE_TETS", "deformedMesh.abc");
	}

	Eigen::Vector3f initialVelocity;
	initialVelocity.x() = 0; initialVelocity.y() = 0; initialVelocity.z() = 0;

	calculateLambdaAndMu();

	settings.youngsModulus = youngsModulus;
	settings.poissonRatio = poissonRatio;
	settings.deltaT = timeStep;
	settings.gravity = -9.8;
	settings.lambda = lambda;
	settings.mu = mu;
	settings.numConstraintIts = numConstraintIts;
	settings.w = 1.0;
	settings.printStrainEnergy = false;
	settings.printStrainEnergyToFile = printStrainEnergyToFile;
	settings.useSOR = useSOR;
	settings.correctStrongForcesWithSubteps = false;
	settings.numTetrahedraIterations = 1;
	settings.print();
	
	std::string nodes("barout.node");
	std::string tets("barout.ele");
	TetGenIO::readNodes(nodes, *particles, invM, initialVelocity);
	TetGenIO::readTetrahedra(tets, tetrahedra, particles);


	settings.numTetrahedra = tetrahedra.size();

	std::vector<int> vertexConstraintIndices;
	ConstraintsIO::readMayaVertexConstraints(vertexConstraintIndices, "barLowVertexConstraints.txt");

	for (int i = 0; i < vertexConstraintIndices.size(); ++i)
	{
		(*particles)[vertexConstraintIndices[i]].inverseMass() = 0.0;
	}

	std::cout << "Finished Reading Data From Disk, starting simulation ... " << std::endl;
	std::cout << "Num Tets: " << tetrahedra.size() << "; Num Nodes: " << particles->size() << std::endl;

	numConstraintInfluences.resize(particles->size());
	currentPositions.resize(particles->size());

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
	rotation[0] = 0.0;
	rotation[1] = 128.0;
	rotation[2] = 0.0;
	zoom = 1.0 / 12.0f;

	if (useFEMSolver)
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

	if (writeToAlembic)
	{
		smHandler->initTopology(*particles, tetrahedra);
		std::cout << "Initialised Topology for Alembic Output!" << std::endl;
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

	TwAddVarRW(solverSettings, "YoungsModulus", TW_TYPE_FLOAT, &youngsModulus,
		" label='Youngs Modulus' min=0.0 max=1000.0 step=0.01 keyIncr=s keyDecr=S help='Stiffness' ");

	TwAddVarRW(solverSettings, "PoissonRatio", TW_TYPE_FLOAT, &poissonRatio,
		" label='Poisson Ratio' min=0.0 max=0.5 step=0.01 keyIncr=s keyDecr=S help='Poisson Ratio' ");

	TwAddVarRW(solverSettings, "rotationX", TW_TYPE_FLOAT, &rotation[0],
		" label='Cam Rotation X' min=0.0 max=360.0 step=1 keyIncr=s keyDecr=S help='Rotation about X' ");
	TwAddVarRW(solverSettings, "rotationY", TW_TYPE_FLOAT, &rotation[1],
		" label='Cam Rotation Y' min=0.0 max=360.0 step=1 keyIncr=s keyDecr=S help='Rotation about X' ");
	TwAddVarRW(solverSettings, "rotationZ", TW_TYPE_FLOAT, &rotation[2],
		" label='Cam Rotation Z' min=0.0 max=360.0 step=1 keyIncr=s keyDecr=S help='Rotation about X' ");

	TwAddVarRW(solverSettings, "zoom", TW_TYPE_FLOAT, &zoom,
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