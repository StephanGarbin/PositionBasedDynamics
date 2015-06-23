#include <iostream>
#include <chrono>
#include <thread>

#include <tbb\tick_count.h>

#include "TetGenIO.h"

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "PBDSolver.h"
#include "PBDSolverSettings.h"

#include "GLUTHelper.h"
#include "AntTweakBar.h"


std::vector<PBDTetrahedra3d> tetrahedra;
std::shared_ptr<std::vector<PBDParticle>> particles = std::make_shared<std::vector<PBDParticle>>();
std::vector<Eigen::Vector3f> temporaryPositions;
std::vector<int> numConstraintInfluences;
PBDSolver solver;

PBDSolverSettings settings;

int numMilliseconds = 0;

double sumExecutionTime;
int timingPrintInterval = 100;
int currentFrame = 1;

float baryCentre[3];
float radius;

void mainLoop();

float youngsModulus;
float poissonRatio;

float lambda;
float mu;

void calculateLambdaAndMu()
{
	mu = youngsModulus / (2.0 * (1.0 + poissonRatio));
	lambda = (youngsModulus * poissonRatio) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio));

	settings.lambda = lambda;
	settings.mu = mu;
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

	radius = radiusTemp / 13;

	std::cout << "Barycentre: " << std::endl;
	std::cout << baryCentreTemp << std::endl;
	std::cout << "Radius: " << std::endl;
	std::cout << radius << std::endl << std::endl;
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

	//glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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
	solver.advanceSystem(tetrahedra, particles, settings, temporaryPositions, numConstraintInfluences);
	tbb::tick_count end = tbb::tick_count::now();
	sumExecutionTime += (end - start).seconds();
	if (currentFrame % timingPrintInterval == 0)
	{
		std::cout << "Average simulation Time: " << sumExecutionTime / currentFrame << "s." << std::endl;
	}

	TwDraw();
	glutSwapBuffers();
}

// Callback function called by GLUT when window size changes
void Reshape(int width, int height)
{
	// Set OpenGL viewport and camera
	//glViewport(0, 0, width, height);
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

	if (argc > 1)
	{
		//Young's modulus
		k = std::stof(std::string(argv[1]));

		//Poisson ratio
		v = std::stof(std::string(argv[2]));

		numConstraintIts = std::stoi(std::string(argv[3]));

		invM = std::stof(std::string(argv[4]));
	}
	else
	{
		k = 1.0;
		v = 0.4333;
		numConstraintIts = 5;
		invM = 1;
	}


	Eigen::Vector3f initialVelocity;
	initialVelocity.x() = 0; initialVelocity.y() = 0; initialVelocity.z() = 0;

	calculateLambdaAndMu();

	youngsModulus = 1.0;
	poissonRatio = 0.3;

	settings.deltaT = 0.005;
	settings.gravity = -9.8;
	settings.lambda = lambda;
	settings.mu = mu;
	settings.numConstraintIts = numConstraintIts;


	//Test mesh

	//particles->emplace_back(Eigen::Vector3d(0.029, 5.173, -0.061), Eigen::Vector3d(0.0, 0.0, 0.0), 0);
	//particles->emplace_back(Eigen::Vector3d(-0.739, 3.002, 1.269), Eigen::Vector3d(0.0, 0.0, 0.0), invM);
	//particles->emplace_back(Eigen::Vector3d(1.564, 3.002, -0.061), Eigen::Vector3d(0.0, 0.0, 0.0), invM);
	//particles->emplace_back(Eigen::Vector3d(-0.739, 3.002, -1.390), Eigen::Vector3d(0.0, 0.0, 0.0), invM);
	//particles->emplace_back(Eigen::Vector3d(0.029, 0.826, -0.061), Eigen::Vector3d(0.0, 0.0, 0.0), invM);

	//std::vector<int> temp = { 1, 0, 2, 3 };
	//tetrahedra.emplace_back(std::move(temp), particles);
	//std::vector<int> temp2 = { 4, 3, 2, 1 };
	//tetrahedra.emplace_back(std::move(temp2), particles);
	
	std::string nodes("barout.node");
	std::string tets("barout.ele");
	TetGenIO::readNodes(nodes, *particles, invM, initialVelocity);
	TetGenIO::readTetrahedra(tets, tetrahedra, particles);


	settings.numTetrahedra = tetrahedra.size();

	for (int i = 0; i < 4; ++i)
	{
		(*particles)[i].inverseMass() = 0;
	}

	std::cout << "Finished Reading Data From Disk, starting simulation ... " << std::endl;
	std::cout << "Num Tets: " << tetrahedra.size() << "; Num Nodes: " << particles->size() << std::endl;

	numConstraintInfluences.resize(particles->size());
	temporaryPositions.resize(particles->size());

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
	helper.initCamera(-(baryCentre[0] + radius), baryCentre[1] + radius, -(baryCentre[2] + radius),
		baryCentre[0], baryCentre[1], baryCentre[2]);

	//TweakBar Interface
	TwInit(TW_OPENGL, NULL);
	TwWindowSize(glutSettings.height, glutSettings.width);
	TwBar* solverSettings;
	solverSettings = TwNewBar("Solver Settings");

	TwDefine(" GLOBAL help='FEM based PBD Solver Demo.' ");
	TwAddVarRW(solverSettings, "stepSize", TW_TYPE_FLOAT, &settings.deltaT,
		" label='Step Size' min=0.0001 max=10 step=0.01 keyIncr=s keyDecr=S help='Internal Solver Step Size (0.005 is stable)' ");

	TwAddVarRW(solverSettings, "constraintIts", TW_TYPE_INT32, &settings.numConstraintIts,
		" label='Constraint Iterations' min=1 max=100 step=1 keyIncr=s keyDecr=S help='Internal Solver Constraint Iterations (5 is stable)' ");


	TwBar* materialSettings;
	materialSettings = TwNewBar("Material Settings");

	TwAddVarRW(materialSettings, "YoungsModulus", TW_TYPE_FLOAT, &youngsModulus,
		" label='Youngs Modulus' min=0.0 max=100.0 step=0.01 keyIncr=s keyDecr=S help='Stiffness' ");

	TwAddVarRW(materialSettings, "PoissonRatio", TW_TYPE_FLOAT, &poissonRatio,
		" label='Poisson Ratio' min=0.0 max=0.5 step=0.01 keyIncr=s keyDecr=S help='Poisson Ratio' ");


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