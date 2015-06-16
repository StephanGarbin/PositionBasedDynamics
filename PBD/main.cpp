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


std::vector<PBDTetrahedra3d> tetrahedra;
std::shared_ptr<std::vector<PBDParticle>> particles = std::make_shared<std::vector<PBDParticle>>();
PBDSolverSettings settings;

int numMilliseconds = 0;

double sumExecutionTime;
int timingPrintInterval = 100;
int currentFrame = 1;

void mainLoopGLUT(void)
{
	//glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);

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
	PBDSolver::advanceSystem(tetrahedra, particles, settings);
	tbb::tick_count end = tbb::tick_count::now();
	sumExecutionTime += (end - start).seconds();
	if (currentFrame % timingPrintInterval == 0)
	{
		std::cout << "Average simulation Time: " << sumExecutionTime / currentFrame << "s." << std::endl;
	}
	glutSwapBuffers();
	++currentFrame;
}

void idleLoopGLUT(void)
{
	//glEnable(GL_DEPTH_TEST);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);

	//Render Tets
	for (int t = 0; t < tetrahedra.size(); ++t)
	{
		tetrahedra[t].glRender(0.5, 0.5, 0.5);
	}

	glEnable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset(-1, -1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int t = 0; t < tetrahedra.size(); ++t)
	{
		tetrahedra[t].glRender(1.0, 1.0, 1.0);
	}
	glDisable(GL_POLYGON_OFFSET_LINE);


	//Advance Solver
	tbb::tick_count start = tbb::tick_count::now();
	PBDSolver::advanceSystem(tetrahedra, particles, settings);
	tbb::tick_count end = tbb::tick_count::now();
	sumExecutionTime += (end - start).seconds();
	if (currentFrame % timingPrintInterval == 0)
	{
		std::cout << "Average simulation Time: " << sumExecutionTime / currentFrame << "s." << std::endl;
	}
	glutSwapBuffers();
	++currentFrame;

	std::this_thread::sleep_for(std::chrono::milliseconds(numMilliseconds));
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

	float mu = k / (2 * (1 + v));
	float lambda = (k * v) / ((1 + v) * (1 - 2 * v));

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


	GLUTSettings glutSsettings;
	glutSsettings.height = 1024;
	glutSsettings.width = 1024;
	glutSsettings.windowName = "PBD Test 1";
	glutSsettings.GLVersionMajor = 3;
	glutSsettings.GLVersionMinor = 0;
	glutSsettings.positionX = 100;
	glutSsettings.positionY = 100;

	GLUTHelper helper;
	helper.initWindow(argc, argv, glutSsettings);
	helper.initCamera(10, 10, 10, 0, -6, 0);

	helper.setIdleFunc(&idleLoopGLUT);

	helper.enterDisplayLoop(&mainLoopGLUT);
}