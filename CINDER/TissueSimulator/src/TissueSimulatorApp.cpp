#include <vector>
#include <memory>
#include <string>

#include <windows.h>
#include <stdio.h>
#include <fcntl.h>
#include <io.h>
#include <iostream>
#include <fstream>

#include "cinder\app\AppNative.h"
#include "cinder\gl\gl.h"
#include "cinder\params\Params.h"
#include "cinder\Camera.h"
#include "cinder\Matrix.h"
#include "cinder\Utilities.h"

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"
#include "Parameters.h"

#include "ConstraintsIO.h"
#include "TetGenIO.h"

#include "PBDGPU_Solver.h"
#include "CUDAMemoryOptimiser.h"

#include "CUDA_GLOBALS.h"

using namespace ci;
using namespace ci::app;
using namespace std;

//! Simple function from the canonical example (http://www.halcyon.com/~ast/dload/guicon.htm)
void redirectCOUT2Console()
{
	int hConHandle;
	long lStdHandle;

	CONSOLE_SCREEN_BUFFER_INFO coninfo;

	FILE *fp;
	// allocate a console for this app
	AllocConsole();

	// set the screen buffer to be big enough to let us scroll text

	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE),

		&coninfo);

	coninfo.dwSize.Y = 2000;

	SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE),

		coninfo.dwSize);

	// redirect unbuffered STDOUT to the console

	lStdHandle = (long)GetStdHandle(STD_OUTPUT_HANDLE);

	hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);

	fp = _fdopen(hConHandle, "w");

	*stdout = *fp;

	setvbuf(stdout, NULL, _IONBF, 0);

	// redirect unbuffered STDIN to the console

	lStdHandle = (long)GetStdHandle(STD_INPUT_HANDLE);

	hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);

	fp = _fdopen(hConHandle, "r");

	*stdin = *fp;

	setvbuf(stdin, NULL, _IONBF, 0);

	// redirect unbuffered STDERR to the console

	lStdHandle = (long)GetStdHandle(STD_ERROR_HANDLE);

	hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);

	fp = _fdopen(hConHandle, "w");

	*stderr = *fp;

	setvbuf(stderr, NULL, _IONBF, 0);

	// make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog

	// point to console as well

	ios::sync_with_stdio();
}

struct simulatorData
{
	simulatorData()
	{
		particles = std::make_shared<std::vector<PBDParticle>>();
	}

	//ACCESSORS
	std::vector<PBDTetrahedra3d>& getTets()
	{
		return tetrahedra;
	}

	std::vector<PBDParticle>& getParticles()
	{
		return *particles;
	}

	std::vector<int>& getPositionConstraints()
	{
		return positionConstrainedParticles;
	}
	
	//DATA
	std::vector<PBDTetrahedra3d> tetrahedra;
	std::shared_ptr<std::vector<PBDParticle>> particles;
	std::vector<int> positionConstrainedParticles;
};

class TissueSimulatorApp : public AppNative {
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();

private:
	//INIT FUNCTIONS
	void setupInterface();
	void handleIO();

	//SIMULATION DATA
	simulatorData m_data;

	//SOLVERs
	PBDGPU_Solver m_solverGPU;

	//INTERFACE
	params::InterfaceGlRef m_interfaceParams;
	Parameters m_params;

	//RENDERING
	CameraPersp m_camera;
};

void TissueSimulatorApp::setup()
{
	//Create Console & Redirect Output Streams
	redirectCOUT2Console();

	//0. Setup Interface
	setupInterface();

	//1. Do IO
	handleIO();

	//2. Setup GPU Solver
	m_solverGPU.setup(m_data.getTets(), m_data.particles);

	//3. Setup Camera
	m_camera.lookAt(Vec3f(10.0, 10.0, 10.0), Vec3f(0.0, 0.0, 0.0), Vec3f(0.0, 1.0, 0.0));

	//Finally, draw one
	draw();
}

void TissueSimulatorApp::mouseDown( MouseEvent event )
{
}

void TissueSimulatorApp::update()
{
	m_solverGPU.advanceSystem(m_data.particles, m_params);
	//sleep(1000);
}

void TissueSimulatorApp::draw()
{
	gl::enableDepthRead();
	gl::enableDepthWrite();

	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) ); 

	gl::setMatrices(m_camera);

	gl::enableWireframe();

	for (int i = 0; i < m_data.getTets().size(); ++i)
	{
		m_data.getTets()[i].glRender(1.0, 1.0, 1.0);
	}

	m_interfaceParams->draw();
}


void
TissueSimulatorApp::setupInterface()
{
	m_params.initialiseToDefaults();

	m_interfaceParams = params::InterfaceGl::create(getWindow(), "App parameters", toPixels(Vec2i(200, 300)));

	m_interfaceParams->addParam("Material", m_params.materialModelNames, &m_params.materialModel);

	//Elasticity
	m_interfaceParams->addSeparator("Elasticity");
	m_interfaceParams->addParam("Young's Modulus", &m_params.youngsModulus).min(0.0f).max(100000.0f).precision(4).step(0.25f);
	m_interfaceParams->addParam("Poisson's Ratio", &m_params.poissonRatio).min(-0.5f).max(0.5f).precision(4).step(0.1f);

	//Ansisotropy
	m_interfaceParams->addSeparator("Anisotropy");
	m_interfaceParams->addParam("Anisotropic Strength", &m_params.anisotropyStrength).min(0.0f).max(10.0f).precision(4).step(0.1f);
	m_interfaceParams->addParam("Anisotropy Direction", &m_params.anisotropyDirection);
	m_params.normaliseAnisotropyDirection();

	//Viscosity
	m_interfaceParams->addSeparator("Viscosity");
	m_interfaceParams->addParam("Alpha", &m_params.alpha);
	m_interfaceParams->addParam("Rho", &m_params.rho);

	m_interfaceParams->addParam("Inverse Mass", &m_params.inverseMass).min(0.0f).max(100000.0f).precision(4).step(1.0f);
	m_interfaceParams->addSeparator();
	m_interfaceParams->addParam("Timestep", &m_params.timeStep).min(0.000001f).max(1.0f).precision(4).step(0.05f);
	m_interfaceParams->addParam("Num Constraint Iterations", &m_params.numConstraintIterations).min(2).max(1000).precision(1).step(1);
	m_interfaceParams->addParam("Num GPU Block Iterations", &m_params.numGPUBlockIterations).min(1).max(1000).precision(1).step(1);
}

void TissueSimulatorApp::handleIO()
{
	Eigen::Vector3f initialVelocity;
	initialVelocity.setZero();

	TetGenIO::readNodes(m_params.tetGenNodeFile, m_data.getParticles(), m_params.inverseMass, initialVelocity);
	TetGenIO::readTetrahedra(m_params.tetGenElementFile, m_data.getTets(), m_data.particles);

	ConstraintsIO::readMayaVertexConstraints(m_data.getPositionConstraints(), m_params.positionConstraintFile);

	//Set inverse mass of constrained nodes to 0
	for (int i = 0; i < m_data.getPositionConstraints().size(); ++i)
	{
		m_data.getParticles()[m_data.getPositionConstraints()[i]].inverseMass() = 0.0;
	}

	//OPTIMISE MEMORY LAYOUT FOR CUDA
	//CUDAMemoryOptimiser::optimiseTetrahedralIndexingBasedOnNodeMemory(m_data.getParticles(), m_data.getTets());
	//CUDAMemoryOptimiser::optimiseTetrahedralIndexingBasedOnMemoryChunks(m_data.getParticles(), m_data.getTets(), NUM_THREADS_PER_BLOCK);
}

CINDER_APP_NATIVE( TissueSimulatorApp, RendererGl )
