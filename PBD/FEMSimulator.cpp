#include "FEMSimulator.h"

#include <iostream>

#include "volumetricMeshLoader.h"
#include "homogeneousNeoHookeanIsotropicMaterial.h"
#include "isotropicHyperelasticFEM.h"
#include "isotropicHyperelasticFEMForceModel.h"
#include "generateMassMatrix.h"
#include "implicitBackwardEulerSparse.h"

struct FEMSimulatorData
{
	std::vector<int>* m_pointConstraints;
	TetMesh* m_mesh;
	HomogeneousNeoHookeanIsotropicMaterial* m_material;
	IsotropicHyperelasticFEM* m_model;
	ForceModel* m_forceModel;
	SparseMatrix* m_massMatrix;
	ImplicitBackwardEulerSparse* m_implicitIntegrator;

	int m_DOF;
	int m_positiveDefiniteSolver = 0;
	double m_timeStep;
	double m_dampingMassCoeff = 0.0;
	double m_dampingStiffnessCoeff = 0.01;

	~FEMSimulatorData()
	{
		delete m_pointConstraints;
	}
};



FEMSimulator::FEMSimulator()
{
	m_data = new FEMSimulatorData();
}


FEMSimulator::~FEMSimulator()
{
	delete m_data;
}

void
FEMSimulator::initSolver(const std::string& meshfileName,
std::vector<int>& pointConstraints,
double youngsModulus, double poissonRatio, double timeStep)
{
	m_meshfileName = meshfileName;

	m_data->m_pointConstraints = new std::vector<int>();
	for (int i = 0; i < pointConstraints.size(); ++i)
	{
		m_data->m_pointConstraints->push_back(pointConstraints[i] * 3);
		m_data->m_pointConstraints->push_back(pointConstraints[i] * 3 + 1);
		m_data->m_pointConstraints->push_back(pointConstraints[i] * 3 + 2);
	}

	m_data->m_mesh = static_cast<TetMesh*>(VolumetricMeshLoader::load(m_meshfileName.c_str()));

	m_data->m_mesh->orient();

	m_data->m_material = new HomogeneousNeoHookeanIsotropicMaterial(youngsModulus, poissonRatio);

	m_data->m_model = new IsotropicHyperelasticFEM(m_data->m_mesh, m_data->m_material,-1.7e+308, true);
	//m_data->m_model = new IsotropicHyperelasticFEM(m_data->m_mesh, m_data->m_material);

	m_data->m_forceModel = new IsotropicHyperelasticFEMForceModel(m_data->m_model);

	m_data->m_DOF = 3 * m_data->m_mesh->getNumVertices();

	m_data->m_timeStep = timeStep;

	GenerateMassMatrix::computeMassMatrix(m_data->m_mesh, &m_data->m_massMatrix, true);

	m_data->m_positiveDefiniteSolver = 0;

	m_data->m_dampingMassCoeff = 0.0;
	m_data->m_dampingStiffnessCoeff = 0.0;

	//std::cout << "Initialising Integrator..." << std::endl;

	m_data->m_implicitIntegrator = new ImplicitBackwardEulerSparse(m_data->m_DOF,
		m_data->m_timeStep, m_data->m_massMatrix, m_data->m_forceModel,
		m_data->m_positiveDefiniteSolver, m_data->m_pointConstraints->size(),&(*m_data->m_pointConstraints)[0] , m_data->m_dampingMassCoeff,
		m_data->m_dampingStiffnessCoeff);
	
	//m_data->m_implicitIntegrator = new ImplicitBackwardEulerSparse(m_data->m_DOF,
	//	m_data->m_timeStep, m_data->m_massMatrix, m_data->m_forceModel,
	//	m_data->m_positiveDefiniteSolver);


	//finally, initialise buffers
	m_displacements.resize(m_data->m_DOF);
	m_externalForces.resize(m_data->m_DOF);

	m_data->m_implicitIntegrator->ResetToRest();

	//and start at iteration 0
	m_numIts = 0;
}


std::vector<double>&
FEMSimulator::getCurrentDisplacements()
{
	m_data->m_implicitIntegrator->GetqState(&m_displacements[0]);
	return m_displacements;
}

std::vector<double>&
FEMSimulator::getExternalForces()
{
	return m_externalForces;
}

void
FEMSimulator::doTimeStep(bool applyExternalForces)
{
	if (applyExternalForces)
	{
		m_data->m_implicitIntegrator->SetExternalForcesToZero();
		//m_data->m_implicitIntegrator->SetExternalForces(&m_externalForces[0]);
	}

	//std::cout << "System Time: " << m_data->m_implicitIntegrator->GetSystemSolveTime() << std::endl;
	m_data->m_implicitIntegrator->DoTimestep();
}
