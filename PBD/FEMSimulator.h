#pragma once

#include <string>
#include <vector>

struct FEMSimulatorData;

class FEMSimulator
{
public:
	FEMSimulator();

	~FEMSimulator();

	void initSolver(const std::string& meshfileName, std::vector<int>& pointConstraints,
		double youngsModulus, double poissonRatio, double timeStep);

	void doTimeStep(bool applyExternalForces);

	std::vector<double>& getCurrentDisplacements();

	std::vector<double>& getExternalForces();

	int howManyTimeSteps() { return m_numIts; }

private:
	std::string m_meshfileName;
	FEMSimulatorData* m_data;
	std::vector<double> m_displacements;
	std::vector<double> m_externalForces;
	int m_numIts;
};

