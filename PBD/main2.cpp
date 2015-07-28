#include <vector>
#include <iostream>

#include "CUDAPBD_SolverSettings.h"
#include "CUDA_WRAPPER.h"

int main()
{
	std::vector<int> indices(4260);
	std::vector<float> positions(426);
	std::vector<float> inverseMasses(426);
	std::vector<float> refShapeMatrixInverses(9585);
	std::vector<float> volumes(1065);
	std::vector<float> positions_result(426);

	CUDAPBD_SolverSettings settings;
	settings.lambda = 11;
	settings.mu = 13;
	settings.numBlocks = 65;
	settings.numThreadsPerBlock = 64;
	settings.trueNumberOfConstraints = 1065;

	queryCUDADevices();

	CUDA_projectConstraints(indices, positions, inverseMasses, refShapeMatrixInverses,
		volumes, settings);

	std::cout << "done!" << std::endl;
}
