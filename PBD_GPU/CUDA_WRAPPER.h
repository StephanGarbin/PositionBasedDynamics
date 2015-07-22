#include <vector>

#include "CUDAPBD_SolverSettings.h"

int CUDA_projectConstraints(std::vector<int> indices,
	std::vector<float> positions,
	std::vector<float> inverseMasses,
	std::vector<float> refShapeMatrixInverses,
	std::vector<float> volumes,
	const CUDAPBD_SolverSettings& settings);