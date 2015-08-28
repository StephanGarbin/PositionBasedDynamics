#include <vector>

#include "Parameters.h"

int CUDA_projectConstraints(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	float* device_F, float* device_U, float* device_V,
	const Parameters& settings);

void queryCUDADevices();


bool CUDA_allocateBuffers(int** device_indices, float** device_positions,
	float** device_inverseMasses, float** device_refShapeMatrixInverses,
	float** device_volumes,
	float** device_F, float** device_U, float** device_V,
	std::vector<int>& indices,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes,
	std::vector<float>& F,
	std::vector<float>& U,
	std::vector<float>& V);

bool CUDA_destroyBuffers(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	float* device_F, float* device_U, float* device_V);

bool CUDA_updateBuffers(float* device_positions, std::vector<float>& positions);

bool CUDA_getBuffers(float* device_positions, std::vector<float>& positions);