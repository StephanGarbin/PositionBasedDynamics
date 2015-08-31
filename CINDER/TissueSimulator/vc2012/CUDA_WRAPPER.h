#include <vector>

#include "Parameters.h"

int CUDA_projectConstraints(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	float* device_F, float* device_U, float* device_V,
	float* device_anisotropyDirection, float* device_viscosityMultiplier,
	const Parameters& settings);

void queryCUDADevices();


bool CUDA_allocateBuffers(int** device_indices, float** device_positions,
	float** device_inverseMasses, float** device_refShapeMatrixInverses,
	float** device_volumes,
	float** device_F, float** device_U, float** device_V,
	float** device_anisotropyDirection, float** device_viscosityMultiplier,
	std::vector<int>& indices,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes,
	std::vector<float>& F,
	std::vector<float>& U,
	std::vector<float>& V,
	std::vector<float>& anisotropyDirection,
	std::vector<float>& viscosityMultiplier);

bool CUDA_destroyBuffers(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	float* device_F, float* device_U, float* device_V,
	float* device_anisotropyDirection, float* device_viscosityMultiplier);

bool CUDA_updateBuffers(float* device_positions, std::vector<float>& positions,
	float* device_anisotropyDirection, std::vector<float>& anisotropyDirection, bool updateAnisotropyDirection);

bool CUDA_getBuffers(float* device_positions, std::vector<float>& positions);