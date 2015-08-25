#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "CUDA_WRAPPER.h"
#include "CUDA_GLOBALS.h"

#define IDX (blockIdx.x * blockDim.x + threadIdx.x)
#define LOCAL_IDX threadIdx.x

__shared__ float F[NUM_THREADS_PER_BLOCK][3][3];

__device__ float sqr(float x)
{
	return x * x;
}

__device__ float determinantF(int idx)
{
	return F[idx][0][0]
		* (F[idx][1][1] * F[idx][2][2] - F[idx][1][2] * F[idx][2][1])
		- F[idx][0][1]
		* (F[idx][1][0] * F[idx][2][2] - F[idx][1][2] * F[idx][2][0])
		+ F[idx][0][2]
		* (F[idx][1][0] * F[idx][2][1] - F[idx][1][1] * F[idx][2][0]);
}

__device__ void calculateF(int globalIdx, int idx, float* positions, float* refShapeMatrixInverse,
	int* indices, int trueNumConstraints, int numParticles)
{
	int localIndices[4];
	localIndices[0] = indices[globalIdx + trueNumConstraints * 0] * 3;
	localIndices[1] = indices[globalIdx + trueNumConstraints * 1] * 3;
	localIndices[2] = indices[globalIdx + trueNumConstraints * 2] * 3;
	localIndices[3] = indices[globalIdx + trueNumConstraints * 3] * 3;
	//localIndices[0] = indices[globalIdx + trueNumConstraints * 0];
	//localIndices[1] = indices[globalIdx + trueNumConstraints * 1];
	//localIndices[2] = indices[globalIdx + trueNumConstraints * 2];
	//localIndices[3] = indices[globalIdx + trueNumConstraints * 3];

	float temp[3][3];

	//float node3[3];
	//node3[0] = positions[localIndices[3] + 0];
	//node3[1] = positions[localIndices[3] + 1];
	//node3[2] = positions[localIndices[3] + 2];

	//1. Calculate Deformed Shape Matrix
	temp[0][0] = positions[localIndices[0] + 0] - positions[localIndices[3] + 0];
	temp[1][0] = positions[localIndices[0] + 1] - positions[localIndices[3] + 1];
	temp[2][0] = positions[localIndices[0] + 2] - positions[localIndices[3] + 2];
	temp[0][1] = positions[localIndices[1] + 0] - positions[localIndices[3] + 0];
	temp[1][1] = positions[localIndices[1] + 1] - positions[localIndices[3] + 1];
	temp[2][1] = positions[localIndices[1] + 2] - positions[localIndices[3] + 2];
	temp[0][2] = positions[localIndices[2] + 0] - positions[localIndices[3] + 0];
	temp[1][2] = positions[localIndices[2] + 1] - positions[localIndices[3] + 1];
	temp[2][2] = positions[localIndices[2] + 2] - positions[localIndices[3] + 2];

	//2. Multiply 
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				sum += temp[row][i] * refShapeMatrixInverse[(i * 3 + col) * trueNumConstraints + globalIdx];
			}

			F[idx][row][col] = sum;
		}
	}
}


__device__ __forceinline__ float calculateStrainEnergy_NEO_HOOKEAN(float volume, float lambda, float mu, float I1, float I3)
{
	return volume * (0.5f * mu * (I1 - log(I3) - 3.0f) + (lambda / 8.0f) * (log(I3) * log(I3)));
}

__device__ void updatePositions_recomputeGradients(int globalIdx, int idx, float lagrangeMultiplier, float* positions,
	float* masses, int* indices, int trueNumConstraints, int numParticles, float volume, float* refShapeMatrixInverse,
	float I3, float lambda, float mu, float det)
{
	//1. Copy refShapeMatrixInverse from global memory
	//for (int row = 0; row < 3; ++row)
	//{
	//	for (int col = 0; col < 3; ++col)
	//	{
	//		//We need the TRANSPOSE of the reference shape matrix inverse
	//		Gradient[idx][row][col] = refShapeMatrixInverse[(col * 3 + row) * trueNumConstraints + globalIdx];
	//	}
	//}

	float temp[3][3];

	//3. Multiply with First Piola-Kirchoff Stress tensor
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				float FInverseTEntry;

				if (row == 0)
				{
					if (i == 0)
					{
						FInverseTEntry = F[idx][1][1] * F[idx][2][2] - F[idx][2][1] * F[idx][1][2];
					}
					else if (i == 1)
					{
						FInverseTEntry = -(F[idx][1][0] * F[idx][2][2] - F[idx][2][0] * F[idx][1][2]);
					}
					else if (i == 2)
					{
						FInverseTEntry = F[idx][1][0] * F[idx][2][1] - F[idx][2][0] * F[idx][1][1];
					}
				}
				else if (row == 1)
				{
					if (i == 0)
					{
						FInverseTEntry = -(F[idx][0][1] * F[idx][2][2] - F[idx][2][1] * F[idx][0][2]);
					}
					else if (i == 1)
					{
						FInverseTEntry = F[idx][0][0] * F[idx][2][2] - F[idx][2][0] * F[idx][0][2];
					}
					else if (i == 2)
					{
						FInverseTEntry = -(F[idx][0][0] * F[idx][2][1] - F[idx][2][0] * F[idx][0][1]);
					}
				}
				else if (row == 2)
				{
					if (i == 0)
					{
						FInverseTEntry = F[idx][0][1] * F[idx][1][2] - F[idx][1][1] * F[idx][0][2];
					}
					else if (i == 1)
					{
						FInverseTEntry = -(F[idx][0][0] * F[idx][1][2] - F[idx][1][0] * F[idx][0][2]);
					}
					else if (i == 2)
					{
						FInverseTEntry = F[idx][0][0] * F[idx][1][1] - F[idx][1][0] * F[idx][0][1];
					}
				}

				FInverseTEntry /= det;

				float PFEntry = F[idx][row][i] * mu - (FInverseTEntry * mu) + FInverseTEntry * ((lambda * log(I3)) / 2.0f);

				sum += PFEntry * refShapeMatrixInverse[(col * 3 + i) * trueNumConstraints + globalIdx];
			}

			temp[row][col] = sum;
		}
	}

	int localIndices[4];
	localIndices[0] = indices[globalIdx + trueNumConstraints * 0] * 3;
	localIndices[1] = indices[globalIdx + trueNumConstraints * 1] * 3;
	localIndices[2] = indices[globalIdx + trueNumConstraints * 2] * 3;
	localIndices[3] = indices[globalIdx + trueNumConstraints * 3] * 3;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			atomicAdd(&positions[localIndices[i] + j], masses[globalIdx + i * trueNumConstraints] * lagrangeMultiplier * temp[j][i] * volume);
		}
	}

	for (int j = 0; j < 3; ++j)
	{
		float sum = 0.0f;
		for (int col = 0; col < 3; ++col)
		{
			sum += temp[j][col] * volume;
		}

		atomicAdd(&positions[localIndices[3] + j], masses[globalIdx + 3 * trueNumConstraints] * lagrangeMultiplier * -sum);
	}
}


__device__ float calculateTraceFTransposeF_INPLACE()
{
	float trace = 0.0f;
	for (int diagIdx = 0; diagIdx < 3; ++diagIdx)
	{
		for (int i = 0; i < 3; ++i)
		{
			trace += F[threadIdx.x][i][diagIdx] * F[threadIdx.x][i][diagIdx];
		}
	}

	return trace;
}

__device__ __forceinline__ float getFtFEntry(int row, int col)
{
	return F[threadIdx.x][0][row] * F[threadIdx.x][0][col]
		+ F[threadIdx.x][1][row] * F[threadIdx.x][1][col]
		+ F[threadIdx.x][2][row] * F[threadIdx.x][2][col];
}

__device__ __forceinline__ float calculatedeterminantFTransposeF_INPLACE()
{
	return (getFtFEntry(0, 0)
		* (getFtFEntry(1, 1) * getFtFEntry(2, 2) - getFtFEntry(1, 2) * getFtFEntry(2, 1)))
		- (getFtFEntry(0, 1)
		* (getFtFEntry(1, 0) * getFtFEntry(2, 2) - getFtFEntry(1, 2) * getFtFEntry(2, 0)))
		+ (getFtFEntry(0, 2)
		* (getFtFEntry(1, 0) * getFtFEntry(2, 1) - getFtFEntry(1, 1) * getFtFEntry(2, 0)));
}

__device__ void calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE(int globalIdx, int idx, float volume, float* refShapeMatrixInverse, int trueNumConstraints, float mu, float lambda, float I3,
	float& snGr0, float& snGr1, float& snGr2, float& snGr3, float det)
{
	//1. Copy refShapeMatrixInverse from global memory
	//for (int row = 0; row < 3; ++row)
	//{
	//	for (int col = 0; col < 3; ++col)
	//	{
	//		//We need the TRANSPOSE of the reference shape matrix inverse
	//		Gradient[idx][row][col] = refShapeMatrixInverse[(col * 3 + row) * trueNumConstraints + globalIdx];
	//	}
	//}

	float temp[3][3];

	//3. Multiply with First Piola-Kirchoff Stress tensor
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				float FInverseTEntry;

				if (row == 0)
				{
					if (i == 0)
					{
						FInverseTEntry = F[idx][1][1] * F[idx][2][2] - F[idx][2][1] * F[idx][1][2];
					}
					else if (i == 1)
					{
						FInverseTEntry = -(F[idx][1][0] * F[idx][2][2] - F[idx][2][0] * F[idx][1][2]);
					}
					else if (i == 2)
					{
						FInverseTEntry = F[idx][1][0] * F[idx][2][1] - F[idx][2][0] * F[idx][1][1];
					}
				}
				else if (row == 1)
				{
					if (i == 0)
					{
						FInverseTEntry = -(F[idx][0][1] * F[idx][2][2] - F[idx][2][1] * F[idx][0][2]);
					}
					else if (i == 1)
					{
						FInverseTEntry = F[idx][0][0] * F[idx][2][2] - F[idx][2][0] * F[idx][0][2];
					}
					else if (i == 2)
					{
						FInverseTEntry = -(F[idx][0][0] * F[idx][2][1] - F[idx][2][0] * F[idx][0][1]);
					}
				}
				else if (row == 2)
				{
					if (i == 0)
					{
						FInverseTEntry = F[idx][0][1] * F[idx][1][2] - F[idx][1][1] * F[idx][0][2];
					}
					else if (i == 1)
					{
						FInverseTEntry = -(F[idx][0][0] * F[idx][1][2] - F[idx][1][0] * F[idx][0][2]);
					}
					else if (i == 2)
					{
						FInverseTEntry = F[idx][0][0] * F[idx][1][1] - F[idx][1][0] * F[idx][0][1];
					}
				}

				FInverseTEntry /= det;

				float PFEntry = F[idx][row][i] * mu - (FInverseTEntry * mu) + FInverseTEntry * ((lambda * log(I3)) / 2.0f);

				sum += PFEntry * refShapeMatrixInverse[(col * 3 + i) * trueNumConstraints + globalIdx];
			}

			temp[row][col] = sum;
		}
	}

	//4. Copy back
	snGr0 = 0.0f;
	for (int i = 0; i < 3; ++i)
	{
		snGr0 += sqr(temp[i][0] * volume);
	}
	snGr0 = sqrtf(snGr0);

	snGr1 = 0.0f;
	for (int i = 0; i < 3; ++i)
	{
		snGr1 += sqr(temp[i][1] * volume);
	}
	snGr1 = sqrtf(snGr1);

	snGr2 = 0.0f;
	for (int i = 0; i < 3; ++i)
	{
		snGr2 += sqr(temp[i][2] * volume);
	}
	snGr2 = sqrtf(snGr2);

	//4. Calculate last column
	snGr3 = 0.0f;
	for (int i = 0; i < 3; ++i)
	{
		float sum = 0.0f;
		for (int col = 0; col < 3; ++col)
		{
			sum += temp[i][col] * volume;
		}
		snGr3 += sqr(sum);
	}
	snGr3 = sqrtf(snGr3);
}


__global__ void solveFEMConstraint(float* positions, int* indices, float* inverseMass, float* volume, float* refShapeMatrixInverse,
	float lambda, float mu, int trueNumConstraints, int numParticles)
{
	if (IDX > trueNumConstraints)
	{
		return;
	}

	//1. Calculate Deformation Gradient F
	calculateF(IDX, threadIdx.x, positions, refShapeMatrixInverse, indices, trueNumConstraints, numParticles);

	//3. Compute Invariants

	float I3 = calculatedeterminantFTransposeF_INPLACE();
	float I1 = calculateTraceFTransposeF_INPLACE();

	//6. Calculate Strain Energy Gradient
	float snGr0;
	float snGr1;
	float snGr2;
	float snGr3;

	float det = determinantF(threadIdx.x);

	calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE(IDX, threadIdx.x, volume[IDX], refShapeMatrixInverse, trueNumConstraints, mu, lambda, I3,
		snGr0, snGr1, snGr2, snGr3, det);

	//7. Calculate Lagrange Multiplier
	float denominator = 0.0f;
	denominator += inverseMass[IDX + 0 * trueNumConstraints] * snGr0;
	denominator += inverseMass[IDX + 1 * trueNumConstraints] * snGr1;
	denominator += inverseMass[IDX + 2 * trueNumConstraints] * snGr2;
	denominator += inverseMass[IDX + 3 * trueNumConstraints] * snGr3;

	float lagrangeMultiplier = -(calculateStrainEnergy_NEO_HOOKEAN(volume[IDX], lambda, mu, I1, I3) / denominator);

	if (isnan(lagrangeMultiplier))
	{
		return;
	}

	//8. Update Positions
	updatePositions_recomputeGradients(IDX, threadIdx.x, lagrangeMultiplier, positions, inverseMass, indices, trueNumConstraints, numParticles,
		volume[IDX], refShapeMatrixInverse, I3, lambda, mu, det);
}

__global__ void computeDiagonalF(float* positions, int* indices, float* F, float* U, float* V, float* refShapeMatrixInverse, int trueNumConstraints)
{
	if (IDX > trueNumConstraints)
	{
		return;
	}

	float F_functionLevel[3][3];

	//1. COMPUTE F--------------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------------------------------
	int localIndices[4];
	localIndices[0] = indices[IDX + trueNumConstraints * 0] * 3;
	localIndices[1] = indices[IDX + trueNumConstraints * 1] * 3;
	localIndices[2] = indices[IDX + trueNumConstraints * 2] * 3;
	localIndices[3] = indices[IDX + trueNumConstraints * 3] * 3;

	float temp[3][3];

	//1. Calculate Deformed Shape Matrix
	temp[0][0] = positions[localIndices[0] + 0] - positions[localIndices[3] + 0];
	temp[1][0] = positions[localIndices[0] + 1] - positions[localIndices[3] + 1];
	temp[2][0] = positions[localIndices[0] + 2] - positions[localIndices[3] + 2];
	temp[0][1] = positions[localIndices[1] + 0] - positions[localIndices[3] + 0];
	temp[1][1] = positions[localIndices[1] + 1] - positions[localIndices[3] + 1];
	temp[2][1] = positions[localIndices[1] + 2] - positions[localIndices[3] + 2];
	temp[0][2] = positions[localIndices[2] + 0] - positions[localIndices[3] + 0];
	temp[1][2] = positions[localIndices[2] + 1] - positions[localIndices[3] + 1];
	temp[2][2] = positions[localIndices[2] + 2] - positions[localIndices[3] + 2];

	//2. Multiply 
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				sum += temp[row][i] * refShapeMatrixInverse[(i * 3 + col) * trueNumConstraints + IDX];
			}

			F_functionLevel[row][col] = sum;
		}
	}

	//1. DIAGONALISE F--------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------------

}


cudaError_t projectConstraints(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	const Parameters& settings);

int CUDA_projectConstraints(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	const Parameters& settings)
{
	//GPU
	cudaError_t cudaStatus = projectConstraints(device_indices, device_positions,
		device_inverseMasses, device_refShapeMatrixInverses, device_volumes, settings);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Critical Error, aborting...");
		return 1;
	}

	return 0;
}

cudaError_t cudaErrorWrapper(cudaError_t status)
{
	if (status != cudaSuccess)
	{
		std::cerr << "ERROR: " << cudaGetErrorString(status) << std::endl;
	}

	return status;
}

bool checkCudaErrorStatus(cudaError_t status)
{
	if (status != cudaSuccess)
	{
		std::cerr << "ERROR: " << cudaGetErrorString(status) << std::endl;
		return false;
	}
	return true;
}

void getCudaDeviceProperties(int device)
{
	cudaDeviceProp properties;
	cudaGetDeviceProperties(&properties, device);

	std::cout << "Compute Capabilities for " << properties.name << " : " << std::endl;
	std::cout << "Major: " << properties.major << ", Minor: " << properties.minor << std::endl;
	std::cout << "Details: " << std::endl;
	std::cout << "	Num of SM    : " << properties.multiProcessorCount << std::endl;
	std::cout << "	Mem per Block: " << properties.sharedMemPerBlock << std::endl;
	std::cout << "	Mem per SM   : " << properties.sharedMemPerMultiprocessor << std::endl;
}

void queryCUDADevices()
{
	cudaError_t deviceStatus;

	int deviceCount = 0;
	deviceStatus = cudaGetDeviceCount(&deviceCount);

	std::cout << "Num CUDA Devices Found: " << deviceCount << std::endl;
	deviceStatus = cudaErrorWrapper(cudaSetDevice(0));
	checkCudaErrorStatus(deviceStatus);
	getCudaDeviceProperties(0);
}

cudaError_t projectConstraints(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	const Parameters& settings)
{
	float* dev_positions = 0;
	float* dev_inverseMasses = 0;
	int* dev_indices = 0;
	float* dev_refShapeMatrixInverses = 0;
	float* dev_volumes = 0;

	cudaError_t deviceStatus;

	//Execute Kernel
	//std::cout << "Executing Kernel..." << settings.numConstraintIterations<< std::endl;
	//std::cout << settings.numBlocks * settings.numThreadsPerBlock << "threads..." << std::endl;
	//std::cout << settings.trueNumberOfConstraints << "true num of constraints..." << std::endl;

	dim3 numBlocks;
	dim3 numThreads;
	numBlocks.x = settings.numBlocks;
	//numBlocks.x = 1;
	numBlocks.y = 1;
	numBlocks.z = 1;

	numThreads.x = settings.numThreadsPerBlock;
	numThreads.y = 1;
	numThreads.z = 1;

	cudaEvent_t start;
	cudaEvent_t end;
	cudaEventCreate(&start);
	cudaEventCreate(&end);

	cudaEventRecord(start);
	for (int it = 0; it < settings.numConstraintIterations; ++it)
	{
		solveFEMConstraint << <numBlocks, numThreads >> >(
			device_positions, device_indices, device_inverseMasses,
			device_volumes, device_refShapeMatrixInverses,
			settings.lambda, settings.mu, settings.trueNumberOfConstraints,
			settings.numParticles);

		cudaErrorWrapper(cudaDeviceSynchronize());
	}
	cudaErrorWrapper(cudaDeviceSynchronize());

	cudaEventRecord(end);
	cudaEventSynchronize(end);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, end);
	std::cout << "Execution Time: " << milliseconds / 1000.0 << "s." << std::endl;

	deviceStatus = cudaGetLastError();
	checkCudaErrorStatus(deviceStatus);

	return deviceStatus;
}

bool CUDA_allocateBuffers(int** device_indices, float** device_positions,
	float** device_inverseMasses, float** device_refShapeMatrixInverses,
	float** device_volumes,
	std::vector<int>& indices,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes)
{
	std::cout << "Allocating CUDA Buffers" << std::endl;
	cudaError_t deviceStatus;

	deviceStatus = cudaSetDevice(0);
	deviceStatus = cudaGetLastError();
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}

	deviceStatus = cudaMalloc((void**)device_indices, indices.size() * sizeof(int));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_positions, positions.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_inverseMasses, inverseMasses.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_refShapeMatrixInverses, refShapeMatrixInverses.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_volumes, volumes.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}

	std::cout << "Copying Memory to CUDA device..." << std::endl;
	deviceStatus = cudaMemcpy(*device_indices, &indices[0], indices.size() * sizeof(int), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_positions, &positions[0], positions.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_inverseMasses, &inverseMasses[0], inverseMasses.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_refShapeMatrixInverses, &refShapeMatrixInverses[0], refShapeMatrixInverses.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_volumes, &volumes[0], volumes.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}


	//Print some statistics
	double totalNumbytes = indices.size() * sizeof(int)
		+positions.size() * sizeof(float)
		+volumes.size() * sizeof(float)
		+inverseMasses.size() * sizeof(float)
		+refShapeMatrixInverses.size() * sizeof(float);

	std::cout << "Memory Usage: " << std::endl;
	std::cout << "	Indices (int32)             : " << indices.size() * sizeof(int) << " bytes" << std::endl;
	std::cout << "	Positions (float)           : " << positions.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Volumes (float)             : " << volumes.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Masses (float)              : " << inverseMasses.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Ref. Shape Matrices: (float): " << refShapeMatrixInverses.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	--------------------------------------------------------" << std::endl;
	std::cout << "	Total                       : " << totalNumbytes << " bytes (" << totalNumbytes / 1000.0 << " kb)" << std::endl;
}

bool CUDA_destroyBuffers(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes)
{
	cudaError_t deviceStatus;

	cudaFree(device_positions);
	cudaFree(device_inverseMasses);
	cudaFree(device_indices);
	cudaFree(device_refShapeMatrixInverses);
	cudaFree(device_volumes);

	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	deviceStatus = cudaDeviceReset();
	if (deviceStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}


	deviceStatus = cudaGetLastError();
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool CUDA_updateBuffers(float* device_positions, std::vector<float>& positions)
{
	cudaError_t deviceStatus;

	deviceStatus = cudaMemcpy(device_positions, &positions[0], positions.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}

	return true;
}

bool CUDA_getBuffers(float* device_positions, std::vector<float>& positions)
{
	cudaError_t deviceStatus;

	deviceStatus = cudaErrorWrapper(cudaMemcpy(&positions[0], device_positions, positions.size() * sizeof(float), cudaMemcpyDeviceToHost));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}

	return true;
}