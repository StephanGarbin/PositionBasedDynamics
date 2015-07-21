
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

const int NUM_THREADS_PER_BLOCK_SINGLE = 8;
const int NUM_THREADS_PER_BLOCK = NUM_THREADS_PER_BLOCK_SINGLE * NUM_THREADS_PER_BLOCK_SINGLE;

__shared__ float F[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float FTransposeF[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float FInverseTranspose[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float FirstPiolaKirchoffTensor[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float Gradient[NUM_THREADS_PER_BLOCK][3][4];
__shared__ int LocalIndices[NUM_THREADS_PER_BLOCK][4];


__device__ float sqr(float x)
{
	return x * x;
}

__device__ float traceFTransposeF(int idx)
{
	return FTransposeF[idx][0][0] + FTransposeF[idx][1][1] + FTransposeF[idx][2][2];
}

__device__ float determinantFTransposeF(int idx)
{
	return FTransposeF[idx][0][0]
		* (FTransposeF[idx][1][1] * FTransposeF[idx][2][2] - FTransposeF[idx][1][2] * FTransposeF[idx][2][1])
		- FTransposeF[idx][0][1]
		* (FTransposeF[idx][1][0] * FTransposeF[idx][2][2] - FTransposeF[idx][1][2] * FTransposeF[idx][2][0])
		+ FTransposeF[idx][0][2]
		* (FTransposeF[idx][1][0] * FTransposeF[idx][2][1] - FTransposeF[idx][1][1] * FTransposeF[idx][2][0]);
}

__device__ float determinantFInverseTranspose(int idx)
{
	return FInverseTranspose[idx][0][0]
		* (FInverseTranspose[idx][1][1] * FInverseTranspose[idx][2][2] - FInverseTranspose[idx][1][2] * FInverseTranspose[idx][2][1])
		- FInverseTranspose[idx][0][1]
		* (FInverseTranspose[idx][1][0] * FInverseTranspose[idx][2][2] - FInverseTranspose[idx][1][2] * FInverseTranspose[idx][2][0])
		+ FInverseTranspose[idx][0][2]
		* (FInverseTranspose[idx][1][0] * FInverseTranspose[idx][2][1] - FInverseTranspose[idx][1][1] * FInverseTranspose[idx][2][0]);
}

__device__ void calculateF(int idx, float* positions, float* refShapeMatrixInverse)
{
	//1. Calculate Deformed Shape Matrix
	FirstPiolaKirchoffTensor[idx][0][0] = positions[LocalIndices[idx][0] * 3 + 0] - positions[LocalIndices[idx][3] * 3 + 0];
	FirstPiolaKirchoffTensor[idx][1][0] = positions[LocalIndices[idx][0] * 3 + 1] - positions[LocalIndices[idx][3] * 3 + 1];
	FirstPiolaKirchoffTensor[idx][2][0] = positions[LocalIndices[idx][0] * 3 + 2] - positions[LocalIndices[idx][3] * 3 + 2];

	FirstPiolaKirchoffTensor[idx][0][1] = positions[LocalIndices[idx][1] * 3 + 0] - positions[LocalIndices[idx][3] * 3 + 0];
	FirstPiolaKirchoffTensor[idx][1][1] = positions[LocalIndices[idx][1] * 3 + 1] - positions[LocalIndices[idx][3] * 3 + 1];
	FirstPiolaKirchoffTensor[idx][2][1] = positions[LocalIndices[idx][1] * 3 + 2] - positions[LocalIndices[idx][3] * 3 + 2];

	FirstPiolaKirchoffTensor[idx][0][2] = positions[LocalIndices[idx][2] * 3 + 0] - positions[LocalIndices[idx][3] * 3 + 0];
	FirstPiolaKirchoffTensor[idx][1][2] = positions[LocalIndices[idx][2] * 3 + 1] - positions[LocalIndices[idx][3] * 3 + 1];
	FirstPiolaKirchoffTensor[idx][2][2] = positions[LocalIndices[idx][2] * 3 + 2] - positions[LocalIndices[idx][3] * 3 + 2];

	//2. Multiply 
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				sum += FirstPiolaKirchoffTensor[idx][row][i] * refShapeMatrixInverse[idx * 3 * 3 + row * 3 + col];
			}

			F[idx][row][col] = sum;
		}
	}

}

__device__ void calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN(int idx, float mu, float lambda, float I3)
{
	//1. Copy over F
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FirstPiolaKirchoffTensor[idx][row][col] = F[idx][row][col];
		}
	}

	//2. Multiply with mu
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FirstPiolaKirchoffTensor[idx][row][col] *= mu;
		}
	}

	//3. Subtract mu times FInverseTranspose
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FirstPiolaKirchoffTensor[idx][row][col] -= mu * FInverseTranspose[idx][row][col];
		}
	}

	//4. Add (ambda * logI3) / 2.0 * FInverseTranspose
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FirstPiolaKirchoffTensor[idx][row][col] += ((lambda * log(I3)) / 2.0f) * FInverseTranspose[idx][row][col];
		}
	}
}

__device__ float calculateStrainEnergy_NEO_HOOKEAN(float volume, float lambda, float mu, float I1, float I3)
{
	return volume * (0.5f * mu * (I1 - log(I3) - 3.0f) + (lambda / 8.0f) * (log(I3) * log(I3)));
}

__device__ void calculateStrainEnergyGradient_NEO_HOOKEAN(int idx, float volume, float* refShapeMatrixInverse)
{
	//1. Copy refShapeMatrixInverse from global memory
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			Gradient[idx][row][col] = refShapeMatrixInverse[idx * 3 + row * 3 + col];
		}
	}

	//2. Multiply by volume
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			Gradient[idx][row][col] *= volume;
		}
	}

	//3. Multiply with First Piola-Kirchoff Stress tensor
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				sum += Gradient[idx][row][i] * FirstPiolaKirchoffTensor[idx][i][row];
			}

			FTransposeF[idx][row][col] = sum;
		}
	}

	//4. Copy back
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			Gradient[idx][row][col] = FTransposeF[idx][row][col];
		}
	}


	//4. Calculate last column
	for (int row = 0; row < 3; ++row)
	{
		float sum = 0.0f;
		for (int col = 0; col < 3; ++col)
		{
			sum += Gradient[idx][row][col];
		}
		Gradient[idx][row][3] = -sum;
	}
}

__device__ void calculateFTransposeF(int idx)
{
	//1. Copy over F (1&2 could be combined in future!)
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FTransposeF[idx][row][col] = F[idx][row][col];
		}
	}

	//2. Transpose F
	float temp;
	temp = FTransposeF[idx][0][1];
	FTransposeF[idx][0][1] = FTransposeF[idx][1][0];
	FTransposeF[idx][1][0] = temp;

	temp = FTransposeF[idx][0][2];
	FTransposeF[idx][0][2] = FTransposeF[idx][2][0];
	FTransposeF[idx][2][0] = temp;

	temp = FTransposeF[idx][1][2];
	FTransposeF[idx][1][2] = FTransposeF[idx][2][1];
	FTransposeF[idx][2][1] = temp;

	//3. Multiply with F
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				sum += FTransposeF[idx][row][i] * F[idx][i][row];
			}

			FirstPiolaKirchoffTensor[idx][row][col] = sum;
		}
	}

	//Copy back
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FTransposeF[idx][row][col] = FirstPiolaKirchoffTensor[idx][row][col];
		}
	}
}

__device__ void calculateFInverseTranspose(int idx)
{
	//1. Calculate cofactors
	FInverseTranspose[idx][0][0] = F[idx][1][1] * F[idx][2][2] - F[idx][2][1] * F[idx][1][2];
	FInverseTranspose[idx][0][1] = F[idx][1][0] * F[idx][2][2] - F[idx][2][0] * F[idx][1][2];
	FInverseTranspose[idx][0][2] = F[idx][1][0] * F[idx][2][1] - F[idx][2][0] * F[idx][1][1];

	FInverseTranspose[idx][1][0] = F[idx][0][1] * F[idx][2][2] - F[idx][2][1] * F[idx][0][2];
	FInverseTranspose[idx][1][1] = F[idx][0][0] * F[idx][2][2] - F[idx][2][0] * F[idx][0][2];
	FInverseTranspose[idx][1][2] = F[idx][0][0] * F[idx][2][1] - F[idx][2][0] * F[idx][0][1];

	FInverseTranspose[idx][2][0] = F[idx][1][1] * F[idx][1][2] - F[idx][2][1] * F[idx][0][2];
	FInverseTranspose[idx][2][1] = F[idx][1][0] * F[idx][1][2] - F[idx][2][0] * F[idx][0][2];
	FInverseTranspose[idx][2][2] = F[idx][1][0] * F[idx][1][1] - F[idx][2][0] * F[idx][0][1];

	//2. Transpose
	float temp;

	temp = FInverseTranspose[idx][0][1];
	FInverseTranspose[idx][0][1] = FInverseTranspose[idx][1][0];
	FInverseTranspose[idx][1][0] = temp;

	temp = FInverseTranspose[idx][0][2];
	FInverseTranspose[idx][0][2] = FInverseTranspose[idx][2][0];
	FInverseTranspose[idx][2][0] = temp;

	temp = FInverseTranspose[idx][1][2];
	FInverseTranspose[idx][1][2] = FInverseTranspose[idx][2][1];
	FInverseTranspose[idx][2][1] = temp;

	//3. Calculate the determinant
	float determinant = determinantFInverseTranspose(idx);

	//4. Multiply
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FInverseTranspose[idx][row][col] /= determinant;
		}
	}
}

__device__ float squaredNormGradient(int idx, int particleIdx)
{
	return sqrtf(sqr(Gradient[idx][0][particleIdx])
		+ sqr(Gradient[idx][1][particleIdx])
		+ sqr(Gradient[idx][2][particleIdx]));
}

__device__ float calculateLagrangeMultiplierDenominator(int idx, float* inverseMass)
{
	float denominator = 0.0f;
	for (int i = 0; i < 4; ++i)
	{
		denominator += inverseMass[idx + i] * squaredNormGradient(idx, i);
	}

	return denominator;
}

__device__ void updatePositions(int idx, float lagrangeMultiplier, float* positions, float* inverseMass)
{
	for (int i = 0; i < 4; ++i)
	{
		/*atomicAdd(&positions[LocalIndices[idx][i] * 3 + 0], inverseMass[idx] * lagrangeMultiplier * Gradient[idx][0][i]);
		atomicAdd(&positions[LocalIndices[idx][i] * 3 + 1], inverseMass[idx] * lagrangeMultiplier * Gradient[idx][1][i]);
		atomicAdd(&positions[LocalIndices[idx][i] * 3 + 2], inverseMass[idx] * lagrangeMultiplier * Gradient[idx][2][i]);*/
	}
}

__device__ void getIndices(int idx, int* indices)
{
	for (int i = 0; i < 4; ++i)
	{
		LocalIndices[idx][i] = indices[idx * 4 + i];
	}
}

__global__ void solveFEMConstraint(float* positions, int* indices, float* inverseMass, float* volume, float* refShapeMatrixInverse,
	float lambda, float mu)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//1. Calculate Deformation Gradient F
	calculateF(idx, positions, refShapeMatrixInverse);

	//2. Compute Cauchy Tensors
	calculateFInverseTranspose(idx);
	calculateFTransposeF(idx);

	//3. Compute Invariants
	float I1 = traceFTransposeF(idx);
	float I3 = determinantFTransposeF(idx);

	//4. Calculate First Piola-Kirchoff Stress Tensor
	calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN(idx, mu, lambda, I3);

	//5. Calculate StrainEnergy
	float strainEnergy = calculateStrainEnergy_NEO_HOOKEAN(volume[idx], lambda, mu, I1, I3);

	//6. Calculate Strain Energy Gradient
	calculateStrainEnergyGradient_NEO_HOOKEAN(idx, volume[idx], refShapeMatrixInverse);

	//7. Calculate Lagrange Multiplier
	float lagrangeMultiplier = - (strainEnergy / calculateLagrangeMultiplierDenominator(idx, inverseMass));

	//8. Update Positions
	updatePositions(idx, lagrangeMultiplier, positions, inverseMass);
}



cudaError_t projectConstraints(std::vector<int>& indices,
	std::vector<float>& originalPositions,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes,
	std::vector<float>& positions_result,
	float lambda, float mu);

void projectConstraintsHOST(std::vector<int>& indices,
	std::vector<float>& originalPositions,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes,
	std::vector<float>& positions_result,
	float lambda, float mu);


void setUpSystem(std::vector<int>& indices, std::vector<float>& originalPositions,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes,
	float gravity, float deltaT)
{
	originalPositions.push_back(0.0f); originalPositions.push_back(0.0f); originalPositions.push_back(0.0f);
	originalPositions.push_back(-0.946f); originalPositions.push_back(0.0f); originalPositions.push_back(-1.114f);
	originalPositions.push_back(0.689f); originalPositions.push_back(0.515f); originalPositions.push_back(-1.114f);
	originalPositions.push_back(0.689f); originalPositions.push_back(-0.757f); originalPositions.push_back(-1.114f);
	originalPositions.push_back(0.0f); originalPositions.push_back(0.0f); originalPositions.push_back(-2.576f);

	indices.push_back(3); indices.push_back(0); indices.push_back(2); indices.push_back(1);
	indices.push_back(3); indices.push_back(4); indices.push_back(1); indices.push_back(2);

	for (int i = 0; i < 5; ++i)
	{
		inverseMasses.push_back(1.0f);
	}
	inverseMasses[0] = 0.0f;

	for (int i = 0; i < originalPositions.size(); ++i)
	{
		positions.push_back(originalPositions[i]);
	}

	//apply one time step of deformations
	for (int i = 0; i < 5; ++i)
	{
		positions[i * 3 + 1] += inverseMasses[i] * gravity * deltaT;
	}

	//FROM MATLAB
	volumes.push_back(0.38613f);
	volumes.push_back(0.50676f);


}


int main()
{
	std::vector<int> indices;
	std::vector<float> originalPositions;
	std::vector<float> positions;
	std::vector<float> inverseMasses;
	std::vector<float> refShapeMatrixInverses;
	std::vector<float> volumes;


	float deltaT = 0.005f;
	float gravity = -9.8f;

	float lambda = 0.769231f;
	float mu = 1.15385f;

	setUpSystem(indices, originalPositions, positions, inverseMasses, refShapeMatrixInverses, volumes, gravity, deltaT);

	std::vector<float> positionsResultDevice(positions.size());
	std::vector<float> positionsResultHost(positions.size());

	//CPU
	projectConstraintsHOST(indices, originalPositions, positions,
		inverseMasses, refShapeMatrixInverses, volumes, positionsResultHost, lambda, mu);
	
	//GPU
	cudaError_t cudaStatus = projectConstraints(indices, originalPositions, positions,
		inverseMasses, refShapeMatrixInverses, volumes, positionsResultDevice, lambda, mu);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Critical Error, aborting...");
		return 1;
	}

	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}

	//Print Some Results

	std::cout << "AFTER PROJECION HOST: " << std::endl;
	for (int row = 0; row < 5; +row)
	{
		for (int col = 0; col < 3; ++col)
		{
			std::cout << positionsResultHost[row * 3 + col] << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;

	std::cout << "AFTER PROJECION DEVICE: " << std::endl;
	for (int row = 0; row < 5; +row)
	{
		for (int col = 0; col < 3; ++col)
		{
			std::cout << positionsResultDevice[row * 3 + col] << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;

	return 0;
}

cudaError_t cudaErrorWrapper(cudaError_t status)
{
	if (status != cudaSuccess) {
		fprintf(stderr, "Critical Error occured!");
		std::cout << "ERROR Details: " << cudaGetErrorString(status) << std::endl;
	}

	return status;
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

cudaError_t projectConstraints(std::vector<int>& indices, std::vector<float>& originalPositions,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes,
	std::vector<float>& positions_result,
	float lambda, float mu)
{
	float* dev_positions;
	float* dev_inverseMasses;
	int* dev_indices;

	float* dev_refShapeMatrixInverses;
	float* dev_volumes;

	cudaError_t deviceStatus;

	//Allocate memory
	int deviceCount = 0;
	deviceStatus == cudaGetDeviceCount(&deviceCount);

	std::cout << "Num CUDA Devices Found: " << deviceCount << std::endl;
	deviceStatus = cudaErrorWrapper(cudaSetDevice(0));
	getCudaDeviceProperties(0);

	std::cout << "Calling CUDA Malloc..." << std::endl;

	deviceStatus = cudaErrorWrapper(cudaMalloc((void**)&dev_indices, indices.size() * sizeof(int)));
	deviceStatus = cudaErrorWrapper(cudaMalloc((void**)&dev_positions, positions.size() * sizeof(float)));
	deviceStatus = cudaErrorWrapper(cudaMalloc((void**)&dev_inverseMasses, inverseMasses.size() * sizeof(float)));
	deviceStatus = cudaErrorWrapper(cudaMalloc((void**)&dev_refShapeMatrixInverses, refShapeMatrixInverses.size() * sizeof(float)));
	deviceStatus = cudaErrorWrapper(cudaMalloc((void**)&dev_volumes, volumes.size() * sizeof(float)));

	//Cpy memory
	deviceStatus = cudaErrorWrapper(cudaMemcpy(dev_indices, &indices[0], indices.size() * sizeof(int), cudaMemcpyHostToDevice));
	deviceStatus = cudaErrorWrapper(cudaMemcpy(dev_positions, &positions[0], positions.size() * sizeof(float), cudaMemcpyHostToDevice));
	deviceStatus = cudaErrorWrapper(cudaMemcpy(dev_inverseMasses, &inverseMasses[0], inverseMasses.size() * sizeof(float), cudaMemcpyHostToDevice));
	deviceStatus = cudaErrorWrapper(cudaMemcpy(dev_refShapeMatrixInverses, &refShapeMatrixInverses[0], refShapeMatrixInverses.size() * sizeof(float), cudaMemcpyHostToDevice));
	deviceStatus = cudaErrorWrapper(cudaMemcpy(dev_volumes, &volumes[0], volumes.size() * sizeof(float), cudaMemcpyHostToDevice));

	//Execute Kernel
	//solveFEMConstraint<<<1, 2>>>(dev_positions, dev_indices, dev_inverseMasses, dev_volumes, dev_refShapeMatrixInverses, lambda, mu);

	cudaDeviceSynchronize();

	//Cpy memory back
	positions_result.resize(positions.size());
	deviceStatus = cudaErrorWrapper(cudaMemcpy(&positions_result[0], dev_positions, positions_result.size() * sizeof(float), cudaMemcpyDeviceToHost));

	//Free memory
	cudaFree(dev_positions);
	cudaFree(dev_inverseMasses);
	cudaFree(dev_indices);
	cudaFree(dev_refShapeMatrixInverses);
	cudaFree(dev_volumes);
	return deviceStatus;
}


void projectConstraintsHOST(std::vector<int>& indices,
	std::vector<float>& originalPositions,
	std::vector<float>& positions,
	std::vector<float>& inverseMasses,
	std::vector<float>& refShapeMatrixInverses,
	std::vector<float>& volumes,
	std::vector<float>& positions_result,
	float lambda, float mu)
{
	positions_result.push_back(0.000000000000000f);
	positions_result.push_back(-4.900000000000000f);
	positions_result.push_back(0.000000000000000f);
	positions_result.push_back(-0.946000000000000f);
	positions_result.push_back(-4.900000000000000f);
	positions_result.push_back(-1.114000000000000f);
	positions_result.push_back(0.689000000000000f);
	positions_result.push_back(-4.385000000000001f);
	positions_result.push_back(-1.114000000000000f);
	positions_result.push_back(0.689000000000000f);
	positions_result.push_back(-5.657000000000000f);
	positions_result.push_back(-1.114000000000000f);
	positions_result.push_back(0.000000000000000f);
	positions_result.push_back(-4.900000000000000f);
	positions_result.push_back(-2.576000000000000f);
}