
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "CUDA_WRAPPER.h"

const int NUM_THREADS_PER_BLOCK = 256;

__shared__ float F[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float TEMP1[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float FirstPiolaKirchoffTensor[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float Gradient[NUM_THREADS_PER_BLOCK][3][4];
__shared__ int LocalIndices[NUM_THREADS_PER_BLOCK][4];
__shared__ float LocalMasses[NUM_THREADS_PER_BLOCK][4];


__device__ float sqr(float x)
{
	return x * x;
}

__device__ float traceFTransposeF(int idx)
{
	return TEMP1[idx][0][0] + TEMP1[idx][1][1] + TEMP1[idx][2][2];
}

__device__ float determinantFTransposeF(int idx)
{
	return TEMP1[idx][0][0]
		* (TEMP1[idx][1][1] * TEMP1[idx][2][2] - TEMP1[idx][1][2] * TEMP1[idx][2][1])
		- TEMP1[idx][0][1]
		* (TEMP1[idx][1][0] * TEMP1[idx][2][2] - TEMP1[idx][1][2] * TEMP1[idx][2][0])
		+ TEMP1[idx][0][2]
		* (TEMP1[idx][1][0] * TEMP1[idx][2][1] - TEMP1[idx][1][1] * TEMP1[idx][2][0]);
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

__device__ void calculateF(int globalIdx, int idx, float* positions, float* refShapeMatrixInverse)
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
				sum += FirstPiolaKirchoffTensor[idx][row][i] * refShapeMatrixInverse[globalIdx * 9 + i * 3 + col];
			}

			F[idx][row][col] = sum;
		}
	}
}

__device__ void calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN(int idx, float mu, float lambda, float I3)
{
	//1. Copy over F multiplied with mu
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FirstPiolaKirchoffTensor[idx][row][col] = F[idx][row][col] * mu;
		}
	}

	//3. Subtract mu times FInverseTranspose
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FirstPiolaKirchoffTensor[idx][row][col] -= TEMP1[idx][row][col] * mu;
		}
	}

	//4. Add (lambda * logI3) / 2.0 * FInverseTranspose
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			FirstPiolaKirchoffTensor[idx][row][col] += TEMP1[idx][row][col] * ((lambda * log(I3)) / 2.0f);
		}
	}
}

__device__ float calculateStrainEnergy_NEO_HOOKEAN(float volume, float lambda, float mu, float I1, float I3)
{
	return volume * (0.5f * mu * (I1 - log(I3) - 3.0f) + (lambda / 8.0f) * (log(I3) * log(I3)));
}

__device__ void calculateStrainEnergyGradient_NEO_HOOKEAN(int globalIdx, int idx, float volume, float* refShapeMatrixInverse)
{
	//1. Copy refShapeMatrixInverse from global memory
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			//We need the TRANSPOSE of the reference shape matrix inverse
			Gradient[idx][row][col] = refShapeMatrixInverse[globalIdx * 9 + col * 3 + row];
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
				sum += FirstPiolaKirchoffTensor[idx][row][i] * Gradient[idx][i][col];
			}

			TEMP1[idx][row][col] = sum;
		}
	}

	//4. Copy back
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			Gradient[idx][row][col] = TEMP1[idx][row][col] * volume;
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
	//Combine all into one loop in future!

	//1. Copy over F
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			TEMP1[idx][row][col] = F[idx][row][col];
		}
	}

	//3. Multiply with F (TEMP1 transposed)
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				sum += TEMP1[idx][i][row] * F[idx][i][col];
			}

			FirstPiolaKirchoffTensor[idx][row][col] = sum;
		}
	}

	//Copy back
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			TEMP1[idx][row][col] = FirstPiolaKirchoffTensor[idx][row][col];
		}
	}
}

__device__ void calculateFInverseTranspose(int idx)
{
	//1. Calculate cofactors
	TEMP1[idx][0][0] = F[idx][1][1] * F[idx][2][2] - F[idx][2][1] * F[idx][1][2];
	TEMP1[idx][0][1] = -(F[idx][1][0] * F[idx][2][2] - F[idx][2][0] * F[idx][1][2]);
	TEMP1[idx][0][2] = F[idx][1][0] * F[idx][2][1] - F[idx][2][0] * F[idx][1][1];

	TEMP1[idx][1][0] = -(F[idx][0][1] * F[idx][2][2] - F[idx][2][1] * F[idx][0][2]);
	TEMP1[idx][1][1] = F[idx][0][0] * F[idx][2][2] - F[idx][2][0] * F[idx][0][2];
	TEMP1[idx][1][2] = -(F[idx][0][0] * F[idx][2][1] - F[idx][2][0] * F[idx][0][1]);

	TEMP1[idx][2][0] = F[idx][0][1] * F[idx][1][2] - F[idx][1][1] * F[idx][0][2];
	TEMP1[idx][2][1] = -(F[idx][0][0] * F[idx][1][2] - F[idx][1][0] * F[idx][0][2]);
	TEMP1[idx][2][2] = F[idx][0][0] * F[idx][1][1] - F[idx][1][0] * F[idx][0][1];

	//3. Calculate the determinant
	float determinant = determinantF(idx);
	//printf("Determinant of F: %4.8f \n", determinant);

	//4. Multiply
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			TEMP1[idx][row][col] /= determinant;
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
		denominator += LocalMasses[idx][i] * squaredNormGradient(idx, i);
		//if (idx == 47)
		//{
		//	printf("[%d]: mass: %4.8f, gradientNorm: %4.8f \n", idx, inverseMass[LocalIndices[idx][i]], squaredNormGradient(idx, i));
		//}
	}
	//printf("Denominator: %4.8f \n", denominator);
	return denominator;
}

__device__ void updatePositions(int idx, float lagrangeMultiplier, float* positions, float* inverseMass)
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			atomicAdd(&positions[LocalIndices[idx][i] * 3 + j], LocalMasses[idx][i] * lagrangeMultiplier * Gradient[idx][j][i]);
			//atomicAdd(&positions[LocalIndices[threadIdx.x][i] * 3 + j], 0.0001f);
			//printf("%d, ", LocalIndices[threadIdx.x][i] * 3 + j);
			//printf("Position Update %4.8f \n", LocalMasses[idx][i] * lagrangeMultiplier * Gradient[idx][j][i]);
		}
		//printf("\n");
	}
}

__device__ void getIndices(int idx, int* indices)
{
	for (int i = 0; i < 4; ++i)
	{
		LocalIndices[threadIdx.x][i] = indices[idx * 4 + i];
	}
}

__device__ void getMasses(int idx, float* masses)
{
	for (int i = 0; i < 4; ++i)
	{
		//printf("%d; ", LocalIndices[idx][i]);
		LocalMasses[threadIdx.x][i] = masses[LocalIndices[threadIdx.x][i]];
	}
}

__device__ bool isFIdentity()
{
	return(abs(F[threadIdx.x][0][0] - 1.0f) < 1e-8f
		&& abs(F[threadIdx.x][0][1] - 0.0f) < 1e-8f
		&& abs(F[threadIdx.x][0][2] - 0.0f) < 1e-8f
		&& abs(F[threadIdx.x][1][0] - 0.0f) < 1e-8f
		&& abs(F[threadIdx.x][1][1] - 1.0f) < 1e-8f
		&& abs(F[threadIdx.x][1][2] - 0.0f) < 1e-8f
		&& abs(F[threadIdx.x][2][0] - 0.0f) < 1e-8f
		&& abs(F[threadIdx.x][2][1] - 0.0f) < 1e-8f
		&& abs(F[threadIdx.x][2][2] - 1.0f) < 1e-8f);
}

__global__ void solveFEMConstraint(float* positions, int* indices, float* inverseMass, float* volume, float* refShapeMatrixInverse,
	float lambda, float mu, int trueNumConstraints)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) % trueNumConstraints;

	getIndices(idx, indices);

	getMasses(idx, inverseMass);

	//1. Calculate Deformation Gradient F
	calculateF(idx, threadIdx.x, positions, refShapeMatrixInverse);

	//if (idx == 47)
	//{
	//	printf("F [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", idx,
	//		F[threadIdx.x][0][0], F[threadIdx.x][0][1], F[threadIdx.x][0][2],
	//		F[threadIdx.x][1][0], F[threadIdx.x][1][1], F[threadIdx.x][1][2],
	//		F[threadIdx.x][2][0], F[threadIdx.x][2][1], F[threadIdx.x][2][2]);
	//	printf("RefShapeInverse [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", idx,
	//	refShapeMatrixInverse[idx * 9 + 0 * 3 + 0], refShapeMatrixInverse[idx * 9 + 0 * 3 + 1], refShapeMatrixInverse[idx * 9 + 0 * 3 + 2],
	//	refShapeMatrixInverse[idx * 9 + 1 * 3 + 0], refShapeMatrixInverse[idx * 9 + 1 * 3 + 1], refShapeMatrixInverse[idx * 9 + 1 * 3 + 2],
	//	refShapeMatrixInverse[idx * 9 + 2 * 3 + 0], refShapeMatrixInverse[idx * 9 + 2 * 3 + 1], refShapeMatrixInverse[idx * 9 + 2 * 3 + 2]);
	//	printf("Deformed Shape [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", threadIdx.x,
	//		FirstPiolaKirchoffTensor[threadIdx.x][0][0], FirstPiolaKirchoffTensor[threadIdx.x][0][1], FirstPiolaKirchoffTensor[threadIdx.x][0][2],
	//		FirstPiolaKirchoffTensor[threadIdx.x][1][0], FirstPiolaKirchoffTensor[threadIdx.x][1][1], FirstPiolaKirchoffTensor[threadIdx.x][1][2],
	//		FirstPiolaKirchoffTensor[threadIdx.x][2][0], FirstPiolaKirchoffTensor[threadIdx.x][2][1], FirstPiolaKirchoffTensor[threadIdx.x][2][2]);
	//}

	//printf("F [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", idx,
	//	F[threadIdx.x][0][0], F[threadIdx.x][0][1], F[threadIdx.x][0][2],
	//	F[threadIdx.x][1][0], F[threadIdx.x][1][1], F[threadIdx.x][1][2],
	//	F[threadIdx.x][2][0], F[threadIdx.x][2][1], F[threadIdx.x][2][2]);

	/*printf("Deformed Shape [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", threadIdx.x,
		FirstPiolaKirchoffTensor[threadIdx.x][0][0], FirstPiolaKirchoffTensor[threadIdx.x][0][1], FirstPiolaKirchoffTensor[threadIdx.x][0][2],
		FirstPiolaKirchoffTensor[threadIdx.x][1][0], FirstPiolaKirchoffTensor[threadIdx.x][1][1], FirstPiolaKirchoffTensor[threadIdx.x][1][2],
		FirstPiolaKirchoffTensor[threadIdx.x][2][0], FirstPiolaKirchoffTensor[threadIdx.x][2][1], FirstPiolaKirchoffTensor[threadIdx.x][2][2]);*/
	/*printf("RefShapeInverse [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", idx,
	refShapeMatrixInverse[idx * 9 + 0 * 3 + 0], refShapeMatrixInverse[idx * 9 + 0 * 3 + 1], refShapeMatrixInverse[idx * 9 + 0 * 3 + 2],
	refShapeMatrixInverse[idx * 9 + 1 * 3 + 0], refShapeMatrixInverse[idx * 9 + 1 * 3 + 1], refShapeMatrixInverse[idx * 9 + 1 * 3 + 2],
	refShapeMatrixInverse[idx * 9 + 2 * 3 + 0], refShapeMatrixInverse[idx * 9 + 2 * 3 + 1], refShapeMatrixInverse[idx * 9 + 2 * 3 + 2]);*/

	//if (isFIdentity())
	//{
	//	return;
	//}

	//printf("Deformed Shape [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", threadIdx.x,
	//	FirstPiolaKirchoffTensor[threadIdx.x][0][0], FirstPiolaKirchoffTensor[threadIdx.x][0][1], FirstPiolaKirchoffTensor[threadIdx.x][0][2],
	//	FirstPiolaKirchoffTensor[threadIdx.x][1][0], FirstPiolaKirchoffTensor[threadIdx.x][1][1], FirstPiolaKirchoffTensor[threadIdx.x][1][2],
	//	FirstPiolaKirchoffTensor[threadIdx.x][2][0], FirstPiolaKirchoffTensor[threadIdx.x][2][1], FirstPiolaKirchoffTensor[threadIdx.x][2][2]);

	//printf("F [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", threadIdx.x,
	//	F[threadIdx.x][0][0], F[threadIdx.x][0][1], F[threadIdx.x][0][2],
	//	F[threadIdx.x][1][0], F[threadIdx.x][1][1], F[threadIdx.x][1][2],
	//	F[threadIdx.x][2][0], F[threadIdx.x][2][1], F[threadIdx.x][2][2]);

	/*printf("RefShapeInverse [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n,%4.8f, %4.8f, %4.8f \n", threadIdx.x,
		refShapeMatrixInverse[idx * 9 + 0 * 3 + 0], refShapeMatrixInverse[idx * 9 + 0 * 3 + 1], refShapeMatrixInverse[idx * 9 + 0 * 3 + 2],
		refShapeMatrixInverse[idx * 9 + 1 * 3 + 0], refShapeMatrixInverse[idx * 9 + 1 * 3 + 1], refShapeMatrixInverse[idx * 9 + 1 * 3 + 2],
		refShapeMatrixInverse[idx * 9 + 2 * 3 + 0], refShapeMatrixInverse[idx * 9 + 2 * 3 + 1], refShapeMatrixInverse[idx * 9 + 2 * 3 + 2]);*/


	//2. Compute Cauchy Tensors
	calculateFTransposeF(threadIdx.x);

	//printf("FTransposeF: \n");
	//for (int row = 0; row < 3; ++row)
	//{
	//	for (int col = 0; col < 3; ++col)
	//	{
	//		printf("%4.8f,", FTransposeF[idx][row][col]);
	//	}
	//	printf("\n");
	//}
	//printf("\n \n");

	//3. Compute Invariants
	float I1 = traceFTransposeF(threadIdx.x);
	float I3 = determinantFTransposeF(threadIdx.x);

	//printf("I1 = %4.8f \n", I1);
	//printf("I3 = %4.8f \n", I3);

	calculateFInverseTranspose(threadIdx.x);

	//4. Calculate First Piola-Kirchoff Stress Tensor
	calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN(threadIdx.x, mu, lambda, I3);

	//printf("PF: \n");
	//for (int row = 0; row < 3; ++row)
	//{
	//	for (int col = 0; col < 3; ++col)
	//	{
	//		printf("%4.8f,", FirstPiolaKirchoffTensor[idx][row][col]);
	//	}
	//	printf("\n");
	//}
	//printf("\n \n");

	//5. Calculate StrainEnergy
	float strainEnergy = calculateStrainEnergy_NEO_HOOKEAN(volume[idx], lambda, mu, I1, I3);

	//printf("StrainEnergy = %4.8f \n", strainEnergy);

	//6. Calculate Strain Energy Gradient
	calculateStrainEnergyGradient_NEO_HOOKEAN(idx, threadIdx.x, volume[idx], refShapeMatrixInverse);

	//printf("Strain Energy Gradient: \n");
	//for (int row = 0; row < 3; ++row)
	//{
	//	for (int col = 0; col < 4; ++col)
	//	{
	//		printf("%4.8f,", Gradient[idx][row][col]);
	//	}
	//	printf("\n");
	//}
	//printf("\n \n");

	//7. Calculate Lagrange Multiplier
	float denominator = calculateLagrangeMultiplierDenominator(threadIdx.x, inverseMass);

	//if (denominator == 0.0f)
	//{
	//	return;
	//}

	float lagrangeMultiplier = -(strainEnergy / denominator);
	//if (idx == 47)
	//{
	//	printf("[%d]: I1 = %8.16f \n", idx, I1);
	//	printf("[%d]: I3 = %8.16f \n", idx, I3);
	//	printf("[%d]: lagrangeMultiplier = %8.16f \n", idx, lagrangeMultiplier);
	//	printf("[%d]: strainEnergy = %8.16f \n", idx, strainEnergy);
	//	printf("[%d]: denominator = %8.16f \n", idx, denominator);
	//	printf("PF [%d]: \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f \n", idx,
	//		FirstPiolaKirchoffTensor[threadIdx.x][0][0], FirstPiolaKirchoffTensor[threadIdx.x][0][1], FirstPiolaKirchoffTensor[threadIdx.x][0][2],
	//		FirstPiolaKirchoffTensor[threadIdx.x][1][0], FirstPiolaKirchoffTensor[threadIdx.x][1][1], FirstPiolaKirchoffTensor[threadIdx.x][1][2],
	//		FirstPiolaKirchoffTensor[threadIdx.x][2][0], FirstPiolaKirchoffTensor[threadIdx.x][2][1], FirstPiolaKirchoffTensor[threadIdx.x][2][2]);
	//	printf("Gradient [%d]: \n %4.8f, %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f, %4.8f \n %4.8f, %4.8f, %4.8f, %4.8f \n", idx,
	//		Gradient[threadIdx.x][0][0], Gradient[threadIdx.x][0][1], Gradient[threadIdx.x][0][2], Gradient[threadIdx.x][0][3],
	//		Gradient[threadIdx.x][1][0], Gradient[threadIdx.x][1][1], Gradient[threadIdx.x][1][2], Gradient[threadIdx.x][1][3],
	//		Gradient[threadIdx.x][2][0], Gradient[threadIdx.x][2][1], Gradient[threadIdx.x][2][2], Gradient[threadIdx.x][2][3]);
	//}

	if (isnan(lagrangeMultiplier))
	{
		return;
	}

	//8. Update Positions
	updatePositions(threadIdx.x, lagrangeMultiplier, positions, inverseMass);
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
	std::cout << "Executing Kernel..." << settings.numConstraintIterations<< std::endl;
	std::cout << settings.numBlocks * settings.numThreadsPerBlock << "threads..." << std::endl;
	std::cout << settings.trueNumberOfConstraints << "true num of constraints..." << std::endl;

	dim3 numBlocks;
	dim3 numThreads;
	numBlocks.x = settings.numBlocks;
	//numBlocks.x = 1;
	numBlocks.y = 1;
	numBlocks.z = 1;

	numThreads.x = settings.numThreadsPerBlock;
	numThreads.y = 1;
	numThreads.z = 1;

	for (int it = 0; it < settings.numConstraintIterations; ++it)
	{
		solveFEMConstraint << <numBlocks, numThreads >> >(
			device_positions, device_indices, device_inverseMasses,
			device_volumes, device_refShapeMatrixInverses,
			settings.lambda, settings.mu, settings.trueNumberOfConstraints);

		cudaErrorWrapper(cudaDeviceSynchronize());
	}

	cudaErrorWrapper(cudaDeviceSynchronize());

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

	return true;
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