
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "CUDA_WRAPPER.h"
#include "CUDA_GLOBALS.h"

#define IDX (blockIdx.x * blockDim.x + threadIdx.x)

__shared__ float F[NUM_THREADS_PER_BLOCK][3][3];
//__shared__ float FirstPiolaKirchoffTensor[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float Gradient[NUM_THREADS_PER_BLOCK][3][4];
//__shared__ int LocalIndices[NUM_THREADS_PER_BLOCK][4];


__device__ float sqr(float x)
{
	return x * x;
}

//__device__ float traceFTransposeF(int idx)
//{
//	return TEMP1[idx][0][0] + TEMP1[idx][1][1] + TEMP1[idx][2][2];
//}

//__device__ float determinantFTransposeF(int idx)
//{
//	return TEMP1[idx][0][0]
//		* (TEMP1[idx][1][1] * TEMP1[idx][2][2] - TEMP1[idx][1][2] * TEMP1[idx][2][1])
//		- TEMP1[idx][0][1]
//		* (TEMP1[idx][1][0] * TEMP1[idx][2][2] - TEMP1[idx][1][2] * TEMP1[idx][2][0])
//		+ TEMP1[idx][0][2]
//		* (TEMP1[idx][1][0] * TEMP1[idx][2][1] - TEMP1[idx][1][1] * TEMP1[idx][2][0]);
//}

__device__ float determinantF(int idx)
{
	return F[idx][0][0]
		* (F[idx][1][1] * F[idx][2][2] - F[idx][1][2] * F[idx][2][1])
		- F[idx][0][1]
		* (F[idx][1][0] * F[idx][2][2] - F[idx][1][2] * F[idx][2][0])
		+ F[idx][0][2]
		* (F[idx][1][0] * F[idx][2][1] - F[idx][1][1] * F[idx][2][0]);
}

__device__ void calculateF(int globalIdx, int idx, float* positions, float* refShapeMatrixInverse, int* indices, int trueNumConstraints)
{
	int localIndices[4];
	localIndices[0] = indices[globalIdx + trueNumConstraints * 0] * 3;
	localIndices[1] = indices[globalIdx + trueNumConstraints * 1] * 3;
	localIndices[2] = indices[globalIdx + trueNumConstraints * 2] * 3;
	localIndices[3] = indices[globalIdx + trueNumConstraints * 3] * 3;

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
				sum += temp[row][i] * refShapeMatrixInverse[(i * 3 + col) * trueNumConstraints + globalIdx];
			}

			F[idx][row][col] = sum;
		}
	}
}

//__device__ void calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN(int idx, float mu, float lambda, float I3)
//{
//	//1. Copy over F multiplied with mu
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			FirstPiolaKirchoffTensor[idx][row][col] = F[idx][row][col] * mu;
//		}
//	}
//
//	//3. Subtract mu times FInverseTranspose
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			FirstPiolaKirchoffTensor[idx][row][col] -= TEMP1[idx][row][col] * mu;
//		}
//	}
//
//	//4. Add (lambda * logI3) / 2.0 * FInverseTranspose
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			FirstPiolaKirchoffTensor[idx][row][col] += TEMP1[idx][row][col] * ((lambda * log(I3)) / 2.0f);
//		}
//	}
//}

__device__ float calculateStrainEnergy_NEO_HOOKEAN(float volume, float lambda, float mu, float I1, float I3)
{
	return volume * (0.5f * mu * (I1 - log(I3) - 3.0f) + (lambda / 8.0f) * (log(I3) * log(I3)));
}

//__device__  void calculateStrainEnergyGradient_NEO_HOOKEAN(int globalIdx, int idx, float volume, float* refShapeMatrixInverse, int trueNumConstraints)
//{
//	//1. Copy refShapeMatrixInverse from global memory
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			//We need the TRANSPOSE of the reference shape matrix inverse
//			Gradient[idx][row][col] = refShapeMatrixInverse[(col * 3 + row) * trueNumConstraints + globalIdx];
//		}
//	}
//
//	//3. Multiply with First Piola-Kirchoff Stress tensor
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			float sum = 0.0f;
//
//			for (int i = 0; i < 3; ++i)
//			{
//				sum += FirstPiolaKirchoffTensor[idx][row][i] * Gradient[idx][i][col];
//			}
//
//			F[idx][row][col] = sum;
//		}
//	}
//
//	//4. Copy back
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			Gradient[idx][row][col] = F[idx][row][col] * volume;
//		}
//	}
//
//
//	//4. Calculate last column
//	for (int row = 0; row < 3; ++row)
//	{
//		float sum = 0.0f;
//		for (int col = 0; col < 3; ++col)
//		{
//			sum += Gradient[idx][row][col];
//		}
//		Gradient[idx][row][3] = -sum;
//	}
//}

//__device__ void calculateFTransposeF(int idx)
//{
//	//1. Multiply with F (TEMP1 transposed)
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			float sum = 0.0f;
//
//			for (int i = 0; i < 3; ++i)
//			{
//				sum += F[idx][i][row] * F[idx][i][col];
//			}
//
//			FirstPiolaKirchoffTensor[idx][row][col] = sum;
//		}
//	}
//
//	//Copy back
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			TEMP1[idx][row][col] = FirstPiolaKirchoffTensor[idx][row][col];
//		}
//	}
//}
//
//__device__ void calculateFInverseTranspose(int idx)
//{
//	//1. Calculate cofactors
//	TEMP1[idx][0][0] = F[idx][1][1] * F[idx][2][2] - F[idx][2][1] * F[idx][1][2];
//	TEMP1[idx][0][1] = -(F[idx][1][0] * F[idx][2][2] - F[idx][2][0] * F[idx][1][2]);
//	TEMP1[idx][0][2] = F[idx][1][0] * F[idx][2][1] - F[idx][2][0] * F[idx][1][1];
//
//	TEMP1[idx][1][0] = -(F[idx][0][1] * F[idx][2][2] - F[idx][2][1] * F[idx][0][2]);
//	TEMP1[idx][1][1] = F[idx][0][0] * F[idx][2][2] - F[idx][2][0] * F[idx][0][2];
//	TEMP1[idx][1][2] = -(F[idx][0][0] * F[idx][2][1] - F[idx][2][0] * F[idx][0][1]);
//
//	TEMP1[idx][2][0] = F[idx][0][1] * F[idx][1][2] - F[idx][1][1] * F[idx][0][2];
//	TEMP1[idx][2][1] = -(F[idx][0][0] * F[idx][1][2] - F[idx][1][0] * F[idx][0][2]);
//	TEMP1[idx][2][2] = F[idx][0][0] * F[idx][1][1] - F[idx][1][0] * F[idx][0][1];
//
//	//3. Calculate the determinant
//	//float determinant = determinantF(idx);
//	//printf("Determinant of F: %4.8f \n", determinant);
//
//	//4. Multiply
//	//for (int row = 0; row < 3; ++row)
//	//{
//	//	for (int col = 0; col < 3; ++col)
//	//	{
//	//		TEMP1[idx][row][col] /= determinant;
//	//	}
//	//}
//}

__device__ __forceinline__ float squaredNormGradient(int idx, int particleIdx)
{
	return sqrtf(sqr(Gradient[idx][0][particleIdx])
		+ sqr(Gradient[idx][1][particleIdx])
		+ sqr(Gradient[idx][2][particleIdx]));
}

__device__ float calculateLagrangeMultiplierDenominator(int globalIdx, int idx, float* masses, int trueNumConstraints)
{
	float denominator = 0.0f;
	for (int i = 0; i < 4; ++i)
	{
		denominator += masses[globalIdx + i * trueNumConstraints] * squaredNormGradient(idx, i);
		//if (idx == 47)
		//{
		//	printf("[%d]: mass: %4.8f, gradientNorm: %4.8f \n", idx, inverseMass[LocalIndices[idx][i]], squaredNormGradient(idx, i));
		//}
	}
	//printf("Denominator: %4.8f \n", denominator);
	return denominator;
}

__device__ void updatePositions(int globalIdx, int idx, float lagrangeMultiplier, float* positions, float* masses, int* indices, int trueNumConstraints)
{
	int localIndices[4];
	localIndices[0] = indices[globalIdx + trueNumConstraints * 0] * 3;
	localIndices[1] = indices[globalIdx + trueNumConstraints * 1] * 3;
	localIndices[2] = indices[globalIdx + trueNumConstraints * 2] * 3;
	localIndices[3] = indices[globalIdx + trueNumConstraints * 3] * 3;

	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			atomicAdd(&positions[localIndices[i] + j], masses[globalIdx + i * trueNumConstraints] * lagrangeMultiplier * Gradient[idx][j][i]);
			//atomicAdd(&positions[LocalIndices[threadIdx.x][i] * 3 + j], 0.0001f);
			//printf("%d, ", LocalIndices[threadIdx.x][i] * 3 + j);
			//printf("Position Update %4.8f \n", LocalMasses[idx][i] * lagrangeMultiplier * Gradient[idx][j][i]);
		}
		//printf("\n");
	}
}

//__device__ __forceinline__ void getIndices(int idx, int* indices, int trueNumConstraints)
//{
//	//for (int i = 0; i < 4; ++i)
//	//{
//	//	LocalIndices[threadIdx.x][i] = indices[idx * 4 + i];
//	//}
//	for (int i = 0; i < 4; ++i)
//	{
//		LocalIndices[threadIdx.x][i] = indices[idx + i * trueNumConstraints];
//	}
//}

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
	//float result = 0.0f;
	//for (int i = 0; i < 3; ++i)
	//{
	//	result += F[threadIdx.x][i][row] * F[threadIdx.x][i][col];
	//}

	//return result;
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

//__device__ void calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN_INPLACE(int idx, float mu, float lambda, float I3)
//{
//	float det = determinantF(idx);
//
//	//3. Subtract mu times FInverseTranspose
//	for (int row = 0; row < 3; ++row)
//	{
//		for (int col = 0; col < 3; ++col)
//		{
//			float FInverseTEntry;
//
//			if (row == 0)
//			{
//				if (col == 0)
//				{
//					FInverseTEntry = F[idx][1][1] * F[idx][2][2] - F[idx][2][1] * F[idx][1][2];
//				}
//				else if (col == 1)
//				{
//					FInverseTEntry = -(F[idx][1][0] * F[idx][2][2] - F[idx][2][0] * F[idx][1][2]);
//				}
//				else if (col == 2)
//				{
//					FInverseTEntry = F[idx][1][0] * F[idx][2][1] - F[idx][2][0] * F[idx][1][1];
//				}
//			}
//			else if (row == 1)
//			{
//				if (col == 0)
//				{
//					FInverseTEntry = -(F[idx][0][1] * F[idx][2][2] - F[idx][2][1] * F[idx][0][2]);
//				}
//				else if (col == 1)
//				{
//					FInverseTEntry = F[idx][0][0] * F[idx][2][2] - F[idx][2][0] * F[idx][0][2];
//				}
//				else if (col == 2)
//				{
//					FInverseTEntry = -(F[idx][0][0] * F[idx][2][1] - F[idx][2][0] * F[idx][0][1]);
//				}
//			}
//			else if (row == 2)
//			{
//				if (col == 0)
//				{
//					FInverseTEntry = F[idx][0][1] * F[idx][1][2] - F[idx][1][1] * F[idx][0][2];
//				}
//				else if (col == 1)
//				{
//					FInverseTEntry = -(F[idx][0][0] * F[idx][1][2] - F[idx][1][0] * F[idx][0][2]);
//				}
//				else if (col == 2)
//				{
//					FInverseTEntry = F[idx][0][0] * F[idx][1][1] - F[idx][1][0] * F[idx][0][1];
//				}
//			}
//
//			FInverseTEntry /= det;
//
//			FirstPiolaKirchoffTensor[idx][row][col] = F[idx][row][col] * mu - (FInverseTEntry * mu) + FInverseTEntry * ((lambda * log(I3)) / 2.0f);
//		}
//	}
//}

__device__ void calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE(int globalIdx, int idx, float volume, float* refShapeMatrixInverse, int trueNumConstraints, float mu, float lambda, float I3)
{
	float det = determinantF(idx);

	//1. Copy refShapeMatrixInverse from global memory
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			//We need the TRANSPOSE of the reference shape matrix inverse
			Gradient[idx][row][col] = refShapeMatrixInverse[(col * 3 + row) * trueNumConstraints + globalIdx];
		}
	}

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
				
				sum += PFEntry * Gradient[idx][i][col];
			}

			temp[row][col] = sum;
		}
	}

	//4. Copy back
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			Gradient[idx][row][col] = temp[row][col] * volume;
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


//__device__ void getMasses(int idx, float* masses, int trueNumConstraints)
//{
//	//for (int i = 0; i < 4; ++i)
//	//{
//	//	//printf("%d; ", LocalIndices[idx][i]);
//	//	LocalMasses[threadIdx.x][i] = masses[LocalIndices[threadIdx.x][i]];
//	//}
//	for (int i = 0; i < 4; ++i)
//	{
//		//printf("%d; ", LocalIndices[idx][i]);
//		LocalMasses[threadIdx.x][i] = masses[idx + i * trueNumConstraints];
//	}
//}

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
	if (IDX > trueNumConstraints)
	{
		return;
	}

	int idx = IDX;

	//getIndices(idx, indices, trueNumConstraints);
	//getMasses(idx, inverseMass, trueNumConstraints);

	//1. Calculate Deformation Gradient F
	calculateF(idx, threadIdx.x, positions, refShapeMatrixInverse, indices, trueNumConstraints);

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
	//calculateFTransposeF(threadIdx.x);

	//3. Compute Invariants
	//float I1 = traceFTransposeF(threadIdx.x);
	//float I3 = determinantFTransposeF(threadIdx.x);
	float I3 = calculatedeterminantFTransposeF_INPLACE();
	float I1 = calculateTraceFTransposeF_INPLACE();

	//printf("I1 = %4.8f \n", I1);
	//printf("I3 = %4.8f \n", I3);

	//calculateFInverseTranspose(threadIdx.x);

	//4. Calculate First Piola-Kirchoff Stress Tensor
	//calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN(threadIdx.x, mu, lambda, I3);
	//calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN_INPLACE(threadIdx.x, mu, lambda, I3);

	//5. Calculate StrainEnergy
	float strainEnergy = calculateStrainEnergy_NEO_HOOKEAN(volume[idx], lambda, mu, I1, I3);

	//printf("StrainEnergy = %4.8f \n", strainEnergy);

	//6. Calculate Strain Energy Gradient
	//calculateStrainEnergyGradient_NEO_HOOKEAN(idx, threadIdx.x, volume[idx], refShapeMatrixInverse, trueNumConstraints);
	calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE(idx, threadIdx.x, volume[idx], refShapeMatrixInverse, trueNumConstraints, mu, lambda, I3);

	//7. Calculate Lagrange Multiplier
	float denominator = calculateLagrangeMultiplierDenominator(idx, threadIdx.x, inverseMass, trueNumConstraints);

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
	updatePositions(idx, threadIdx.x, lagrangeMultiplier, positions, inverseMass, indices, trueNumConstraints);
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
			settings.lambda, settings.mu, settings.trueNumberOfConstraints);

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