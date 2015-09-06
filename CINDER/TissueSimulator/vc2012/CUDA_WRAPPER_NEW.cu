//NOTE: As outlined in thesis, code for diagonalisation is largely based on Kopp's work, see http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "CUDA_WRAPPER.h"
#include "CUDA_GLOBALS.h"

#define IDX (blockIdx.x * blockDim.x + threadIdx.x)
#define LOCAL_IDX threadIdx.x

#define SHIFT_1 1
#define SHIFT_2 0
//Indexing for F in global memory
#define GLOBAL_F_IDX_0 ((0 * (trueNumConstraints + SHIFT_1)) + IDX)
#define GLOBAL_F_IDX_1 ((1 * (trueNumConstraints + SHIFT_1)) + IDX)
#define GLOBAL_F_IDX_2 ((2 * (trueNumConstraints + SHIFT_1)) + IDX)

#define GLOBAL_U_IDX ((row * 3 + col) * (trueNumConstraints + SHIFT_2) + IDX)

#define GLOBAL_Vt_MatrixMult_IDX ((col * 3 + i) * (trueNumConstraints + SHIFT_2) + IDX)
#define GLOBAL_Vt_IDX ((row * 3 + col) * (trueNumConstraints + SHIFT_2) + IDX)


#define JACOBI_LIMIT_1 1.0e-22f


__shared__ float F[NUM_THREADS_PER_BLOCK][3];

#define sqr(x) (x * x)

__device__ __forceinline__ float calculateStrainEnergy_NEO_HOOKEAN(float volume, float lambda, float mu, float I1, float I3)
{
	return volume * (0.5f * mu * (I1 - log(I3) - 3.0f) + (lambda / 8.0f) * (log(I3) * log(I3)));
}

__device__ __forceinline__ float calculateStrainEnergy_NEO_HOOKEAN_ANISOTROPIC(float volume, float lambda, float mu, float I1, float I3,
	float anisotropyStrength, float stretch)
{
	float result = (volume * (0.5f * mu * (I1 - log(I3) - 3.0f) + (lambda / 8.0f) * (log(I3) * log(I3))));

	if (stretch - 1.0f > JACOBI_LIMIT_1)
	{
		result += ((anisotropyStrength / 2.0f) * powf(stretch - 1.0f, 2.0f));
	}

	return result;
}

__device__ void updatePositions_recomputeGradients(float lagrangeMultiplier, float* positions,
	float* masses, int* indices, int trueNumConstraints, int numParticles, float volume, float* refShapeMatrixInverse,
	float I3, float lambda, float mu, float det,
	float* globalU, float* globalV)
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

	float temp0[3][3];
	float temp[3][3];
	// Load U
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			temp0[row][col] = globalU[GLOBAL_U_IDX];
		}
	}

	//printf("%.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n",
	//	temp0[0][0], temp0[0][1], temp0[0][2],
	//	temp0[1][0], temp0[1][1], temp0[1][2],
	//	temp0[2][0], temp0[2][1], temp0[2][2]);

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			float FInverseTEntry = 0.0f;

			FInverseTEntry = 1.0f / F[LOCAL_IDX][col];

			sum += temp0[row][col] * (F[LOCAL_IDX][col] * mu - (FInverseTEntry * mu) + FInverseTEntry * ((lambda * log(I3)) / 2.0f));
			temp[row][col] = sum;
		}
	}

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			for (int i = 0; i < 3; ++i)
			{
				sum += temp[row][i] * globalV[GLOBAL_Vt_MatrixMult_IDX];
			}
			temp0[row][col] = sum;
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
				sum += temp0[row][i] * refShapeMatrixInverse[(col * 3 + i) * trueNumConstraints + IDX];
			}

			temp[row][col] = sum;
		}
	}

	int localIndices[4];
	localIndices[0] = indices[IDX + trueNumConstraints * 0] * 3;
	localIndices[1] = indices[IDX + trueNumConstraints * 1] * 3;
	localIndices[2] = indices[IDX + trueNumConstraints * 2] * 3;
	localIndices[3] = indices[IDX + trueNumConstraints * 3] * 3;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			atomicAdd(&positions[localIndices[i] + j], masses[IDX + i * trueNumConstraints] * lagrangeMultiplier * temp[j][i] * volume);
		}
	}

	for (int j = 0; j < 3; ++j)
	{
		float sum = 0.0f;
		for (int col = 0; col < 3; ++col)
		{
			sum += temp[j][col] * volume;
		}

		atomicAdd(&positions[localIndices[3] + j], masses[IDX + 3 * trueNumConstraints] * lagrangeMultiplier * -sum);
	}
}

__device__ __forceinline__ void updatePositions_recomputeGradients_ANISOTROPIC(float lagrangeMultiplier, float* positions,
	float* masses, int* indices, int trueNumConstraints, int numParticles, float volume, float* refShapeMatrixInverse,
	float I3, float lambda, float mu, float det,
	float* globalU, float* globalV,
	float anisotropyStrength, float stretch,
	float rotA0, float rotA1, float rotA2)
{
	float temp0[3][3];
	float temp[3][3];
	// Load U
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			temp0[row][col] = globalU[GLOBAL_U_IDX];
		}
	}
	
	float rotatedA[3];
	rotatedA[0] = rotA0;
	rotatedA[1] = rotA1;
	rotatedA[2] = rotA2;
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			float FInverseTEntry = 0.0f;

			FInverseTEntry = 1.0f / F[LOCAL_IDX][col];

			sum += (F[LOCAL_IDX][col] * mu - (FInverseTEntry * mu) + FInverseTEntry * ((lambda * log(I3)) / 2.0f));

			if (stretch - 1.0f > JACOBI_LIMIT_1)
			{
				sum += ((
					powf(F[LOCAL_IDX][0] * F[LOCAL_IDX][1] * F[LOCAL_IDX][2], -2.0f / 3.0f)
					* (anisotropyStrength * (stretch - 1.0f))
					* ((rotatedA[row] * rotatedA[row]) + (stretch / 3.0f) * (1.0f / std::powf(F[LOCAL_IDX][col], 2.0f)))
					)
					);
			}
			temp[row][col] = temp0[row][col] * sum;
		}
	}

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			for (int i = 0; i < 3; ++i)
			{
				sum += temp[row][i] * globalV[GLOBAL_Vt_MatrixMult_IDX];
			}
			temp0[row][col] = sum;
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
				sum += temp0[row][i] * refShapeMatrixInverse[(col * 3 + i) * trueNumConstraints + IDX];
			}

			temp[row][col] = sum;
		}
	}

	int localIndices[4];
	localIndices[0] = indices[IDX + trueNumConstraints * 0] * 3;
	localIndices[1] = indices[IDX + trueNumConstraints * 1] * 3;
	localIndices[2] = indices[IDX + trueNumConstraints * 2] * 3;
	localIndices[3] = indices[IDX + trueNumConstraints * 3] * 3;

	if (fabs(lagrangeMultiplier) > 0.5f)
	{
		//printf("%.4f, ", lagrangeMultiplier);
		lagrangeMultiplier = 0.0f;
	}
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			atomicAdd(&positions[localIndices[i] + j], masses[IDX + i * trueNumConstraints] * lagrangeMultiplier * temp[j][i] * volume);
		}
	}

	for (int j = 0; j < 3; ++j)
	{
		float sum = 0.0f;
		for (int col = 0; col < 3; ++col)
		{
			sum += temp[j][col] * volume;
		}

		atomicAdd(&positions[localIndices[3] + j], masses[IDX + 3 * trueNumConstraints] * lagrangeMultiplier * -sum);
	}
}

__device__ __forceinline__ float calculateTraceFTransposeF_INPLACE()
{
	float trace = 0.0f;
	for (int diagIdx = 0; diagIdx < 3; ++diagIdx)
	{
		trace += F[LOCAL_IDX][diagIdx] * F[LOCAL_IDX][diagIdx];
	}

	return trace;
}

__device__ __forceinline__ float calculatedeterminantFTransposeF_INPLACE()
{
	return sqr(F[LOCAL_IDX][0]) * sqr(F[LOCAL_IDX][1]) * sqr(F[LOCAL_IDX][2]);
}

__device__ __forceinline__ void calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE(float volume, float* refShapeMatrixInverse, int trueNumConstraints, float mu, float lambda, float I3,
	float& snGr0, float& snGr1, float& snGr2, float& snGr3, float det,
	float* globalU, float* globalV)
{
	float temp0[3][3];
	float temp[3][3];
	// Load U
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			temp0[row][col] = globalU[GLOBAL_U_IDX];
		}
	}

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			float FInverseTEntry = 0.0f;

			FInverseTEntry = 1.0f / F[LOCAL_IDX][col];

			sum += temp0[row][col] * (F[LOCAL_IDX][col] * mu - (FInverseTEntry * mu) + FInverseTEntry * ((lambda * log(I3)) / 2.0f));
			temp[row][col] = sum;
		}
	}

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			for (int i = 0; i < 3; ++i)
			{
				sum += temp[row][i] * globalV[GLOBAL_Vt_MatrixMult_IDX];
			}
			temp0[row][col] = sum;
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
				sum += temp0[row][i] * refShapeMatrixInverse[(col * 3 + i) * trueNumConstraints + IDX];
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


__device__ __forceinline__ void calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE_ANISOTROPIC(float volume, float* refShapeMatrixInverse, int trueNumConstraints, float mu, float lambda, float I3,
	float& snGr0, float& snGr1, float& snGr2, float& snGr3, float det,
	float* globalU, float* globalV,
	float anisotropyStrength, float stretch,
	float rotA0, float rotA1, float rotA2)
{
	float temp0[3][3];
	float temp[3][3];
	// Load U
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			temp0[row][col] = globalU[GLOBAL_U_IDX];
		}
	}


	float rotatedA[3];
	rotatedA[0] = rotA0;
	rotatedA[1] = rotA1;
	rotatedA[2] = rotA2;
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			float FInverseTEntry = 0.0f;

			FInverseTEntry = 1.0f / F[LOCAL_IDX][col];

			sum += (F[LOCAL_IDX][col] * mu - (FInverseTEntry * mu) + FInverseTEntry * ((lambda * log(I3)) / 2.0f));

			if (stretch - 1.0f > JACOBI_LIMIT_1)
			{
				sum += ((
					powf(F[LOCAL_IDX][0] * F[LOCAL_IDX][1] * F[LOCAL_IDX][2], -2.0f / 3.0f)
					* (anisotropyStrength * (stretch - 1.0f))
					* ((rotatedA[row] * rotatedA[row]) + (stretch / 3.0f) * (1.0f / std::powf(F[LOCAL_IDX][col], 2.0f)))
					)
					);
			}
			temp[row][col] = temp0[row][col] * sum;
		}
	}

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;
			for (int i = 0; i < 3; ++i)
			{
				sum += temp[row][i] * globalV[GLOBAL_Vt_MatrixMult_IDX];
			}
			temp0[row][col] = sum;
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
				sum += temp0[row][i] * refShapeMatrixInverse[(col * 3 + i) * trueNumConstraints + IDX];
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
	float lambda, float mu, int trueNumConstraints, int numParticles, float* globalF, float* globalU, float* globalV)
{
	if (IDX > trueNumConstraints)
	{
		return;
	}

	//1. Calculate Deformation Gradient F
	//calculateF(IDX, threadIdx.x, positions, refShapeMatrixInverse, indices, trueNumConstraints, numParticles);

	//1. Load Deformation Gradient

	F[LOCAL_IDX][0] = globalF[GLOBAL_F_IDX_0];
	F[LOCAL_IDX][1] = globalF[GLOBAL_F_IDX_1];
	F[LOCAL_IDX][2] = globalF[GLOBAL_F_IDX_2];

	//printf("%.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n",
	//	F[LOCAL_IDX][0][0], F[LOCAL_IDX][0][1], F[LOCAL_IDX][0][2],
	//	F[LOCAL_IDX][1][0], F[LOCAL_IDX][1][1], F[LOCAL_IDX][1][2],
	//	F[LOCAL_IDX][2][0], F[LOCAL_IDX][2][1], F[LOCAL_IDX][2][2]);

	//3. Compute Invariants

	float I3 = calculatedeterminantFTransposeF_INPLACE();
	float I1 = calculateTraceFTransposeF_INPLACE();

	//6. Calculate Strain Energy Gradient
	float snGr0;
	float snGr1;
	float snGr2;
	float snGr3;

	float det = 0.0f;

	calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE(volume[IDX], refShapeMatrixInverse, trueNumConstraints, mu, lambda, I3,
		snGr0, snGr1, snGr2, snGr3, det,
		globalU, globalV);

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
	updatePositions_recomputeGradients(lagrangeMultiplier, positions, inverseMass, indices, trueNumConstraints, numParticles,
		volume[IDX], refShapeMatrixInverse, I3, lambda, mu, det,
		globalU, globalV);
}


__global__ void solveFEMConstraint_ANISOTROPIC(float* positions, int* indices, float* inverseMass, float* volume, float* refShapeMatrixInverse,
	float lambda, float mu, int trueNumConstraints, int numParticles, float* globalF, float* globalU, float* globalV,
	float* anisotropyDirection, float anisotropyStrength)
{
	if (IDX > trueNumConstraints)
	{
		return;
	}

	//anisotropyStrength = 0.0f;

	//1. Load Deformation Gradient
	F[LOCAL_IDX][0] = globalF[GLOBAL_F_IDX_0];
	F[LOCAL_IDX][1] = globalF[GLOBAL_F_IDX_1];
	F[LOCAL_IDX][2] = globalF[GLOBAL_F_IDX_2];

	//3. Compute Invariants

	float I3 = calculatedeterminantFTransposeF_INPLACE();
	float I1 = calculateTraceFTransposeF_INPLACE();

	//Compute rotated direction vector
	float rotatedDirection[3];
	float temp = anisotropyDirection[0 * trueNumConstraints + IDX];
	for (int i = 0; i < 3; ++i)
	{
		rotatedDirection[i] = globalV[(i * 3 + 0) * (trueNumConstraints + 1) + IDX] * temp;
		//rotatedDirection[i] = globalV[(0 * 3 + i) * (trueNumConstraints + 1) + IDX] * temp;
	}

	temp = anisotropyDirection[1 * trueNumConstraints + IDX];
	for (int i = 0; i < 3; ++i)
	{
		rotatedDirection[i] += globalV[(i * 3 + 1) * (trueNumConstraints + 1) + IDX] * temp;
		//rotatedDirection[i] += globalV[(1 * 3 + i) * (trueNumConstraints + 1) + IDX] * temp;
	}

	temp = anisotropyDirection[2 * trueNumConstraints + IDX];
	for (int i = 0; i < 3; ++i)
	{
		rotatedDirection[i] += globalV[(i * 3 + 2) * (trueNumConstraints + 1) + IDX] * temp;
		//rotatedDirection[i] += globalV[(2 * 3 + i) * (trueNumConstraints + 1) + IDX] * temp;
	}

	//rotatedDirection[0] = 0.0f;
	//rotatedDirection[1] = 0.0f;
	//rotatedDirection[2] = 0.0f;

	//Compute stretch
	float stretch = sqrtf((sqr(F[LOCAL_IDX][0]) * sqr(rotatedDirection[0]))
		+ (sqr(F[LOCAL_IDX][1]) * sqr(rotatedDirection[1]))
		+ (sqr(F[LOCAL_IDX][2]) * sqr(rotatedDirection[2])));

	//printf("%.4f, ", sqr(F[LOCAL_IDX][0]) * sqr(rotatedDirection[0])
	//	+ sqr(F[LOCAL_IDX][1]) * sqr(rotatedDirection[1])
	//	+ sqr(F[LOCAL_IDX][2]) * sqr(rotatedDirection[2]));

	//stretch = 0.0f;
	//anisotropyStrength = 0.0f;

	//printf("[%.4f, %.4f, %.4f]; ", anisotropyDirection[0 * trueNumConstraints + IDX], anisotropyDirection[1 * trueNumConstraints + IDX], anisotropyDirection[2 * trueNumConstraints + IDX]);
	//printf("%.4f, ", stretch);
	//anisotropyStrength = 0.0f;
	//6. Calculate Strain Energy Gradient
	float snGr0;
	float snGr1;
	float snGr2;
	float snGr3;

	float det = 0.0f;

	calculateStrainEnergyGradient_NEO_HOOKEAN_INPLACE_ANISOTROPIC(volume[IDX], refShapeMatrixInverse, trueNumConstraints, mu, lambda, I3,
		snGr0, snGr1, snGr2, snGr3, det,
		globalU, globalV,
		anisotropyStrength, stretch,
		rotatedDirection[0], rotatedDirection[1], rotatedDirection[2]);

	//7. Calculate Lagrange Multiplier
	float denominator = 0.0f;
	denominator += inverseMass[IDX + 0 * trueNumConstraints] * snGr0;
	denominator += inverseMass[IDX + 1 * trueNumConstraints] * snGr1;
	denominator += inverseMass[IDX + 2 * trueNumConstraints] * snGr2;
	denominator += inverseMass[IDX + 3 * trueNumConstraints] * snGr3;

	float lagrangeMultiplier = -(calculateStrainEnergy_NEO_HOOKEAN_ANISOTROPIC(volume[IDX], lambda, mu, I1, I3, anisotropyStrength, stretch) / denominator);

	if (isnan(lagrangeMultiplier))
	{
		return;
	}

	//8. Update Positions
	updatePositions_recomputeGradients_ANISOTROPIC(lagrangeMultiplier, positions, inverseMass, indices, trueNumConstraints, numParticles,
		volume[IDX], refShapeMatrixInverse, I3, lambda, mu, det,
		globalU, globalV,
		anisotropyStrength, stretch,
		rotatedDirection[0], rotatedDirection[1], rotatedDirection[2]);
}


__global__ void computeDiagonalF(float* positions, int* indices, float* globalF, float* globalU, float* globalV, float* refShapeMatrixInverse, int trueNumConstraints)
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

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			float sum = 0.0f;

			for (int i = 0; i < 3; ++i)
			{
				sum += F_functionLevel[i][row] * F_functionLevel[i][col];
			}

			temp[row][col] = sum;
		}
	}

	//printf("%.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n",
	//	F_functionLevel[0][0], F_functionLevel[0][1], F_functionLevel[0][2],
	//	F_functionLevel[1][0], F_functionLevel[1][1], F_functionLevel[1][2],
	//	F_functionLevel[2][0], F_functionLevel[2][1], F_functionLevel[2][2]);

	//1. DIAGONALISE F--------------------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------------------------------------
	float Q[3][3];
	float w[3];
	// Sums of diagonal resp. off-diagonal elements
	float s = 0.0f;
	float c = 0.0f;
	float t = 0.0f;                 // sin(phi), cos(phi), tan(phi) and temporary storage
	float h = 0.0f;

	// Initialize Q to the identitity matrix
	Q[0][0] = 1.0f;
	Q[0][1] = 0.0f;
	Q[0][2] = 0.0f;

	Q[1][0] = 0.0f;
	Q[1][1] = 1.0f;
	Q[1][2] = 0.0f;

	Q[2][0] = 0.0f;
	Q[2][1] = 0.0f;
	Q[2][2] = 1.0f;

	// Initialize w to diag(temp)
	w[0] = temp[0][0];
	w[1] = temp[1][1];
	w[2] = temp[2][2];

	// Main iteration loop
	for (int nIter = 0; nIter < 4; ++nIter)
	{
		// Do sweep
		for (int p = 0; p < 3; p++)
		{
			for (int q = p + 1; q < 3; q++)
			{
				// Calculate Jacobi transformation
				h = (w[q] - w[p]);
				if (fabs(h) + 100.0f * fabs(temp[p][q]) == fabs(h))
				{
					if(h > JACOBI_LIMIT_1)
					{
						t = temp[p][q] / h;
					}
					else
					{
						t = 0.0f;
					}
					//printf("%.4f, h= %.4f, t= %.4f; ", t, h, temp[p][q]);
				}
				else
				{
					//theta = (0.5f * h / temp[p][q]);
					if ((0.5f * h / temp[p][q]) < 0.0f)
					{
						t = -1.0f / (sqrtf(1.0f + sqr((0.5f * h / temp[p][q]))) - (0.5f * h / temp[p][q]));
					}
					else
					{
						t = 1.0f / (sqrtf(1.0f + sqr((0.5f * h / temp[p][q]))) + (0.5f * h / temp[p][q]));
					}
				}
				c = 1.0f / sqrtf(1.0f + sqr(t));
				s = (t * c);

				// Apply Jacobi transformation
				w[p] -= (t * temp[p][q]);
				w[q] += (t * temp[p][q]);
				temp[p][q] = 0.0f;
				for (int r = 0; r < p; r++)
				{
					t = temp[r][p];
					temp[r][p] = c*t - s*temp[r][q];
					temp[r][q] = s*t + c*temp[r][q];
				}
				for (int r = p + 1; r < q; r++)
				{
					t = temp[p][r];
					temp[p][r] = c*t - s*temp[r][q];
					temp[r][q] = s*t + c*temp[r][q];
				}
				for (int r = q + 1; r < 3; r++)
				{
					t = temp[p][r];
					temp[p][r] = c*t - s*temp[q][r];
					temp[q][r] = s*t + c*temp[q][r];
				}

				// Update eigenvectors        
				Q[0][p] = c * Q[0][p] - s * Q[0][q];
				Q[0][q] = s * Q[0][p] + c * Q[0][q];
				Q[1][p] = c * Q[1][p] - s * Q[1][q];
				Q[1][q] = s * Q[1][p] + c * Q[1][q];
				Q[2][p] = c * Q[2][p] - s * Q[2][q];
				Q[2][q] = s * Q[2][p] + c * Q[2][q];
			}
		}
	}

	for (int nIter = 4; nIter < 100; ++nIter)
	{
		// Test for convergence 
		if (fabs(temp[0][1]) + fabs(temp[0][2])
			+ fabs(temp[1][2])
			+ fabs(temp[2][2]) == 0.0f)
		{
			break;
		}

		// Do sweep
		for (int p = 0; p < 3; p++)
		{
			for (int q = p + 1; q < 3; q++)
			{
				if (fabs(w[p]) + 100.0f * fabs(temp[p][q]) == fabs(w[p])
					&& fabs(w[q]) + 100.0f * fabs(temp[p][q]) == fabs(w[q]))
				{
					temp[p][q] = 0.0f;
				}
				else if (fabs(temp[p][q]) > 0.0f)
				{
					// Calculate Jacobi transformation
					h = w[q] - w[p];
					if (fabs(h) + 100.0f * fabs(temp[p][q]) == fabs(h))
					{
						if (h > JACOBI_LIMIT_1)
						{
							t = temp[p][q] / h;
						}
						else
						{
							t = 0.0f;
						}
						//printf("%.4f, h= %.4f, t= %.4f; ", t, h, temp[p][q]);
						//t = temp[p][q] / h;
					}
					else
					{
						//theta = 0.5f * h / temp[p][q];
						if ((0.5f * h / temp[p][q]) < 0.0f)
						{
							t = -1.0f / (sqrtf(1.0f + sqr((0.5f * h / temp[p][q]))) - (0.5f * h / temp[p][q]));
						}
						else
						{
							t = 1.0f / (sqrtf(1.0f + sqr((0.5f * h / temp[p][q]))) + (0.5f * h / temp[p][q]));
						}
					}
					c = 1.0f / sqrtf(1.0f + sqr(t));
					s = t * c;

					// Apply Jacobi transformation
					w[p] -= (t * temp[p][q]);
					w[q] += (t * temp[p][q]);
					temp[p][q] = 0.0f;
					for (int r = 0; r < p; r++)
					{
						t = temp[r][p];
						temp[r][p] = c*t - s*temp[r][q];
						temp[r][q] = s*t + c*temp[r][q];
					}
					for (int r = p + 1; r < q; r++)
					{
						t = temp[p][r];
						temp[p][r] = c*t - s*temp[r][q];
						temp[r][q] = s*t + c*temp[r][q];
					}
					for (int r = q + 1; r < 3; r++)
					{
						t = temp[p][r];
						temp[p][r] = c*t - s*temp[q][r];
						temp[q][r] = s*t + c*temp[q][r];
					}

					// Update eigenvectors        
					Q[0][p] = c * Q[0][p] - s * Q[0][q];
					Q[0][q] = s * Q[0][p] + c * Q[0][q];
					Q[1][p] = c * Q[1][p] - s * Q[1][q];
					Q[1][q] = s * Q[1][p] + c * Q[1][q];
					Q[2][p] = c * Q[2][p] - s * Q[2][q];
					Q[2][q] = s * Q[2][p] + c * Q[2][q];
				}
			}
		}
	}

	//printf("%.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n",
	//	Q[0][0], Q[0][1], Q[0][2],
	//	Q[1][0], Q[1][1], Q[1][2],
	//	Q[2][0], Q[2][1], Q[2][2]);

	//printf("%.4f, %.4f, %.4f", w[0], w[1], w[2]);

	//Correct diagonalised F
	//for (int i = 0; i < 3; ++i)
	//{
	//	w[i] = fmaxf(w[i], 0.0f);
	//}

	float determinantQ = Q[0][0]
		* (Q[1][1] * Q[2][2] - Q[1][2] * Q[2][1])
		- Q[0][1]
		* (Q[1][0] * Q[2][2] - Q[1][2] * Q[2][0])
		+ Q[0][2]
		* (Q[1][0] * Q[2][1] - Q[1][1] * Q[2][0]);

	//printf("%.4f, ", determinantQ);

	//remove reflection from V if necessary
	if (determinantQ < 0.0f)
	{
		float minElementValue = 1.5e+38f;
		int minElementIdx = 111;
		for (int i = 0; i < 3; ++i)
		{
			if (w[i] < minElementValue)
			{
				minElementValue = w[i];
				minElementIdx = i;
			}
		}
		for (int row = 0; row < 3; ++row)
		{
			Q[row][minElementIdx] *= -1.0f;
		}
	}

	//printf("%.4f, %.4f, %.4f", w[0], w[1], w[2]);

	//determine entries of F
	for (int i = 0; i < 3; ++i)
	{
		w[i] = sqrtf(w[i]);
	}

	//printf("%.4f, %.4f, %.4f", w[0], w[1], w[2]);

	//temp = F_orig * V * F.inverse();
	float sum;
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			sum = 0.0f;
			for (int i = 0; i < 3; i++)
			{
				sum += F_functionLevel[row][i] * Q[i][col];
			}

			temp[row][col] = sum;
		}
	}

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			globalV[GLOBAL_Vt_IDX] = Q[row][col];
			//printf("%.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n %.4f, %.4f, %.4f \n",
			//	Q[0][0], Q[0][1], Q[0][2],
			//	Q[1][0], Q[1][1], Q[1][2],
			//	Q[2][0], Q[2][1], Q[2][2]);
		}
	}

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			temp[row][col] *= (1.0f / w[col]);
		}
	}

	int numEntriesNearZero = 0;
	int idx = 0;
	for (int i = 0; i < 3; ++i)
	{
		if (fabs(w[i]) < 1.0e-4f)
		{
			++numEntriesNearZero;
			idx = i;
		}
	}

	if (numEntriesNearZero > 0)
	{
		if (numEntriesNearZero == 1)
		{
			if (idx == 0)
			{
				//temp.col(0) = temp.col(1).cross(temp.col(2)).normalized();
				temp[0][0] = temp[1][1] * temp[2][2] - temp[2][1] * temp[1][2];
				temp[1][0] = temp[2][1] * temp[0][2] - temp[0][1] * temp[2][2];
				temp[2][0] = temp[0][1] * temp[1][2] - temp[1][1] * temp[0][2];
				temp[0][0] /= sqrtf((temp[0][0] * temp[0][0]) + (temp[1][0] * temp[1][0]) + (temp[2][0] * temp[2][0]));
				temp[1][0] /= sqrtf((temp[0][0] * temp[0][0]) + (temp[1][0] * temp[1][0]) + (temp[2][0] * temp[2][0]));
				temp[2][0] /= sqrtf((temp[0][0] * temp[0][0]) + (temp[1][0] * temp[1][0]) + (temp[2][0] * temp[2][0]));
			}
			else if (idx == 1)
			{
				//temp.col(1) = temp.col(0).cross(temp.col(2)).normalized();
				temp[0][1] = temp[1][0] * temp[2][2] - temp[2][0] * temp[1][2];
				temp[1][1] = temp[2][0] * temp[0][2] - temp[0][0] * temp[2][2];
				temp[2][1] = temp[0][0] * temp[1][2] - temp[1][0] * temp[0][2];
				temp[0][1] /= sqrtf((temp[0][1] * temp[0][1]) + (temp[1][1] * temp[1][1]) + (temp[2][1] * temp[2][1]));
				temp[1][1] /= sqrtf((temp[0][1] * temp[0][1]) + (temp[1][1] * temp[1][1]) + (temp[2][1] * temp[2][1]));
				temp[2][1] /= sqrtf((temp[0][1] * temp[0][1]) + (temp[1][1] * temp[1][1]) + (temp[2][1] * temp[2][1]));
			}
			else
			{
				//temp.col(2) = temp.col(0).cross(temp.col(1)).normalized();
				temp[0][2] = temp[1][0] * temp[2][1] - temp[2][0] * temp[1][1];
				temp[1][2] = temp[2][0] * temp[0][1] - temp[0][0] * temp[2][1];
				temp[2][2] = temp[0][0] * temp[1][1] - temp[1][0] * temp[0][1];
				temp[0][2] /= sqrtf((temp[0][2] * temp[0][2]) + (temp[1][2] * temp[1][2]) + (temp[2][2] * temp[2][2]));
				temp[1][2] /= sqrtf((temp[0][2] * temp[0][2]) + (temp[1][2] * temp[1][2]) + (temp[2][2] * temp[2][2]));
				temp[2][2] /= sqrtf((temp[0][2] * temp[0][2]) + (temp[1][2] * temp[1][2]) + (temp[2][2] * temp[2][2]));
			}
		}
		else
		{
			//set temp to identity
			temp[0][0] = 1.0f;
			temp[1][1] = 1.0f;
			temp[2][2] = 1.0f;

			temp[0][1] = 0.0f;
			temp[0][2] = 0.0f;

			temp[1][0] = 0.0f;
			temp[1][2] = 0.0f;

			temp[2][0] = 0.0f;
			temp[2][1] = 0.0f;
		}
	}

	float determinantU = temp[0][0]
		* (temp[1][1] * temp[2][2] - temp[1][2] * temp[2][1])
		- temp[0][1]
		* (temp[1][0] * temp[2][2] - temp[1][2] * temp[2][0])
		+ temp[0][2]
		* (temp[1][0] * temp[2][1] - temp[1][1] * temp[2][0]);

	//remove reflection from temp if necessary
	if (determinantU < 0.0f)
	{
		float minElementValue = 1.5e+38f;
		int minElementIdx = 111;
		for (int i = 0; i < 3; ++i)
		{
			if (w[i] < minElementValue)
			{
				minElementValue = w[i];
				minElementIdx = i;
			}
		}

		w[minElementIdx] *= -1.0f;
		for (int row = 0; row < 3; ++row)
		{
			temp[row][minElementIdx] *= -1.0f;
		}
	}

	for (int i = 0; i < 3; i++)
	{
		w[i] = fmaxf(w[i], 0.577f);
	}

	//Store temp, V & F
	globalF[GLOBAL_F_IDX_0] = w[0];
	globalF[GLOBAL_F_IDX_1] = w[1];
	globalF[GLOBAL_F_IDX_2] = w[2];

	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			globalU[GLOBAL_U_IDX] = temp[row][col];
		}
	}
}


cudaError_t projectConstraints(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	float* device_F, float* device_U, float* device_V,
	float* device_anisotropyDirection, float* device_viscosityMultiplier,
	const Parameters& settings);

int CUDA_projectConstraints(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	float* device_F, float* device_U, float* device_V,
	float* device_anisotropyDirection, float* device_viscosityMultiplier,
	const Parameters& settings)
{
	//GPU
	cudaError_t cudaStatus = projectConstraints(device_indices, device_positions,
		device_inverseMasses, device_refShapeMatrixInverses, device_volumes, device_F, device_U, device_V,
		device_anisotropyDirection, device_viscosityMultiplier,
		settings);

	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Critical Error, aborting...");
		return 1;
	}

	return 0;
}

cudaError_t cudaErrorWrapper(const cudaError_t& status)
{
	if (status != cudaSuccess)
	{
		std::cerr << "ERROR: " << cudaGetErrorString(status) << std::endl;
	}

	return status;
}

bool checkCudaErrorStatus(const cudaError_t& status)
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
	float* device_F, float* device_U, float* device_V,
	float* device_anisotropyDirection, float* device_viscosityMultiplier,
	const Parameters& settings)
{
	cudaError_t deviceStatus;

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
		computeDiagonalF << <numBlocks, numThreads >> >(
			device_positions, device_indices,
			device_F, device_U, device_V,
			device_refShapeMatrixInverses,
			settings.trueNumberOfConstraints);

		cudaErrorWrapper(cudaDeviceSynchronize());
		//std::cout << "Projection-------------------------------------------------------------------------" << std::endl;
		switch (settings.materialModel)
		{
			case Parameters::CONSTITUTIVE_MODEL::NEO_HOOKEAN:
			{
				solveFEMConstraint << <numBlocks, numThreads >> >(
					device_positions, device_indices, device_inverseMasses,
					device_volumes, device_refShapeMatrixInverses,
					settings.lambda, settings.mu, settings.trueNumberOfConstraints,
					settings.numParticles,
					device_F, device_U, device_V);
				break;
			}
			case Parameters::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER:
			{
				solveFEMConstraint_ANISOTROPIC << <numBlocks, numThreads >> >(
					device_positions, device_indices, device_inverseMasses,
					device_volumes, device_refShapeMatrixInverses,
					settings.lambda, settings.mu, settings.trueNumberOfConstraints,
					settings.numParticles,
					device_F, device_U, device_V,
					device_anisotropyDirection, settings.anisotropyStrength);
				break;
			}
			case Parameters::CONSTITUTIVE_MODEL::NEO_HOOKEAN_FIBER_VISCOELASTIC:
			{

			}
		}

		cudaErrorWrapper(cudaDeviceSynchronize());
		//std::cout << "Projection done-------------------------------------------------------------------------" << std::endl;
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
	std::vector<float>& viscosityMultiplier)
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
	deviceStatus = cudaMalloc((void**)device_F, F.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_U, U.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_V, V.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_anisotropyDirection, anisotropyDirection.size() * sizeof(float));
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMalloc((void**)device_viscosityMultiplier, viscosityMultiplier.size() * sizeof(float));
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
	deviceStatus = cudaMemcpy(*device_F, &F[0], F.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_U, &U[0], U.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_V, &V[0], V.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_anisotropyDirection, &anisotropyDirection[0], anisotropyDirection.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}
	deviceStatus = cudaMemcpy(*device_viscosityMultiplier, &viscosityMultiplier[0], viscosityMultiplier.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}


	//Print some statistics
	double totalNumbytes = indices.size() * sizeof(int)
		+positions.size() * sizeof(float)
		+volumes.size() * sizeof(float)
		+inverseMasses.size() * sizeof(float)
		+refShapeMatrixInverses.size() * sizeof(float)
		+F.size() * sizeof(float)
		+U.size() * sizeof(float)
		+V.size() * sizeof(float);
		+anisotropyDirection.size() * sizeof(float);
		+viscosityMultiplier.size() * sizeof(float);

	std::cout << "Memory Usage: " << std::endl;
	std::cout << "	Indices (int32)             : " << indices.size() * sizeof(int) << " bytes" << std::endl;
	std::cout << "	Positions (float)           : " << positions.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Volumes (float)             : " << volumes.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Masses (float)              : " << inverseMasses.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Ref. Shape Matrices: (float): " << refShapeMatrixInverses.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	F                           : " << F.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	U                           : " << U.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	V                           : " << V.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Anisotropy Directions       : " << anisotropyDirection.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	Viscosity Multipliers       : " << viscosityMultiplier.size() * sizeof(float) << " bytes" << std::endl;
	std::cout << "	--------------------------------------------------------" << std::endl;
	std::cout << "	Total                       : " << totalNumbytes << " bytes (" << totalNumbytes / 1000.0 << " kb)" << std::endl;
}

bool CUDA_destroyBuffers(int* device_indices, float* device_positions,
	float* device_inverseMasses, float* device_refShapeMatrixInverses,
	float* device_volumes,
	float* device_F, float* device_U, float* device_V,
	float* device_anisotropyDirection, float* device_viscosityMultiplier)
{
	cudaError_t deviceStatus;

	cudaFree(device_positions);
	cudaFree(device_inverseMasses);
	cudaFree(device_indices);
	cudaFree(device_refShapeMatrixInverses);
	cudaFree(device_volumes);
	cudaFree(device_F);
	cudaFree(device_U);
	cudaFree(device_V);
	cudaFree(device_anisotropyDirection);
	cudaFree(device_viscosityMultiplier);

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

bool CUDA_updateBuffers(float* device_positions, std::vector<float>& positions,
	float* device_anisotropyDirection, std::vector<float>& anisotropyDirection, bool updateAnisotropyDirection)
{
	cudaError_t deviceStatus;

	deviceStatus = cudaMemcpy(device_positions, &positions[0], positions.size() * sizeof(float), cudaMemcpyHostToDevice);
	if (!checkCudaErrorStatus(deviceStatus))
	{
		return false;
	}

	if (updateAnisotropyDirection)
	{
		deviceStatus = cudaMemcpy(device_anisotropyDirection, &anisotropyDirection[0], anisotropyDirection.size() * sizeof(float), cudaMemcpyHostToDevice);
		if (!checkCudaErrorStatus(deviceStatus))
		{
			return false;
		}
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