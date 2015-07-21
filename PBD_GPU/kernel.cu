
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>

const int NUM_THREADS_PER_BLOCK_SINGLE = 8;
const int NUM_THREADS_PER_BLOCK = NUM_THREADS_PER_BLOCK_SINGLE * NUM_THREADS_PER_BLOCK_SINGLE;

__shared__ float F[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float FTransposeF[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float FInverseTranspose[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float FirstPiolaKirchoffTensor[NUM_THREADS_PER_BLOCK][3][3];
__shared__ float Gradient[NUM_THREADS_PER_BLOCK][3][4];
__shared__ float DeltaX[NUM_THREADS_PER_BLOCK][3];

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

__device__ void calculateF(int idx)
{

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
			//Gradient[idx][row][col] = F[idx][row][col];
		}
	}

	//2. Multiply by volume

	//3. Multiply with First Piola-Kirchoff Stress tensor

	//4. Calculate last column

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

			FTransposeF[idx][row][col] = sum;
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

__global__ void solveFEMConstraint(float* positions, float* volume, float* refShapeMatrixInverse,
	float lambda, float mu)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//1. Calculate Deformation Gradient F
	calculateF(idx);

	//2. Compute Cauchy Tensors
	calculateFInverseTranspose(idx);
	calculateFTransposeF(idx);

	//3. Compute Invariants
	float I1 = traceFTransposeF(idx);
	float I3 = determinantFTransposeF(idx);

	//4. Calculate First Piola-Kirchoff Stress Tensor
	calculateFirstPiolaKirchoffTensor_NEO_HOOKEAN(idx, mu, lambda, I3);

}



cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
	int i = threadIdx.x;
	c[i] = a[i] + b[i];
}




int main()
{
	const int arraySize = 5;
	const int a[arraySize] = { 1, 2, 3, 4, 5 };
	const int b[arraySize] = { 10, 20, 30, 40, 50 };
	int c[arraySize] = { 0 };

	// Add vectors in parallel.
	cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addWithCuda failed!");
		return 1;
	}

	printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
		c[0], c[1], c[2], c[3], c[4]);

	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}

	return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
	int *dev_a = 0;
	int *dev_b = 0;
	int *dev_c = 0;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	// Launch a kernel on the GPU with one thread for each element.
	addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

Error:
	cudaFree(dev_c);
	cudaFree(dev_a);
	cudaFree(dev_b);
	
	return cudaStatus;
}
