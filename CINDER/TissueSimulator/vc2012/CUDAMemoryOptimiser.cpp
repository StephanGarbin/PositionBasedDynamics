#include "CUDAMemoryOptimiser.h"

#include <list>
#include <iostream>
#include <algorithm>

CUDAMemoryOptimiser::CUDAMemoryOptimiser()
{

}


CUDAMemoryOptimiser::~CUDAMemoryOptimiser()
{

}


void
CUDAMemoryOptimiser::optimiseTetrahedralIndexingBasedOnNodeMemory(std::vector<PBDParticle>& nodes,
std::vector<PBDTetrahedra3d>& tetrahedra)
{
	std::cout << "MEM OPTIMISER: Initial Score: " << calculateMemoryAlignmentScore(tetrahedra) << std::endl;

	std::vector<PBDTetrahedra3d> newTetOrder;

	std::list<PBDTetrahedra3d> unusedNodes;
	for (int i = 1; i < tetrahedra.size(); ++i)
	{
		unusedNodes.push_back(tetrahedra[i]);
	}

	newTetOrder.push_back(tetrahedra[0]);

	//now push tets back in order
	int currentTet = 0;
	while (!unusedNodes.empty())
	{
		//1. Find 'best fitting' tetrahedra in currently unused collection
		auto currentListElement = unusedNodes.begin();
		int currentBestScore = calculateMemoryAlignmentScore(newTetOrder[currentTet], *currentListElement);

		for (auto e = unusedNodes.begin(); e != unusedNodes.end(); ++e)
		{
			int localScore = calculateMemoryAlignmentScore(newTetOrder[currentTet], *e);
			if (localScore < currentBestScore)
			{
				currentBestScore = localScore;
				currentListElement = e;
			}
		}

		newTetOrder.push_back(*currentListElement);
		unusedNodes.erase(currentListElement);
		
		++currentTet;
	}

	tetrahedra = newTetOrder;
	std::cout << "MEM OPTIMISER: Corrected Score: " << calculateMemoryAlignmentScore(tetrahedra) << std::endl;
}

int
CUDAMemoryOptimiser::calculateMemoryAlignmentScore(const PBDTetrahedra3d& reference,
const PBDTetrahedra3d& target)
{
	int score = 0;
	
	//calculate score based on whether nodes are sequential
	for (int i = 0; i < 4; ++i)
	{
		score += std::abs(reference.getVertexIndices()[i] - target.getVertexIndices()[i] + 1);
	}

	return score;
}

int
CUDAMemoryOptimiser::calculateMemoryAlignmentScore(const std::vector<PBDTetrahedra3d>& order)
{
	int totalScore = 0;
	for (int i = 0; i < order.size() - 1; ++i)
	{
		totalScore += calculateMemoryAlignmentScore(order[i], order[i + 1]);
	}

	return totalScore;
}

void
CUDAMemoryOptimiser::optimiseTetrahedralIndexingBasedOnMemoryChunks(std::vector<PBDParticle>& nodes,
std::vector<PBDTetrahedra3d>& tetrahedra, int numThreadsPerBlock)
{
	for (int i = 0; i < tetrahedra.size() / numThreadsPerBlock; ++i)
	{
		int minParticleIdx = 100000;
		int maxParticleIdx = -100000;
		for (int t = i * numThreadsPerBlock; t < i * numThreadsPerBlock + numThreadsPerBlock; ++t)
		{
			for (int n = 0; n < 4; ++n)
			{
				int currentIdx = tetrahedra[t].getVertexIndices()[n];

				if (currentIdx < minParticleIdx)
				{
					minParticleIdx = currentIdx;
				}
				if (currentIdx > maxParticleIdx)
				{
					maxParticleIdx = currentIdx;
				}
			}
		}

		std::cout << "BLOCK [" << i << "]: " << minParticleIdx << " - " << maxParticleIdx << "( " << maxParticleIdx - minParticleIdx << " )" << std::endl;
	}

	for (int i = 0; i < nodes.size() / numThreadsPerBlock; ++i)
	{
		int count = 0;
		for (int t = 0; t < tetrahedra.size(); ++t)
		{
			bool fits = true;
			for (int n = 0; n < 4; ++n)
			{
				int currentIdx = tetrahedra[t].getVertexIndices()[n];
				if (currentIdx < i * numThreadsPerBlock || currentIdx > i * numThreadsPerBlock + numThreadsPerBlock)
				{
					fits = false;
				}
			}

			if (fits)
			{
				++count;
			}
		}
		std::cout << "NODE BLOCK [" << i << "]: " << count << std::endl;
	}

	//NUM REQUIRED INDICES PER BLOCK OF TETS
	std::cout << "Num Required Indices Per Block: " << std::endl;

	std::vector<std::list<int>> blockIndices;

	for (int i = 0; i < tetrahedra.size() / numThreadsPerBlock; ++i)
	{
		std::list<int> requiredIndices;

		for (int t = i * numThreadsPerBlock; t < i * numThreadsPerBlock + numThreadsPerBlock; ++t)
		{
			for (int n = 0; n < 4; ++n)
			{
				int currentIdx = tetrahedra[t].getVertexIndices()[n];

				if (std::find(requiredIndices.begin(), requiredIndices.end(), currentIdx) == requiredIndices.end())
				{
					requiredIndices.push_back(currentIdx);
				}
			}
		}

		std::cout << "NUM REQUIRED INDIICES BLOCK [" << i << "]: " << requiredIndices.size() << std::endl;
		blockIndices.push_back(requiredIndices);
	}

	//Now determine how much overlap there is
	for (int b = 0; b < blockIndices.size(); ++b)
	{
		std::cout << "FOR [" << b << "]:" << std::endl;
		for (int otherB = (b + 1) % blockIndices.size(); otherB != b; otherB = (otherB + 1) % blockIndices.size())
		{
			int numOverlappingValues = 0;

			for (auto i = blockIndices[b].begin(); i != blockIndices[b].end(); ++i)
			{
				if (std::find(blockIndices[otherB].begin(), blockIndices[otherB].end(), *i) != blockIndices[otherB].end())
				{
					++numOverlappingValues;
				}
			}

			float percentage = 0.0f;
			
			if (blockIndices[b].size() != 0)
			{
				percentage = (float)numOverlappingValues / (float)blockIndices[b].size();
			}
			std::cout << "	WITH [" << otherB << "]: " << percentage << "%" << std::endl;
		} 
	}
}