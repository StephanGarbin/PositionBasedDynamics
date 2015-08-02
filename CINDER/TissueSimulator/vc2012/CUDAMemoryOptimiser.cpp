#include "CUDAMemoryOptimiser.h"

#include <list>
#include <iostream>

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