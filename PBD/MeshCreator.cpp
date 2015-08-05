#include "MeshCreator.h"

#include <iostream>

MeshCreator::MeshCreator()
{
}


MeshCreator::~MeshCreator()
{
}


void
MeshCreator::generateTetBar(std::shared_ptr<std::vector<PBDParticle>>& particles, std::vector<PBDTetrahedra3d>& tets,
int width, int height, int depth)
{
	std::vector<Eigen::Vector3f> points(width*height*depth);

	for (unsigned int i = 0; i < width; i++)
	{
		for (unsigned int j = 0; j < height; j++)
		{
			for (unsigned int k = 0; k < depth; k++)
			{
				points[i*height*depth + j*depth + k] = 0.3f*Eigen::Vector3f((float)i, (float)j, (float)k);
			}
		}
	}

	std::vector<int> indices;
	for (unsigned int i = 0; i < width - 1; i++)
	{
		for (unsigned int j = 0; j < height - 1; j++)
		{
			for (unsigned int k = 0; k < depth - 1; k++)
			{
				// For each block, the 8 corners are numerated as:
				//     4*-----*7
				//     /|    /|
				//    / |   / |
				//  5*-----*6 |
				//   | 0*--|--*3
				//   | /   | /
				//   |/    |/
				//  1*-----*2
				unsigned int p0 = i*height*depth + j*depth + k;
				unsigned int p1 = p0 + 1;
				unsigned int p3 = (i + 1)*height*depth + j*depth + k;
				unsigned int p2 = p3 + 1;
				unsigned int p7 = (i + 1)*height*depth + (j + 1)*depth + k;
				unsigned int p6 = p7 + 1;
				unsigned int p4 = i*height*depth + (j + 1)*depth + k;
				unsigned int p5 = p4 + 1;

				// Ensure that neighboring tetras are sharing faces
				if ((i + j + k) % 2 == 1)
				{
					indices.push_back(p2); indices.push_back(p1); indices.push_back(p6); indices.push_back(p3);
					indices.push_back(p6); indices.push_back(p3); indices.push_back(p4); indices.push_back(p7);
					indices.push_back(p4); indices.push_back(p1); indices.push_back(p6); indices.push_back(p5);
					indices.push_back(p3); indices.push_back(p1); indices.push_back(p4); indices.push_back(p0);
					indices.push_back(p6); indices.push_back(p1); indices.push_back(p4); indices.push_back(p3);
				}
				else
				{
					indices.push_back(p0); indices.push_back(p2); indices.push_back(p5); indices.push_back(p1);
					indices.push_back(p7); indices.push_back(p2); indices.push_back(p0); indices.push_back(p3);
					indices.push_back(p5); indices.push_back(p2); indices.push_back(p7); indices.push_back(p6);
					indices.push_back(p7); indices.push_back(p0); indices.push_back(p5); indices.push_back(p4);
					indices.push_back(p0); indices.push_back(p2); indices.push_back(p7); indices.push_back(p5);
				}
			}
		}
	}

	int numVertices = width*height*depth;
	int numTets = indices.size() / 4;

	std::cout << "Number of tets: " << numTets << "\n";
	std::cout << "Number of vertices: " << numVertices << "\n";

	//Set Masses
	std::vector<float> InverseMasses(numVertices);

	for (unsigned int i = 0; i < numVertices; i++)
	{
		InverseMasses[i] = 1.0;
	}
	for (unsigned int i = 0; i < 1; i++)
	{
		for (unsigned int j = 0; j < height; j++)
		{
			for (unsigned int k = 0; k < depth; k++)
				InverseMasses[i*height*depth + j*depth + k] = 0.0f;
		}
	}


	//Generate Instances

	//Particles
	particles = std::make_shared<std::vector<PBDParticle>>();

	Eigen::Vector3f velocity;
	velocity.setZero();

	for (int i = 0; i < numVertices; ++i)
	{
		(*particles).push_back(PBDParticle(points[i], velocity, InverseMasses[i]));
	}

	//Tets
	for (int i = 0; i < numTets; ++i)
	{
		std::vector<int> localIndices(4);
		localIndices[0] = indices[i * 4 + 0];
		localIndices[1] = indices[i * 4 + 1];
		localIndices[2] = indices[i * 4 + 2];
		localIndices[3] = indices[i * 4 + 3];

		tets.push_back(PBDTetrahedra3d(localIndices, particles));
	}
}

