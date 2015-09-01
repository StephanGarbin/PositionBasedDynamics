//NOTE: The core code (i.e. what is found int the function 'generateTetBar') is DIRECTLY from Bender et al. (the 'Positions-Based Simulation of Continuous Materials' paper)
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
				points[i*height*depth + j*depth + k] = 0.3f * Eigen::Vector3f((float)i, (float)j, (float)k);
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

void
MeshCreator::generateTetBarToFit(std::shared_ptr<std::vector<PBDParticle>>& particles, std::vector<PBDTetrahedra3d>& tets,
int numDivsWidth, int numDivsHeight, int numDivsDepth,
Eigen::Vector2f bottomLeft_above, Eigen::Vector2f topLeft_above, Eigen::Vector2f topRight_above,
float thickness)
{
	//Eigen::Vector3f originBottomLeft_above = Eigen::Vector3f(bottomLeft_above[0], 0.0f, bottomLeft_above[1]);
	Eigen::Vector3f originBottomLeft_above = Eigen::Vector3f(0.0f, 0.0f, 0.0f);

	float width = std::abs(topRight_above[0] - topLeft_above[0]);
	float height = std::abs(topLeft_above[1] - bottomLeft_above[1]);
	float depth = thickness;

	float scalingFactor = 0.01;
	width *= scalingFactor;
	height *= scalingFactor;
	depth *= scalingFactor;

	//std::cout << width << ", " << height << ", " << depth << std::endl;

	//float width = 10;
	//float height = 10;
	//float depth = 10;

	std::vector<Eigen::Vector3f> points(numDivsWidth * numDivsDepth * numDivsHeight);

	for (unsigned int i = 0; i < numDivsWidth; i++)
	{
		for (unsigned int j = 0; j < numDivsHeight; j++)
		{
			for (unsigned int k = 0; k < numDivsDepth; k++)
			{
				Eigen::Vector3f point;
				point.setZero();
				point += originBottomLeft_above;
				point[0] += i * (width / (float)numDivsWidth);
				point[1] += k * (depth / (float)numDivsDepth);
				point[2] += j * (height / (float)numDivsHeight);

				points[i * numDivsHeight * numDivsDepth + j * numDivsDepth + k] =
					point;

				//points[i*height*depth + j*depth + k] = 0.3f * Eigen::Vector3f((float)i, (float)j, (float)k);
			}
		}
	}

	std::vector<int> indices;
	for (unsigned int i = 0; i < numDivsWidth - 1; i++)
	{
		for (unsigned int j = 0; j < numDivsHeight - 1; j++)
		{
			for (unsigned int k = 0; k < numDivsDepth - 1; k++)
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
				unsigned int p0 = i * numDivsHeight * numDivsDepth + j * numDivsDepth + k;
				unsigned int p1 = p0 + 1;
				unsigned int p3 = (i + 1) * numDivsHeight * numDivsDepth + j * numDivsDepth + k;
				unsigned int p2 = p3 + 1;
				unsigned int p7 = (i + 1) * numDivsHeight * numDivsDepth + (j + 1) * numDivsDepth + k;
				unsigned int p6 = p7 + 1;
				unsigned int p4 = i * numDivsHeight * numDivsDepth + (j + 1) * numDivsDepth + k;
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

	int numVertices = numDivsWidth * numDivsHeight * numDivsDepth;
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
		for (unsigned int j = 0; j < numDivsHeight; j++)
		{
			for (unsigned int k = 0; k < numDivsDepth; k++)
				InverseMasses[i * numDivsHeight * numDivsDepth + j * numDivsDepth + k] = 0.0f;
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

void
MeshCreator::generateSingleTet(std::shared_ptr<std::vector<PBDParticle>>& particles, std::vector<PBDTetrahedra3d>& tets,
int width, int height, int depth)
{
	//1. Generate 4 Nodes
	Eigen::Vector3f position;
	Eigen::Vector3f velocity;
	velocity.setZero();
	float invMass = 1.0f;

	position[0] = 1.0f;
	position[1] = 0.0f;
	position[2] = 0.0f;
	(*particles).push_back(PBDParticle(position, velocity, invMass));

	position[0] = 0.0f;
	position[1] = 1.0f;
	position[2] = 0.0f;
	(*particles).push_back(PBDParticle(position, velocity, invMass));

	position[0] = 0.0f;
	position[1] = 0.0f;
	position[2] = 1.0f;
	(*particles).push_back(PBDParticle(position, velocity, invMass));

	position[0] = 0.0f;
	position[1] = 0.0f;
	position[2] = 0.0f;
	(*particles).push_back(PBDParticle(position, velocity, invMass));

	//2. Constrain one face
	(*particles)[0].inverseMass() = 0.0f;
	(*particles)[1].inverseMass() = 0.0f;
	(*particles)[2].inverseMass() = 0.0f;
	(*particles)[3].inverseMass() = 0.0f;

	//3. Generate 1 Tet Element
	std::vector<int> indices = { 0, 1, 2, 3 };
	tets.push_back(PBDTetrahedra3d(std::move(indices), particles));
}