#include "FiberMesh.h"

#include <iostream>


FiberMesh::FiberMesh(std::shared_ptr<std::vector<PBDParticle>> particles, std::vector<PBDTetrahedra3d>* tetrahedra)
{
}


FiberMesh::~FiberMesh()
{
}

void
FiberMesh::drawOpenGL()
{

}

void
FiberMesh::generateTriangleMesh(float tubeThickness, float* vertices, int* faces, int numVertices, int numFaces)
{

}

void
FiberMesh::solveActiveFiberConstraints(float stiffness)
{

}

void
FiberMesh::generateFibersToFillCube(const Eigen::Vector3f& origin,
const Eigen::Vector3f& rotation,
const Eigen::Vector2f& dimension,
int numFibersX, int numFibersY, float randomPerturbationStrength, bool addRandomPerturbation)
{
	//generate origin points in xz plane
	float stepSizeX = dimension.x() / (numFibersX + 2);
	float stepSizeY = dimension.y() / (numFibersY + 2);

	std::vector<Eigen::Vector3f> originPoints(numFibersX * numFibersY);

	//generate regular grid of points ([[ADD RANDOM LATER!]])
	for (int row = 1; row < numFibersX; ++row)
	{
		for (int col = 1; col < numFibersY; ++col)
		{
			originPoints[row * numFibersX + numFibersY] = Eigen::Vector3f(stepSizeX * row, 0.0f, stepSizeY * col);
		}
	}

	//rotate points into position
	Eigen::Matrix3f rotateX;
	rotateX.setZero();
	rotateX(0, 0) = 1.0f;
	rotateX(1, 1) = cos(rotation.x());
	rotateX(1, 2) = sin(rotation.x());
	rotateX(2, 1) = -sin(rotation.x());
	rotateX(2, 2) = cos(rotation.x());
	
	Eigen::Matrix3f rotateY;
	rotateY.setZero();
	rotateY(1, 1) = 1.0f;
	rotateY(0, 0) = cos(rotation.y());
	rotateY(0, 2) = -sin(rotation.y());
	rotateY(2, 0) = sin(rotation.y());
	rotateY(2, 2) = cos(rotation.y());

	Eigen::Matrix3f rotateZ;
	rotateZ.setZero();
	rotateZ(2, 2) = 1.0f;
	rotateZ(0, 0) = cos(rotation.z());
	rotateZ(0, 1) = sin(rotation.z());
	rotateZ(1, 0) = -sin(rotation.z());
	rotateZ(1, 1) = cos(rotation.z());
#
	for (int i = 0; i < numFibersX * numFibersY; ++i)
	{
		originPoints[i] *= rotateX * rotateY * rotateZ;
	}

	Eigen::Vector3f V_01;
	Eigen::Vector3f V_02;

	Eigen::Vector3f n;

	Eigen::Vector3f direction = Eigen::Vector3f(1.0f, 0.0f, 0.0f).cross(Eigen::Vector3f(0.0f, 0.0f, 1.0f)) * rotateX * rotateY * rotateZ;


	//intersect with tets
	for (int i = 0; i < numFibersX * numFibersY; ++i)
	{
		//for each tet
		for (int t = 0; t < m_tetrahedra->size(); ++t)
		{
			int numIntersectionsInThisTetrahedron = 0;

			//for each face
			for (int f = 0; f < 4; ++f)
			{
				//Check PLANE Intersection
				V_01 = (*m_tetrahedra)[t].getFaceVertex(f, 0) - (*m_tetrahedra)[t].getFaceVertex(f, 1);
				V_02 = (*m_tetrahedra)[t].getFaceVertex(f, 0) - (*m_tetrahedra)[t].getFaceVertex(f, 2);

				n = V_01.cross(V_02).normalized();

				float d = (-(*m_tetrahedra)[t].getFaceVertex(f, 0)).dot(n);

				float denominator = n.dot(direction);
				if (denominator == 0.0f)
				{
					continue;
				}
				
				float t = - (d + n.dot(originPoints[i])) / denominator;

				//Check POLYGON Intersection
				Eigen::Vector3f P = originPoints[i] + direction * t;
				float u_0 = P[1] - (*m_tetrahedra)[t].getFaceVertex(f, 0)[1];
				float v_0 = P[2] = (*m_tetrahedra)[t].getFaceVertex(f, 0)[2];
				float u_1 = (*m_tetrahedra)[t].getFaceVertex(f, 1)[1] - (*m_tetrahedra)[t].getFaceVertex(f, 0)[1];
				float u_2 = (*m_tetrahedra)[t].getFaceVertex(f, 2)[1] - (*m_tetrahedra)[t].getFaceVertex(f, 0)[1];
				float v_1 = (*m_tetrahedra)[t].getFaceVertex(f, 1)[2] - (*m_tetrahedra)[t].getFaceVertex(f, 0)[2];
				float v_2 = (*m_tetrahedra)[t].getFaceVertex(f, 2)[2] - (*m_tetrahedra)[t].getFaceVertex(f, 0)[2];
				float alpha;
				float beta;
				if (u_1 == 0.0f)
				{
					beta = u_0 / u_2;
					if (beta >= 0.0f && beta <= 1.0f)
					{
						alpha = (v_0 - beta * v_2) / v_1;
					}
				}
				else
				{
					beta = (v_0 * u_1 - u_0 * v_1) / (v_2 * u_1 - u_2 * v_1);
					if (beta >= 0.0f && beta <= 1.0f)
					{
						alpha = (u_0 - beta * u_2) / u_1;
					}
				}

				bool intersects = alpha >= 0.0f && beta >= 0.0f && (alpha + beta) <= 1.0f;

				if (intersects)
				{
					++numIntersectionsInThisTetrahedron;
				}
			}

			std::cout << "Num Intersections: " << numIntersectionsInThisTetrahedron << std::endl;
		}
	}

	//save intersection points
}