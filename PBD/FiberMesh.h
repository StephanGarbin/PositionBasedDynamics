#pragma once

#include <vector>
#include <string>
#include <memory>

#include <Eigen\Dense>

#include "PBDTetrahedra3d.h"
#include "PBDParticle.h"

class FiberMesh
{
public:
	FiberMesh(std::shared_ptr<std::vector<PBDParticle>> particles, std::vector<PBDTetrahedra3d>* tetrahedra);
	~FiberMesh();


	void drawOpenGL();

	void generateTriangleMesh(float tubeThickness, float* vertices, int* faces, int numVertices, int numFaces);

	void solveActiveFiberConstraints(float stiffness);

	void generateFibersToFillCube(const Eigen::Vector3f& origin,
		const Eigen::Vector3f& rotation,
		const Eigen::Vector2f& dimension,
		int numFibersX, int numFibersY, float randomPerturbationStrength, bool addRandomPerturbation);

private:
	std::vector<float> m_vertices;
	std::vector<int> faces;

	std::shared_ptr<std::vector<PBDParticle>> m_particles;
	std::vector<PBDTetrahedra3d>* m_tetrahedra;
};

