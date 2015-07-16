#pragma once

#include <string>
#include <vector>
#include <memory>

#include <Eigen\Dense>
#include <Alembic\Abc\All.h>

#include "AbcWriter.h"

#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"


class SurfaceMeshHandler
{
public:
	SurfaceMeshHandler(const std::string& surfaceMeshFile, const std::string& abcFile);
	~SurfaceMeshHandler();

	void initTopology(std::vector<PBDParticle>& particles, std::vector<PBDTetrahedra3d>& tets);

	void setSample(const std::vector<Eigen::Vector3f>& positions);

private:
	int m_numVerticesInMesh;
	std::string m_surfaceMeshFileName;
	std::shared_ptr<AbcWriter> m_AbcWriter;
	std::vector<Alembic::Abc::V3f> m_vertices;
	std::vector<int> m_faceIndices;
	std::vector<int> m_faceCounts;
};

