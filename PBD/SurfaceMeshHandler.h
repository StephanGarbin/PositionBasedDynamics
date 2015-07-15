#pragma once

#include <string>
#include <vector>
#include <memory>

#include <Eigen\Dense>
#include <Alembic\Abc\All.h>

#include "AbcWriter.h"


class SurfaceMeshHandler
{
public:
	SurfaceMeshHandler(const std::string& surfaceMeshFile);
	~SurfaceMeshHandler();

	void setSample(const std::vector<Eigen::Vector3f>& positions);

private:
	int m_numVerticesInMesh;
	std::string m_surfaceMeshFileName;
	std::shared_ptr<AbcWriter> m_AbcWriter;
	std::vector<Alembic::Abc::V3f> m_faceIndices;
	std::vector<int> m_faceCounts;
};

