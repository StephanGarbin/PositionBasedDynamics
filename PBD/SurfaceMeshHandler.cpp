#include "SurfaceMeshHandler.h"


SurfaceMeshHandler::SurfaceMeshHandler(const std::string& surfaceMeshFile, const std::string& abcFile) : m_surfaceMeshFileName(surfaceMeshFile)
{
	if (m_surfaceMeshFileName == "WRITE_TETS")
	{
		std::cout << "Writing Tets as if surface!" << std::endl;
	}

	m_AbcWriter = std::make_shared<AbcWriter>(abcFile, "deformedMesh");
}


SurfaceMeshHandler::~SurfaceMeshHandler()
{
}


void
SurfaceMeshHandler::setSample(const std::vector<Eigen::Vector3f>& positions)
{
	for (int i = 0; i < positions.size(); ++i)
	{
		for (int c = 0; c < 3; ++c)
		{
			m_vertices[i][c] = positions[i][c];
		}
	}

	m_AbcWriter->addSample(m_vertices, m_faceIndices, m_faceCounts);
}

void
SurfaceMeshHandler::initTopology(std::vector<PBDParticle>& particles, std::vector<PBDTetrahedra3d>& tets)
{
	m_vertices.resize(particles.size());

	int numElements = tets.size();
	
	m_faceCounts.resize(numElements);

	for (int i = 0; i < numElements; ++i)
	{
		m_faceCounts[i] = 4;
	}

	for (int i = 0; i < numElements; ++i)
	{
		for (int v = 0; v < 4; ++v)
		{
			m_faceIndices.push_back(tets[i].getVertexIndices()[v]);
		}
	}
}