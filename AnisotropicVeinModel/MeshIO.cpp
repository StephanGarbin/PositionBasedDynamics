#include "MeshIO.h"

void loadMesh(const std::string& fileName, meshType& mesh)
{

	OpenMesh::IO::Options opt;
	//opt += OpenMesh::IO::Options::VertexColor;
	opt += OpenMesh::IO::Options::VertexNormal;
	//opt += OpenMesh::IO::Options::VertexTexCoord;
	//opt += OpenMesh::IO::Options::FaceColor;
	//opt += OpenMesh::IO::Options::FaceNormal;
	//opt += OpenMesh::IO::Options::FaceTexCoord;

	if (!OpenMesh::IO::read_mesh(mesh, fileName, opt))
	{
		std::cout << "ERROR: Cannot read: " << fileName << std::endl;
	}
}

void saveMesh(const std::string& fileName, meshType& mesh)
{
	OpenMesh::IO::Options opt;
	//opt += OpenMesh::IO::Options::VertexColor;
	opt += OpenMesh::IO::Options::VertexNormal;

	if (!OpenMesh::IO::write_mesh(mesh, fileName, opt))
	{
		std::cout << "ERROR: Cannot write: " << fileName << std::endl;
	}
}