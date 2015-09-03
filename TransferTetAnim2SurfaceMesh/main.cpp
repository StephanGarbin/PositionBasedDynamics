#include <string>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "AbcReader.h"
#include "AbcWriter.h"
#include "MeshIO.h"

#include "OMDefinitions.h"

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		std::cout << "Please provide: surface mesh name, alembic archive name, target alembic archive name!" << std::endl;
		return 0;
	}

	std::string surfaceMeshName(argv[1]);
	std::string tetMeshName(argv[2]);
	std::string targetMeshName(argv[3]);

	//1. Read Surface Mesh
	meshType surfaceMesh;
	loadMesh(surfaceMeshName, surfaceMesh);

	//2. Read Alembic Archive
	AbcReader tetMeshSequence;
	tetMeshSequence.openArchive(tetMeshName);

	//4. Write Alembic Archive
	AbcWriter targetMeshSequence(targetMeshName, "animatedMesh");
	std::vector<Imath::V3f> vertexPositions;
	std::vector<int> faceIndices;
	std::vector<int> faceCounts;

	for (meshType::FaceIter f_it = surfaceMesh.faces_begin(); f_it != surfaceMesh.faces_end(); ++f_it)
	{
		int count = 0;
		for (meshType::FaceVertexIter v_it = surfaceMesh.fv_begin(*f_it); v_it != surfaceMesh.fv_end(*f_it); ++v_it)
		{
			faceIndices.push_back(v_it->idx());
			++count;
		}
		faceCounts.push_back(count);
	}

	std::reverse(faceIndices.begin(), faceIndices.end());
	
	std::cout << faceIndices.size() << " Face Indices Recorded." << std::endl;
	std::cout << faceCounts.size() << " Face Counts Recorded." << std::endl;
	std::cout << surfaceMesh.n_vertices() << " Vertices Recorded." << std::endl;

	vertexPositions.resize(surfaceMesh.n_vertices());

	std::cout << "Finished IO, converting sequence..." << std::endl;


	for (int s = 0; s < tetMeshSequence.getNumSamples(); ++s)
	{
		tetMeshSequence.sampleSpecific(s);
		for (int i = 0; i < tetMeshSequence.getPositions().size(); ++i)
		{
			if (std::isnan(vertexPositions[i][0]) || std::isinf(vertexPositions[i][0])
				|| std::isnan(vertexPositions[i][1]) || std::isinf(vertexPositions[i][1])
				|| std::isnan(vertexPositions[i][2]) || std::isinf(vertexPositions[i][2]))
			{
				//std::cout << "Warning: Nan/Inf detected at position " << i << " (Sample: " << s << ")" << std::endl;
			}
			vertexPositions[i][0] = tetMeshSequence.getPositions()[i][0];
			vertexPositions[i][1] = tetMeshSequence.getPositions()[i][1];
			vertexPositions[i][2] = tetMeshSequence.getPositions()[i][2];

			//std::cout << vertexPositions[i][0] << ", " << vertexPositions[i][1] << ", " << vertexPositions[i][2] << std::endl;
		}
		targetMeshSequence.addSample(vertexPositions, faceIndices, faceCounts);
	}

	std::cout << "Tidying up, ..." << std::endl;
}