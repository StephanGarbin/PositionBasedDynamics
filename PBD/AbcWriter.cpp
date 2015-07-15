#include "AbcWriter.h"

#include <iostream>

#include <Alembic\Abc\All.h>
#include <Alembic\AbcGeom\All.h>
#include <Alembic\AbcCoreHDF5\All.h>
#include <Alembic\AbcCoreOgawa\All.h>

struct AbcWriterImp
{
	std::shared_ptr<Alembic::Abc::OArchive> archive;
	std::shared_ptr<Alembic::AbcGeom::OPolyMesh> mesh;
};

AbcWriter::AbcWriter(const std::string& fileName, const std::string& objectName) : m_archiveName(fileName), m_objectName(objectName)
{
	m_data = std::make_shared<AbcWriterImp>();
	m_data->archive = std::make_shared<Alembic::Abc::OArchive>(Alembic::AbcCoreOgawa::WriteArchive(), m_archiveName);
	m_data->mesh = std::make_shared<Alembic::AbcGeom::OPolyMesh>(Alembic::Abc::OObject(*m_data->archive, Alembic::Abc::kTop), m_objectName);
}


AbcWriter::~AbcWriter()
{

}


bool
AbcWriter::addSample(std::vector<Eigen::Vector3f>& vertices, std::vector<int>& faceIndices)
{
	// Make sure that our file is open
	if (!m_fileIsOpen)
	{
		std::cout << "ERROR: Alembic Archive [" <<
			m_archiveName << "] is not open! Sample could not be written." << std::endl;
		return false;
	}

	//get schema
	Alembic::AbcGeom::OPolyMeshSchema& schema = m_data->mesh->getSchema();

	//create a sample
	Alembic::AbcGeom::OPolyMeshSchema::Sample sample;

	//POSITION
	//this is horrible, remove in future
	std::vector<Alembic::Abc::V3f> positions;
	positions.resize(vertices.size());
	for (int i = 0; i < vertices.size(); ++i)
	{
		positions[i][0] = vertices[i][0];
		positions[i][1] = vertices[i][1];
		positions[i][2] = vertices[i][2];
	}

	sample.setPositions(Alembic::Abc::P3fArraySample(&positions[0], positions.size()));

	//FACE-INDICES
	sample.setFaceIndices(Alembic::Abc::Int32ArraySample(&faceIndices[0], faceIndices.size()));

	//FACE COUNTS
	std::vector<int> faceCounts;
	faceCounts.resize(faceIndices.size());
	for (int i = 0; faceCounts.size(); ++i)
	{
		faceCounts[i] = 3;
	}
	sample.setFaceCounts(Alembic::Abc::Int32ArraySample(&faceCounts[0], faceCounts.size()));


	//set mesh sample
	schema.set(sample);

	return true;
}
