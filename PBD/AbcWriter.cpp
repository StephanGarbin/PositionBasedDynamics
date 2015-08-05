#include "AbcWriter.h"

#include <iostream>

#include <Alembic\AbcGeom\All.h>
#include <Alembic\AbcCoreHDF5\All.h>
#include <Alembic\AbcCoreOgawa\All.h>
#include <Alembic\AbcCoreFactory\All.h>


struct AbcWriterImp
{
	std::shared_ptr<Alembic::Abc::OArchive> archive;
	//std::shared_ptr<Alembic::Abc::OObject> topObject;
	std::shared_ptr<Alembic::AbcGeom::OPolyMesh> mesh;
	std::shared_ptr<Alembic::AbcGeom::OXform> transform;
};

AbcWriter::AbcWriter(const std::string& fileName, const std::string& objectName) : m_archiveName(fileName), m_objectName(objectName)
{
	m_data = std::make_shared<AbcWriterImp>();
	m_data->archive = std::make_shared<Alembic::Abc::OArchive>(Alembic::AbcCoreOgawa::WriteArchive(), m_archiveName);
	//m_data->topObject = std::make_shared<Alembic::Abc::OObject>(m_data->archive->getTop());
	m_data->transform = std::make_shared<Alembic::AbcGeom::OXform>(Alembic::Abc::OObject(m_data->archive->getTop()), "XForm");
	m_data->mesh = std::make_shared<Alembic::AbcGeom::OPolyMesh>(*m_data->transform, m_objectName);

	if (m_data->archive->valid())
	{
		m_fileIsOpen = true;
	}
}


AbcWriter::~AbcWriter()
{
	if (m_fileIsOpen)
	{
		std::cout << "Writing [ " << m_data->archive->getName() << " ]" << std::endl;
		std::cout << "Num Samples Saved: " << m_data->mesh->getSchema().getNumSamples() << std::endl;
	}
}


bool
AbcWriter::addSample(std::vector<Alembic::Abc::V3f>& vertices, std::vector<int>& faceIndices, std::vector<int>& faceCounts)
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
	sample.setPositions(Alembic::Abc::P3fArraySample(&vertices[0], vertices.size()));

	//FACE-INDICES
	sample.setFaceIndices(Alembic::Abc::Int32ArraySample(&faceIndices[0], faceIndices.size()));

	//FACE COUNTS
	sample.setFaceCounts(Alembic::Abc::Int32ArraySample(&faceCounts[0], faceCounts.size()));

	//set mesh sample
	schema.set(sample);

	return true;
}
