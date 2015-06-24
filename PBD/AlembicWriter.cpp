#include "AlembicWriter.h"


AlembicWriter::AlembicWriter(const std::string& archive)
{
	m_archive = new Alembic::AbcGeom::OArchive(Alembic::AbcCoreHDF5::WriteArchive(), archive);
	m_meshObject = new Alembic::AbcGeom::OPolyMesh(Alembic::AbcGeom::OObject(archive, Alembic::AbcGeom::kTop), "beam");
}


AlembicWriter::~AlembicWriter()
{
	//Destructor writes the archive
	std::cout << "Writing: " << m_archive->getName() << std::endl;
	m_archive->~OArchive();
	delete m_meshObject;
}


void
AlembicWriter::writeSample()
{
	Alembic::AbcGeom::OPolyMeshSchema& meshSchema = m_meshObject->getSchema();

	// Set a mesh sample.
	// We're creating the sample inline here,
	// but we could create a static sample and leave it around,
	// only modifying the parts that have changed.

	// Alembic is strongly typed. P3fArraySample is for an array
	// of 32-bit points, which are the mesh vertices. g_verts etc.
	// are defined in MeshData.cpp.
	Alembic::AbcGeom::OPolyMeshSchema::Sample mesh_samp(
		Alembic::AbcGeom::P3fArraySample((const V3f *)g_verts, g_numVerts));

	// Set the sample twice.
	// Because the data is the same in both samples, Alembic will
	// store only one copy, but note that two samples have been set.
	meshSchema.set(mesh_samp);
}