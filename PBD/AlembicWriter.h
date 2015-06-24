#pragma once

#include <string>
#include <vector>
#include <memory>

#include <Alembic\AbcGeom\All.h>
#include <Alembic\AbcCoreHDF5\All.h>

class AlembicWriter
{
public:
	AlembicWriter(const std::string& archive);
	~AlembicWriter();

	void writeSample();

private:
	Alembic::AbcGeom::OArchive* m_archive;
	Alembic::AbcGeom::OPolyMesh* m_meshObject;
};

