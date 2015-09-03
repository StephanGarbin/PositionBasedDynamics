#pragma once

#include <string>
#include <vector>
#include <memory>

#include <Alembic\Abc\All.h>
#include <Eigen\Dense>

struct AbcWriterImp;

class AbcWriter
{
public:
	AbcWriter(const std::string& fileName, const std::string& objectName);
	~AbcWriter();

	bool addSample(std::vector<Alembic::Abc::V3f>& vertices, std::vector<int>& faceIndices, std::vector<int>& faceCounts);

private:
	std::string m_archiveName;
	std::string m_objectName;

	bool m_fileIsOpen;

	std::shared_ptr<AbcWriterImp> m_data;
};

