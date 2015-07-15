#pragma once

#include <string>
#include <vector>
#include <memory>

#include <Eigen\Dense>

struct AbcWriterImp;

class AbcWriter
{
public:
	AbcWriter(const std::string& fileName, const std::string& objectName);
	~AbcWriter();

	bool addSample(std::vector<Eigen::Vector3f>& vertices, std::vector<int>& faceIndices);

private:
	std::string m_archiveName;
	std::string m_objectName;

	bool m_fileIsOpen;

	std::shared_ptr<AbcWriterImp> m_data;
};

