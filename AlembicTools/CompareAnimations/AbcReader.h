#pragma once

#include <string>
#include <vector>

#include <Eigen\Dense>

#include <Alembic\Abc\All.h>

struct AbcReaderImp;

class AbcReader
{
public:
	AbcReader();
	~AbcReader();

	bool openArchive(const std::string& file);

	void sampleForward();
	void sampleBackward();
	void sampleSpecific(int sample);

	int getNumSamples();

	//Data Accessors
	std::vector<Eigen::Vector3f>& getPositions() { return m_positions; }

	std::vector<Eigen::Vector3f>& getVelocities() { return m_velocities; }

	std::vector<int>& getFaceIndices(int idx) { return m_faceIndices[idx]; }
private:

	void readCurrentSampleIntoMemory();

	std::shared_ptr<AbcReaderImp> m_data;

	std::vector<Eigen::Vector3f> m_positions;
	std::vector<Eigen::Vector3f> m_velocities;
	std::vector<std::vector<int>> m_faceIndices;
};

