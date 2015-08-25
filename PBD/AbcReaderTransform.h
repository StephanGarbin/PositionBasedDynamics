#pragma once

#include <string>
#include <vector>

#include <Eigen\Dense>

#include <Alembic\Abc\All.h>

struct AbcTransformReaderImp;

class AbcReaderTransform
{
public:
	AbcReaderTransform();
	~AbcReaderTransform();

	bool openArchive(const std::string& file, const std::vector<std::string>& transformNames);

	bool sampleForward(int idx);
	bool sampleBackward(int idx);
	bool sampleSpecific(int idx, int sample);

	int getNumSamples(int idx);

	//Data Accessors
	const Eigen::Vector3f& getTranslation(int idx){ return m_translation[idx]; }
	const Eigen::Vector3f& getRotation(int idx){ return m_rotation[idx]; }
	const Eigen::Matrix3f& getScale(int idx){ return m_scale[idx]; }

private:

	void readCurrentSampleIntoMemory(int idx);

	std::shared_ptr<AbcTransformReaderImp> m_data;
	std::vector<Eigen::Vector3f> m_translation;
	std::vector<Eigen::Vector3f> m_rotation;
	std::vector<Eigen::Matrix3f> m_scale;
};

