#include "AbcReaderTransform.h"

#include <sstream>

#include <Alembic\AbcGeom\All.h>
#include <Alembic\AbcCoreHDF5\All.h>
#include <Alembic\AbcCoreOgawa\All.h>
#include <Alembic\AbcCoreFactory\All.h>

struct AbcTransformReaderImp
{
	std::shared_ptr<Alembic::Abc::IArchive> archive;
	std::vector<std::shared_ptr<Alembic::AbcGeom::IXform>> transform;

	std::vector<int> numSamples;
	std::vector<int> currentSample;
};

AbcReaderTransform::AbcReaderTransform()
{
	m_data = std::make_shared<AbcTransformReaderImp>();
}


AbcReaderTransform::~AbcReaderTransform()
{
}


bool
AbcReaderTransform::openArchive(const std::string& file, const std::vector<std::string>& transformNames)
{
	Alembic::AbcCoreFactory::IFactory factory;
	factory.setPolicy(Alembic::Abc::ErrorHandler::kQuietNoopPolicy);
	Alembic::AbcCoreFactory::IFactory::CoreType coreType;

	m_data->archive =
		std::make_shared<Alembic::Abc::IArchive>(factory.getArchive(file, coreType));

	Alembic::Abc::IObject topObject(*m_data->archive, Alembic::Abc::kTop);

	m_data->transform.resize(transformNames.size());

	m_translation.resize(transformNames.size());
	m_rotation.resize(transformNames.size());
	m_scale.resize(transformNames.size());

	for (int i = 0; i < transformNames.size(); ++i)
	{
		std::cout << "Finding node: " << transformNames[i] << std::endl;
		m_data->transform[i] =
			std::make_shared<Alembic::AbcGeom::IXform>(Alembic::AbcGeom::IObject(topObject, transformNames[i]), Alembic::AbcGeom::kWrapExisting);

		Alembic::AbcGeom::IXformSchema schema(m_data->transform[i]->getSchema());

		std::cout << "Num Transform Schema Samples for [ " << transformNames[i] << " ]: " << schema.getNumSamples() << std::endl;

		m_data->numSamples.push_back(schema.getNumSamples());
		m_data->currentSample.push_back(0);

		sampleSpecific(i, 0);
	}

	return true;
}

void
AbcReaderTransform::readCurrentSampleIntoMemory(int idx)
{
	//get the Sample
	Alembic::AbcGeom::ISampleSelector sampleSelector((Alembic::Abc::index_t)m_data->currentSample[idx]);
	Alembic::AbcGeom::XformSample transformSample;
	m_data->transform[idx]->getSchema().get(transformSample, sampleSelector);
	
	m_translation[idx] = Eigen::Vector3f(transformSample.getTranslation()[0],
		transformSample.getTranslation()[1], transformSample.getTranslation()[2]);

	m_rotation[idx] = Eigen::Vector3f(transformSample.getXRotation(),
		transformSample.getYRotation(), transformSample.getZRotation());

	Imath::V3d scale = transformSample.getScale();
	m_scale[idx].setZero();
	m_scale[idx](0, 0) = scale[0];
	m_scale[idx](1, 1) = scale[1];
	m_scale[idx](2, 2) = scale[2];
}

bool
AbcReaderTransform::sampleForward(int idx)
{
	if (m_data->currentSample[idx] < m_data->numSamples[idx])
	{
		m_data->currentSample[idx] += 1;
		readCurrentSampleIntoMemory(idx);
		return true;
	}
	else
	{
		return false;
	}
}

bool
AbcReaderTransform::sampleBackward(int idx)
{
	if (m_data->currentSample[idx] >= 0)
	{
		m_data->currentSample[idx] -= 1;
		readCurrentSampleIntoMemory(idx);
		return true;
	}
	else
	{
		return false;
	}
}

bool
AbcReaderTransform::sampleSpecific(int idx, int sample)
{
	if (sample >= 0 && sample < m_data->numSamples[idx])
	{
		m_data->currentSample[idx] = sample;
		readCurrentSampleIntoMemory(idx);
		return true;
	}
	else
	{
		return false;
	}
}

int
AbcReaderTransform::getNumSamples(int idx)
{
	return m_data->numSamples[idx];
}