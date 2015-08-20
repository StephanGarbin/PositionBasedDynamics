#include "AbcReader.h"

#include <sstream>

#include <Alembic\AbcGeom\All.h>
#include <Alembic\AbcCoreHDF5\All.h>
#include <Alembic\AbcCoreOgawa\All.h>
#include <Alembic\AbcCoreFactory\All.h>

struct AbcReaderImp
{
	std::shared_ptr<Alembic::Abc::IArchive> archive;
	std::shared_ptr<Alembic::AbcGeom::IPolyMesh> mesh;
	std::shared_ptr<Alembic::AbcGeom::IXform> transform;

	int numSamples;
	int currentSample;
};

AbcReader::AbcReader()
{
	m_data = std::make_shared<AbcReaderImp>();
}


AbcReader::~AbcReader()
{
}


bool
AbcReader::openArchive(const std::string& file, const std::string& meshName)
{
	Alembic::AbcCoreFactory::IFactory factory;
	factory.setPolicy(Alembic::Abc::ErrorHandler::kQuietNoopPolicy);
	Alembic::AbcCoreFactory::IFactory::CoreType coreType;

	m_data->archive =
		std::make_shared<Alembic::Abc::IArchive>(factory.getArchive(file, coreType));

	Alembic::Abc::IObject topObject(*m_data->archive, Alembic::Abc::kTop);

	m_data->transform =
		std::make_shared<Alembic::AbcGeom::IXform>(Alembic::AbcGeom::IObject(topObject, meshName), Alembic::AbcGeom::kWrapExisting);

	std::stringstream ss;
	ss << meshName << "Shape";

	m_data->mesh = std::make_shared<Alembic::AbcGeom::IPolyMesh>(Alembic::AbcGeom::IObject(*m_data->transform, ss.str()), Alembic::AbcGeom::kWrapExisting);

	ss.clear();

	Alembic::AbcGeom::IPolyMeshSchema& schema = m_data->mesh->getSchema();

	std::cout << "Num Poly Mesh Schema Samples Read From file: " << schema.getNumSamples() << std::endl;
	m_data->numSamples = schema.getNumSamples();
	sampleSpecific(0);


	return true;
}

void
AbcReader::readCurrentSampleIntoMemory()
{
	//get the Sample
	Alembic::AbcGeom::ISampleSelector sampleSelector((Alembic::Abc::index_t)m_data->currentSample);

	Alembic::AbcGeom::IPolyMeshSchema::Sample sample;
	m_data->mesh->getSchema().get(sample, sampleSelector);

	Alembic::AbcGeom::XformSample transformSample;
	m_data->transform->getSchema().get(transformSample, sampleSelector);

	Eigen::Vector3f translationFromXForm(transformSample.getTranslation()[0],
		transformSample.getTranslation()[1], transformSample.getTranslation()[2]);
	Imath::V3d scale = transformSample.getScale();
	Eigen::Matrix3f scaleMatrix;
	scaleMatrix.setZero();
	scaleMatrix(0, 0) = scale[0];
	scaleMatrix(1, 1) = scale[0];
	scaleMatrix(2, 2) = scale[0];

	//1. Positions
	int numPositions = sample.getPositions()->size();
	m_positions.resize(numPositions);
	for (int i = 0; i < numPositions; ++i)
	{
		m_positions[i][0] = (*sample.getPositions())[i][0];
		m_positions[i][1] = (*sample.getPositions())[i][1];
		m_positions[i][2] = (*sample.getPositions())[i][2];
		//m_positions[i] += translationFromXForm.transpose() * scaleMatrix;
		//m_positions[i] += translationFromXForm.transpose();
		//m_positions[i] = m_positions[i].transpose() * scaleMatrix;

		//std::cout << m_positions[i] << std::endl;
	}

	//2. FaceIndices
	int numFaceIndices = sample.getFaceIndices()->size();
	m_faceIndices.resize(numFaceIndices);
	Alembic::Abc::Int32ArraySamplePtr faceIndices = sample.getFaceIndices();
	for (int i = 0; i < numFaceIndices / 3; ++i)
	{
		m_faceIndices[i].resize(3);
		m_faceIndices[i][0] = (*faceIndices)[i * 3 + 0];
		m_faceIndices[i][1] = (*faceIndices)[i * 3 + 1];
		m_faceIndices[i][2] = (*faceIndices)[i * 3 + 2];
	}

	//3. Velocities
	//int numVelocities = sample.getVelocities()->size();
	//m_velocities.resize(numVelocities);
	//for (int i = 0; i < numVelocities; ++i)
	//{
	//	m_velocities[i][0] = (*sample.getVelocities())[i][0];
	//	m_velocities[i][1] = (*sample.getVelocities())[i][1];
	//	m_velocities[i][2] = (*sample.getVelocities())[i][2];
	//}
}

bool
AbcReader::sampleForward()
{
	if (m_data->currentSample < m_data->numSamples)
	{
		m_data->currentSample += 1;
		readCurrentSampleIntoMemory();
		return true;
	}
	else
	{
		return false;
	}
}

bool
AbcReader::sampleBackward()
{
	if (m_data->currentSample >= 0)
	{
		m_data->currentSample -= 1;
		readCurrentSampleIntoMemory();
		return true;
	}
	else
	{
		return false;
	}
}

bool
AbcReader::sampleSpecific(int sample)
{
	if (sample >= 0 && sample < m_data->numSamples)
	{
		m_data->currentSample = sample;
		readCurrentSampleIntoMemory();
		return true;
	}
	else
	{
		return false;
	}
}

int
AbcReader::getNumSamples()
{
	return m_data->numSamples;
}