#include "AbcReader.h"

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
AbcReader::openArchive(const std::string& file)
{
	Alembic::AbcCoreFactory::IFactory factory;
	factory.setPolicy(Alembic::Abc::ErrorHandler::kQuietNoopPolicy);
	Alembic::AbcCoreFactory::IFactory::CoreType coreType;

	m_data->archive =
		std::make_shared<Alembic::Abc::IArchive>(factory.getArchive(file, coreType));

	Alembic::Abc::IObject topObject(*m_data->archive, Alembic::Abc::kTop);

	m_data->transform =
		std::make_shared<Alembic::AbcGeom::IXform>(Alembic::AbcGeom::IObject(topObject, "XForm"), Alembic::AbcGeom::kWrapExisting);

	
	m_data->mesh = std::make_shared<Alembic::AbcGeom::IPolyMesh>(Alembic::AbcGeom::IObject(*m_data->transform, "deformedMesh"), Alembic::AbcGeom::kWrapExisting);
	

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

	//1. Positions
	int numPositions = sample.getPositions()->size();
	m_positions.resize(numPositions);
	for (int i = 0; i < numPositions; ++i)
	{
		m_positions[i][0] = (*sample.getPositions())[i][0];
		m_positions[i][1] = (*sample.getPositions())[i][1];
		m_positions[i][2] = (*sample.getPositions())[i][2];
	}

	//2. FaceIndices
	//int numFaceIndices = sample.getFaceIndices()->size();
	//m_positions.resize(numFaceIndices);
	//for (int i = 0; i < numFaceIndices; ++i)
	//{
	//	m_positions[i][0] = (*sample.getPositions())[i][0];
	//	m_positions[i][1] = (*sample.getPositions())[i][1];
	//	m_positions[i][2] = (*sample.getPositions())[i][2];
	//}

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