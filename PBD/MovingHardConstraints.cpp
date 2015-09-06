#include "MovingHardConstraints.h"

#include "ConstraintsIO.h"

MovingHardConstraints::MovingHardConstraints()
{
	m_speed = 1.0f;
	m_reader = std::make_shared<AbcReaderTransform>();
}


MovingHardConstraints::~MovingHardConstraints()
{
}
void
MovingHardConstraints::readConstraintFiles(const std::vector<std::string>& files)
{
	m_constraintIndices.resize(files.size());

	for (int i = 0; i < files.size(); ++i)
	{
		ConstraintsIO::readMayaVertexConstraints(m_constraintIndices[i], files[i]);
	}

	m_previousPosition.resize(files.size());
	for (int i = 0; i < files.size(); ++i)
	{
		m_previousPosition[i].setZero();
	}
}

void
MovingHardConstraints::readFromAbc(const std::string& fileName, const std::vector<std::string>& topBottomTransformNames)
{
	m_reader->openArchive(fileName, topBottomTransformNames);
}

void
MovingHardConstraints::initialisePositionMasses(std::vector<PBDParticle>& positions)
{
	for (int i = 0; i < m_constraintIndices.size(); ++i)
	{
		for (int c = 0; c < m_constraintIndices[i].size(); ++c)
		{
			positions[m_constraintIndices[i][c]].inverseMass() = 0.0f;
		}
	}
}

void
MovingHardConstraints::updatePositions(std::vector<PBDParticle>& positions, int systemFrame, float timeStep, int locatorIdx)
{
	timeStep *= m_speed;

	Eigen::Vector3f position;

	int frame1 = std::floorf((float)systemFrame * timeStep);
	int frame2 = std::ceilf((float)systemFrame * timeStep);

	float weightFrame1 = ((float)systemFrame * timeStep) - (float)frame1;
	float weightFrame2 = 1.0f - weightFrame1;

	Eigen::Vector3f frame1Position;
	Eigen::Vector3f frame2Position;
	m_reader->sampleSpecific(locatorIdx, frame1);
	frame1Position = m_reader->getTranslation(locatorIdx);
	m_reader->sampleSpecific(locatorIdx, frame2);
	frame2Position = m_reader->getTranslation(locatorIdx);

	position = weightFrame2 * frame1Position + weightFrame1 * frame2Position;

	if (m_previousPosition[locatorIdx].squaredNorm() == 0.0f)
	{
		m_previousPosition[locatorIdx] = position;
	}
	//std::cout << position - m_previousPosition << std::endl;
	//for (int i = 0; i < m_constraintIndices.size(); ++i)
	//{

		for (int c = 0; c < m_constraintIndices[locatorIdx].size(); ++c)
		{
			positions[m_constraintIndices[locatorIdx][c]].position() += position - m_previousPosition[locatorIdx];
			positions[m_constraintIndices[locatorIdx][c]].previousPosition() = positions[m_constraintIndices[locatorIdx][c]].position();
		}
	//}

	m_previousPosition[locatorIdx] = position;
}