#include "CollisionRod.h"

#include <GL\glew.h>
#include <gl\GL.h>
#include <GL\glut.h>

CollisionRod::CollisionRod()
{
	m_reader = std::make_shared<AbcReaderTransform>();
}


CollisionRod::~CollisionRod()
{
}

void
CollisionRod::glRender(int systemFrame, float timeStep, int numSpheres, float sphereRadius)
{
	timeStep *= 8.0f;

	//get start & end point
	Eigen::Vector3f top;
	Eigen::Vector3f bottom;

	int frame1 = std::floorf((float)systemFrame * timeStep);
	int frame2 = std::ceilf((float)systemFrame * timeStep);

	float weightFrame1 = ((float)systemFrame * timeStep) - (float)frame1;
	float weightFrame2 = 1.0f - weightFrame1;

	Eigen::Vector3f frame1Top;
	Eigen::Vector3f frame2Top;
	m_reader->sampleSpecific(0, frame1);
	frame1Top = m_reader->getTranslation(0);
	m_reader->sampleSpecific(0, frame2);
	frame2Top = m_reader->getTranslation(0);

	Eigen::Vector3f frame1Bottom;
	Eigen::Vector3f frame2Bottom;
	m_reader->sampleSpecific(1, frame1);
	frame1Bottom = m_reader->getTranslation(1);
	m_reader->sampleSpecific(1, frame2);
	frame2Bottom = m_reader->getTranslation(1);

	top = frame1Top * weightFrame2 + frame2Top * weightFrame1;
	bottom = frame1Bottom * weightFrame2 + frame2Bottom * weightFrame1;

	//now enforce collision spheres
	float stepSize = 1.0f / (float)numSpheres;
	for (int s = 0; s < numSpheres; ++s)
	{
		Eigen::Vector3f sphereCentre = (s * stepSize) * top + (1.0f - (s * stepSize)) * bottom;
		
		glPushMatrix();
		glTranslated(sphereCentre[0], sphereCentre[1], sphereCentre[2]);
		glColor3f(1.0, 0.0, 0.0);
		glutSolidSphere(sphereRadius, 50, 50);
		glColor3f(1.0, 1.0, 1.0);
		glPopMatrix();
	}
}

void
CollisionRod::readFromAbc(const std::string& fileName, const std::vector<std::string>& topBottomTransformNames)
{
	m_reader->openArchive(fileName, topBottomTransformNames);
}

void
CollisionRod::computeDeltaXPositionConstraint(float w1, float w2, float restDistance,
const Eigen::Vector3f& x1, const Eigen::Vector3f& x2, Eigen::Vector3f& temp, Eigen::Vector3f& deltaX)
{
	temp = x1 - x2;

	float squaredNorm = temp.squaredNorm();

	if (w1 + w2 == 0.0f || squaredNorm == 0.0f)
	{
		deltaX.setZero();
	}
	else
	{
		deltaX = (squaredNorm - restDistance) * (temp.normalized()) / (w1 + w2);
	}
}


void
CollisionRod::resolveParticleCollisions(std::vector<PBDParticle>& particles, int systemFrame, float timeStep,
int numSpheres, float sphereRadius)
{
	timeStep *= 8.0f;

	//get start & end point
	Eigen::Vector3f top;
	Eigen::Vector3f bottom;

	int frame1 = std::floorf((float)systemFrame * timeStep);
	int frame2 = std::ceilf((float)systemFrame * timeStep);

	float weightFrame1 = ((float)systemFrame * timeStep) - (float)frame1;
	float weightFrame2 = 1.0f - weightFrame1;

	Eigen::Vector3f frame1Top;
	Eigen::Vector3f frame2Top;
	m_reader->sampleSpecific(0, frame1);
	frame1Top = m_reader->getTranslation(0);
	m_reader->sampleSpecific(0, frame2);
	frame2Top = m_reader->getTranslation(0);

	Eigen::Vector3f frame1Bottom;
	Eigen::Vector3f frame2Bottom;
	m_reader->sampleSpecific(1, frame1);
	frame1Bottom = m_reader->getTranslation(0);
	m_reader->sampleSpecific(1, frame2);
	frame2Bottom = m_reader->getTranslation(1);

	top = frame1Top * weightFrame2 + frame2Top * weightFrame1;
	bottom = frame1Bottom * weightFrame2 + frame2Bottom * weightFrame1;

	//now enforce collision spheres
	float stepSize = 1.0f / (float)numSpheres;
	for (int s = 0; s < numSpheres; ++s)
	{
		Eigen::Vector3f sphereCentre = (s * stepSize) * top + (1.0f - (s * stepSize)) * bottom;
		Eigen::Vector3f temp;
		for (int p = 0; p < particles.size(); ++p)
		{
			if ((particles[p].position() - sphereCentre).squaredNorm() < sphereRadius)
			{
				//enforce distant constraint
				float w2 = particles[p].inverseMass();
				Eigen::Vector3f x2 = particles[p].position();
				Eigen::Vector3f deltaX;
				//computeDeltaXPositionConstraint(0.0f, w2, sphereRadius + 1e-15, sphereCentre, x2, temp, deltaX);
				computeDeltaXPositionConstraint(w2, 0.0f, sphereRadius, x2, sphereCentre, temp, deltaX);
				particles[p].position() += deltaX;
			}
		}
	}
}
