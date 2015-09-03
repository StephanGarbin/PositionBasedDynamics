#include "CollisionSphere.h"

#include <GL\glew.h>
#include <gl\GL.h>
#include <GL\glut.h>

#include "commonMath.h"

CollisionSphere::CollisionSphere()
{
	m_reader = std::make_shared<AbcReaderTransform>();
	m_frameLimit == -1;
	m_translation.setZero();
	m_previousCollisionSphereCentre = Eigen::Vector3f(-1.119f, -1.119f, -1.771f);
	m_lastProcessedFrame = -1;
}


CollisionSphere::~CollisionSphere()
{
}


void
CollisionSphere::readFromAbc(const std::string& fileName, const std::string& transformName)
{
	std::vector<std::string> temp;
	temp.push_back(transformName);
	m_reader->openArchive(fileName, temp);
}

void
CollisionSphere::resolveParticleCollisions(std::vector<PBDParticle>& particles, int systemFrame, float timeStep,
float sphereRadius)
{
	timeStep *= 4.0f;
	//get start & end point
	Eigen::Vector3f top;
	Eigen::Vector3f bottom;

	int frame1 = std::floorf((float)systemFrame * timeStep);
	int frame2 = std::ceilf((float)systemFrame * timeStep);

	if (m_frameLimit > 0)
	{
		frame1 = (frame1 > m_frameLimit) ? m_frameLimit : frame1;
		frame2 = (frame2 > m_frameLimit) ? m_frameLimit : frame2;
	}

	float weightFrame1 = ((float)systemFrame * timeStep) - (float)frame1;
	float weightFrame2 = 1.0f - weightFrame1;

	Eigen::Vector3f frame1Top;
	Eigen::Vector3f frame2Top;
	m_reader->sampleSpecific(0, frame1);
	frame1Top = m_reader->getTranslation(0);
	m_reader->sampleSpecific(0, frame2);
	frame2Top = m_reader->getTranslation(0);

	top = frame1Top * weightFrame2 + frame2Top * weightFrame1;

	//now enforce collision sphere
	Eigen::Vector3f sphereCentre = top;
	//std::cout << sphereCentre << std::endl;
	for (int p = 0; p < particles.size(); ++p)
	{;
		if (std::sqrtf((sphereCentre - particles[p].position()).squaredNorm()) < sphereRadius)
		{
			Eigen::Vector3f sphereMotion;
			Eigen::Vector3f particleMotion;
			if (m_previousCollisionSphereCentre[0] == -1.119f &&m_previousCollisionSphereCentre[1] == -1.119f && m_collisionSphereCentre[2] == -1.771f)
			{
				sphereMotion.setZero();
			}
			else
			{
				sphereMotion = sphereCentre - m_previousCollisionSphereCentre;
			}

			particleMotion = particles[p].position() - particles[p].previousPosition();

			float penetrationAmount = sphereRadius - std::sqrtf((sphereCentre - particles[p].position()).squaredNorm());

			//std::cout << penetrationAmount << std::endl;

			//if (sphereMotion.squaredNorm() == 0.0f) //interpenetration caused by particles
			//{
			//	particles[p].position() -= particleMotion.normalized() * penetrationAmount;
			//}
			//else //interpenetration caused by sphere
			//{
			//	std::vector<Eigen::Vector3f> intersectionPoints(2);
			//	int numIntersections;

			//	particles[p].position() += sphereMotion.normalized() * penetrationAmount;
			//}
			particles[p].position() -= (sphereCentre - particles[p].position()).normalized() * penetrationAmount;

			//std::cout << "eC";
			//enforce distant constraint
			//float w2 = particles[p].inverseMass();
			//Eigen::Vector3f x2 = particles[p].position();
			//Eigen::Vector3f deltaX;
			//Eigen::Vector3f temp;
			//computeDeltaXPositionConstraint(0.0f, w2, sphereRadius + 1e-15, sphereCentre, x2, temp, deltaX);
			//particles[p].position() += deltaX;

			m_previousCollisionSphereCentre = sphereCentre;
		}
	}
}

void
CollisionSphere::calculateNewSphereCentre(int systemFrame, float timeStep)
{
	//Don't process frames twice
	if (systemFrame == m_lastProcessedFrame)
	{
		return;
	}

	timeStep *= 1.0f;
	//get start & end point
	Eigen::Vector3f top;
	Eigen::Vector3f bottom;

	int frame1 = std::floorf((float)systemFrame * timeStep);
	int frame2 = std::ceilf((float)systemFrame * timeStep);

	if (m_frameLimit > 0)
	{
		frame1 = (frame1 > m_frameLimit) ? m_frameLimit : frame1;
		frame2 = (frame2 > m_frameLimit) ? m_frameLimit : frame2;
	}

	float weightFrame1 = ((float)systemFrame * timeStep) - (float)frame1;
	float weightFrame2 = 1.0f - weightFrame1;

	Eigen::Vector3f frame1Top;
	Eigen::Vector3f frame2Top;
	m_reader->sampleSpecific(0, frame1);
	frame1Top = m_reader->getTranslation(0);
	m_reader->sampleSpecific(0, frame2);
	frame2Top = m_reader->getTranslation(0);

	m_previousCollisionSphereCentre = m_collisionSphereCentre;
	m_collisionSphereCentre = frame1Top * weightFrame2 + frame2Top * weightFrame1;

	m_lastProcessedFrame = systemFrame;
}

void
CollisionSphere::resolveParticleCollisions_SAFE(std::vector<PBDParticle>& particles, int systemFrame, float timeStep,
	float sphereRadius, int start, int end)
{
	Eigen::Vector3f sphereCentre = m_collisionSphereCentre;
	for (int p = start; p != end; ++p)
	{
		if (std::sqrtf((sphereCentre - particles[p].position()).squaredNorm()) < sphereRadius)
		{
			Eigen::Vector3f sphereMotion;
			Eigen::Vector3f particleMotion;
			if (m_previousCollisionSphereCentre[0] == -1.119f &&m_previousCollisionSphereCentre[1] == -1.119f && m_collisionSphereCentre[2] == -1.771f)
			{
				sphereMotion.setZero();
			}
			else
			{
				sphereMotion = sphereCentre - m_previousCollisionSphereCentre;
			}

			particleMotion = particles[p].position() - particles[p].previousPosition();

			float penetrationAmount = sphereRadius - std::sqrtf((sphereCentre - particles[p].position()).squaredNorm());
			
			Eigen::Vector3f correction = (sphereCentre - particles[p].position()).normalized() * penetrationAmount;

			if (std::isnan(correction[0]) || std::isinf(correction[0])
				|| std::isnan(correction[1]) || std::isinf(correction[1])
				|| std::isnan(correction[2]) || std::isinf(correction[2]))
			{
				std::cout << "ERROR in spher collision: NaN result!" << std::endl;
			}

			particles[p].position() -= correction;
			//for (int n = 0; n < particles[p].getContainingTetIdxs().size(); ++n)
			//{
			//	
			//}
			continue;
		}
		//sphereRadius += sphereRadius / 16.0f;
		//if (std::sqrtf((sphereCentre - particles[p].position()).squaredNorm()) < sphereRadius)
		//{
		//	Eigen::Vector3f sphereMotion;
		//	Eigen::Vector3f particleMotion;
		//	if (m_previousCollisionSphereCentre[0] == -1.119f &&m_previousCollisionSphereCentre[1] == -1.119f && m_collisionSphereCentre[2] == -1.771f)
		//	{
		//		sphereMotion.setZero();
		//	}
		//	else
		//	{
		//		sphereMotion = sphereCentre - m_previousCollisionSphereCentre;
		//	}

		//	particleMotion = particles[p].position() - particles[p].previousPosition();

		//	float penetrationAmount = (sphereRadius - std::sqrtf((sphereCentre - particles[p].position()).squaredNorm())) / 2;

		//	particles[p].position() -= (sphereCentre - particles[p].position()).normalized() * penetrationAmount;
		//}
	}
}


void
CollisionSphere::checkForSinglePointIntersection_SAFE(const Eigen::Vector3f& point, float& penetrationDistance, bool& penetrates, float sphereRadius)
{
	if (std::sqrtf((m_collisionSphereCentre - point).squaredNorm()) < sphereRadius)
	{
		penetrates = true;
		penetrationDistance = sphereRadius - std::sqrtf((m_collisionSphereCentre - point).squaredNorm());
	}
	else
	{
		penetrates = false;
		penetrationDistance = 0.0f;
	}

}

void
CollisionSphere::checkForSinglePointIntersection_normalReflection_SAFE(const Eigen::Vector3f& previousPoint, const Eigen::Vector3f& point, float& penetrationDistance,
	bool& penetrates, Eigen::Vector3f& normal, Eigen::Vector3f& correction, float sphereRadius)
{
	if (std::sqrtf((m_collisionSphereCentre - point).squaredNorm()) < sphereRadius)
	{
		penetrates = true;
		penetrationDistance = sphereRadius - std::sqrtf((m_collisionSphereCentre - point).squaredNorm());
		Eigen::Vector3f intersectionPoint = (point - previousPoint) - (point - previousPoint).normalized() * penetrationDistance;

		normal = (intersectionPoint - m_collisionSphereCentre).normalized();

		correction = intersectionPoint + (intersectionPoint - 2.0f * (intersectionPoint.dot(normal)) * normal).normalized() * penetrationDistance;
	}
	else
	{
		penetrates = false;
		penetrationDistance = 0.0f;
		normal.setZero();
		correction.setZero();
	}

}

void
CollisionSphere::glRender(int systemFrame, float timeStep, float sphereRadius)
{
	timeStep *= 1.0f;
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

	top = frame1Top * weightFrame2 + frame2Top * weightFrame1;

	//now enforce collision sphere
	Eigen::Vector3f sphereCentre = top;

	glPushMatrix();
	glTranslated(sphereCentre[0], sphereCentre[1], sphereCentre[2]);
	glColor3f(1.0, 0.0, 0.0);
	glutSolidSphere(sphereRadius, 50, 50);
	glColor3f(1.0, 1.0, 1.0);
	glPopMatrix();
}


void
CollisionSphere::computeDeltaXPositionConstraint(float w1, float w2, float restDistance,
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
