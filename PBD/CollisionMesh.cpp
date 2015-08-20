#include "CollisionMesh.h"


CollisionMesh::CollisionMesh()
{
	m_translation.setZero();
}


CollisionMesh::~CollisionMesh()
{
}


void
CollisionMesh::readFromAbc(const std::string& fileName)
{
	m_reader = std::make_shared<AbcReader>();
	m_reader->openArchive(fileName, "collisionGeometry");
}

//adapted from http://www.blackpawn.com/texts/pointinpoly/
bool SameSide(Eigen::Vector3f& p1, Eigen::Vector3f& p2, Eigen::Vector3f& a, Eigen::Vector3f& b)
{
	Eigen::Vector3f cp1 = (b - a).cross(p1 - a);
	Eigen::Vector3f cp2 = (b - a).cross(p2 - a);

	return cp1.dot(cp2) >= 0.0f;
}

//adapted from http://www.blackpawn.com/texts/pointinpoly/
bool PointInTriangle(Eigen::Vector3f& p, Eigen::Vector3f& a, Eigen::Vector3f& b, Eigen::Vector3f& c)
{
	return (SameSide(p, a, b, c) && SameSide(p, b, a, c) && SameSide(p, c, a, b));
}
void
CollisionMesh::resolveParticleCollisions(std::vector<PBDParticle>& particles)
{
	//for all particles
	for (int p = 0; p < particles.size(); ++p)
	{
		//test whether they intersect the collision shape
		int numTriangles = m_reader->getNumFaces();

		for (int t = 0; t < numTriangles; ++t)
		{
			Eigen::Vector3f p0 = m_reader->getPositions()[m_reader->getFaceIndices(t)[0]] + m_translation;
			Eigen::Vector3f p1 = m_reader->getPositions()[m_reader->getFaceIndices(t)[1]] + m_translation;
			Eigen::Vector3f p2 = m_reader->getPositions()[m_reader->getFaceIndices(t)[2]] + m_translation;


			Eigen::Vector3f particle = particles[p].position();

			p0[1] = 0.0f;
			p1[1] = 0.0f;
			p2[1] = 0.0f;
			particle[1] = 0.0f;

			if (PointInTriangle(particle, p0, p1, p2) && particles[p].position()[1] > (m_reader->getPositions()[m_reader->getFaceIndices(t)[0]] + m_translation)[1])
			{
				particles[p].position()[1] = (m_reader->getPositions()[m_reader->getFaceIndices(t)[0]] + m_translation)[1] - 1e-12f;
			}
		}
	}
}