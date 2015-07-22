#include "PBDTetrahedra3d.h"

#include <iostream>

#define pbdX1 (*m_particles)[m_vertexIndices[0]].previousPosition()
#define pbdX2 (*m_particles)[m_vertexIndices[1]].previousPosition()
#define pbdX3 (*m_particles)[m_vertexIndices[2]].previousPosition()
#define pbdX4 (*m_particles)[m_vertexIndices[3]].previousPosition()

#define pbdx1 (*m_particles)[m_vertexIndices[0]].position() 
#define pbdx2 (*m_particles)[m_vertexIndices[1]].position() 
#define pbdx3 (*m_particles)[m_vertexIndices[2]].position() 
#define pbdx4 (*m_particles)[m_vertexIndices[3]].position() 

#define pbdX1_raw (*m_particles)[m_vertexIndices[0]]
#define pbdX2_raw (*m_particles)[m_vertexIndices[1]]
#define pbdX3_raw (*m_particles)[m_vertexIndices[2]]
#define pbdX4_raw (*m_particles)[m_vertexIndices[3]]

#define pbdx1_raw (*m_particles)[m_vertexIndices[0]] 
#define pbdx2_raw (*m_particles)[m_vertexIndices[1]] 
#define pbdx3_raw (*m_particles)[m_vertexIndices[2]] 
#define pbdx4_raw (*m_particles)[m_vertexIndices[3]] 

#include <GL\glew.h>
#include <gl\GL.h>

PBDTetrahedra3d::PBDTetrahedra3d(std::vector<int>&& vertexIndices, const std::shared_ptr<std::vector<PBDParticle>>& particles)
{
	m_vertexIndices = std::move(vertexIndices);
	m_particles = particles;
	calculateReferenceShapeMatrix();
	calculateReferenceShapeMatrixInverseTranspose();
	calculateUndeformedVolume();
	calculateUndeformedSideLengths();
}

PBDTetrahedra3d::PBDTetrahedra3d(std::vector<int>& vertexIndices, const std::shared_ptr<std::vector<PBDParticle>>& particles)
{
	m_vertexIndices = vertexIndices;
	m_particles = particles;
	calculateReferenceShapeMatrix();
	calculateReferenceShapeMatrixInverseTranspose();
	calculateUndeformedVolume();
	calculateUndeformedSideLengths();
}


PBDTetrahedra3d::~PBDTetrahedra3d()
{
}


const Eigen::Matrix3f&
PBDTetrahedra3d::getDeformedShapeMatrix()
{
	calculateDeformedShapeMatrix();
	return m_deformedShapeMatrix;
}


void
PBDTetrahedra3d::calculateReferenceShapeMatrix()
{
	m_referenceShapeMatrix.col(0) = pbdX1 - pbdX4;
	m_referenceShapeMatrix.col(1) = pbdX2 - pbdX4;
	m_referenceShapeMatrix.col(2) = pbdX3 - pbdX4;
}

void PBDTetrahedra3d::calculateReferenceShapeMatrixInverseTranspose()
{
	m_referenceShapeMatrixInverse = m_referenceShapeMatrix.inverse();
	m_referenceShapeMatrixInverseTranspose = m_referenceShapeMatrix.inverse().transpose();
}

void
PBDTetrahedra3d::calculateDeformedShapeMatrix()
{
	m_deformedShapeMatrix.col(0) = pbdx1 - pbdx4;
	m_deformedShapeMatrix.col(1) = pbdx2 - pbdx4;
	m_deformedShapeMatrix.col(2) = pbdx3 - pbdx4;
}


Eigen::Matrix3f
PBDTetrahedra3d::getDeformationGradient()
{
	calculateDeformedShapeMatrix();
	//std::cout << "Reference Shape Matrix: " << std::endl;
	//std::cout << m_referenceShapeMatrix << std::endl;

	//std::cout << "Deformed Shape Matrix: " << std::endl;
	//std::cout << m_deformedShapeMatrix << std::endl;


	return m_deformedShapeMatrix * m_referenceShapeMatrixInverse;
}

float
PBDTetrahedra3d::getVolume()
{
	//http://mathworld.wolfram.com/Tetrahedron.html
	Eigen::Matrix4d temp;
	temp(0, 0) = pbdx1.x();
	temp(0, 1) = pbdx1.y();
	temp(0, 2) = pbdx1.z();
	temp(0, 3) = 1.0;

	temp(1, 0) = pbdx2.x();
	temp(1, 1) = pbdx2.y();
	temp(1, 2) = pbdx2.z();
	temp(1, 3) = 1.0;

	temp(2, 0) = pbdx3.x();
	temp(2, 1) = pbdx3.y();
	temp(2, 2) = pbdx3.z();
	temp(2, 3) = 1.0;

	temp(3, 0) = pbdx4.x();
	temp(3, 1) = pbdx4.y();
	temp(3, 2) = pbdx4.z();
	temp(3, 3) = 1.0;

	return (1.0 / 6.0) * std::abs(temp.determinant());
}

//float
//PBDTetrahedra3d::getUndeformedVolume()
//{
//	return m_undeformedVolume;
//}

void
PBDTetrahedra3d::calculateUndeformedVolume()
{
	//http://mathworld.wolfram.com/Tetrahedron.html
	Eigen::Matrix4f temp;
	temp(0, 0) = pbdX1.x();
	temp(0, 1) = pbdX1.y();
	temp(0, 2) = pbdX1.z();
	temp(0, 3) = 1.0;

	temp(1, 0) = pbdX2.x();
	temp(1, 1) = pbdX2.y();
	temp(1, 2) = pbdX2.z();
	temp(1, 3) = 1.0;

	temp(2, 0) = pbdX3.x();
	temp(2, 1) = pbdX3.y();
	temp(2, 2) = pbdX3.z();
	temp(2, 3) = 1.0;

	temp(3, 0) = pbdX4.x();
	temp(3, 1) = pbdX4.y();
	temp(3, 2) = pbdX4.z();
	temp(3, 3) = 1.0;

	m_undeformedVolume = (1.0 / 6.0) * std::abs(temp.determinant());
}

PBDParticle&
PBDTetrahedra3d::get_x(int index)
{
	switch (index)
	{
	case 0:
		return pbdx1_raw;
		break;
	case 1:
		return pbdx2_raw;
		break;
	case 2:
		return pbdx3_raw;
		break;
	case 3:
		return pbdx4_raw;
		break;
	default:
		return PBDParticle();
		break;
	}
}


PBDParticle&
PBDTetrahedra3d::get_X(int index)
{
	switch (index)
	{
	case 0:
		return pbdX1_raw;
		break;
	case 1:
		return pbdX2_raw;
		break;
	case 2:
		return pbdX3_raw;
		break;
	case 3:
		return pbdX4_raw;
		break;
	default:
		return PBDParticle();
		break;
	}
}

void
PBDTetrahedra3d::glRender(double r, double g, double b)
{
	glBegin(GL_TRIANGLES);
		glColor3d(r, g, b);
		glVertex3d(pbdx1.x(), pbdx1.y(), pbdx1.z());
		glVertex3d(pbdx2.x(), pbdx2.y(), pbdx2.z());
		glVertex3d(pbdx3.x(), pbdx3.y(), pbdx3.z());

		glVertex3d(pbdx1.x(), pbdx1.y(), pbdx1.z());
		glVertex3d(pbdx3.x(), pbdx3.y(), pbdx3.z());
		glVertex3d(pbdx4.x(), pbdx4.y(), pbdx4.z());

		glVertex3d(pbdx1.x(), pbdx1.y(), pbdx1.z());
		glVertex3d(pbdx3.x(), pbdx3.y(), pbdx3.z());
		glVertex3d(pbdx2.x(), pbdx2.y(), pbdx2.z());

		glVertex3d(pbdx2.x(), pbdx2.y(), pbdx2.z());
		glVertex3d(pbdx3.x(), pbdx3.y(), pbdx3.z());
		glVertex3d(pbdx4.x(), pbdx4.y(), pbdx4.z());
	glEnd();
}

void
PBDTetrahedra3d::calculateUndeformedSideLengths()
{
	m_undeformedSideLengths.push_back((pbdX1 - pbdX3).squaredNorm());
	m_undeformedSideLengths.push_back((pbdX1 - pbdX4).squaredNorm());
	m_undeformedSideLengths.push_back((pbdX1 - pbdX2).squaredNorm());
	m_undeformedSideLengths.push_back((pbdX3 - pbdX4).squaredNorm());
	m_undeformedSideLengths.push_back((pbdX3 - pbdX2).squaredNorm());
	m_undeformedSideLengths.push_back((pbdX4 - pbdX2).squaredNorm());
}

float
PBDTetrahedra3d::getUndeformedSideLength(int idx)
{
	return m_undeformedSideLengths[idx];
}