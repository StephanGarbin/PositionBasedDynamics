#pragma once

#include <vector>
#include <memory>

#include <Eigen/Dense>

#include "PBDParticle.h"

class PBDTetrahedra3d
{
public:
	PBDTetrahedra3d(std::vector<int>&& vertexIndices, const std::shared_ptr<std::vector<PBDParticle>>& particles);
	PBDTetrahedra3d(std::vector<int>& vertexIndices, const std::shared_ptr<std::vector<PBDParticle>>& particles);
	~PBDTetrahedra3d();

	const std::vector<int>& getVertexIndices() const { return m_vertexIndices; }

	const Eigen::Matrix3f& getReferenceShapeMatrix() const { return m_referenceShapeMatrix; }
	const Eigen::Matrix3f& getReferenceShapeMatrixInverseTranspose() const { return m_referenceShapeMatrixInverseTranspose; }
	const Eigen::Matrix3f& getDeformedShapeMatrix();

	Eigen::Matrix3f getDeformationGradient();

	Eigen::Matrix3f getVelocityGradient();

	float getVolume();

	float getUndeformedVolume()
	{
		return m_undeformedVolume;
	}

	float getUndeformedVolumeAlternative()
	{
		return m_undeformedVolumeAlternative;
	}

	float getUndeformedSideLength(int idx);

	PBDParticle& get_x(int index);
	PBDParticle& get_X(int index);

	Eigen::Matrix3f& getUpsilon()
	{
		return m_upsilon;
	}

	Eigen::Matrix3f& getFullUpsilon(int idx)
	{
		return m_upsilonFull[idx];
	}

	void setFullUpsilonCount(int num)
	{
		m_upsilonFull.resize(num);
		for (int i = 0; i < num; ++i)
		{
			m_upsilonFull[i].setZero();
		}
	}

	Eigen::Vector3f& getFaceVertex(int face, int vertex);

	Eigen::Matrix3f& getDistortionDissipative() { return m_distortionDissipative; }
	Eigen::Matrix3f& getDistortionElastic() { return m_distortionElastic; }

	void setPreviousDeformedShapeMatrix()
	{
		calculateDeformedShapeMatrix();
		m_deformedShapeMatrix_previousPosition = m_deformedShapeMatrix;
	}

	Eigen::Matrix3f getRelativeDeformationGradient()
	{
		calculateDeformedShapeMatrix();
		return m_deformedShapeMatrix * m_deformedShapeMatrix_previousPosition.inverse();
	}

	void calculateReferenceVelocityMatrix();

	Eigen::Matrix3f calculateRelativeDeformationGradientVelocity();

	Eigen::Matrix3f calculateRelativeDeformationGradientPosition();

	void glRender(double r, double g, double b);
private:

	void initialise(std::vector<int>& vertexIndices, const std::shared_ptr<std::vector<PBDParticle>>& particles);

	void calculateUndeformedVolume();
	void calculateReferenceShapeMatrix();
	void calculateReferenceShapeMatrixInverseTranspose();
	void calculateDeformedShapeMatrix();

	void calculateUndeformedSideLengths();

	Eigen::Matrix3f m_referenceShapeMatrix;
	Eigen::Matrix3f m_referenceShapeMatrixInverseTranspose;
	Eigen::Matrix3f m_referenceShapeMatrixInverse;
	Eigen::Matrix3f m_deformedShapeMatrix;

	Eigen::Matrix3f m_deformedShapeMatrix_previousVelocity;
	Eigen::Matrix3f m_deformedShapeMatrix_previousPosition;

	Eigen::Matrix3f m_velocityMatrix;

	std::vector<int> m_vertexIndices;
	std::shared_ptr<std::vector<PBDParticle>> m_particles;
	float m_undeformedVolume;

	std::vector<float> m_undeformedSideLengths;

	float m_undeformedVolumeAlternative;


	//viscoelasticity
	Eigen::Matrix3f m_upsilon;
	std::vector<Eigen::Matrix3f> m_upsilonFull;

	//Rubin-Bodner stuff
	Eigen::Matrix3f m_distortionElastic;
	Eigen::Matrix3f m_distortionDissipative;

};

