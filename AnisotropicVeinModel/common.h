/*
Please note that the first three functions are the same as in the other courseworks.
*/
#pragma once

#include "OMDefinitions.h"

#include <OpenMesh\Core\Geometry\VectorT.hh>
#include <Eigen\Dense>

Eigen::Vector3d OMVec3_2_Eigen(const OpenMesh::Vec3d& x);

Eigen::Vector3d OMVec3f_2_Eigend(const OpenMesh::Vec3f& x);

Eigen::Vector3f OMVec3_2_Eigen(const OpenMesh::Vec3f& x);

//! Computes differences in Laplacians based on given rotations
double computeDeformationSurfaceEnergy(meshType& oldMesh, meshType& newMesh,
	const std::vector<Eigen::Matrix3f>& rotations);


//! Computes differences in edge length
double computeEdgeLengthErrors(meshType& oldMesh, meshType& newMesh);

//! Computes differences in face area
double computeFaceAreaErrors(meshType& oldMesh, meshType& newMesh);