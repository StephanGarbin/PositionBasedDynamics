/*
Please note that this is the same file as in all the other courseworks.
*/
#pragma once

#include <OpenMesh\Core\Mesh\Traits.hh>
#include <OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>

struct myMeshTraits : public OpenMesh::DefaultTraits
{
	VertexAttributes(OpenMesh::Attributes::Normal |
	OpenMesh::Attributes::Color);
};

typedef OpenMesh::TriMesh_ArrayKernelT<myMeshTraits> meshType;

typedef meshType::ConstVertexIter vertIt_t;
typedef meshType::VertexVertexIter oneRingIt_t;
typedef meshType::VertexHandle vert_h;

typedef meshType::VertexOHalfedgeIter oneRingItOHalfEdge_t;
typedef meshType::VertexIHalfedgeIter oneRingItIHalfEdge_t;


typedef meshType::VertexEdgeIter oneRingItEdge_t;
typedef meshType::FaceVertexIter faceVertIt_t;

