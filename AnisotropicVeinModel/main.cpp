#include <string>
#include <vector>

#include <Eigen\Dense>

#include "MeshIO.h"
#include "OMDefinitions.h"
#include "common.h"

void computePrincipalSurfaceCurvatureDirections(meshType& inputMesh, meshType& outputMesh);

int main(int argc, char* argv[])
{
	std::string inputMeshFile = "C:/Users/Stephan/Desktop/veinMesh_smoothed.ply";
	std::string outputMeshFile = "C:/Users/Stephan/Desktop/veinMesh_curv.ply";
	
	meshType inputMesh;
	loadMesh(inputMeshFile, inputMesh);
	meshType outputMesh = inputMesh;

	computePrincipalSurfaceCurvatureDirections(inputMesh, outputMesh);

	saveMesh(outputMeshFile, outputMesh);
}

//Code based on the implementation by Alliez http://www-sop.inria.fr/members/Pierre.Alliez/demos/curvature/

float arc_sinus(float sinus)
{
	if (sinus >= 1.0f)
	{
		return 1.57079632679489661f;
	}
	if (sinus <= -1.0f)
	{
		return -1.57079632679489661f;
	}
	return std::asin(sinus);
}

void computePrincipalSurfaceCurvatureDirections(meshType& inputMesh, meshType& outputMesh)
{
	outputMesh.request_vertex_normals();

	if (!outputMesh.has_vertex_normals())
	{
		std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
	}

	for (int i = 0; i < inputMesh.n_vertices(); ++i)
	{
		OpenMesh::VertexHandle currentVertex(i);

		Eigen::Matrix3f curvatureTensor;
		curvatureTensor.setZero();

		float minLength = 100000000.0f;
		for (oneRingItOHalfEdge_t heIt(inputMesh.voh_begin(currentVertex)); heIt.is_valid(); ++heIt)
		{
			OpenMesh::FaceHandle currentFace(inputMesh.face_handle(*heIt));
			OpenMesh::FaceHandle adjacentFace(inputMesh.opposite_face_handle(*heIt));

			Eigen::Vector3f p1 = OMVec3_2_Eigen(OpenMesh::Vec3f(inputMesh.point(inputMesh.to_vertex_handle(*heIt))));
			Eigen::Vector3f p2 = OMVec3_2_Eigen(OpenMesh::Vec3f(inputMesh.point(inputMesh.to_vertex_handle(inputMesh.opposite_halfedge_handle(*heIt)))));

			Eigen::Vector3f edge = p1 - p2;
			Eigen::Vector3f edge_norm = edge.normalized();

			float edge_length = edge.squaredNorm();

			if (minLength > edge_length && minLength > 0.0f)
			{
				minLength = edge_length;
			}
		}

		for (oneRingItOHalfEdge_t heIt(inputMesh.voh_begin(currentVertex)); heIt.is_valid(); ++heIt)
		{
			OpenMesh::FaceHandle currentFace(inputMesh.face_handle(*heIt));
			OpenMesh::FaceHandle adjacentFace(inputMesh.opposite_face_handle(*heIt));

			Eigen::Vector3f p1 = OMVec3_2_Eigen(OpenMesh::Vec3f(inputMesh.point(inputMesh.to_vertex_handle(*heIt))));
			Eigen::Vector3f p2 = OMVec3_2_Eigen(OpenMesh::Vec3f(inputMesh.point(inputMesh.to_vertex_handle(inputMesh.opposite_halfedge_handle(*heIt)))));

			Eigen::Vector3f edge = p1 - p2;
			Eigen::Vector3f edge_norm = edge.normalized();

			float edge_length = edge.squaredNorm();

			Eigen::Vector3f normal1 = OMVec3_2_Eigen(inputMesh.calc_face_normal(currentFace));
			Eigen::Vector3f normal2 = OMVec3_2_Eigen(inputMesh.calc_face_normal(adjacentFace));

			float sinus = normal1.cross(normal2).dot(edge_norm);

			float beta = arc_sinus(sinus);

			curvatureTensor += (beta * minLength) * edge_norm * edge_norm.transpose();
		}
		//std::cout << curvatureTensor << std::endl;
		//Diagonalise
		Eigen::EigenSolver<Eigen::Matrix3f> eigenSolver(curvatureTensor);
		Eigen::Matrix3f S = eigenSolver.pseudoEigenvalueMatrix(); //squared eigenvalues of F
		Eigen::Matrix3f V = eigenSolver.pseudoEigenvectors(); //eigenvectors

		Eigen::Matrix3f S_sorted;
		Eigen::Matrix3f V_sorted;
		S_sorted.setZero();
		V_sorted.setZero();

		//sort
		if (std::abs(S(0, 0)) > std::abs(S(1, 1)) && std::abs(S(0, 0)) > std::abs(S(2, 2)))
		{
			S_sorted(2, 2) = S(0, 0);
			V_sorted.col(2) = V.col(0);
			if (std::abs(S(1, 1)) > std::abs(S(2, 2)))
			{
				S_sorted(1, 1) = S(1, 1);
				V_sorted.col(1) = V.col(1);

				S_sorted(0, 0) = S(2, 2);
				V_sorted.col(0) = V.col(2);
			}
			else
			{
				S_sorted(1, 1) = S(2, 2);
				V_sorted.col(1) = V.col(2);

				S_sorted(0, 0) = S(1, 1);
				V_sorted.col(0) = V.col(1);
			}
		}
		else if (std::abs(S(1, 1)) > std::abs(S(0, 0)) && std::abs(S(1, 1)) > std::abs(S(2, 2)))
		{
			S_sorted(2, 2) = S(1, 1);
			V_sorted.col(2) = V.col(1);
			if (std::abs(S(0, 0)) > std::abs(S(2, 2)))
			{
				S_sorted(1, 1) = S(0, 0);
				V_sorted.col(1) = V.col(0);

				S_sorted(0, 0) = S(2, 2);
				V_sorted.col(0) = V.col(2);
			}
			else
			{
				S_sorted(1, 1) = S(2, 2);
				V_sorted.col(1) = V.col(1);

				S_sorted(0, 0) = S(0, 0);
				V_sorted.col(0) = V.col(0);
			}
		}
		else
		{
			S_sorted(2, 2) = S(2, 2);
			V_sorted.col(2) = V.col(2);
			if (S(0, 0) > S(1, 1))
			{
				S_sorted(1, 1) = S(0, 0);
				V_sorted.col(1) = V.col(0);

				S_sorted(0, 0) = S(1, 1);
				V_sorted.col(0) = V.col(1);
			}
			else
			{
				S_sorted(1, 1) = S(1, 1);
				V_sorted.col(1) = V.col(1);

				S_sorted(0, 0) = S(0, 0);
				V_sorted.col(0) = V.col(0);
			}
		}

		Eigen::Vector3f direction = V_sorted.col(2);

		if (direction.squaredNorm() > 0.0f)
		{
			direction.normalize();
		}

		//std::cout << direction << std::endl;

		outputMesh.set_normal(currentVertex, OpenMesh::Vec3f(direction[0], direction[1], direction[2]));
	}
}