#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <memory>

#include <boost\algorithm\string.hpp>
#include <boost\algorithm\string\trim_all.hpp>

#include <Eigen\Dense>

#include <ANN\ANN.h>

#include "MeshIO.h"
#include "OMDefinitions.h"
#include "common.h"

#include "TetGenIO.h"
#include "PBDParticle.h"
#include "PBDTetrahedra3d.h"

#include "CustomTetAttributeIO.h"


void computePrincipalSurfaceCurvatureDirections(meshType& inputMesh, meshType& outputMesh);

void computeStiffnessAndAnisotropyFromParticles(std::vector<Eigen::Vector3f>& particles, std::vector<Eigen::Vector3f>& directionsPerParticle,
	std::vector<PBDParticle>& nodes, std::vector<PBDTetrahedra3d>& tets,
	std::vector<float>& youngsModulus,
	std::vector<float>& anisotropyStrength, std::vector<Eigen::Vector3f>& anisotropyDirections);

bool isInTet(PBDTetrahedra3d& tet, Eigen::Vector3f p);

void readParticlesFromRib(const std::string& file, std::vector<Eigen::Vector3f>& particles);

void assignDirection2Particles(const std::vector<Eigen::Vector3f>& particles, meshType& mesh,
	std::vector<Eigen::Vector3f>& particleDirections);

int main(int argc, char* argv[])
{
	if (argc < 7)
	{
		std::cout << "Please provide: input mesh file, output mesh file, particle RIB file, element file, node file, target attribute file." << std::endl;

		return 1;
	}
	std::string inputMeshFile(argv[1]);
	std::string outputMeshFile(argv[2]);

	std::string particleRIBFile(argv[3]);
	std::string elementFile(argv[4]);
	std::string nodeFile(argv[5]);
	std::string attributeFileName(argv[6]);

	std::cout << "Reading meshes..." << std::endl;
	meshType inputMesh;
	loadMesh(inputMeshFile, inputMesh);
	meshType outputMesh = inputMesh;

	computePrincipalSurfaceCurvatureDirections(inputMesh, outputMesh);
	saveMesh(outputMeshFile, outputMesh);
	std::cout << "Finished Computing principal curvature directions..." << std::endl;

	std::vector<Eigen::Vector3f> particles;
	readParticlesFromRib(particleRIBFile, particles);
	std::cout << "Finished Reading particles..." << std::endl;

	std::vector<Eigen::Vector3f> particleDirections;
	assignDirection2Particles(particles, outputMesh, particleDirections);
	std::cout << "Finished setting particle directions..." << std::endl;

	std::vector<float> youngsModulus;
	std::vector<float> anisotropyStrength;
	std::vector<Eigen::Vector3f> anisotropyDirections;
	std::vector<PBDTetrahedra3d> tets;
	std::shared_ptr<std::vector<PBDParticle>> nodes = std::make_shared<std::vector<PBDParticle>>();
	Eigen::Vector3f initialVelocity; initialVelocity.setZero();
	TetGenIO::readNodes(nodeFile, *nodes, 1.0f, initialVelocity);
	TetGenIO::readTetrahedra(elementFile, tets, nodes);
	std::cout << "Read Tetrahedra from disk..." << std::endl;

	computeStiffnessAndAnisotropyFromParticles(particles, particleDirections, *nodes, tets, youngsModulus, anisotropyStrength, anisotropyDirections);
	std::cout << "Finished computing per-tet values..." << std::endl;

	writeCustomTetAttributes(youngsModulus, anisotropyStrength, anisotropyDirections, attributeFileName);
	std::cout << "All done." << std::endl;
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

void computeStiffnessAndAnisotropyFromParticles(std::vector<Eigen::Vector3f>& particles, std::vector<Eigen::Vector3f>& directionsPerParticle,
	std::vector<PBDParticle>& nodes, std::vector<PBDTetrahedra3d>& tets,
	std::vector<float>& youngsModulus,
	std::vector<float>& anisotropyStrength, std::vector<Eigen::Vector3f>& anisotropyDirections)
{
	std::vector<float> numParticlesPerTet(tets.size());

	std::vector<std::vector<Eigen::Vector3f>> directionsPerTet(tets.size());

	for (int p = 0; p < particles.size(); ++p)
	{
		int tetIdx;
		for (int t = 0; t < tets.size(); ++t)
		{
			if (isInTet(tets[t], particles[p]))
			{
				numParticlesPerTet[t] += 1.0f;
				directionsPerTet[t].push_back(directionsPerParticle[p]);
			}
		}
	}

	float maxParticles = 0;
	for (int t = 0; t < tets.size(); ++t)
	{
		maxParticles = (numParticlesPerTet[t] > maxParticles) ? numParticlesPerTet[t] : maxParticles;
	}

	//assign mean directions
	anisotropyDirections.resize(directionsPerTet.size());
	for (int i = 0; i < anisotropyDirections.size(); ++i)
	{
		Eigen::Vector3f sum; sum.setZero();
		for (int d = 0; d < directionsPerTet[i].size(); ++d)
		{
			sum += directionsPerTet[i][d];
		}

		if (sum.squaredNorm() != 0.0f)
		{
			anisotropyDirections[i] = sum / (float)directionsPerTet[i].size();
		}
		else
		{
			anisotropyDirections[i].setZero();
		}
	}

	youngsModulus.resize(anisotropyDirections.size());
	anisotropyStrength.resize(anisotropyDirections.size());
	for (int i = 0; i < anisotropyDirections.size(); ++i)
	{
		if (maxParticles > 0.0f)
		{
			youngsModulus[i] = numParticlesPerTet[i] / maxParticles;
			anisotropyStrength[i] = numParticlesPerTet[i] / maxParticles;
		}
		else
		{
			youngsModulus[i] = 0.0f;
			anisotropyStrength[i] = 0.0f;
		}
	}
}

void readParticlesFromRib(const std::string& fileName, std::vector<Eigen::Vector3f>& particles)
{
	std::ifstream file;
	file.open(fileName);
	
	bool readingPoints = false;
	std::string line;
	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			boost::algorithm::trim_all(line);

			if (readingPoints)
			{
				if (line.find("]", 0) != std::string::npos)
				{
					readingPoints = false;
					break;
				}
				else
				{
					std::vector<std::string> floatVector;
					boost::split(floatVector, line, boost::is_any_of(" "));

					for (int i = 0; i < floatVector.size() / 3; ++i)
					{
						particles.push_back(Eigen::Vector3f(std::stof(floatVector[i * 3 + 0]),
							std::stof(floatVector[i * 3 + 1]),
							std::stof(floatVector[i * 3 + 2])));
					}
				}
			}
			else
			{
				if (line[0] == 'P' && line[1] == 'o' && line[2] == 'i' && line[3] == 'n' && line[4] == 't' && line[5] == 's')
				{
					readingPoints = true;
					size_t start = line.find("[", 0);
					line.erase(0, start + 1);
					std::vector<std::string> floatVector;
					boost::split(floatVector, line, boost::is_any_of(" "));

					for (int i = 0; i < floatVector.size() / 3; ++i)
					{
						particles.push_back(Eigen::Vector3f(std::stof(floatVector[i * 3 + 0]),
							std::stof(floatVector[i * 3 + 1]),
							std::stof(floatVector[i * 3 + 2])));
					}
				}

			}
		}
		file.close();
	}

	std::cout << "Read " << particles.size() << " particles from RIB." << std::endl;
}

bool isInTet(PBDTetrahedra3d& tet, Eigen::Vector3f p)
{
	Eigen::Vector4f point;
	point[0] = p[0];
	point[1] = p[1];
	point[2] = p[2];
	point[3] = 1.0f;

	Eigen::Vector4f x1;
	x1[0] = tet.get_x(0).position()[0];
	x1[1] = tet.get_x(0).position()[1];
	x1[2] = tet.get_x(0).position()[2];
	x1[3] = 1.0f;

	Eigen::Vector4f x2;
	x2[0] = tet.get_x(1).position()[0];
	x2[1] = tet.get_x(1).position()[1];
	x2[2] = tet.get_x(1).position()[2];
	x2[3] = 1.0f;

	Eigen::Vector4f x3;
	x3[0] = tet.get_x(2).position()[0];
	x3[1] = tet.get_x(2).position()[1];
	x3[2] = tet.get_x(2).position()[2];
	x3[3] = 1.0f;

	Eigen::Vector4f x4;
	x4[0] = tet.get_x(3).position()[0];
	x4[1] = tet.get_x(3).position()[1];
	x4[2] = tet.get_x(3).position()[2];
	x4[3] = 1.0f;

	bool isInTet = false;
	float lbd1;
	float lbd2;
	float lbd3;
	float lbd4;

	Eigen::Matrix4f detMatrix_0;
	detMatrix_0.row(0) = x1; detMatrix_0.row(1) = x2; detMatrix_0.row(2) = x3; detMatrix_0.row(3) = x4;

	Eigen::Matrix4f detMatrix_1;
	detMatrix_1.row(0) = point; detMatrix_1.row(1) = x2; detMatrix_1.row(2) = x3; detMatrix_1.row(3) = x4;

	Eigen::Matrix4f detMatrix_2;
	detMatrix_2.row(0) = x1; detMatrix_2.row(1) = point; detMatrix_2.row(2) = x3; detMatrix_2.row(3) = x4;

	Eigen::Matrix4f detMatrix_3;
	detMatrix_3.row(0) = x1; detMatrix_3.row(1) = x2; detMatrix_3.row(2) = point; detMatrix_3.row(3) = x4;

	Eigen::Matrix4f detMatrix_4;
	detMatrix_4.row(0) = x1; detMatrix_4.row(1) = x2; detMatrix_4.row(2) = x3; detMatrix_4.row(3) = point;

	if (detMatrix_0.determinant() == 0.0f)
	{
		std::cout << "Degenerate Tetrahedron detected! " << std::endl;
		return false;
	}

	lbd1 = detMatrix_1.determinant() / detMatrix_0.determinant();
	lbd2 = detMatrix_2.determinant() / detMatrix_0.determinant();
	lbd3 = detMatrix_3.determinant() / detMatrix_0.determinant();
	lbd4 = detMatrix_4.determinant() / detMatrix_0.determinant();

	if ((lbd1 < 0.0f && lbd2 < 0.0f && lbd3 < 0.0f && lbd4 < 0.0f)
		|| (lbd1 > 0.0f && lbd2 > 0.0f && lbd3 > 0.0f && lbd4 > 0.0f))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void assignDirection2Particles(const std::vector<Eigen::Vector3f>& particles, meshType& mesh,
	std::vector<Eigen::Vector3f>& particleDirections)
{
	std::cout << "Building kd-tree..." << std::endl;

	ANNpointArray p1 = annAllocPts(mesh.n_vertices(), 3);

	ANNkd_tree* kdTree = new ANNkd_tree(p1, mesh.n_vertices(), 3);

	std::vector<Eigen::Vector3f> meshNormals(mesh.n_vertices());
	size_t p = 0;
	for (vertIt_t it(mesh.vertices_begin()); it != mesh.vertices_end(); ++it)
	{
		p1[p][0] = mesh.point(*it)[0];
		p1[p][1] = mesh.point(*it)[1];
		p1[p][2] = mesh.point(*it)[2];

		meshNormals[p][0] = mesh.normal(*it)[0];
		meshNormals[p][1] = mesh.normal(*it)[1];
		meshNormals[p][2] = mesh.normal(*it)[2];
		++p;
	}

	int knn = 1;
	ANNpoint point = annAllocPt(3);
	ANNidxArray nnIndices = new ANNidx[knn];
	ANNdistArray nnDistances = new ANNdist[knn];

	particleDirections.resize(particles.size());
	for (int i = 0; i < particles.size(); ++i)
	{
		point[0] = particles[i][0];
		point[1] = particles[i][1];
		point[2] = particles[i][2];
		//std::cout << "3" << std::endl;
		//std::cout << point[0] << ", " << point[1] << ", " << point[2] << std::endl;
		kdTree->annkSearch(point, knn, nnIndices, nnDistances);
		//std::cout << "4" << std::endl;
		int ptIdx = nnIndices[0];
		particleDirections[i] = meshNormals[ptIdx];
	}
	std::cout << "5" << std::endl;
	annClose();

	annDeallocPt(point);
	delete[] nnIndices;
	delete[] nnDistances;

	annDeallocPts(p1);
	delete kdTree;

	std::cout << "Finished particle-direction assingment." << std::endl;
}