#include "common.h"
#include "commonMath.h"

Eigen::Vector3d OMVec3_2_Eigen(const OpenMesh::Vec3d& x)
{
	return Eigen::Vector3d(x[0], x[1], x[2]);
}

Eigen::Vector3d OMVec3f_2_Eigend(const OpenMesh::Vec3f& x)
{
	return Eigen::Vector3d(x[0], x[1], x[2]);
}

Eigen::Vector3f OMVec3_2_Eigen(const OpenMesh::Vec3f& x)
{
	return Eigen::Vector3f(x[0], x[1], x[2]);
}


void computeLaplacianSums(std::vector<std::vector<Eigen::Vector3f>>& components, meshType& mesh)
{
	//Add entries
	for (int i = 0; i < mesh.n_vertices(); ++i)
	{
		OpenMesh::Vec3f centrePoint = mesh.point(vert_h(i));

		std::vector<Eigen::Vector3f> neighbourPoints;
		std::vector<int> neighbourPointsIndices;
		std::vector<float> coefficients;


		for (oneRingIt_t vvIt(mesh.vv_begin(vert_h(i))); vvIt != mesh.vv_end(vert_h(i)); ++vvIt)
		{
			neighbourPoints.push_back(OMVec3_2_Eigen(mesh.point(*vvIt)));
			neighbourPointsIndices.push_back(vvIt->idx());
		}

		//two boolen vectors (incoming, outgoing)
		std::vector<bool> incoming;
		std::vector<bool> outgoing;

		std::vector<bool> isboundary;

		for (oneRingIt_t vvIt(mesh.vv_begin(vert_h(i))); vvIt != mesh.vv_end(vert_h(i)); ++vvIt)
		{
			isboundary.push_back(mesh.is_boundary(*vvIt));
		}

		//std::cout << "Valence at " << vert_h(i).idx() << " : " << neighbourPointsIndices.size() << std::endl;

		float sumArea = 0.0f;
		float coeffSum = 0.0f;

		Eigen::Vector3f centralPoint(OMVec3_2_Eigen(centrePoint));
		for (int iu = 0; iu < neighbourPoints.size(); ++iu)
		{
			Eigen::Vector3f previousPoint;
			Eigen::Vector3f nextPoint;

			bool previousPointIsBoundary;
			bool nextPointIsBoundary;
			bool currentPointIsBoundary;

			if (iu == 0)
			{
				nextPoint = neighbourPoints[iu + 1];
				previousPoint = neighbourPoints[neighbourPoints.size() - 1];

				nextPointIsBoundary = isboundary[iu + 1];
				previousPointIsBoundary = isboundary[neighbourPoints.size() - 1];
			}
			else if (iu == neighbourPoints.size() - 1)
			{
				nextPoint = neighbourPoints[0];
				previousPoint = neighbourPoints[iu - 1];

				nextPointIsBoundary = isboundary[0];
				previousPointIsBoundary = isboundary[iu - 1];
			}
			else
			{
				nextPoint = neighbourPoints[iu + 1];
				previousPoint = neighbourPoints[iu - 1];

				nextPointIsBoundary = isboundary[iu + 1];
				previousPointIsBoundary = isboundary[iu - 1];
			}

			currentPointIsBoundary = isboundary[iu];

			Eigen::Vector3f a(centralPoint - nextPoint);
			Eigen::Vector3f b(neighbourPoints[iu] - nextPoint);

			float beta = std::atan2f((a.cross(b) / (a.norm() * b.norm())).norm(), a.dot(b) / (a.norm() * b.norm()));

			Eigen::Vector3f c(centralPoint - previousPoint);
			Eigen::Vector3f d(neighbourPoints[iu] - previousPoint);

			float alpha = std::atan2f((c.cross(d) / (c.norm() * d.norm())).norm(), c.dot(d) / (c.norm() * d.norm()));

			float coeff;

			if (currentPointIsBoundary)
			{
				if (nextPointIsBoundary && previousPointIsBoundary)
				{
					if (std::isinf(cotan(alpha)))
					{
						coeff = 1.0 * (cotan(beta));
					}
					else
					{
						coeff = 1.0 * (cotan(alpha));
					}
				}
				else if (nextPointIsBoundary)
				{
					coeff = 1.0 * (cotan(alpha));
				}
				else if (previousPointIsBoundary)
				{
					coeff = 1.0 * (cotan(beta));
				}
			}
			else
			{
				coeff = 0.5 * (cotan(alpha) + cotan(beta));
			}

			Eigen::Vector3f r = previousPoint - OMVec3_2_Eigen(centrePoint);
			Eigen::Vector3f s = neighbourPoints[iu] - OMVec3_2_Eigen(centrePoint);

			float theta = std::acosf((r.dot(s)) / (r.norm() * s.norm()));

			if (std::isnan(coeff))
			{
				std::cout << "NAN detected! " << std::endl;
			}

			if (std::isinf(coeff))
			{
				std::cout << "INF detected! " <<
					"Cotan alpha: " << cotan(alpha) <<
					"; Cotan beta: " << cotan(beta) <<
					"; Idx: " << neighbourPointsIndices[iu] << std::endl;
			}

			coeffSum += coeff;

			coefficients.push_back(coeff);


			if (neighbourPointsIndices.size() == 2 && iu == 0)
			{
				Eigen::Vector3f x = neighbourPoints[0] - centralPoint;
				Eigen::Vector3f y = neighbourPoints[1] - centralPoint;
				sumArea += (x.cross(y).norm() / 2.0f) / 3.0f;
				continue;
			}

			if (neighbourPointsIndices.size() == 2 && iu == 1)
			{
				continue; //we've already compute this for 0
			}

			//Area
			if (previousPointIsBoundary && currentPointIsBoundary && !nextPointIsBoundary)
			{
				//do nothing
			}
			else
			{
				//else calculate area regularly
				Eigen::Vector3f x = nextPoint - centralPoint;
				Eigen::Vector3f y = neighbourPoints[iu] - centralPoint;
				sumArea += (x.cross(y).norm() / 2.0f) / 3.0f;
			}
		}

		if (sumArea == 0)
		{
			std::cout << "ERROR: zero area detected" << std::endl;
		}

		std::vector<Eigen::Vector3f> currentComponent;

		for (int iy = 0; iy < neighbourPointsIndices.size(); ++iy)
		{
			currentComponent.push_back(coefficients[iy] * (centralPoint - neighbourPoints[iy]));
		}

		components.emplace_back(std::move(currentComponent));
	}

}

double computeDeformationSurfaceEnergy(meshType& oldMesh, meshType& newMesh,
	const std::vector<Eigen::Matrix3f>& rotations)
{
	std::vector<std::vector<Eigen::Vector3f>> oldLaplacianComponents;
	std::vector<std::vector<Eigen::Vector3f>> newLaplacianComponents;
	computeLaplacianSums(oldLaplacianComponents, oldMesh);
	computeLaplacianSums(newLaplacianComponents, newMesh);

	double totalEnergy = 0.0;

	for (int i = 0; i < oldLaplacianComponents.size(); ++i)
	{
		for (int n = 0; n < oldLaplacianComponents[i].size(); ++n)
		{
			totalEnergy += (oldLaplacianComponents[i][n] - rotations[i] * newLaplacianComponents[i][n]).norm();
		}
	}

	std::cout << "Total Energy: " << totalEnergy << std::endl;
	std::cout << "Energy Mean : " << totalEnergy / oldLaplacianComponents.size() << std::endl;

	return totalEnergy;
}


double computeEdgeLengthErrors(meshType& oldMesh, meshType& newMesh)
{
	double error = 0.0;

	for (int i = 0; i < oldMesh.n_edges(); ++i)
	{
		double oldEdgeLength = oldMesh.calc_edge_sqr_length(meshType::EdgeHandle(i));

		double newEdgeLength = newMesh.calc_edge_sqr_length(meshType::EdgeHandle(i));

		//error += std::pow(oldEdgeLength - newEdgeLength, 2.0);
		error += std::abs(oldEdgeLength - newEdgeLength);
	}

	std::cout << "Squared Edge Length Error: " << error << std::endl;
	std::cout << "Square Mean Edge Error   : " << error / oldMesh.n_edges() << std::endl;

	return error;
}

double computeFaceAreaErrors(meshType& oldMesh, meshType& newMesh)
{
	double error = 0.0;



	for (int i = 0; i < oldMesh.n_faces(); ++i)
	{
		std::vector<Eigen::Vector3f> verticesOld;
		for (meshType::FaceEdgeIter fEIt = oldMesh.fe_begin(meshType::FaceHandle(i)); fEIt != oldMesh.fe_end(meshType::FaceHandle(i)); ++fEIt)
		{
			verticesOld.push_back(OMVec3_2_Eigen(oldMesh.calc_edge_vector(*fEIt)));
		}

		std::vector<Eigen::Vector3f> verticesNew;
		for (meshType::FaceEdgeIter fEIt = newMesh.fe_begin(meshType::FaceHandle(i)); fEIt != newMesh.fe_end(meshType::FaceHandle(i)); ++fEIt)
		{
			verticesNew.push_back(OMVec3_2_Eigen(newMesh.calc_edge_vector(*fEIt)));
		}

		double oldArea = verticesOld[0].cross(verticesOld[1]).norm() / 2.0;

		double newArea = verticesNew[0].cross(verticesNew[1]).norm() / 2.0;

		//error += std::pow(oldEdgeLength - newEdgeLength, 2.0);
		error += std::abs(oldArea - newArea);
	}

	std::cout << "Squared Edge Area Error: " << error << std::endl;
	std::cout << "Square Mean Area Error   : " << error / oldMesh.n_edges() << std::endl;

	return error;
}