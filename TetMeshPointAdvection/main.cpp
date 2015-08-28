#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "AbcReader.h"

int main(int argc, char* argv[])
{
	//Open Alembic Archive
	std::string input1;
	input1 = std::string(argv[1]);
	std::cout << "Abc File: " << input1 << std::endl;

	//Read the archives
	AbcReader reader1;
	if (!reader1.openArchive(input1))
	{
		std::cout << "Error Opening Abc archive, aborting..." << std::endl;
		return 1;
	}

	std::cout << "Done reading tet mesh..." << std::endl;

	//Open Point Data Files
	std::vector<int> sampleFrames = { 7, 11, 15, 19, 23, 28 };
	std::vector<std::vector<Eigen::Vector3f>> samplePointsFrames;
	std::vector<std::string> samplePointsFramesFileNames =
	{ "MarkerDisp_Liver7.txt",
	"MarkerDisp_Liver11.txt",
	"MarkerDisp_Liver15.txt",
	"MarkerDisp_Liver19.txt",
	"MarkerDisp_Liver23.txt",
	"MarkerDisp_Liver28.txt" };

	std::string fileLocation = "C:/Users/Stephan/Desktop/NEW_LIVER_DATA/LiverPhantomIndentation/MarkerDisplacement/";

	samplePointsFrames.resize(samplePointsFramesFileNames.size());

	for (int i = 0; i < samplePointsFramesFileNames.size(); ++i)
	{
		std::ifstream file;
		file.open(fileLocation + samplePointsFramesFileNames[i]);

		if (!file.is_open())
		{
			std::cout << "ERROR: Could not open file: [ " << fileLocation + samplePointsFramesFileNames[i] << " ]" << std::endl;
			continue;
		}

		std::string currentLine;
		while (std::getline(file, currentLine))
		{
			std::vector<std::string> inputs;
			boost::split(inputs, currentLine, boost::is_any_of(" "));

			samplePointsFrames[i].push_back(Eigen::Vector3f(std::stof(inputs[0]), std::stof(inputs[1]), std::stof(inputs[2])));
		}
	}

	std::cout << "Done reading points..." << std::endl;

	//Find Barycentric Coordinates-------------------------------------------------------------------------------------------------

	int numTrackedPoints = samplePointsFrames[0].size();
	int numTets = reader1.getNumFaces();
	int numExcludedPoints = 0;

	std::cout << "Num Tetrahedra: " << numTets << std::endl;
	std::cout << "Num Tracked Points: " << numTrackedPoints << std::endl;

	std::vector<Eigen::Vector4f> baryCentricCoords(numTrackedPoints);
	std::vector<int> baryCentricCoordsTetIdx(numTrackedPoints);

	//determine barycentric coordinates for each initial tracked internal point
	for (int p = 0; p < numTrackedPoints; ++p)
	{
		int numMatchingTets = 0;
		for (int t = 0; t < numTets; ++t)
		{
			//std::cout << t << ", ";
			//roughly follows the method from http://dennis2society.de/main/painless-tetrahedral-barycentric-mapping

			Eigen::Vector4f point;
			point[0] = samplePointsFrames[0][p][0];
			point[1] = samplePointsFrames[0][p][1];
			point[2] = samplePointsFrames[0][p][2];
			point[3] = 1.0f;

			Eigen::Vector4f x1;
			x1[0] = reader1.getPositions()[reader1.getFaceIndices(t)[0]][0];
			x1[1] = reader1.getPositions()[reader1.getFaceIndices(t)[0]][1];
			x1[2] = reader1.getPositions()[reader1.getFaceIndices(t)[0]][2];
			x1[3] = 1.0f;

			Eigen::Vector4f x2;
			x2[0] = reader1.getPositions()[reader1.getFaceIndices(t)[1]][0];
			x2[1] = reader1.getPositions()[reader1.getFaceIndices(t)[1]][1];
			x2[2] = reader1.getPositions()[reader1.getFaceIndices(t)[1]][2];
			x2[3] = 1.0f;

			Eigen::Vector4f x3;
			x3[0] = reader1.getPositions()[reader1.getFaceIndices(t)[2]][0];
			x3[1] = reader1.getPositions()[reader1.getFaceIndices(t)[2]][1];
			x3[2] = reader1.getPositions()[reader1.getFaceIndices(t)[2]][2];
			x3[3] = 1.0f;

			Eigen::Vector4f x4;
			x4[0] = reader1.getPositions()[reader1.getFaceIndices(t)[3]][0];
			x4[1] = reader1.getPositions()[reader1.getFaceIndices(t)[3]][1];
			x4[2] = reader1.getPositions()[reader1.getFaceIndices(t)[3]][2];
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
				continue;
			}

			lbd1 = detMatrix_1.determinant() / detMatrix_0.determinant();
			lbd2 = detMatrix_2.determinant() / detMatrix_0.determinant();
			lbd3 = detMatrix_3.determinant() / detMatrix_0.determinant();
			lbd4 = detMatrix_4.determinant() / detMatrix_0.determinant();

			if ((lbd1 < 0.0f && lbd2 < 0.0f && lbd3 < 0.0f && lbd4 < 0.0f)
				|| (lbd1 > 0.0f && lbd2 > 0.0f && lbd3 > 0.0f && lbd4 > 0.0f))
			{
				isInTet = true;
			}

			if (isInTet)
			{
				//save tet idx
				baryCentricCoordsTetIdx[p] = t;

				//save barycentric coords
				baryCentricCoords[p] = Eigen::Vector4f(lbd1, lbd2, lbd3, lbd4);
				//std::cout << baryCentricCoords[p] << std::endl;

				++numMatchingTets;
			}

			//std::cout << x1 << std::endl;
		}

		if (numMatchingTets != 1)
		{
			std::cout << "ERROR: More than one tet match found! ( " << numMatchingTets << " )" << std::endl;
			baryCentricCoords[p] = Eigen::Vector4f::Zero();
			baryCentricCoordsTetIdx[p] = -1;
			numExcludedPoints += 1;
		}
		//else
		//{
		//	std::cout << "point: " << p << ": " << baryCentricCoordsTetIdx[p] << std::endl;
		//}
	}
	

	std::cout << "Finished determining barycentric coordinates, analysing mesh..." << std::endl;
	//-----------------------------------------------------------------------------------------------------------------------------
	//CHECK POINT ALIGNMENT

	std::vector<float> sumofSquaredPositionDifferences;;
	std::vector<Eigen::Vector3f> meanPositionDifferences;

	sumofSquaredPositionDifferences.resize(sampleFrames.size());
	meanPositionDifferences.resize(sampleFrames.size());
	for (int i = 0; i < sumofSquaredPositionDifferences.size(); ++i)
	{
		sumofSquaredPositionDifferences[i] = 0.0f;
		meanPositionDifferences[i].setZero();
	}

	for (int i = 0; i < sampleFrames.size(); ++i)
	{
		int currentFrame = sampleFrames[i];
		int fileFrame = (int)((float)currentFrame / 0.005f);

		if (!reader1.sampleSpecific(fileFrame))
		{
			std::cout << "ERROR: Could not sample [ " << currentFrame << " => " << fileFrame << " ]." << std::endl;
			continue;
		}
		
		float sum = 0.0f;
		Eigen::Vector3f pointSum; pointSum.setZero();

		for (int p = 0; p < numTrackedPoints; ++p)
		{
			//skip points for which we have no match
			if (baryCentricCoordsTetIdx[p] < 0)
			{
				continue;
			}

			Eigen::Vector3f groundTruthPointPosition = samplePointsFrames[i][p];

			Eigen::Vector3f simulatedPoint =
				reader1.getPositions()[reader1.getFaceIndices(baryCentricCoordsTetIdx[p])[0]] * baryCentricCoords[p][0]
				+ reader1.getPositions()[reader1.getFaceIndices(baryCentricCoordsTetIdx[p])[1]] * baryCentricCoords[p][1]
				+ reader1.getPositions()[reader1.getFaceIndices(baryCentricCoordsTetIdx[p])[2]] * baryCentricCoords[p][2]
				+ reader1.getPositions()[reader1.getFaceIndices(baryCentricCoordsTetIdx[p])[3]] * baryCentricCoords[p][3];

			sum += std::sqrtf((groundTruthPointPosition - simulatedPoint).squaredNorm());
			pointSum += groundTruthPointPosition - simulatedPoint;
		}

		sumofSquaredPositionDifferences[i] = (sum / (float)(numTrackedPoints - numExcludedPoints));
		meanPositionDifferences[i] = pointSum / (float)(numTrackedPoints - numExcludedPoints);
	}

	std::cout << "Writing results..." << std::endl;

	//Write Results for Matlab
	std::ofstream fileStream;
	fileStream.open("PhantomDeformationResults.txt");

	std::stringstream stream;

	stream << "sumSquaredPositionDifferences = [ ";
	for (int i = 0; i < sumofSquaredPositionDifferences.size(); ++i)
	{
		if (i < sumofSquaredPositionDifferences.size() - 1)
		{
			stream << sumofSquaredPositionDifferences[i] << ", ";
		}
		else
		{
			stream << sumofSquaredPositionDifferences[i] << "]; " << std::endl;
		}
	}
	stream << std::endl;
	stream << "meanSumSquaredPositionDifferencesX = [ ";
	for (int i = 0; i < meanPositionDifferences.size(); ++i)
	{
		if (i < meanPositionDifferences.size() - 1)
		{
			stream << meanPositionDifferences[i][0] << ", ";
		}
		else
		{
			stream << meanPositionDifferences[i][0] << "]; " << std::endl;
		}
	}

	stream << std::endl;
	stream << "meanSumSquaredPositionDifferencesY = [ ";
	for (int i = 0; i < meanPositionDifferences.size(); ++i)
	{
		if (i < meanPositionDifferences.size() - 1)
		{
			stream << meanPositionDifferences[i][1] << ", ";
		}
		else
		{
			stream << meanPositionDifferences[i][1] << "]; " << std::endl;
		}
	}

	stream << std::endl;
	stream << "meanSumSquaredPositionDifferencesZ = [ ";
	for (int i = 0; i < meanPositionDifferences.size(); ++i)
	{
		if (i < meanPositionDifferences.size() - 1)
		{
			stream << meanPositionDifferences[i][2] << ", ";
		}
		else
		{
			stream << meanPositionDifferences[i][2] << "]; " << std::endl;
		}
	}

	fileStream << stream.str();

	fileStream.close();
	fileStream.clear();

	std::cout << "---------------------------------------------------" << std::endl;
	std::cout << "RESULTS:" << std::endl;
	std::cout << stream.str() << std::endl;

	stream.clear();

	return 0;
}