#include "AbcReader.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
	std::string input1;
	std::string input2;

	if (argc >= 2)
	{
		input1 = std::string(argv[1]);
		input2 = std::string(argv[2]);
	}
	else
	{
		input1 = "C:/Users/Stephan/Desktop/deformedMeshFEM.abc";
		input2 = "C:/Users/Stephan/Desktop/deformedMeshFEM.abc";
	}

	std::cout << "FILE 1: " << input1 << std::endl;
	std::cout << "FILE 2: " << input2 << std::endl;

	//Read the archives
	AbcReader reader1;
	if (!reader1.openArchive(input1))
	{
		std::cout << "Error Opening first file, aborting..." << std::endl;
		return 1;
	}

	AbcReader reader2;
	if (!reader2.openArchive(input2))
	{
		std::cout << "Error Opening first file, aborting..." << std::endl;
		return 1;
	}

	if (reader1.getNumSamples() != reader2.getNumSamples())
	{
		std::cout << "WARNING: Meshes contain a different Number of Samples" << std::endl;
	}

	std::vector<float> sumofSquaredPositionDifferences;;

	std::vector<float> meanPositionDifferences;

	int meanLength = 200;

	//Compare Mesh Samples
	while (reader1.sampleForward() & reader2.sampleForward())
	{
		std::vector<Eigen::Vector3f>& positions1 = reader1.getPositions();
		std::vector<Eigen::Vector3f>& positions2 = reader2.getPositions();

		float sosd = 0.0f;

		for (int i = 0; i < positions1.size(); ++i)
		{
			sosd += (positions1[i] - positions2[i]).squaredNorm();
		}

		sumofSquaredPositionDifferences.push_back(sosd);
		meanPositionDifferences.push_back(sosd / (float)positions1.size());
	}

	//Write Results for Matlab
	std::ofstream stream;
	stream.open("ComparisonResults.txt");
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
	stream << "meanSumSquaredPositionDifferences = [ ";
	for (int i = 0; i < meanPositionDifferences.size(); ++i)
	{
		if (i < meanPositionDifferences.size() - 1)
		{
			stream << meanPositionDifferences[i] << ", ";
		}
		else
		{
			stream << meanPositionDifferences[i] << "]; " << std::endl;
		}
	}

	stream.close();
	stream.clear();

	return 0;
}
